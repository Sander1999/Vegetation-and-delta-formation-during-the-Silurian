#!/usr/bin/env python3
"""Generate UGRID and documented Delft3D FM input files; does not run the engine."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import netcdf_file


def bed_elevation(x, y, c):
    """Idealised sloping floodplain/channel grading into a level flooded basin."""
    mouth = c["basin_start_m"]
    floodplain = (c["outlet_stage_m"] + c["floodplain_mouth_elevation_m"]
                  + np.tan(np.deg2rad(c["tilt_deg"])) * (mouth - x))
    channel = c["channel_depth_m"] * np.exp(
        -(2 * (y - c["width_m"] / 2) / c["channel_width_m"]) ** 4)
    t = np.clip((x - mouth) / c["basin_transition_m"], 0, 1)
    blend = t * t * (3 - 2 * t)
    return (1 - blend) * (floodplain - channel) + blend * (
        c["outlet_stage_m"] - c["basin_depth_m"])


def vegetation_mask(x, y, case, c):
    delta = (x >= c["basin_start_m"] - 0.25) & (x <= c["length_m"] - 0.25)
    if case == "bare":
        return np.zeros_like(x, dtype=bool)
    if case == "half_delta":
        return delta & (y > c["width_m"] / 2)
    if case == "half_delta_mirrored":
        return delta & (y < c["width_m"] / 2)
    if case in ("river_and_delta", "river_margins"):
        offset = np.abs(y - c["width_m"] / 2)
        margins = (x < c["basin_start_m"]) & (offset > c["channel_width_m"] / 2)
        margins &= offset < 1.5 * c["channel_width_m"]
        return margins if case == "river_margins" else margins | delta
    raise ValueError(f"Unknown case {case}")


def validate(c):
    for key in ("length_m", "width_m", "grid_spacing_m", "duration_s",
                "grain_density_kg_m3", "grain_diameter_m", "max_timestep_s",
                "basin_transition_m", "channel_width_m", "initial_sediment_thickness_m"):
        if c[key] <= 0:
            raise ValueError(f"{key} must be positive")
    if not 0 < c["bed_porosity"] < 1:
        raise ValueError("bed_porosity must be between zero and one")
    if not 0 < c["basin_start_m"] < c["length_m"] - c["basin_transition_m"]:
        raise ValueError("Basin transition must lie within the domain")
    if c.get("hydrodynamic_spinup_s", 0) < 0:
        raise ValueError("hydrodynamic_spinup_s must be nonnegative")
    for length in (c["length_m"], c["width_m"]):
        if not np.isclose(length / c["grid_spacing_m"], round(length / c["grid_spacing_m"])):
            raise ValueError("Grid spacing must divide both box dimensions")


def mesh_arrays(c):
    nx, ny = [round(c[k] / c["grid_spacing_m"]) for k in ("length_m", "width_m")]
    x, y = np.meshgrid(np.linspace(0, c["length_m"], nx + 1),
                       np.linspace(0, c["width_m"], ny + 1))
    nodes = np.arange((ny + 1) * (nx + 1)).reshape(ny + 1, nx + 1)
    faces = np.stack([nodes[:-1, :-1], nodes[:-1, 1:], nodes[1:, 1:], nodes[1:, :-1]], axis=-1).reshape(-1, 4)
    horizontal = np.stack([nodes[:, :-1], nodes[:, 1:]], axis=-1).reshape(-1, 2)
    vertical = np.stack([nodes[:-1], nodes[1:]], axis=-1).reshape(-1, 2)
    edges = np.concatenate([horizontal, vertical])
    return x.ravel(), y.ravel(), faces, edges


def write_mesh(path, c):
    x, y, faces, edges = mesh_arrays(c)
    with netcdf_file(path, "w", version=2) as f:
        f.Conventions = "CF-1.8 UGRID-1.0"
        f.title = "Idealised EM4-scale vegetation delta; Cartesian metres"
        for name, size in (("mesh2d_nNodes", len(x)), ("mesh2d_nFaces", len(faces)),
                           ("mesh2d_nEdges", len(edges)), ("mesh2d_nMax_face_nodes", 4), ("Two", 2)):
            f.createDimension(name, size)
        mesh = f.createVariable("mesh2d", "i", ())
        mesh[...] = 0
        mesh.cf_role = "mesh_topology"
        mesh.topology_dimension = 2
        mesh.node_coordinates = "mesh2d_node_x mesh2d_node_y"
        mesh.face_node_connectivity = "mesh2d_face_nodes"
        mesh.edge_node_connectivity = "mesh2d_edge_nodes"
        mesh.face_coordinates = "mesh2d_face_x mesh2d_face_y"
        for name, data, dim in (("node_x", x, "mesh2d_nNodes"), ("node_y", y, "mesh2d_nNodes"),
                                ("face_x", x[faces].mean(axis=1), "mesh2d_nFaces"),
                                ("face_y", y[faces].mean(axis=1), "mesh2d_nFaces")):
            v = f.createVariable("mesh2d_" + name, "d", (dim,))
            v[:] = data
            v.units = "m"
            v.standard_name = "projection_x_coordinate" if name.endswith("x") else "projection_y_coordinate"
        for name, data, dims, role in (
            ("face_nodes", faces, ("mesh2d_nFaces", "mesh2d_nMax_face_nodes"), "face_node_connectivity"),
            ("edge_nodes", edges, ("mesh2d_nEdges", "Two"), "edge_node_connectivity")):
            v = f.createVariable("mesh2d_" + name, "i", dims)
            v[:] = data + 1
            v.start_index = 1
            v.cf_role = role
        z = f.createVariable("mesh2d_node_z", "d", ("mesh2d_nNodes",))
        z[:] = bed_elevation(x, y, c)
        z.standard_name = "altitude"
        z.units = "m"
        z.mesh = "mesh2d"
        z.location = "node"
    return x, y, faces


def text(path, content):
    path.write_text(content.strip() + "\n", encoding="utf-8")


def sample_block(quantity, filename):
    # Static FM fields use method5=triangulation. Method4 means inside-polygon
    # in timespaceinitialfield, despite differing generic forcing method tables.
    return f"QUANTITY={quantity}\nFILENAME={filename}\nFILETYPE=7\nMETHOD=5\nOPERAND=O\n"


def build_case(directory, case, c):
    directory.mkdir(parents=True, exist_ok=True)
    x, y, faces = write_mesh(directory / "flume_net.nc", c)
    # Nodes include the boundary of the interpolation hull. Face centres are added
    # so a cell-centred bed can be supplied without smoothing its channel shape.
    sx = np.concatenate([x, x[faces].mean(axis=1)])
    sy = np.concatenate([y, y[faces].mean(axis=1)])
    bed = bed_elevation(sx, sy, c)
    mask = vegetation_mask(sx, sy, case, c)
    fields = {"bedlevel": bed,
              "stemheight": np.full_like(sx, c["stem_height_m"]),
              "stemdiameter": np.full_like(sx, c["stem_diameter_m"]),
              "stemdensity": mask.astype(float) * c["stem_density_m2"]}
    oldext = []
    for name, values in fields.items():
        np.savetxt(directory / f"{name}.xyz", np.column_stack([sx, sy, values]), fmt="%.10g")
        oldext.append(sample_block(name, f"{name}.xyz"))
    np.savetxt(directory / "massbalancearea.xyz", np.column_stack([sx, sy, np.ones_like(sx)]), fmt="%.10g")
    oldext.append(sample_block("massbalanceareaflume", "massbalancearea.xyz"))
    text(directory / "spatial.ext", "\n".join(oldext))
    half = c["channel_width_m"] / 2
    boundaries = [("inlet", -1e-6, c["width_m"] / 2 - half, c["width_m"] / 2 + half,
                   "dischargebnd", c["discharge_m3_s"], "m3/s"),
                  ("outlet", c["length_m"] + 1e-6, 0, c["width_m"],
                   "waterlevelbnd", c["outlet_stage_m"], "m")]
    ext = "[General]\nfileVersion=2.01\nfileType=extForce\n"
    bc = "[General]\nfileVersion=1.01\nfileType=boundConds\n"
    for name, bx, y1, y2, quantity, value, unit in boundaries:
        text(directory / f"{name}.pli", f"{name}\n2 2\n{bx} {y1}\n{bx} {y2}")
        ext += f"\n[Boundary]\nquantity={quantity}\nlocationFile={name}.pli\nforcingFile=boundaries.bc\n"
        for point in (1, 2):
            bc += (f"\n[forcing]\nName={name}_{point:04d}\nFunction=constant\n"
                   f"Quantity={quantity}\nUnit={unit}\n{value:.12g}\n")
        # Clear water at either boundary if inflow occurs. Outflow transports the
        # interior concentration; the zero value is an inflow condition.
        ext += f"\n[Boundary]\nquantity=sedfracbndplastic\nlocationFile={name}.pli\nforcingFile=sediment_boundary.bc\n"
    text(directory / "boundaries.ext", ext)
    text(directory / "boundaries.bc", bc)
    sediment_bc = "[General]\nfileVersion=1.01\nfileType=boundConds\n"
    for name in ("inlet", "outlet"):
        for point in (1, 2):
            sediment_bc += (f"\n[forcing]\nName={name}_{point:04d}\nFunction=constant\n"
                            "Quantity=sedfracbndplastic\nUnit=kg/m3\n0\n")
    text(directory / "sediment_boundary.bc", sediment_bc)
    # Native morphological table format: prescribed zero solid-volume bedload
    # at the upstream boundary. Native applies this when that boundary is inflow.
    text(directory / "bedload_boundary.bcm", f"""
location 'inlet'
time-function 'non-equidistant'
reference-time 20220101
time-unit 'seconds'
interpolation 'linear'
parameter 'time' unit '[s]'
parameter 'transport excl pores plastic' unit '[m3/s/m]'
0 0
{max(86400, c['duration_s'] * 10)} 0
""")
    text(directory / "sediment.sed", f"""
[SedimentFileInformation]
FileVersion=02.00
[SedimentOverall]
Cref={c['grain_density_kg_m3'] * (1-c['bed_porosity'])}
[Sediment]
Name=plastic
SedTyp=sand
RhoSol={c['grain_density_kg_m3']}
SedDia={c['grain_diameter_m']}
CDryB={c['grain_density_kg_m3'] * (1-c['bed_porosity'])}
IniSedThick={c['initial_sediment_thickness_m']}
FacDSS=1.0
""")
    text(directory / "morphology.mor", f"""
[MorphologyFileInformation]
FileVersion=02.00
[Morphology]
MorFac=1
MorStt={c.get('hydrodynamic_spinup_s', 0)}
SedTransStt={c.get('hydrodynamic_spinup_s', 0)}
CmpUpdStt={c.get('hydrodynamic_spinup_s', 0)}
BedUpd=true
CmpUpd=true
DensIn=false
NeuBcSand=false
BcFil=bedload_boundary.bcm
SedThr={c['sediment_depth_threshold_m']}
Thresh=0.001
Sus=1
Bed=1
SusW=0
BedW=0
[Boundary]
Name=inlet
IBedCond=5
[Underlayer]
IUnderLyr=1
[Output]
BedLayerSedimentMass=true
BedLayerThickness=true
Concentration=true
TotalTransport=true
Taub=true
TranspType=0
""")
    text(directory / "flow.mdu", f"""
[General]
Program=D-Flow FM
Version=1.2
FileVersion=1.02
[geometry]
NetFile=flume_net.nc
BedlevType=1
# BedlevMode native default is1; HYDROLIB1.0.1 does not expose this key.
Conveyance2D=-1
Kmx=0
AngLat=0
WaterLevIni={c['outlet_stage_m']}
[physics]
UnifFrictType=1
UnifFrictCoef={c['manning_n']}
Vicouv={c['horizontal_eddy_viscosity_m2_s']}
Dicouv={c['horizontal_eddy_diffusivity_m2_s']}
Salinity=0
Temperature=0
[numerics]
Epshu={c['wet_depth_threshold_m']}
[time]
RefDate=20220101
Tunit=S
TStart=0
TStop={c['duration_s'] + c.get('hydrodynamic_spinup_s', 0)}
DtUser=1
DtInit=0.01
DtMax={c['max_timestep_s']}
[external forcing]
ExtForceFile=spatial.ext
ExtForceFileNew=boundaries.ext
[sediment]
Sedimentmodelnr=4
SedFile=sediment.sed
MorFile=morphology.mor
[veg]
Vegetationmodelnr=1
Cdveg={c['vegetation_drag_coefficient']}
[output]
OutputDir=output
MapInterval={c['output_interval_s']}
HisInterval={c['output_interval_s']}
MbaInterval={c['output_interval_s']}
MbaWriteNetCDF=1
MbaWriteCsv=1
RstInterval={c['duration_s'] + c.get('hydrodynamic_spinup_s', 0)}
# WriMap_sediment native default is1; HYDROLIB1.0.1 does not expose this key.
""")
    return {"name": case, "directory": case, "mdu": f"{case}/flow.mdu",
            "duration_s": c["duration_s"] + c.get("hydrodynamic_spinup_s", 0),
            "morphological_duration_s": c["duration_s"],
            "hydrodynamic_spinup_s": c.get("hydrodynamic_spinup_s", 0), "initial_state": "same_fresh_idealised_bed",
            "vegetated_area_fraction_nodes": float(mask[:len(x)].mean()),
            "status": "inputs_generated_not_executed"}


def build(config_path, output, spacing=None, duration=None):
    c = json.loads(Path(config_path).read_text())
    if spacing is not None:
        c["grid_spacing_m"] = spacing
    if duration is not None:
        c["duration_s"] = duration
    validate(c)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    cases = [build_case(output / case, case, c) for case in c["cases"]]
    manifest = {"configuration": c, "cases": cases,
                "execution_status": "not_executed_by_generator",
                "fresh_comparison": "All generated cases start from the same bed; half-delta masks occupy a prescribed downstream region, not an observed existing delta.",
                "historical_sequence_not_automatically_executed": [
                    "First establish flow on the fixed bed, then run the bare morphological stage and branch its full restart into bare and half-delta continuations",
                    "Fresh river-margin vegetation with flow spin-up, followed by its morphological stage and then a restart adding delta vegetation; no repeated spin-up"],
                "limitations": ["Hydraulic Baptist vegetation only; no mechanical root reinforcement",
                                "Single noncohesive plastic fraction; not calibrated to the 2022 mixture",
                                "Default D-Morphology Van Rijn transport closure is not calibrated to angular plastic",
                                "No imposed external sediment feed; finite erodible bed is the sediment source",
                                "Native engine acceptance, conservation and resolution checks are required"]}
    text(output / "manifest.json", json.dumps(manifest, indent=2))
    return manifest


def main():
    root = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=root / "experiment.json")
    parser.add_argument("--output", type=Path, default=root / "cases")
    parser.add_argument("--spacing", type=float)
    parser.add_argument("--duration", type=float, help="Active morphological duration in seconds; initial hydrodynamic spin-up is additional")
    args = parser.parse_args()
    result = build(args.config, args.output, args.spacing, args.duration)
    print(f"Generated {len(result['cases'])} cases in {args.output}; no simulation executed.")


if __name__ == "__main__":
    main()
