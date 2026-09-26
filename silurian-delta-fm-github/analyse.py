#!/usr/bin/env python3
"""Plot declared input geometry or analyse actual serial Delft3D FM map output."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def _plotting():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def preview(config_path, output):
    """Inputs only: never synthesises a simulation outcome."""
    from build_cases import mesh_arrays, bed_elevation, vegetation_mask, validate
    c = json.loads(Path(config_path).read_text())
    validate(c)
    x, y, faces, _ = mesh_arrays(c)
    xc, yc = x[faces].mean(axis=1), y[faces].mean(axis=1)
    z = bed_elevation(xc, yc, c)
    nx = round(c["length_m"] / c["grid_spacing_m"])
    ny = round(c["width_m"] / c["grid_spacing_m"])
    plt = _plotting()
    fig = plt.figure(figsize=(13, 8))
    for i, case in enumerate(c["cases"]):
        ax = fig.add_subplot(2, 2, i + 1, projection="3d")
        colors = plt.get_cmap("copper")((z-z.min())/(z.max()-z.min()))
        colors[vegetation_mask(xc, yc, case, c)] = [0.08, 0.40, 0.18, 1]
        ax.plot_surface(xc.reshape(ny, nx), yc.reshape(ny, nx), z.reshape(ny, nx),
                        facecolors=colors.reshape(ny, nx, 4), rstride=1, cstride=1,
                        linewidth=0, antialiased=False)
        ax.set(xlabel="x (m)", ylabel="y (m)", zlabel="bed elevation (m)",
               title=case.replace("_", " ") + " — INPUT")
        ax.set_box_aspect((3.7, 1.2, 0.8))
        ax.set_yticks([0, c["width_m"] / 2, c["width_m"]])
        ax.set_zticks([z.min(), 0, z.max()])
        ax.zaxis.set_major_formatter(plt.FuncFormatter(lambda value, pos: f"{value:.2f}"))
        ax.view_init(25, -120)
    fig.suptitle("Prescribed initial geometry and vegetation mask\nGreen = prescribed vegetation; vertical scale exaggerated; no simulated result")
    fig.tight_layout()
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    path = output / "initial_geometry_INPUT.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    return path


def _array(variable):
    data = np.ma.asarray(variable[:])
    if np.any(np.ma.getmaskarray(data)):
        raise ValueError(f"{variable.name}: masked values require explicit domain handling")
    data = np.asarray(data)
    if not np.isfinite(data).all():
        raise ValueError(f"{variable.name}: non-finite values")
    return data


def _select(ds, candidates, purpose):
    names = [name for name in candidates if name in ds.variables]
    if not names:
        raise ValueError(f"Missing {purpose}. Supported names: {candidates}; present: {list(ds.variables)}")
    return ds.variables[names[0]]


def read_map(path, end_time=None):
    """Read documented FM fields, enforcing cell/time alignment and units."""
    from netCDF4 import Dataset
    with Dataset(path) as ds:
        bedvar = _select(ds, ("mesh2d_mor_bl", "mesh2d_flowelem_bedlevel_bl"), "time-varying bed elevation")
        # mesh2d_flowelem_bl is static initial geometry: never use it as a final bed.
        if getattr(bedvar, "units", None) != "m":
            raise ValueError("Bed elevation must have units m")
        if "time" not in bedvar.dimensions or bedvar.ndim != 2:
            raise ValueError(f"Expected time × face bed array, got {bedvar.dimensions}")
        timevar = ds.variables.get("time")
        if timevar is None or not str(getattr(timevar, "units", "")).startswith("seconds since"):
            raise ValueError("Expected native FM time in seconds since reference date")
        times = _array(timevar)
        if len(times) < 2 or np.any(np.diff(times) <= 0):
            raise ValueError("At least two strictly increasing map times are required")
        beds = np.moveaxis(_array(bedvar), bedvar.dimensions.index("time"), 0)
        if end_time is not None:
            match = np.flatnonzero(np.isclose(times, end_time, rtol=0, atol=1e-6))
            if len(match) != 1 or match[0] == 0:
                raise ValueError("Requested end time must match a stored map after the initial state")
            stop = int(match[0]) + 1
            times, beds = times[:stop], beds[:stop]
        face_dim = next(d for d in bedvar.dimensions if d != "time")
        area = _select(ds, ("mesh2d_flowelem_ba",), "flow-cell area")
        if area.dimensions != (face_dim,) or getattr(area, "units", "") not in ("m2", "m^2", "m²"):
            raise ValueError(f"Cell area dimensions/units do not match bed: {area.dimensions}")
        areas = _array(area)
        if np.any(areas <= 0):
            raise ValueError("Non-positive face areas")
        topology = ds.variables.get("mesh2d")
        if topology is None or getattr(topology, "cf_role", "") != "mesh_topology":
            raise ValueError("Expected mesh2d UGRID topology")
        coord_names = str(getattr(topology, "face_coordinates", "")).split()
        if len(coord_names) != 2:
            raise ValueError("Expected two face_coordinates declared by mesh2d")
        coords = [ds.variables[n] for n in coord_names]
        if any(v.dimensions != (face_dim,) for v in coords):
            raise ValueError("Face coordinates do not match bed dimension")
        x, y = [_array(v) for v in coords]
        return {"x": x, "y": y, "area": areas, "initial": beds[0], "final": beds[-1],
                "times_s": times, "bed_variable": bedvar.name,
                "time_units": timevar.units, "face_dimension": face_dim}


def compute_metrics(data, c, threshold):
    if threshold <= 0:
        raise ValueError("Bed-change footprint threshold must be positive")
    dz = data["final"] - data["initial"]
    area = data["area"]
    gain = np.maximum(dz, 0) * area
    loss = np.maximum(-dz, 0) * area
    downstream = data["x"] >= c["basin_start_m"]
    weights = gain * downstream
    deposited = float(weights.sum())
    footprint = downstream & (dz >= threshold)
    return {"map_start_time_s": float(data["times_s"][0]),
            "map_end_time_s": float(data["times_s"][-1]),
            "map_duration_s": float(data["times_s"][-1] - data["times_s"][0]),
            "bed_variable": data["bed_variable"],
            "positive_bed_change_bulk_volume_m3": float(gain.sum()),
            "negative_bed_change_bulk_volume_m3": float(loss.sum()),
            "net_bed_bulk_volume_change_m3": float(np.sum(dz * area)),
            "net_bed_solid_volume_change_m3": float((1 - c["bed_porosity"]) * np.sum(dz * area)),
            "maximum_aggradation_m": float(max(0, dz.max())),
            "maximum_erosion_m": float(max(0, -dz.min())),
            "downstream_deposition_footprint_m2": float(area[footprint].sum()),
            "footprint_threshold_m": threshold,
            "downstream_deposition_centroid_x_m": float(np.sum(data["x"] * weights) / deposited) if deposited > 0 else None,
            "downstream_deposition_centroid_y_m": float(np.sum(data["y"] * weights) / deposited) if deposited > 0 else None,
            "interpretation": "Positive/negative net bed changes between first and last stored maps; not gross sediment deposition/erosion histories, not a conservation residual. Footprint includes submerged deposition and is not emergent land area."}


def water_balance(history):
    """Read the native water ledger; do not infer a balance from map snapshots."""
    from netCDF4 import Dataset
    with Dataset(history) as ds:
        names = ('water_balance_total_volume', 'water_balance_volume_error')
        if any(n not in ds.variables for n in names):
            return {'status': 'unavailable', 'reason': 'Native water ledger missing'}
        variables = [ds[n] for n in names]
        if any(getattr(v, 'units', '') not in ('m3', 'm^3', 'm³') for v in variables):
            raise ValueError('Native water ledger must use cubic metres')
        if variables[0].dimensions != variables[1].dimensions:
            raise ValueError('Water volume and error dimensions do not match')
        volume, error = [_array(ds[n]) for n in names]
        if not volume.size or volume.shape != error.shape:
            raise ValueError('Empty or inconsistent native water ledger')
        scale = float(np.max(np.abs(volume)))
        if scale == 0:
            return {'status': 'unavailable', 'reason': 'Zero water-volume scale; relative diagnostic undefined',
                    'max_absolute_volume_error_m3': float(np.max(np.abs(error)))}
        relative = float(np.max(np.abs(error)) / scale)
        return {'status': 'read_from_native_history', 'total_volume_scale_m3': scale,
                'max_absolute_volume_error_m3': float(np.max(np.abs(error))),
                'final_signed_volume_error_m3': float(error.ravel()[-1]),
                'max_error_relative_to_volume': relative,
                'below_diagnostic_threshold': relative < 1e-5,
                'diagnostic_threshold': 1e-5,
                'scope': 'Raw native water ledger. The morphology update moves the bed while retaining wet-cell surface elevation; the ledger omits this geometric displacement. This is not a complete moving-bed water-conservation test or a calibration test.'}


def analyse(map_path, config_path, output, scenario=None, threshold=0.001):
    map_path = Path(map_path)
    if not map_path.is_file():
        raise FileNotFoundError(f"No native FM map output: {map_path}")
    c = json.loads(Path(config_path).read_text())
    data = read_map(map_path)
    metrics = compute_metrics(data, c, threshold)
    metrics["source_map"] = str(map_path.resolve())
    metrics["scenario"] = scenario or map_path.parent.parent.name
    metrics["config_used_for_interpretation"] = c
    history = list(map_path.parent.glob('*_his.nc'))
    metrics['water_balance'] = (water_balance(history[0]) if len(history) == 1 else
                                {'status': 'unavailable', 'reason': 'No unique serial history file'})
    metrics['water_balance']['net_bed_bulk_volume_change_m3'] = metrics['net_bed_bulk_volume_change_m3']
    metrics['water_balance']['moving_bed_note'] = 'Bed-volume change is shown separately to diagnose operator-split geometry. Adding it to the raw ledger must not be presented as proof of physical water conservation.'
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    (output / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    plt = _plotting()
    x, y, initial, final = [data[k] for k in ("x", "y", "initial", "final")]
    dz = final - initial
    fig, axes = plt.subplots(3, 1, figsize=(11, 9), constrained_layout=True)
    limits = (min(initial.min(), final.min()), max(initial.max(), final.max()))
    for ax, z, label in zip(axes, [initial, final, dz], ["First stored bed", "Last stored bed", "Net bed change"]):
        options = {"cmap": "terrain", "vmin": limits[0], "vmax": limits[1]}
        if label == "Net bed change":
            limit = max(float(np.max(np.abs(dz))), 1e-12)
            options = {"cmap": "RdBu_r", "vmin": -limit, "vmax": limit}
        artist = ax.scatter(x, y, c=z, s=9, marker="s", linewidths=0, **options)
        ax.set(xlabel="x (m)", ylabel="y (m)", title=label, aspect="equal")
        fig.colorbar(artist, ax=ax, label="m")
    fig.suptitle(f"{metrics['scenario']} — native FM output\nStored interval {metrics['map_start_time_s']:g}–{metrics['map_end_time_s']:g} s")
    fig.savefig(output / "bed_change.png", dpi=170)
    plt.close(fig)
    fig = plt.figure(figsize=(11, 6))
    ax = fig.add_subplot(projection="3d")
    surf = ax.plot_trisurf(x, y, final, cmap="terrain", linewidth=0, antialiased=False)
    ax.set(xlabel="x (m)", ylabel="y (m)", zlabel="Bed elevation (m)",
           title=f"{metrics['scenario']} — final stored FM bed\nTriangulated display; vertical scale exaggerated")
    ax.set_box_aspect((3.7, 1.2, 0.8))
    ax.set_yticks([0, c["width_m"] / 2, c["width_m"]])
    ax.view_init(25, -120)
    fig.colorbar(surf, ax=ax, shrink=0.65, label="Bed elevation (m)")
    fig.savefig(output / "bed_3d.png", dpi=170)
    plt.close(fig)
    return metrics


def main():
    root = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, default=root / "experiment.json")
    p.add_argument("--map", type=Path)
    p.add_argument("--output", type=Path)
    p.add_argument("--scenario")
    p.add_argument("--threshold", type=float, default=0.001, help="Downstream deposition footprint bed-change threshold (m)")
    p.add_argument("--preview", action="store_true")
    a = p.parse_args()
    if a.preview:
        print(preview(a.config, a.output or root / "figures"))
    elif a.map:
        output = a.output or a.map.parent.parent / "outputanalysis"
        print(json.dumps(analyse(a.map, a.config, output, a.scenario, a.threshold), indent=2))
    else:
        p.error("Choose --preview for labelled input geometry or --map for actual engine output")


if __name__ == "__main__":
    main()
