"""Export actual depth-averaged FM UGRID frames to binary VTU/PVD, without smoothing."""
from __future__ import annotations

import argparse
import base64
import json
from pathlib import Path
import struct
import xml.etree.ElementTree as ET

import numpy as np
from netCDF4 import Dataset


def array(parent, name, values, kind="Float64", components=1):
    dtype = {"Float64": "<f8", "Int64": "<i8", "UInt8": "u1"}[kind]
    data = np.ascontiguousarray(values, dtype=dtype).tobytes()
    element = ET.SubElement(parent, "DataArray", type=kind, Name=name,
                            NumberOfComponents=str(components), format="binary")
    element.text = base64.b64encode(struct.pack("<Q", len(data)) + data).decode("ascii")


def read_face(dataset, name, frame, count):
    variable = dataset.variables[name]
    selections = tuple(frame if dim == "time" else slice(None) for dim in variable.dimensions)
    result = np.ma.filled(variable[selections], np.nan).astype(float)
    if result.shape != (count,):
        raise ValueError(f"{name}: expected one value per 2D face, got {result.shape}; 3D layered output is unsupported")
    return result


def export(source: Path, target: Path, stride=1, last=False):
    if target.exists() and any(target.iterdir()):
        raise ValueError(f"Output directory must be new or empty: {target}")
    target.mkdir(parents=True, exist_ok=True)
    with Dataset(source) as dataset:
        faces = dataset.variables["mesh2d_face_nodes"]
        topology = np.ma.asarray(faces[:])
        start = int(getattr(faces, "start_index", 1))
        polygons = [np.asarray(row.compressed(), dtype=int) - start for row in topology]
        counts = np.array([len(row) for row in polygons])
        if np.any(counts < 3):
            raise ValueError("Non-polygon face encountered")
        nodes = np.concatenate(polygons)
        x = np.asarray(dataset.variables["mesh2d_node_x"][:])
        y = np.asarray(dataset.variables["mesh2d_node_y"][:])
        if np.any(nodes < 0) or np.any(nodes >= len(x)):
            raise ValueError("Invalid mesh connectivity")
        nfaces = len(polygons)
        owners = np.repeat(np.arange(nfaces), counts)
        times = np.asarray(dataset.variables["time"][:], dtype=float)
        if not len(times):
            raise ValueError("No saved timesteps")
        time_units = getattr(dataset.variables["time"], "units", "")
        if time_units.split(" since ")[0].strip().lower() not in ("s", "second", "seconds"):
            raise ValueError(f"Time must be recorded in seconds, got {time_units!r}")
        if not np.all(np.isfinite(times)) or np.any(np.diff(times) <= 0):
            raise ValueError("Source times must be finite and strictly increasing")
        bed_name = next((name for name in ("mesh2d_mor_bl", "mesh2d_flowelem_bedlevel_bl")
                         if name in dataset.variables and "time" in dataset.variables[name].dimensions), None)
        if bed_name is None:
            raise ValueError("A saved time-varying bed field is required; static bed coordinates are not an evolving bed")
        if getattr(dataset.variables[bed_name], "units", "").strip() != "m":
            raise ValueError("Bed elevation must have metre units")
        initial_bed = read_face(dataset, bed_name, 0, nfaces)
        frames = [len(times)-1] if last else sorted(set(range(0, len(times), stride)) | {len(times)-1})
        collection = ET.Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
        entries = ET.SubElement(collection, "Collection")
        fields_record = {}
        for frame in frames:
            bed = read_face(dataset, bed_name, frame, nfaces)
            if not np.all(np.isfinite(bed)):
                raise ValueError("Nonfinite bed cannot define the display geometry")
            fields = {"bed_elevation_m": bed, "bed_change_m": bed-initial_bed}
            mapping = {"water_depth_m": "mesh2d_waterdepth", "water_level_m": "mesh2d_s1",
                       "vegetation_density_m-2": "mesh2d_rnveg", "vegetation_height_m": "mesh2d_stemheight",
                       "vegetation_diameter_m": "mesh2d_diaveg"}
            for output, original in mapping.items():
                if original in dataset.variables:
                    fields[output] = read_face(dataset, original, frame, nfaces)
                    fields_record[output] = original
            if "water_level_m" in fields:
                fields["water_level_above_bed_m"] = fields["water_level_m"] - bed
            if all(n in dataset.variables for n in ("mesh2d_ucx", "mesh2d_ucy")):
                u, v = (read_face(dataset, n, frame, nfaces) for n in ("mesh2d_ucx", "mesh2d_ucy"))
                fields["depth_averaged_velocity_m_s"] = np.column_stack((u, v, np.zeros(nfaces)))
                fields["speed_m_s"] = np.hypot(u, v)
            if "mesh2d_sedfrac_concentration" in dataset.variables:
                variable = dataset.variables["mesh2d_sedfrac_concentration"]
                remaining = [d for d in variable.dimensions if d != "time"]
                data = np.ma.filled(variable[tuple(frame if d == "time" else slice(None) for d in variable.dimensions)], np.nan)
                if len(remaining) != 2 or "mesh2d_nFaces" not in remaining:
                    raise ValueError("Only depth-averaged sediment concentrations are supported")
                data = np.moveaxis(data, remaining.index("mesh2d_nFaces"), 0)
                for fraction in range(data.shape[1]):
                    fields[f"SSC_fraction_{fraction+1}_kg_m3"] = data[:, fraction]
                fields["SSC_total_kg_m3"] = np.sum(data, axis=1)
            points = np.column_stack((x[nodes], y[nodes], bed[owners]))
            root = ET.Element("VTKFile", type="UnstructuredGrid", version="1.0", byte_order="LittleEndian", header_type="UInt64")
            grid = ET.SubElement(root, "UnstructuredGrid")
            piece = ET.SubElement(grid, "Piece", NumberOfPoints=str(len(nodes)), NumberOfCells=str(nfaces))
            array(ET.SubElement(piece, "Points"), "Points", points, components=3)
            cells = ET.SubElement(piece, "Cells")
            array(cells, "connectivity", np.arange(len(nodes)), "Int64")
            array(cells, "offsets", np.cumsum(counts), "Int64")
            array(cells, "types", np.full(nfaces, 7), "UInt8")  # VTK_POLYGON
            cell_data = ET.SubElement(piece, "CellData", Scalars="bed_elevation_m")
            for name, values in fields.items():
                array(cell_data, name, values, components=values.shape[1] if values.ndim == 2 else 1)
            filename = f"frame_{frame:05d}.vtu"
            ET.ElementTree(root).write(target / filename, encoding="utf-8", xml_declaration=True)
            ET.SubElement(entries, "DataSet", timestep=str(times[frame]-times[0]), group="", part="0", file=filename)
        ET.ElementTree(collection).write(target / "delta.pvd", encoding="utf-8", xml_declaration=True)
        metadata = {"source_map_filename": source.name, "source_time_units": getattr(dataset.variables["time"], "units", "unknown"),
                    "pvd_time": "seconds relative to first saved source frame", "first_source_time": float(times[0]),
                    "frame_indices": frames, "cells": nfaces, "points": len(nodes), "fields": list(fields),
                    "source_fields": fields_record, "bed_source": bed_name,
                    "bed_change_reference": "first saved source frame, including when --last is selected",
                    "geometry": "Each original face has disconnected copies of its original XY nodes, all at that face's bed elevation. No interpolation or vertical exaggeration in exported coordinates.",
                    "velocity": "Depth-averaged horizontal components; displayed third component is zero, not a computed vertical velocity.",
                    "missing_fields": "Variables absent from source are omitted, never replaced by assumed values."}
        (target / "export_metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    return target / "delta.pvd"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("map_file", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--last", action="store_true")
    args = parser.parse_args()
    if args.stride < 1:
        parser.error("--stride must be positive")
    print(export(args.map_file.resolve(), args.output_dir.resolve(), args.stride, args.last))
