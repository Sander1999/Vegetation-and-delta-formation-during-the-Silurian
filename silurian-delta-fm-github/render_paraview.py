"""Run with ParaView pvpython: render actual VTU/PVD cells and save reloadable state."""
import argparse
import os
from pathlib import Path
import xml.etree.ElementTree as ET

from paraview.simple import (PVDReader, GetActiveViewOrCreate, Show, ColorBy, GetColorTransferFunction,
    GetScalarBar, Text, SaveScreenshot, SaveState, Render, Threshold, CellDatatoPointData, WarpByScalar)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pvd", type=Path)
    parser.add_argument("--output", type=Path, required=True, help="PNG path; PVSM saved beside it")
    parser.add_argument("--field", default="bed_change_m")
    parser.add_argument("--label", default="Delft3D FM", help="Scenario label shown on the image")
    parser.add_argument("--z-scale", type=float, default=3.0, help="Display-only vertical exaggeration")
    parser.add_argument("--time", type=float, help="Nearest saved time in seconds; default last")
    parser.add_argument("--color-limit", type=float, help="Symmetric limits for bed change comparison, metres")
    parser.add_argument("--water", action="store_true", help="Add translucent actual wet-cell water surface")
    args = parser.parse_args()
    if args.z_scale <= 0:
        parser.error("--z-scale must be positive")
    source = args.pvd.resolve(); output = args.output.resolve(); output.parent.mkdir(parents=True, exist_ok=True)
    reader = PVDReader(registrationName="FM cell-centred bed", FileName=str(source))
    times = list(reader.TimestepValues)
    time = times[-1] if args.time is None else min(times, key=lambda t: abs(t-args.time))
    reader.UpdatePipeline(time)
    if args.field not in reader.CellData.keys():
        raise ValueError(f"Field absent: {args.field}; available: {list(reader.CellData.keys())}")
    view = GetActiveViewOrCreate("RenderView"); view.ViewSize = [1800, 1100]; view.ViewTime = time
    view.UseColorPaletteForBackground = 0; view.Background = [0.97, 0.98, 0.99]
    display = Show(reader, view); display.Representation = "Surface"; display.Scale = [1, 1, args.z_scale]
    ColorBy(display, ("CELLS", args.field)); display.SetScalarBarVisibility(view, True)
    lookup = GetColorTransferFunction(args.field)
    low, high = reader.CellData[args.field].GetRange()
    if args.field == "bed_change_m":
        limit = args.color_limit if args.color_limit is not None else max(abs(low), abs(high), 1e-12)
        lookup.ApplyPreset("Cool to Warm (Extended)", True); lookup.RescaleTransferFunction(-limit, limit)
    else:
        upper = high if high > low else low+1e-12
        lookup.RGBPoints = [low, 0.267, 0.005, 0.329,
                            (low+upper)/2, 0.128, 0.567, 0.551,
                            upper, 0.993, 0.906, 0.144]
        lookup.ColorSpace = "RGB"
    bar = GetScalarBar(lookup, view); bar.Title = args.field; bar.ComponentTitle = ""; bar.TitleColor=[0.1,0.1,0.1];bar.LabelColor=[0.1,0.1,0.1]
    bar.TitleFontSize=28;bar.LabelFontSize=26
    if args.water:
        if "water_level_above_bed_m" not in reader.CellData.keys():
            raise ValueError("Source does not include water level")
        wet=Threshold(Input=reader);wet.Scalars=["CELLS","water_depth_m"];wet.LowerThreshold=1e-5;wet.UpperThreshold=1e30
        point=CellDatatoPointData(Input=wet)
        water=WarpByScalar(Input=point);water.Scalars=["POINTS","water_level_above_bed_m"];water.UseNormal=1;water.Normal=[0,0,1]
        wd=Show(water,view);ColorBy(wd,None);wd.DiffuseColor=[0.16,0.52,0.76];wd.Opacity=0.2;wd.Scale=[1,1,args.z_scale]
    bounds=reader.GetDataInformation().GetBounds()
    cx=(bounds[0]+bounds[1])/2;cy=(bounds[2]+bounds[3])/2;cz=args.z_scale*(bounds[4]+bounds[5])/2
    span=max(bounds[1]-bounds[0],bounds[3]-bounds[2],args.z_scale*(bounds[5]-bounds[4]))
    view.CameraPosition=[cx+0.9*span,cy-1.15*span,cz+1.1*span];view.CameraFocalPoint=[cx,cy,cz];view.CameraViewUp=[0,0,1];view.CameraParallelProjection=1;view.CameraParallelScale=0.72*span
    title=Text();title.Text=f"{args.label}\nt = {time:g} s since first saved frame; vertical display scale {args.z_scale:g}x\nCell-centred 2D model shown in 3D; disconnected flat faces"
    td=Show(title,view);td.WindowLocation="Upper Left Corner";td.FontSize=28;td.Color=[0.1,0.1,0.1]
    Render(view);SaveScreenshot(str(output),view,ImageResolution=[1800,1100])
    state=output.with_suffix(".pvsm");SaveState(str(state))
    # Make the saved reader path relative to the state directory; the PVD's
    # frame paths are already relative. Load state from this directory or use
    # ParaView's data-directory relocation dialog after moving the package.
    tree=ET.parse(state)
    for element in tree.iter():
        for key,value in list(element.attrib.items()):
            if value == str(source):element.set(key,os.path.relpath(source,state.parent))
    tree.write(state,encoding="utf-8",xml_declaration=True)
    print(f"Rendered {output}; state {state}; cells {reader.GetDataInformation().GetNumberOfCells()}")


if __name__ == "__main__":main()
