# Actual FM results in ParaView

`export_paraview.py` reads the saved 2D FM UGRID map and writes binary VTU frames plus `delta.pvd`. Use the project's numerical Python profile (NumPy and netCDF4), then ParaView's own `pvpython` for rendering:

```sh
"$HOME/.local/share/python-envs/delft3d-build-py313/bin/python" export_paraview.py /path/to/flow_map.nc /path/to/paraview_export
"/Applications/ParaView-6.1.1.app/Contents/bin/pvpython" render_paraview.py /path/to/paraview_export/delta.pvd --output /path/to/paraview_export/bed_change.png --field bed_change_m --z-scale 3
```

Replace the ParaView executable location on other machines. `--water` adds translucent water only for cells with saved water depth above 0.00001 m. `--field bed_elevation_m` renders bed elevation. `--color-limit` sets a common symmetric bed-change range for comparing cases; use the same limits, camera and vertical scale for comparisons. `--time` selects the nearest saved frame; the default is the last. Export `--stride` subsamples frames while retaining the last, or `--last` exports only the final frame.

The export preserves original face XY polygons and cell values. Every face has its own disconnected vertices at that cell's bed elevation. There is **no smoothing or interpolation between neighbouring cells**. This produces flat terraces and possible gaps between unequal neighbouring elevations; these are a display representation of cell-centred bed levels, not inferred cliffs. Coordinates remain metres with no exaggeration. The renderer applies a labelled vertical display scale only. A depth-averaged 2D simulation shown obliquely is not a 3D hydrodynamic calculation.

Available arrays include bed elevation, change relative to the first saved frame, water level/depth, horizontal depth-averaged velocity, vegetation stem density/height/diameter, and each suspended sediment fraction plus their sum in kg/m³. Missing source fields are omitted. The vector's vertical component is zero for display, not a modelled vertical velocity. Bed change is relative to the first saved map frame, which need not be the original model initial condition. PVD time is seconds since that first saved frame; the source time origin and field mapping are recorded in `export_metadata.json`. Layered 3D fields are rejected rather than silently averaged. Export requires saved time-varying bed elevation in metres and finite, strictly increasing times in seconds; static bed geometry is not substituted for evolving morphology.

For the optional water surface, cells are selected by saved depth and their disconnected vertices are displaced by saved water level minus bed elevation. Cell-to-point conversion here cannot average neighbouring cells because their vertices are disconnected. No artificial water level is introduced.

Open `delta.pvd` in ParaView and press **Apply** to inspect arrays and play saved times. The PNG renderer also saves `.pvsm` with a relative PVD filename. Load the state from its directory; after moving it, use ParaView's state-loading data-directory/file relocation controls to locate `delta.pvd`. Keep the PVD and all its VTU files together. Regenerating the state with the renderer is another portable option. The figures visualise simulation output; they do not independently establish conservation or experimental agreement.

Official documentation: [loading datasets and time series](https://docs.paraview.org/en/latest/UsersGuide/dataIngestion.html), [saving images and state](https://docs.paraview.org/en/latest/UsersGuide/savingResults.html).
