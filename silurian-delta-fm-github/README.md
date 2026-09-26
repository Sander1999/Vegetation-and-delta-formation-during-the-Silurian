# Vegetation and delta formation

This project extends the small second-year EMriver practical in **Effects of the Silurian greening on delta structure and delta formation** using D-Flow FM and D-Morphology. It asks how sparse, small hydraulic obstacles change water movement and the location of erosion and deposition. The inverted plastic trees in the original practical represented vegetation qualitatively. A modest response is a useful result; a large tree-like stabilisation effect is not the intended benchmark.

The essay supplies the 35 ml/s water discharge, 2° table tilt and 30-minute stages. The initial bed, water level, grain properties and vegetation dimensions are declared scenarios. The model is an experimental analogue, not a calibrated reconstruction of the original table or a Silurian landscape. The original essay is not redistributed in this repository.

## Repository contents

This source-only distribution includes the simulation scripts, an unexecuted notebook, scenario configuration, pinned Python dependencies, native Mac setup scripts, ParaView export/rendering tools and compact verification records. The report, photographs, generated inputs, simulation outputs and ParaView frames are not included. Generate them using the commands below.

The original four-case run completed 30 minutes of active morphology and produced roughly 6–7% lower net erosion and aggradation volumes in the vegetation cases. These are scenario results rather than measured Cooksonia effects. Numerical scope and limitations are documented in `docs/EXPERIMENT.md` and `verification/`.

## Run on an Apple Silicon Mac

The native engine is built from the community Mac port of Delft3D. Its pinned source and libraries live outside OneDrive, in `~/.local/share/delft3d-workspace`. The Python environment is shared at `~/.local/share/python-envs/delft3d-build-py313`. No Docker or Windows emulator is involved.

From this folder, using the explicit interpreter:

```sh
~/.local/share/python-envs/delft3d-build-py313/bin/python runtime/setup_mac.py --check
~/.local/share/python-envs/delft3d-build-py313/bin/python build_cases.py --output runs/comparison
~/.local/share/python-envs/delft3d-build-py313/bin/python run_cases.py --cases runs/comparison
~/.local/share/python-envs/delft3d-build-py313/bin/python compare_cases.py --cases runs/comparison --output results/comparison
```

Use a new output directory for each experiment. The runner refuses an existing simulation output directory rather than mixing runs. `--case bare` selects one case. Input hashes, engine hash, requested/final time and execution status are saved per case. A successful process exit alone is not a conservation or calibration check.

`Silurian_delta.ipynb` explains the parameters, equations and plotting workflow. It reads actual engine outputs when present. `analyse.py --preview` draws only the prescribed initial geometry, clearly labelled as input.

For a fresh Mac installation, see `runtime/setup_mac.py --help` and `docs/INSTALL.md`. Keep the documented source/build/dependency locations: the port's installed libraries are not yet independently relocatable.

## ParaView

`export_paraview.py` converts native UGRID results into a VTU/PVD time series. `render_paraview.py` uses the installed ParaView Python to save a 3D figure and reloadable state. Exported coordinates remain in metres; vertical exaggeration changes only the display. See [docs/PARAVIEW.md](docs/PARAVIEW.md) for commands, available fields and the depth-averaged interpretation.

## Comparisons

Each fresh case first establishes flow for 60 seconds on a fixed bed, then computes 30 minutes of active morphology. A full run therefore stops at 1,860 seconds. Continuations do not repeat the initial spin-up.

The fresh-bed comparison holds initial bed, discharge, sediment and duration constant:

| Case | Prescribed vegetation |
| --- | --- |
| `bare` | None |
| `half_delta` | One half of the downstream region |
| `half_delta_mirrored` | Opposite half, checking geometric symmetry |
| `river_and_delta` | River margins and downstream region |

The staged experiment is described in `docs/EXPERIMENT.md`. Run `run_sequence.py --output runs/sequence --execute` with the same explicit Python interpreter to compute its six stages. Without `--execute`, it only prepares the inputs. Compare vegetation additions against a bare continuation from the same restart, not merely against an earlier photograph of a younger delta.

Vegetation acts through the native hydraulic resistance model. Mechanical root reinforcement, biological growth and plant-induced sediment cohesion are not implemented. A finite initial erodible layer supplies sediment; no measured external sediment feed is available. Neither a smaller delta nor more stable channels is imposed as the answer.

## Interpretation and verification

`docs/EXPERIMENT.md` gives the mathematics and numerical checks. `docs/SOURCES.md` distinguishes the essay, manufacturer information and research literature. The compact JSON records in `verification/` document the completed numerical checks.

Bed-change figures and metrics distinguish net aggradation, erosion, downstream deposition footprint and deposition centroid. Positive bed change between two snapshots is not cumulative gross deposition. Bed volume alone is not a full sediment balance; suspended storage and boundary fluxes must also be included.

The flume's millimetre particles and shallow flow make sediment-transport closure, wetting/drying, mesh size and viscosity sensitivities important. This is a research proof of concept until those checks and suitable observations support a stronger interpretation.

## Code and source provenance

Engine: [tdamsma/Delft3D Mac port](https://github.com/tdamsma/Delft3D/tree/01bc2f3e53e2d5319696c5ed2b8a74f5ca904220), based on [Deltares/Delft3D](https://github.com/Deltares/Delft3D). Its source remains subject to the licences included in that repository. This project does not relicense or bundle the engine. The archive hash and source commit are pinned in the setup script.

Generated inputs, simulation outputs, figures and ParaView files are ignored by Git. No environment, engine binary or report is bundled. Nothing has been uploaded to GitHub.

## Checks

Run the existing lightweight project tests with the dedicated interpreter:

```sh
~/.local/share/python-envs/delft3d-build-py313/bin/python -m unittest test_project -v
```

These tests check input generation and analysis logic; they do not launch a new full simulation. Native engine verification is available through `runtime/verify_mac.py --help`.

## Redistribution

The engine and dependencies retain their upstream licences. This folder does not assign a new licence to the project-specific code. Choose an appropriate project licence before offering reuse permissions; keep the upstream source and licence references in `docs/INSTALL.md`.
