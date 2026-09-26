# Native Delft3D FM on Apple Silicon

This setup builds a pinned community macOS port of Delft3D FM using Apple Clang, GNU Fortran, Open MPI and Conan dependencies. It does not require Docker or a Deltares container account. The native build, installation and help probe completed on this Mac. The packaged verification workflow passed 361 upstream unit tests (one additional test disabled upstream) and five analytic hydrodynamic benchmarks. Project-specific sediment and vegetation verification is recorded separately in `../verification/` and `EXPERIMENT.md`.

## Runtime locations

Default durable locations, outside OneDrive:

- Python: `~/.local/share/python-envs/delft3d-build-py313/bin/python`
- Source: `~/.local/share/delft3d-workspace/source`
- Build and installation: `source/build_fm-suite_release/` and its `install/` subdirectory
- Conan packages: `~/.local/share/delft3d-workspace/conan`

**Retain the source, build and Conan trees.** This port's binaries depend on their original build/cache locations. Copying only `install/`, relocating the workspace or deleting the Conan cache can break the runtime. These directories are required runtime storage, not disposable scratch.

## Prepare or resume a build

Use macOS on Apple Silicon with Xcode Command Line Tools and Homebrew. Install native dependencies if absent:

```bash
brew install gcc cmake ninja open-mpi openblas boost libxml2 precice googletest pugixml xerces-c pkgconf libomp
```

The full suite also needs GoogleTest, pugixml and Xerces-C. `pkgconf` supplies CMake package discovery, and `libomp` supplies the Apple Clang OpenMP runtime used by FBC. Conan builds the pinned GDAL, NetCDF and PETSc dependencies; a second Homebrew GDAL installation is not required.

From this project directory, use the dedicated central Python profile:

```bash
"$HOME/.local/share/python-envs/delft3d-build-py313/bin/python" runtime/setup_mac.py
```

The script checks the existing environment and refuses incompatible versions instead of upgrading a shared profile. It pins Conan 2.32.0, dfm-tools 0.47.0, hydrolib-core 1.0.1, meshkernel 8.3.0, xugrid 0.15.3, nbformat 5.11.1, nbclient 0.11.0 and ipykernel 7.3.0. It downloads and checks the source only when the source directory is absent or empty. Existing nonempty source must carry a matching `.delft3d-source.json`; unknown contents are never overwritten. Resume uses upstream's `--keep-build` option.

On another Mac without that profile, consult its shared-environment registry and select an existing Python 3.13 interpreter. Invoke that interpreter with `runtime/setup_mac.py --bootstrap-python /absolute/path/to/python3.13`. Only when the dedicated central profile is absent will the script create it and install the pinned packages. Register that profile in the machine's environment registry. Do not create a project-local `.venv` or use Anaconda base for this specialised stack.

## Open the notebook

Select **Delft3D Python 3.13** when opening `Silurian_delta.ipynb`. Register the kernel for your dedicated central profile with:

```bash
"$HOME/.local/share/python-envs/delft3d-build-py313/bin/python" -m ipykernel install --user --name delft3d-build-py313 --display-name "Delft3D Python 3.13"
```

The first notebook cell prints the actual interpreter. Check that it matches this profile.

## Check the engine and run a model

```bash
"$HOME/.local/share/python-envs/delft3d-build-py313/bin/python" runtime/setup_mac.py --check
"$HOME/.local/share/python-envs/delft3d-build-py313/bin/python" runtime/run_fm.py /absolute/path/to/model.mdu
```

`--check` launches the installed kernel with `--help`; it does not run a simulation. The model launcher uses `dflowfm --autostartstop`, sets the model directory as its working directory, configures the library path and writes `dflowfm.log` there. Paths containing spaces are supported. `--log /absolute/path/to/run.log` selects another log; an existing selected log is replaced. Check final model time, finite outputs and conservation, then timestep/mesh sensitivity before interpreting differences. A zero exit code alone is insufficient.

The bootstrap configures dependencies first, then builds with six parallel jobs and installs the suite. GoogleTest discovery uses CMake's `PRE_TEST` mode: discovery happens when CTest starts instead of competing with compiler jobs during the build. This addresses a startup discovery timeout observed under heavy build load; it does not disable tests or change numerical acceptance criteria.

The bootstrap finishes with the help probe. The source also includes CTest and analytic verification cases. Run CTest serially from `source/build_fm-suite_release` with `ctest -C Release --parallel 1`; `PRE_TEST` still discovers and executes the tests. Their execution must be reported separately from the port author's claimed results. Vegetation and sediment parameters remain scenarios unless compared with observations.

## Source provenance and third-party licences

Source: [tdamsma/Delft3D macOS port](https://github.com/tdamsma/Delft3D/tree/01bc2f3e53e2d5319696c5ed2b8a74f5ca904220), derived from [Deltares/Delft3D](https://github.com/Deltares/Delft3D). Exact commit: `01bc2f3e53e2d5319696c5ed2b8a74f5ca904220`. The archive at `https://api.github.com/repos/tdamsma/Delft3D/tarball/01bc2f3e53e2d5319696c5ed2b8a74f5ca904220` has SHA256 `57dbc0aaf398c267a0f8c32a5b0c660c0b335d9d108b203564d2690eec9eac58`. A changed archive hash stops setup.

This is a community port, not an officially supported Deltares macOS distribution. See its [Mac build notes](https://github.com/tdamsma/Delft3D/blob/01bc2f3e53e2d5319696c5ed2b8a74f5ca904220/doc/building-macos.md) and [GNU compatibility notes](https://github.com/tdamsma/Delft3D/blob/01bc2f3e53e2d5319696c5ed2b8a74f5ca904220/doc/gnu-toolchain-support.md). Limitations include the non-relocatable installation, unavailable WAVE ESMF regridding wrapper, and omitted GNU-incompatible `io_netcdf_dll`; FM uses the ordinary Fortran `io_netcdf` target.

The engine and dependencies retain their existing licences. The source includes `LICENSE`, AGPL/GPL/LGPL texts and component-specific notices; consult these before redistributing binaries or modified source. This project does not relicense Delft3D, the port, Homebrew or Conan dependencies. Retain their notices with redistribution, and keep source identity, configuration, dependency and verification records with scientific results.
