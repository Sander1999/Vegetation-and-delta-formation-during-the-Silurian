"""Launch the native D-Flow FM kernel; a help probe is not a model validation."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import subprocess


DEFAULT_WORKSPACE = Path.home() / ".local/share/delft3d-workspace"


def runtime(workspace: Path) -> tuple[Path, dict[str, str]]:
    prefix = workspace.expanduser().resolve() / "source/build_fm-suite_release/install"
    binary = prefix / "bin/dflowfm"
    if not binary.is_file():
        raise FileNotFoundError(f"Native kernel not installed: {binary}. Run runtime/setup_mac.py first.")
    env = os.environ.copy()
    env["DYLD_LIBRARY_PATH"] = str(prefix / "lib") + (":" + env["DYLD_LIBRARY_PATH"] if env.get("DYLD_LIBRARY_PATH") else "")
    env.setdefault("OMP_NUM_THREADS", "1")
    return binary, env


def probe(workspace: Path) -> int:
    binary, env = runtime(workspace)
    result = subprocess.run([str(binary), "--help"], env=env, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=60)
    print(result.stdout)
    if result.returncode:
        raise RuntimeError(f"Kernel help probe failed with status {result.returncode}")
    print("Kernel help probe passed. This does not execute or validate a simulation.")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mdu", nargs="?", type=Path)
    parser.add_argument("--workspace", type=Path, default=DEFAULT_WORKSPACE)
    parser.add_argument("--check", action="store_true", help="Run the installed kernel's --help only")
    parser.add_argument("--log", type=Path, help="Default: <model directory>/dflowfm.log (replaced)")
    args = parser.parse_args()
    if args.check:
        if args.mdu:
            parser.error("--check does not accept a model")
        return probe(args.workspace)
    if not args.mdu:
        parser.error("provide an MDU file or use --check")
    mdu = args.mdu.expanduser().resolve(strict=True)
    if mdu.suffix.lower() != ".mdu":
        parser.error("the model must be an .mdu file")
    binary, env = runtime(args.workspace)
    log = args.log.expanduser().resolve() if args.log else mdu.parent / "dflowfm.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    command = [str(binary), "--autostartstop", mdu.name]
    print(f"Running {mdu}; log: {log}", flush=True)
    with log.open("w") as output:
        result = subprocess.run(command, cwd=mdu.parent, env=env, stdout=output, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"D-Flow FM failed with status {result.returncode}; inspect {log}")
    print("Kernel exited successfully. Check final model time and conservation before interpreting results.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
