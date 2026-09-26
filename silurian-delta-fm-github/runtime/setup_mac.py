"""Build the pinned community Apple Silicon port in a durable central workspace.

Source: https://github.com/tdamsma/Delft3D/tree/weekend-ai-experiments/macos-port
The build and Conan trees are runtime dependencies: do not remove or relocate them.
Run --check for an actual installed-kernel help probe, not a simulation test.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import urllib.request

from run_fm import DEFAULT_WORKSPACE, probe

COMMIT = "01bc2f3e53e2d5319696c5ed2b8a74f5ca904220"
URL = f"https://api.github.com/repos/tdamsma/Delft3D/tarball/{COMMIT}"
SHA256 = "57dbc0aaf398c267a0f8c32a5b0c660c0b335d9d108b203564d2690eec9eac58"
REQUIREMENTS = ["conan==2.32.0", "dfm-tools==0.47.0", "hydrolib-core==1.0.1", "meshkernel==8.3.0", "xugrid==0.15.3", "nbformat==5.11.1", "nbclient==0.11.0", "ipykernel==7.3.0"]
DEFAULT_ENV = Path.home() / ".local/share/python-envs/delft3d-build-py313"
MARKER = ".delft3d-source.json"


def run(command: list[str], **kwargs) -> None:
    print("Running:", repr(command), flush=True)
    subprocess.run(command, check=True, **kwargs)


def prepare_source(workspace: Path) -> Path:
    source = workspace / "source"
    if source.exists() and any(source.iterdir()):
        marker = source / MARKER
        if not marker.is_file():
            raise RuntimeError(f"Refusing to overwrite or adopt unidentified nonempty source: {source}")
        identity = json.loads(marker.read_text())
        if identity.get("commit") != COMMIT or identity.get("archive_sha256") != SHA256:
            raise RuntimeError(f"Source identity does not match the pinned release: {marker}")
        return source
    workspace.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="delft3d-download-", dir=workspace) as temporary:
        temporary = Path(temporary)
        archive = temporary / "source.tar.gz"
        request = urllib.request.Request(URL, headers={"User-Agent": "silurian-delta-mac-setup"})
        with urllib.request.urlopen(request, timeout=120) as response, archive.open("wb") as output:
            shutil.copyfileobj(response, output)
        with archive.open("rb") as downloaded:
            digest = hashlib.file_digest(downloaded, "sha256").hexdigest()
        if digest != SHA256:
            raise RuntimeError(f"Source checksum mismatch: {digest}; expected {SHA256}")
        extraction = temporary / "unpacked"
        extraction.mkdir()
        with tarfile.open(archive) as contents:
            contents.extractall(extraction, filter="data")
        roots = list(extraction.iterdir())
        if len(roots) != 1 or not (roots[0] / "build.py").is_file():
            raise RuntimeError("Unexpected source archive layout")
        (roots[0] / MARKER).write_text(json.dumps({"commit": COMMIT, "archive_sha256": SHA256, "url": URL}, indent=2) + "\n")
        if source.exists():
            source.rmdir()  # Only the empty, previously inspected directory.
        roots[0].rename(source)
    return source


def prepare_python(environment: Path, bootstrap: Path | None) -> Path:
    python = environment / "bin/python"
    if environment.exists():
        if not python.is_file():
            raise RuntimeError(f"Refusing to replace an existing unrecognised environment: {environment}")
    else:
        if bootstrap is None:
            raise RuntimeError("Central build environment is absent. Supply --bootstrap-python /absolute/path/to/python3.13 after consulting the environment registry.")
        if not bootstrap.is_absolute():
            raise RuntimeError("--bootstrap-python must be an absolute interpreter path")
        version = subprocess.check_output([str(bootstrap), "-c", "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"], text=True).strip()
        if version != "3.13":
            raise RuntimeError("The bootstrap interpreter must be Python 3.13")
        environment.parent.mkdir(parents=True, exist_ok=True)
        run([str(bootstrap), "-m", "venv", str(environment)])
        lock = Path(__file__).resolve().with_name("requirements-lock.txt")
        if not lock.is_file():
            raise RuntimeError(f"Missing shipped dependency lock: {lock}")
        run([str(python), "-m", "pip", "install", "-r", str(lock)])
        print(f"Register {environment} in the shared Python environment registry before using it for other projects.")
    # Existing shared environments are checked, never silently upgraded.
    check = "import importlib.metadata as m; " + "; ".join(
        f"assert m.version({item.split('==')[0]!r}) == {item.split('==')[1]!r}, {('Version mismatch: ' + item)!r}"
        for item in REQUIREMENTS)
    run([str(python), "-c", check])
    run([str(python), "-c", "import conan, dfm_tools, hydrolib.core, meshkernel, xugrid"])
    run([str(python), "-m", "pip", "check"])
    return python


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, default=DEFAULT_WORKSPACE)
    parser.add_argument("--environment", type=Path, default=DEFAULT_ENV)
    parser.add_argument("--bootstrap-python", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    workspace = args.workspace.expanduser().resolve()
    if args.check:
        return probe(workspace)
    if platform.system() != "Darwin" or platform.machine() != "arm64":
        parser.error("This setup targets native Apple Silicon macOS")
    brew = shutil.which("brew") or "/opt/homebrew/bin/brew"
    dependencies = ["gcc", "cmake", "ninja", "open-mpi", "openblas", "boost", "libxml2", "precice", "googletest", "pugixml", "xerces-c", "pkgconf", "libomp"]
    for name in dependencies:
        result = subprocess.run([brew, "list", "--versions", name], capture_output=True, text=True)
        if result.returncode or not result.stdout.strip():
            raise RuntimeError("Install the documented native dependencies first: brew install " + " ".join(dependencies))
    python = prepare_python(args.environment.expanduser().absolute(), args.bootstrap_python)
    source = prepare_source(workspace)
    env = os.environ.copy()
    prefix = subprocess.check_output([brew, "--prefix"], text=True).strip()
    env["PATH"] = os.pathsep.join([str(python.parent), prefix + "/bin", env.get("PATH", "")])
    env["CONAN_HOME"] = str(workspace / "conan")
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    run([str(python), "run_conan.py", "initialize", "external"], cwd=source, env=env)
    run([str(python), "build.py", "--config", "fm-suite", "--build-type", "Release", "--build-dependencies", "--keep-build"], cwd=source, env=env)
    build = source / "build_fm-suite_release"
    # Discover GoogleTests when CTest starts, not while compiler jobs compete
    # for resources. This changes neither the test inventory nor its criteria.
    run(["cmake", "-S", str(source / "src/cmake"), "-B", str(build),
         "-DCMAKE_GTEST_DISCOVER_TESTS_DISCOVERY_MODE=PRE_TEST"], cwd=source, env=env)
    run(["cmake", "--build", str(build), "--config", "Release", "--parallel", "6"], cwd=source, env=env)
    run(["cmake", "--install", str(build), "--config", "Release"], cwd=source, env=env)
    return probe(workspace)


if __name__ == "__main__":
    raise SystemExit(main())
