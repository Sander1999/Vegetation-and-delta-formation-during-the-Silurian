"""Run upstream unit and analytic tests; retain reports, clean analytic scratch."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET

from run_fm import DEFAULT_WORKSPACE, runtime
from setup_mac import COMMIT


def sanitize(value, replacements):
    if isinstance(value, str):
        for original, replacement in replacements:
            value = value.replace(original, replacement)
        return value
    if isinstance(value, dict):
        return {key: sanitize(item, replacements) for key, item in value.items()}
    if isinstance(value, list):
        return [sanitize(item, replacements) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workspace", type=Path, default=DEFAULT_WORKSPACE)
    parser.add_argument("--output-dir", type=Path, required=True, help="New or empty directory for retained reports")
    parser.add_argument("--skip-unit-tests", action="store_true", help="Explicitly omit CTest; report scope remains analytic-only")
    parser.add_argument("--case-timeout", type=float, default=600.0, help="Native engine timeout per analytic case, seconds")
    args = parser.parse_args()
    if args.case_timeout <= 0:
        parser.error("--case-timeout must be positive")
    workspace = args.workspace.expanduser().resolve()
    destination = args.output_dir.expanduser().resolve()
    if destination.exists() and any(destination.iterdir()):
        parser.error("--output-dir must be new or empty; existing results are not overwritten")
    destination.mkdir(parents=True, exist_ok=True)
    source = workspace / "source"
    build = source / "build_fm-suite_release"
    suite = source / "tools/verification_cases/run_all.py"
    report = {"source_commit_expected": COMMIT,
              "started_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "analytic-only" if args.skip_unit_tests else "unit-and-analytic",
              "unit_tests": {"status": "explicitly_skipped" if args.skip_unit_tests else "not_run"},
              "analytic_tests": {"status": "not_run"}, "passed_requested_scope": False}
    with tempfile.TemporaryDirectory(prefix="delft3d-verification-") as temporary:
        scratch = Path(temporary).resolve()
        replacements = [(str(scratch), "<temporary-verification>"),
                        (str(workspace), "<delft3d-workspace>"),
                        (str(destination), "<reports>"), (str(Path.home()), "<home>")]
        try:
            binary, env = runtime(workspace)
            env.update({"PYTHONDONTWRITEBYTECODE": "1", "NUMBA_CACHE_DIR": str(scratch / "numba"),
                        "MPLCONFIGDIR": str(scratch / "matplotlib"), "XDG_CACHE_HOME": str(scratch / "cache")})
            if not suite.is_file():
                raise FileNotFoundError(f"Upstream analytic suite missing: {suite}")
            if not args.skip_unit_tests:
                ctest = shutil.which("ctest") or "/opt/homebrew/bin/ctest"
                unit_log = scratch / "ctest.log"
                junit = scratch / "ctest.xml"
                with unit_log.open("w") as output:
                    result = subprocess.run([ctest, "-C", "Release", "--parallel", "1",
                                             "--output-on-failure", "--output-junit", str(junit)],
                                            cwd=build, env=env, stdout=output, stderr=subprocess.STDOUT)
                unit = {"returncode": result.returncode, "status": "failed", "cases": []}
                if junit.is_file():
                    for case in ET.parse(junit).getroot().iter("testcase"):
                        status = "failed" if case.find("failure") is not None or case.find("error") is not None else "disabled" if case.get("status") == "disabled" else "skipped" if case.find("skipped") is not None or case.get("status") in ("notrun", "skipped") else "passed"
                        unit["cases"].append({"name": case.get("name"), "status": status, "time_s": case.get("time")})
                    if result.returncode == 0 and any(c["status"] == "passed" for c in unit["cases"]) and not any(c["status"] == "failed" for c in unit["cases"]):
                        unit["status"] = "passed"
                report["unit_tests"] = unit
            analytic_log = scratch / "analytic.log"
            analytic_json = scratch / "analytic.json"
            with analytic_log.open("w") as output:
                result = subprocess.run([sys.executable, str(suite), "--runs-root", str(scratch / "cases"),
                                         "--dflowfm-binary", str(binary), "--timeout", str(args.case_timeout),
                                         "--json-out", str(analytic_json)], cwd=source, env=env,
                                        stdout=output, stderr=subprocess.STDOUT)
            report["analytic_tests"] = {"status": "failed", "returncode": result.returncode}
            if analytic_json.is_file():
                analytic = json.loads(analytic_json.read_text())
                (destination / "analytic-results.json").write_text(json.dumps(sanitize(analytic, replacements), indent=2) + "\n")
                if result.returncode == 0 and analytic.get("summary", {}).get("all_passed") is True and len(analytic.get("cases", {})) == 5:
                    report["analytic_tests"]["status"] = "passed"
            report["passed_requested_scope"] = report["analytic_tests"]["status"] == "passed" and (args.skip_unit_tests or report["unit_tests"]["status"] == "passed")
        except Exception as exc:
            report["error"] = repr(exc)
        finally:
            # Preserve diagnostics, including failures, before deleting large model outputs.
            for log in scratch.rglob("*.log"):
                relative = log.relative_to(scratch)
                target = destination / "logs" / relative
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_text(sanitize(log.read_text(errors="replace"), replacements))
            report["finished_utc"] = datetime.now(timezone.utc).isoformat()
            report["scratch_outputs_removed"] = True
            (destination / "verification-summary.json").write_text(json.dumps(sanitize(report, replacements), indent=2) + "\n")
    print(f"Verification reports: {destination}; scope: {report['scope']}; passed: {report['passed_requested_scope']}")
    return 0 if report["passed_requested_scope"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
