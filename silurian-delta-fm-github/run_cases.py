#!/usr/bin/env python3
"""Run generated cases with the native FM executable and check completion."""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def input_hashes(directory):
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(directory.iterdir()) if p.is_file()
            and p.suffix.lower() in {'.mdu', '.ext', '.bc', '.pli', '.xyz', '.sed', '.mor', '.bcm', '.nc'}}


def run_case(directory, engine, duration, timeout=7200):
    directory = Path(directory).resolve()
    engine = Path(engine).resolve(strict=True)
    # Never silently reuse or overwrite a previous simulation's outputs.
    output = directory / 'output'
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f'{output} is not empty. Generate a new run directory.')
    log_path = directory / 'engine.log'
    command = [str(engine), '--autostartstop', 'flow.mdu']
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '1'
    record = {'command': command, 'input_sha256': input_hashes(directory),
              'engine_sha256': hashlib.sha256(engine.read_bytes()).hexdigest(),
              'status': 'started', 'requested_duration_s': duration,
              'omp_num_threads': 1}
    start = time.monotonic()
    try:
        with log_path.open('w') as log:
            result = subprocess.run(command, cwd=directory, env=env, stdout=log,
                                    stderr=subprocess.STDOUT, timeout=timeout)
        record['returncode'] = result.returncode
        if result.returncode:
            raise RuntimeError(f'Engine returned {result.returncode}; see {log_path}')
        maps = list(output.glob('*_map.nc'))
        if len(maps) != 1:
            raise RuntimeError(f'Expected one serial map file; found {len(maps)}')
        from netCDF4 import Dataset
        import numpy as np
        with Dataset(maps[0]) as ds:
            t = ds.variables['time']
            if not str(t.units).startswith('seconds since'):
                raise ValueError(f'Unexpected time units: {t.units}')
            times = np.asarray(t[:])
            if len(times) < 2 or not np.isfinite(times).all():
                raise ValueError('Missing or invalid map time steps')
            record['last_map_time_s'] = float(times[-1])
            if abs(times[-1] - duration) > 1e-6:
                raise RuntimeError(f'Run incomplete: final map at {times[-1]} s, expected {duration}')
        record['map'] = str(maps[0].relative_to(directory))
        record['status'] = 'executed_to_requested_time'
        record['verification_scope'] = 'Execution and output time only; inspect conservation and physics separately.'
    except Exception as exc:
        record['status'] = 'failed'
        record['error'] = str(exc)
        raise
    finally:
        record['elapsed_wall_s'] = time.monotonic() - start
        (directory / 'run_record.json').write_text(json.dumps(record, indent=2) + '\n')
    return record


def main():
    root = Path(__file__).resolve().parent
    default_engine = Path.home() / '.local/share/delft3d-workspace/source/build_fm-suite_release/install/bin/dflowfm'
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cases', type=Path, default=root / 'cases')
    p.add_argument('--engine', type=Path, default=default_engine)
    p.add_argument('--case', action='append', help='Case name; repeat to select several')
    p.add_argument('--timeout', type=float, default=7200)
    a = p.parse_args()
    manifest = json.loads((a.cases / 'manifest.json').read_text())
    names = {case['name'] for case in manifest['cases']}
    if a.case and set(a.case) - names:
        p.error(f'Unknown case(s): {set(a.case) - names}')
    for case in manifest['cases']:
        if a.case and case['name'] not in a.case:
            continue
        print(f"Running {case['name']}...", flush=True)
        record = run_case(a.cases / case['directory'], a.engine,
                          case['duration_s'], a.timeout)
        print(record['status'], flush=True)


if __name__ == '__main__':
    main()
