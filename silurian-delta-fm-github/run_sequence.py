#!/usr/bin/env python3
"""Prepare or run staged vegetation experiments using full native FM restarts.

Default: prepare inputs and a plan only. --execute runs the two independent
first stages, then their continuations. No stage is a vegetation-growth model:
planting is a prescribed intervention between runs on the unchanged grid.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
import re
import shutil

from build_cases import build_case, validate, vegetation_mask
from run_cases import run_case

STAGES = (
    ('bare_first', 'bare', None),
    ('bare_continued', 'bare', 'bare_first'),
    ('half_delta', 'half_delta', 'bare_first'),
    ('half_delta_mirrored', 'half_delta_mirrored', 'bare_first'),
    ('river_margins_first', 'river_margins', None),
    ('river_then_delta', 'river_and_delta', 'river_margins_first'),
)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, data):
    path.write_text(json.dumps(data, indent=2) + '\n')


def set_mdu(path, start, stop, restart):
    text = path.read_text()
    for key, value in [('TStart', start), ('TStop', stop), ('RstInterval', stop-start)]:
        text, count = re.subn(rf'(?im)^{key}\s*=.*$', f'{key}={value:g}', text)
        if count != 1:
            raise ValueError(f'Expected one {key} in {path}')
    if restart:
        text += '\n[restart]\nRestartFile=parent_rst.nc\nRstIgnoreBl=0\n'
    path.write_text(text)


def prepare(config, output, stage_seconds, spacing=None):
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f'{output} is not empty; select a new sequence directory')
    output.mkdir(parents=True, exist_ok=True)
    c = copy.deepcopy(config)
    c['duration_s'] = stage_seconds
    if spacing is not None:
        c['grid_spacing_m'] = spacing
    validate(c)
    spinup = float(c.get('hydrodynamic_spinup_s', 0))
    if spinup < 0 or spinup % c['output_interval_s']:
        raise ValueError('Spin-up must be nonnegative and a multiple of output_interval_s')
    if stage_seconds <= 0 or stage_seconds % c['output_interval_s']:
        raise ValueError('Stage duration must be positive and a multiple of output_interval_s')
    manifest = {
        'status': 'prepared_not_executed', 'stage_duration_s': stage_seconds,
        'hydrodynamic_spinup_s': spinup,
        'time_convention': 'First stages include hydraulic spin-up followed by stage_duration_s of morphology; restart stages have no repeated spin-up',
        'method': 'Full native restart on an identical grid; prescribed vegetation intervention',
        'interpretation': 'Scenario analogue of the essay sequence, not a calibrated reconstruction',
        'restart_source_review': {
            'commit': '01bc2f3e53e2d5319696c5ed2b8a74f5ca904220',
            'initialization': 'flow_flowinit: external spatial fields initialized before load_restart_file',
            'restored': 'unc_read_map_or_rst restores mor_bl, morft, sediment concentrations and bodsed',
            'vegetation': 'Reader does not restore rnveg/diaveg/stemheight; new spatial inputs remain',
            'scope': 'Source inspection; actual restart acceptance and continuity checked separately',
        },
        'config': c, 'stages': [],
    }
    for name, vegetation, parent in STAGES:
        directory = output / name
        build_case(directory, vegetation, c)
        start = spinup + stage_seconds if parent else 0
        stop = start + stage_seconds if parent else spinup + stage_seconds
        set_mdu(directory / 'flow.mdu', start, stop, bool(parent))
        if parent:
            # FM morphology start offsets are relative to TStart, including restarts.
            morpath = directory / 'morphology.mor'
            morphology = morpath.read_text()
            for key in ('MorStt', 'SedTransStt', 'CmpUpdStt'):
                morphology, count = re.subn(rf'(?im)^{key}\s*=.*$', f'{key}=0', morphology)
                if count != 1:
                    raise ValueError(f'Expected one {key} in {morpath}')
            morpath.write_text(morphology)
        manifest['stages'].append({
            'name': name, 'directory': name, 'vegetation': vegetation,
            'parent': parent, 'start_s': start, 'stop_s': stop,
            'hydrodynamic_spinup_s': 0 if parent else spinup,
            'morphological_duration_s': stage_seconds,
            'status': 'waiting_for_parent_restart' if parent else 'inputs_generated_not_executed',
        })
    save(output / 'sequence.json', manifest)
    return manifest


def native_restart(parent, requested_time):
    """Require morphology and hydrodynamic state; never silently use a map file."""
    import numpy as np
    from netCDF4 import Dataset
    matches = []
    for path in sorted((parent / 'output').glob('*_rst.nc')):
        with Dataset(path) as ds:
            if 'time' not in ds.variables:
                continue
            time = ds.variables['time']
            if not str(time.units).startswith('seconds since 2022-01-01'):
                raise ValueError(f'Unexpected restart time origin: {time.units}')
            values = np.asarray(time[:]).ravel()
            if len(values) != 1 or not np.isclose(values[0], requested_time, rtol=0, atol=1e-6):
                continue
            names = set(ds.variables)
            required = ('mor_bl', 'morft', 'bodsed', 's1', 'unorm', 'plastic', 'plastic_bnd')
            missing = [n for n in required if n not in names]
            if missing:
                raise ValueError(f'{path} lacks full native restart variables: {missing}')
            for name in required:
                a = np.ma.asarray(ds.variables[name][:])
                if a.count() == 0 or not np.isfinite(a.compressed()).all():
                    raise ValueError(f'{path}: invalid restart state {name}')
            matches.append(path)
    if len(matches) != 1:
        raise ValueError(f'Expected one full restart at {requested_time}s in {parent}; found {len(matches)}')
    return matches[0]


def continuity(parent, child, start, vegetation, config):
    """Verify physical state continuity and the prescribed new vegetation mask."""
    import numpy as np
    from netCDF4 import Dataset
    pmap = list((parent / 'output').glob('*_map.nc'))
    cmap = list((child / 'output').glob('*_map.nc'))
    if len(pmap) != 1 or len(cmap) != 1:
        raise ValueError('Expected one serial map per stage for continuity check')
    with Dataset(pmap[0]) as p, Dataset(cmap[0]) as c:
        pt, ct = np.asarray(p['time'][:]), np.asarray(c['time'][:])
        if not np.isclose(pt[-1], start, rtol=0, atol=1e-6) or not np.isclose(ct[0], start, rtol=0, atol=1e-6):
            raise ValueError('Maps do not expose the common restart instant')
        name = next((n for n in ('mesh2d_mor_bl', 'mesh2d_flowelem_bedlevel_bl')
                     if n in p.variables and n in c.variables), None)
        if name is None:
            raise ValueError('No dynamic bed field available to verify restart continuity')
        before, after = np.ma.asarray(p[name][-1]), np.ma.asarray(c[name][0])
        if before.shape != after.shape or not np.array_equal(np.ma.getmaskarray(before), np.ma.getmaskarray(after)):
            raise ValueError('Restart bed grid/mask mismatch')
        delta = np.abs(before-after).compressed()
        error = float(delta.max()) if delta.size else float('inf')
        # 0.1 micrometre permits float map output rounding, not a morphodynamic reset.
        if not np.isfinite(error) or error > 1e-7:
            raise ValueError(f'Restart bed discontinuity: {error} m')
        checks = {'bed_max_abs_difference_m': error, 'bed_tolerance_m': 1e-7}
        for field in ('mesh2d_bodsed', 'mesh2d_sedfrac_concentration', 'mesh2d_s1'):
            left, right = np.ma.asarray(p[field][-1]), np.ma.asarray(c[field][0])
            if left.shape != right.shape or not np.array_equal(np.ma.getmaskarray(left), np.ma.getmaskarray(right)):
                raise ValueError(f'Restart shape/mask mismatch in {field}')
            if not np.isfinite(left.compressed()).all() or not np.isfinite(right.compressed()).all():
                raise ValueError(f'Nonfinite restored state in {field}')
            # Same single-precision map representation: allow two rounding units,
            # not a physical adjustment to bed mass, concentration or water level.
            eps = np.finfo(c[field].dtype).eps
            scale = np.maximum(np.abs(left), np.abs(right))
            tolerance = 2 * eps * np.maximum(scale, np.finfo(c[field].dtype).tiny)
            difference = np.abs(left-right)
            if np.any(difference > tolerance):
                raise ValueError(f'Restart state discontinuity in {field}: {difference.max()}')
            checks[field] = {'max_abs_difference': float(difference.max()),
                             'tolerance': 'two output floating-point rounding units'}
        face_coords = c['mesh2d'].face_coordinates.split()
        x, y = np.asarray(c[face_coords[0]][:]), np.asarray(c[face_coords[1]][:])
        expected = vegetation_mask(x, y, vegetation, config) * config['stem_density_m2']
        density = np.asarray(c['mesh2d_rnveg'][0])
        if not np.allclose(density, expected, rtol=1e-6, atol=1e-6):
            raise ValueError('Continuation vegetation density differs from prescribed mask')
        for field, key in [('mesh2d_stemheight', 'stem_height_m'), ('mesh2d_diaveg', 'stem_diameter_m')]:
            # FM sets diameter to zero where no stems are present.
            active = expected > 0 if field == 'mesh2d_diaveg' else np.ones_like(expected, dtype=bool)
            if not np.allclose(np.asarray(c[field][0])[active], config[key], rtol=1e-6, atol=1e-9):
                raise ValueError(f'Continuation vegetation geometry differs: {field}')
        changed = int(np.count_nonzero(np.abs(density-np.asarray(p['mesh2d_rnveg'][-1])) > 1e-6))
        checks['vegetation'] = {'mask_matches': True, 'changed_density_cells': changed,
                               'note': 'Intervention changes stem density footprint; stem height and diameter stay fixed'}
        checks['scope'] = 'State continuity and prescribed vegetation; budgets/calibration checked separately'
        return checks


def execute(output, manifest, engine, timeout):
    output = Path(output)
    manifest['status'] = 'execution_started'
    save(output / 'sequence.json', manifest)
    try:
        for stage in manifest['stages']:
            directory = output / stage['directory']
            if stage['parent']:
                parent = output / stage['parent']
                if sha(parent / 'flume_net.nc') != sha(directory / 'flume_net.nc'):
                    raise ValueError('Restart requires identical grids')
                restart = native_restart(parent, stage['start_s'])
                destination = directory / 'parent_rst.nc'
                shutil.copy2(restart, destination)
                stage['restart'] = {'source': str(restart.relative_to(output)), 'sha256': sha(destination)}
            stage['status'] = 'started'
            save(output / 'sequence.json', manifest)
            print(f"Running {stage['name']} ({stage['start_s']}–{stage['stop_s']} s)", flush=True)
            record = run_case(directory, engine, stage['stop_s'], timeout)
            stage['run_record'] = record
            if stage['parent']:
                stage['restart_continuity'] = continuity(output / stage['parent'], directory, stage['start_s'], stage['vegetation'], manifest['config'])
            stage['status'] = 'executed_to_requested_time'
            save(output / 'sequence.json', manifest)
        manifest['status'] = 'executed_all_stages'
        manifest['verification_scope'] = 'Execution time, bed/sediment/water-level restart continuity and prescribed vegetation mask; no conservation or calibration claim'
    except Exception as exc:
        stage['status'] = 'failed'
        stage['error'] = str(exc)
        manifest['status'] = 'failed'
        raise
    finally:
        save(output / 'sequence.json', manifest)


def main():
    root = Path(__file__).resolve().parent
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--config', type=Path, default=root / 'experiment.json')
    p.add_argument('--output', type=Path, default=root / 'sequence')
    p.add_argument('--stage-seconds', type=float, default=1800)
    p.add_argument('--spacing', type=float)
    p.add_argument('--execute', action='store_true')
    p.add_argument('--engine', type=Path, default=Path.home() / '.local/share/delft3d-workspace/source/build_fm-suite_release/install/bin/dflowfm')
    p.add_argument('--timeout', type=float, default=7200)
    a = p.parse_args()
    manifest_path = a.output / 'sequence.json'
    if a.execute and manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
        if manifest['status'] != 'prepared_not_executed':
            p.error('Existing sequence has started; use a new output directory to avoid overwriting results')
        if a.stage_seconds != manifest['stage_duration_s']:
            p.error('--stage-seconds differs from the prepared sequence')
        if a.spacing is not None and a.spacing != manifest['config']['grid_spacing_m']:
            p.error('--spacing differs from the prepared sequence')
    else:
        manifest = prepare(json.loads(a.config.read_text()), a.output, a.stage_seconds, a.spacing)
    if a.execute:
        execute(a.output, manifest, a.engine, a.timeout)
    else:
        print(f'Prepared {len(STAGES)} stages in {a.output}; no engine run performed.')


if __name__ == '__main__':
    main()
