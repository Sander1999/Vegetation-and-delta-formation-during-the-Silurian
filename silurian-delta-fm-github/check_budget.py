#!/usr/bin/env python3
"""Check native sediment mass accounting, including bed and boundary transport.

Requires D-Flow FM mass-balance-area NetCDF output for the entire domain.
Map bed-volume change alone cannot establish conservation. Native MBA fluxes
are masses integrated over each reporting interval (kg), not kg/s: do not
multiply them by the reporting interval or apply trapezoidal integration.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import re

import numpy as np
from netCDF4 import Dataset


def labels(var):
    a = var[:]
    if a.dtype.kind == 'S' and a.dtype.itemsize == 1:
        a = np.asarray(a).view(f'S{a.shape[-1]}').reshape(a.shape[:-1])
    return [str(x.decode() if isinstance(x, bytes) else x).strip().rstrip('\x00')
            for x in np.asarray(a).ravel()]


def values(var, expected_unit=None):
    if expected_unit and str(getattr(var, 'units', '')).strip() != expected_unit:
        raise ValueError(f'{var.name}: expected {expected_unit}, got {getattr(var, "units", None)}')
    a = np.ma.asarray(var[:])
    if np.ma.getmaskarray(a).any() or not np.isfinite(a).all():
        raise ValueError(f'{var.name} contains missing/nonfinite budget data')
    return np.asarray(a, dtype=float)


def find_budget(case):
    matches = []
    for path in sorted((case / 'output').glob('*.nc')):
        with Dataset(path) as ds:
            if {'area_id', 'flux_dir', 'balarea'} <= set(ds.variables):
                matches.append(path)
    if len(matches) != 1:
        raise ValueError(f'Expected one native mass-balance-area NetCDF, found {len(matches)}. '
                         'Enable MbaInterval, MbaWriteNetCDF and a whole-domain massbalancearea.')
    return matches[0]


def read_area_series(ds, name, area):
    var = ds.variables[name]
    if var.dimensions != ('time', 'area_id'):
        raise ValueError(f'Unexpected storage dimensions {var.dimensions} in {name}')
    return values(var, 'kg')[:, area]


def mapped_redistribution(path, fraction, end_time):
    """Lower bound on moved mass from successive spatial storage snapshots.

    Half the L1 change across bed and suspended cell inventories counts a
    transfer A→B once. Coarse snapshots miss within-interval cycling, so this
    is a resolved lower bound, not gross erosion or total transport history.
    """
    with Dataset(path) as ds:
        times = values(ds['time'])
        if len(times) < 2 or not np.isclose(times[-1], end_time, rtol=0, atol=1e-6):
            raise ValueError('Map snapshots do not cover the full budget end time')
        area = values(ds['mesh2d_flowelem_ba'], 'm2')
        names = labels(ds['sedfrac_name'])
        if names.count(fraction) != 1:
            raise ValueError(f'Cannot identify unique map sediment fraction {fraction}')
        index = names.index(fraction)
        bedvar = ds['mesh2d_bodsed']
        if bedvar.dimensions != ('time', 'mesh2d_nFaces', 'nSedTot'):
            raise ValueError(f'Unsupported bed storage dimensions: {bedvar.dimensions}')
        bed = values(bedvar, 'kg m-2')[:, :, index] * area[None, :]
        bed_l1 = float(np.abs(np.diff(bed, axis=0)).sum())
        suspended_l1 = 0.0
        if 'sussedfrac_name' in ds.variables and fraction in labels(ds['sussedfrac_name']):
            index = labels(ds['sussedfrac_name']).index(fraction)
            concvar = ds['mesh2d_sedfrac_concentration']
            if concvar.dimensions != ('time', 'nSedSus', 'mesh2d_nFaces'):
                raise ValueError('Only depth-averaged concentration maps are supported')
            concentration = values(concvar, 'kg m-3')[:, index, :]
            depth = values(ds['mesh2d_waterdepth'], 'm')
            if np.min(depth) < -1e-10:
                raise ValueError('Negative water depth in sediment storage calculation')
            suspended = concentration * np.maximum(depth, 0) * area[None, :]
            suspended_l1 = float(np.abs(np.diff(suspended, axis=0)).sum())
        return {'lower_bound_moved_mass_kg': 0.5 * (bed_l1+suspended_l1),
                'bed_cell_absolute_changes_kg': bed_l1,
                'suspended_cell_absolute_changes_kg': suspended_l1,
                'method': 'Half summed absolute cell-inventory changes between successive maps; resolves internal redistribution without normalising by inactive bed inventory.',
                'limitation': 'Output-sampled lower bound; misses transfers reversed between snapshots and inherits map precision.'}


def check(case, relative_tolerance=1e-3, absolute_tolerance=1e-8):
    case = Path(case).resolve()
    # MBA suspended storage includes MorFac weighting in native accounting.
    # This checker intentionally supports the physical-time experiment only.
    mor = (case / 'morphology.mor').read_text()
    fac = re.search(r'(?im)^MorFac\s*=\s*([^#\n]+)', mor)
    if not fac or float(fac.group(1)) != 1:
        raise ValueError('This checker requires explicit constant MorFac=1')
    path = find_budget(case)
    with Dataset(path) as ds:
        areas = labels(ds['area_id'])
        candidates = [i for i, a in enumerate(areas) if a.casefold() == 'whole model']
        if len(candidates) != 1:
            raise ValueError(f'Expected unique Whole model area, got {areas}')
        area = candidates[0]
        directions = labels(ds['flux_dir'])
        if set(directions) != {'from', 'to'}:
            raise ValueError(f'Unknown flux direction convention: {directions}')
        incoming, outgoing = directions.index('from'), directions.index('to')
        time = values(ds['time'])
        if len(time) == 0 or np.any(np.diff(time) <= 0):
            raise ValueError('Empty or non-increasing mass-balance times')
        area_value = float(values(ds['balarea'])[area])
        maps = list((case / 'output').glob('*_map.nc'))
        if len(maps) != 1:
            raise ValueError('Expected one serial map to verify whole-domain area')
        with Dataset(maps[0]) as mapds:
            maparea = float(values(mapds['mesh2d_flowelem_ba']).sum())
        if not np.isclose(area_value, maparea, rtol=1e-8, atol=1e-10):
            raise ValueError(f'MBA covers {area_value} m2 but map covers {maparea} m2')
        found = []
        for name, var in ds.variables.items():
            if name.endswith('_flux_values') and getattr(var, 'balance_area', '').strip() == areas[area]:
                match = re.fullmatch(r'area\d+(const\d+)_flux_values', name)
                if match and match.group(1) + '_bed_mass' in ds.variables:
                    found.append((name, match.group(1)))
        if not found:
            raise ValueError('No whole-domain sediment mass balance exists')
        result = {'case': str(case), 'native_budget': str(path.relative_to(case)),
                  'status': 'checked', 'time_end_s': float(time[-1]),
                  'balance_area_m2': area_value, 'fractions': [],
                  'equation': 'residual = change(suspended + bed + native bed-shortage + fluff storage) + outward mass - inward mass',
                  'native_flux_units': 'kg integrated per output interval',
                  'tolerance': {'absolute_kg': absolute_tolerance, 'relative': relative_tolerance,
                    'basis': 'Relative to actual cumulative boundary exchange or the map-resolved lower bound on internal mass redistribution, not net storage change or the large inactive initial bed. 0.1% is a declared numerical screening target; refine timestep/grid for physical convergence. Absolute 1e-8 kg avoids zero-throughput division.'},
                  'scope': 'Native discrete mass accounting, not calibration or validation of sediment closure'}
        for name, prefix in found:
            var = ds[name]
            if var.dimensions[0] != 'time' or var.dimensions[-1] != 'flux_dir':
                raise ValueError(f'Unexpected flux dimensions {var.dimensions}')
            flux = values(var, 'kg')
            names = labels(ds[name.replace('_flux_values', '_fluxes')])
            if flux.shape != (len(time), len(names), 2):
                raise ValueError('Flux data shape does not match coordinates')
            groups = [n.split(',', 1)[0].strip() for n in names]
            storage = np.array([g == 'From/to storage' for g in groups])
            boundary = np.array([g.startswith('Boundary') for g in groups])
            other = ~(storage | boundary)
            if other.any() and np.max(np.abs(flux[:, other, :]), initial=0) > absolute_tolerance:
                raise ValueError(f'Unsupported nonzero terms in finite-bed experiment: {[n for n, m in zip(names, other) if m]}')
            if not storage.any() or not boundary.any():
                raise ValueError('Missing storage or boundary terms; cannot close budget')
            change = (flux[:, storage, outgoing]-flux[:, storage, incoming]).sum(axis=1)
            inward = flux[:, boundary, incoming].sum(axis=1)
            outward = flux[:, boundary, outgoing].sum(axis=1)
            residual = change + outward - inward
            cumulative = np.cumsum(residual)
            redistribution = mapped_redistribution(maps[0], str(var.balance_quantity), float(time[-1]))
            movement = max(redistribution['lower_bound_moved_mass_kg'], float((inward+outward).sum()))
            threshold = absolute_tolerance + relative_tolerance * movement
            native = read_area_series(ds, prefix+'_balance_error', area)
            native_cum = read_area_series(ds, prefix+'_balance_cumerror', area)
            if not np.allclose(residual, native, rtol=1e-8, atol=1e-10) or not np.allclose(cumulative, native_cum, rtol=1e-8, atol=1e-10):
                raise ValueError('Independent flux arithmetic differs from native reported errors')
            components = {}
            storage_sum = np.zeros(len(time))
            for suffix in ('mass', 'bed_mass', 'bedshort_mass', 'fluff_mass'):
                key = prefix+'_'+suffix
                if key in ds.variables:
                    component = read_area_series(ds, key, area)
                    storage_sum += component
                    components[suffix] = {'end_kg': float(component[-1])}
            # MBA records end-of-period storage. Infer initial storage from its
            # first reported change; compare later differences independently.
            if len(time) > 1 and not np.allclose(np.diff(storage_sum), change[1:], rtol=1e-7, atol=1e-8):
                raise ValueError('Recorded state changes do not match the native storage terms')
            initial = float(storage_sum[0]-change[0])
            worst = float(np.max(np.abs(cumulative)))
            passed = worst <= threshold
            result['fractions'].append({'name': str(var.balance_quantity),
                'passed': passed, 'initial_storage_inferred_kg': initial,
                'final_storage_kg': float(storage_sum[-1]), 'storage_components': components,
                'boundary_inward_kg': float(inward.sum()), 'boundary_outward_kg': float(outward.sum()),
                'zero_external_feed_verified': bool(float(inward.sum()) <= absolute_tolerance),
                'storage_change_kg': float(change.sum()), 'residual_kg': float(cumulative[-1]),
                'max_abs_cumulative_residual_kg': worst, 'movement_scale_kg': movement,
                'mapped_redistribution': redistribution,
                'allowed_residual_kg': threshold,
                'relative_to_movement': worst/movement if movement > 0 else None,
                'boundary_terms': {n: {'inward_kg': float(flux[:, i, incoming].sum()),
                                      'outward_kg': float(flux[:, i, outgoing].sum())}
                                   for i, n in enumerate(names) if boundary[i]}})
        result['status'] = 'passed' if all(f['passed'] for f in result['fractions']) else 'failed'
        return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('case', type=Path)
    p.add_argument('--relative-tolerance', type=float, default=1e-3)
    p.add_argument('--absolute-tolerance-kg', type=float, default=1e-8)
    a = p.parse_args()
    if a.relative_tolerance < 0 or a.absolute_tolerance_kg < 0:
        p.error('Tolerances must be nonnegative')
    try:
        result = check(a.case, a.relative_tolerance, a.absolute_tolerance_kg)
    except (ValueError, KeyError, OSError) as exc:
        result = {'status': 'not_verified', 'reason': str(exc)}
    print(json.dumps(result, indent=2))
    if a.case.is_dir():
        (a.case/'sediment_budget.json').write_text(json.dumps(result, indent=2)+'\n')
    raise SystemExit(0 if result['status'] == 'passed' else 2)


if __name__ == '__main__':
    main()
