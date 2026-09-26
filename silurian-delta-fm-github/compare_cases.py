#!/usr/bin/env python3
"""Analyse completed matched FM cases and write a comparison and budget report."""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import numpy as np
from analyse import analyse, read_map, _plotting, _array
from check_budget import check
from netCDF4 import Dataset


def wide_comparison(fields, output):
    """An additional report layout of the same final-minus-initial fields."""
    plt = _plotting()
    limit = max(float(np.max(np.abs(d['final']-d['initial']))) for _, d in fields)
    with plt.rc_context({'font.size': 12}):
        fig, axes = plt.subplots(2, 2, figsize=(11, 6), layout='constrained')
        for ax, (name, data) in zip(axes.flat, fields):
            artist = ax.scatter(data['x'], data['y'], c=data['final']-data['initial'],
                                s=10, marker='s', linewidths=0, cmap='RdBu_r', vmin=-limit, vmax=limit)
            ax.set(title=name.replace('_', ' ').capitalize(), xlabel='x (m)', ylabel='y (m)', aspect='equal')
        fig.colorbar(artist, ax=axes.ravel().tolist(), label='Bed elevation change (m)', shrink=.85)
        fig.suptitle('Matched native FM scenarios — common colour scale', fontsize=15)
        fig.savefig(output, dpi=200)
        plt.close(fig)


def compare(cases, output):
    cases, output = Path(cases), Path(output)
    manifest = json.loads((cases/'manifest.json').read_text())
    output.mkdir(parents=True, exist_ok=True)
    config = output/'configuration.json'
    config.write_text(json.dumps(manifest['configuration'], indent=2)+'\n')
    rows, fields = [], []
    for entry in manifest['cases']:
        directory = cases/entry['directory']
        record = json.loads((directory/'run_record.json').read_text())
        if record['status'] != 'executed_to_requested_time':
            raise RuntimeError(f'{directory}: run incomplete; no matched comparison')
        map_path = directory/record['map']
        metrics = analyse(map_path, config, output/entry['name'], entry['name'])
        budget = check(directory)
        (output/entry['name']/'sediment_budget.json').write_text(json.dumps(budget, indent=2)+'\n')
        fraction = budget['fractions'][0]
        with Dataset(map_path) as ds:
            depth = _array(ds['mesh2d_waterdepth']).astype(float)
            u = _array(ds['mesh2d_ucx']).astype(float)
            v = _array(ds['mesh2d_ucy']).astype(float)
            concentration = _array(ds['mesh2d_sedfrac_concentration']).astype(float)
            if not all(np.isfinite(a).all() for a in (depth,u,v,concentration)):
                raise ValueError(f'{directory}: nonfinite hydrodynamic or concentration output')
            if depth.min() < -1e-8 or concentration.min() < -1e-8:
                raise ValueError(f'{directory}: negative depth or concentration beyond 1e-8 output tolerance')
            wet = depth >= manifest['configuration']['sediment_depth_threshold_m']
            speed = np.hypot(u,v)
            froude = speed[wet]/np.sqrt(9.81*depth[wet])
            times = np.asarray(ds['time'][:])
            spinup = entry.get('hydrodynamic_spinup_s',0)
            inactive = times < spinup if spinup > 0 else times == times[0]
            beds = _array(ds['mesh2d_mor_bl'])
            before = float(np.max(np.abs(beds[inactive]-beds[0])))
            if before > 1e-7:
                raise ValueError(f'{directory}: bed changed during fixed-bed spin-up')
            physical = {'minimum_depth_m':float(depth.min()),
                        'minimum_concentration_kg_m3':float(concentration.min()),
                        'maximum_speed_m_s':float(speed.max()),
                        'maximum_froude_at_sediment_active_depth':float(froude.max()) if froude.size else None,
                        'max_bed_change_before_transport_start_m':before,
                        'scope':'Finite stored fields and nonnegative depth/concentration; extrema are scenario diagnostics, not calibration.'}
            (output/entry['name']/'field_checks.json').write_text(json.dumps(physical,indent=2)+'\n')
        rows.append({
            'scenario':entry['name'],
            'active_morphology_s':entry.get('morphological_duration_s',entry['duration_s']),
            'spinup_s':entry.get('hydrodynamic_spinup_s',0),
            'aggradation_bulk_m3':metrics['positive_bed_change_bulk_volume_m3'],
            'erosion_bulk_m3':metrics['negative_bed_change_bulk_volume_m3'],
            'deposition_footprint_m2':metrics['downstream_deposition_footprint_m2'],
            'centroid_x_m':metrics['downstream_deposition_centroid_x_m'],
            'centroid_y_m':metrics['downstream_deposition_centroid_y_m'],
            'sediment_budget_status':budget['status'],
            'sediment_max_residual_kg':fraction['max_abs_cumulative_residual_kg'],
            'sediment_relative_to_movement':fraction['relative_to_movement'],
            'sediment_inward_kg':fraction['boundary_inward_kg'],
            'zero_external_feed_verified':fraction['zero_external_feed_verified'],
            'raw_water_error_m3':metrics['water_balance'].get('max_absolute_volume_error_m3'),
        })
        fields.append((entry['name'],read_map(map_path)))
    with (output/'comparison.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    plt=_plotting()
    fig,axes=plt.subplots(len(fields),1,figsize=(11,2.7*len(fields)),layout='constrained',squeeze=False)
    limit=max(float(np.max(np.abs(d['final']-d['initial']))) for _,d in fields)
    for ax,(name,d) in zip(axes[:,0],fields):
        artist=ax.scatter(d['x'],d['y'],c=d['final']-d['initial'],s=8,marker='s',linewidths=0,cmap='RdBu_r',vmin=-limit,vmax=limit)
        ax.set(title=name.replace('_',' '),xlabel='x (m)',ylabel='y (m)',aspect='equal')
    fig.colorbar(artist,ax=axes[:,0].tolist(),label='Net bed elevation change (m)',shrink=.7)
    fig.suptitle('Matched native FM scenarios — common colour scale\nHydraulic vegetation effect; assumed flume geometry and plastic sediment')
    fig.savefig(output/'comparison.png',dpi=170);plt.close(fig)
    wide_comparison(fields, output/'comparison_wide.png')
    summary={'cases':rows,'scope':'Numerical scenario comparison; not empirical calibration. Inspect conservation and resolution sensitivity before interpretation.','all_sediment_budgets_passed':all(r['sediment_budget_status']=='passed' for r in rows),'water_limitation':'Morphological bed update retains wet-cell water-surface elevation. Raw water ledger and bed displacement are reported separately; no full moving-bed water-conservation claim.'}
    (output/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cases',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    print(json.dumps(compare(a.cases,a.output),indent=2))
