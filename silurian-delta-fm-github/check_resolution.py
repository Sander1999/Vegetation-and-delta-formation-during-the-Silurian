#!/usr/bin/env python3
"""Compare matching stored times across grid/timestep refinements; report sensitivity."""
import argparse
import json
from pathlib import Path
import numpy as np
from analyse import read_map, compute_metrics


def compare(baseline, refinement, time):
    baseline,refinement=Path(baseline),Path(refinement)
    b=json.loads((baseline/'manifest.json').read_text())
    r=json.loads((refinement/'manifest.json').read_text())
    differences={k:[b['configuration'].get(k),r['configuration'].get(k)] for k in set(b['configuration'])|set(r['configuration']) if b['configuration'].get(k)!=r['configuration'].get(k)}
    allowed={'duration_s','cases','grid_spacing_m','max_timestep_s'}
    if set(differences)-allowed:raise ValueError(f'Unmatched physical configuration: {differences}')
    base_names={e['name'] for e in b['cases']}
    rows=[]
    for entry in r['cases']:
        name=entry['name']
        if name not in base_names:raise ValueError(f'No baseline for {name}')
        records=[]
        for root,c in [(baseline,b['configuration']),(refinement,r['configuration'])]:
            maps=list((root/name/'output').glob('*_map.nc'))
            if len(maps)!=1:raise ValueError(f'Expected unique serial map in {root/name}')
            d=read_map(maps[0],end_time=time)
            records.append(compute_metrics(d,c,.001))
        keys=['positive_bed_change_bulk_volume_m3','negative_bed_change_bulk_volume_m3','downstream_deposition_footprint_m2','downstream_deposition_centroid_x_m','downstream_deposition_centroid_y_m']
        metrics={}
        for key in keys:
            bv,rv=records[0][key],records[1][key]
            metrics[key]={'baseline':bv,'refinement':rv,'difference':rv-bv if bv is not None and rv is not None else None,'relative_change':(rv-bv)/abs(bv) if bv not in (None,0) and rv is not None else None}
        rows.append({'scenario':name,'metrics':metrics})
    contrast = {}
    lookup = {row['scenario']:row['metrics'] for row in rows}
    if {'bare','half_delta'} <= set(lookup):
        for metric in ['positive_bed_change_bulk_volume_m3','negative_bed_change_bulk_volume_m3']:
            contrast[metric] = {run:100*(lookup['half_delta'][metric][run]/lookup['bare'][metric][run]-1) for run in ['baseline','refinement']}
    return {'vegetation_contrast_percent_of_bare':contrast, 'stored_time_s':time,'active_morphology_s':time-b['configuration'].get('hydrodynamic_spinup_s',0),'configuration_differences':differences,'cases':rows,'scope':'Matched early-time sensitivity. Two resolutions do not establish asymptotic convergence or calibration; report effects on both individual metrics and the vegetation contrast.','boundary_discretisation':'The nominal inlet polyline is identical. Its mesh intersection may open a different resolved width on different grids; this contributes to grid sensitivity.'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--baseline',type=Path,required=True)
    p.add_argument('--refinement',type=Path,required=True)
    p.add_argument('--time',type=float,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args(); result=compare(a.baseline,a.refinement,a.time)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
