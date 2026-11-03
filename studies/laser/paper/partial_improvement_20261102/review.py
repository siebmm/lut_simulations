"""Native reconstruction and scoring for the paper's partial-weld figures."""
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[3]; STUDY=ROOT/'studies/laser'
OUT=STUDY/'results/partial_improvement_20261102'; DEST=OUT/'review'
CONFIRMATION=DEST/'confirmation'
sys.path.insert(0,str(HERE.parent/'filter_study_20261101'))
import study as s
s.OUT=OUT
s.LOCAL=STUDY/'build/partial_improvement_20261102'


def execute(jobs,dest):
    s.execute('partial',dict(requests=jobs),dest)
    for job in jobs:s.compact(dest/job['id'])


def locate_confirmation(job_id):
    paths=list(CONFIRMATION.glob('b*/'+job_id+'/arrays.npz'))
    if len(paths)!=1:raise RuntimeError(f'Expected one confirmation for {job_id}: {paths}')
    return paths[0]


def metric(v,points,distance):
    valid=np.isfinite(v);near=valid&(distance<=.5);back=valid&(distance>1)
    n=float(np.max(v[near])) if near.any() else 0
    b=float(np.max(v[back])) if back.any() else 0
    pk=int(np.argmax(np.where(valid,v,-np.inf))) if valid.any() else None
    return dict(cmax_db=float(20*np.log10(n/b)) if n>0 and b>0 else None,
        near_support_pct=float(valid[distance<=.5].mean()*100),support_pct=float(valid.mean()*100),
        near_um_s=n*1e6,background_um_s=b*1e6,peak_distance_mm=float(distance[pk]) if pk is not None else None,
        peak_x_mm=float(points[pk,0]*1000) if pk is not None else None,peak_depth_mm=float(points[pk,1]*1000) if pk is not None else None)
