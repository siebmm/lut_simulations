"""Collect the fixed acquisition test and generate portable publication evidence."""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import numpy as np

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
CORE=HERE.parent/'core_saft_20261101';EXT=HERE.parent/'comsol_chapter_20261101'
OUT=ROOT/'studies/laser/results/partial_acquisition_20261101'
sys.path.insert(0,str(HERE.parent/'filter_study_20261101'))
import study as s


def write(p,v):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(v,indent=2)+'\n')


def fmt(v,d=2):return '--' if v is None else f'{v:.{d}f}'


def main():
    protocol=json.loads((HERE/'protocol.json').read_text())
    extra=json.loads((HERE/'left_confirmation_protocol.json').read_text())
    jobs=[(j,OUT/j['id']) for j in protocol['requests']]+[(j,OUT/'left_confirmation'/j['id']) for j in extra['requests']]
    metrics=json.loads((OUT/'metrics.json').read_text())['rows']+json.loads((OUT/'left_confirmation/metrics.json').read_text())['rows']
    records=[]
    for dest in [CORE,EXT]:(dest/'data/joint_acquisition').mkdir(parents=True,exist_ok=True)
    for j,folder in jobs:
        with np.load(folder/'arrays.npz') as a:
            points=a['points'].reshape(-1,2);x=np.unique(points[:,0])*1000;z=np.unique(points[:,1])*1000
            xx,zz=np.meshgrid(x,z);np.testing.assert_allclose(np.c_[xx.ravel(),zz.ravel()],points*1000,atol=1e-10)
            arrays=dict(x_mm=x,depth_mm=z)
            for case in ['lof','reference','difference']:
                for field in ['amplitude','count','cf','retained','screened']:
                    arrays[case+'_'+field]=a[case+'_2_'+field].reshape(len(z),len(x))
                row=next(r for r in metrics if r['id']==j['id'] and r['case']==case)
                recomputed=s.metrics(a[case+'_2_amplitude'],points,'partial')
                for key in ['cmax_db','near_support_pct','support_pct']:
                    if row[key] is None:assert recomputed[key] is None
                    else:np.testing.assert_allclose(row[key],recomputed[key],atol=1e-10,rtol=0)
            np.savez_compressed(EXT/'data/joint_acquisition'/(j['id']+'.npz'),**arrays)
            core={k:v for k,v in arrays.items() if not k.startswith(('reference_','difference_'))}
            np.savez_compressed(CORE/'data/joint_acquisition'/(j['id']+'.npz'),**core)
        records.append(dict(id=j['id'],selection=j['selection'],metrics=[r for r in metrics if r['id']==j['id']],
                            original_array_sha256=s.sha(folder/'arrays.npz')))
    for dest in [CORE,EXT]:
        data=dest/'data/joint_acquisition'
        records_out=records if dest==EXT else [dict(r,metrics=[m for m in r['metrics'] if m['case']=='lof']) for r in records]
        write(data/'index.json',records_out)
        write(data/'protocol.json',protocol);write(data/'left_confirmation_protocol.json',extra)
        write(data/'verification.json',dict(status='passed',scientific_configurations=25,reproduction_controls=1,
              metric_cases_checked=78,all_processing='raw',dataset_fingerprint=protocol['dataset_fingerprint']))
        write(data/'sha256.json',{p.name:s.sha(p) for p in data.iterdir() if p.is_file() and p.name!='sha256.json'})
    def metric(id,case='lof'):return next(r for r in metrics if r['id']==id and r['case']==case)
    names=[('control__all','All sources / all receivers'),('fine__both_sides','Weld sources / both side surfaces'),
           ('fine__left_only','Weld sources / left side surface'),('fine__right_only','Weld sources / right side surface')]
    lines=[]
    for id,name in names:
        r=metric(id);common=r['common_support_joint']['cmax_db'] if id!='control__all' else metric('fine__both_sides')['common_support_all']['cmax_db']
        lines.append(name+' & '+fmt(r['cmax_db'])+' & '+fmt(r['near_support_pct'],1)+' & '+fmt(common)+' & '+fmt(r['peak']['distance_mm'],3)+r' \\')
    for dest,folder in [(CORE,'core_saft'),(EXT,'comsol')]:
        (dest/f'tables/{folder}/partial_joint.tex').write_text('\n'.join(lines)+'\n')
    lines=[]
    for id,name in names[1:]+[(j['id'],j['id'].replace('__',' / ').replace('_',' ')) for j in extra['requests']]:
        rs=[metric(id,c) for c in ['lof','reference','difference']]
        lines.append(name+' & '+' & '.join(fmt(r['cmax_db']) for r in rs)+' & '+fmt(rs[0]['near_support_pct'],1)+r' \\')
    (EXT/'tables/comsol/partial_joint_controls.tex').write_text('\n'.join(lines)+'\n')
    lines=[]
    for j in protocol['requests'][:18]:
        id=j['id'];path=j['selection']['path'];n=int(path.split('-')[0]);m=path.split('_')[1];label='$'+m[:n]+r'\mid '+m[n:]+'$'
        rs=[metric(id,c) for c in ['lof','reference','difference']]
        lines.append(label+' & '+' & '.join(fmt(r['cmax_db']) for r in rs)+' & '+fmt(rs[0]['near_support_pct'],1)+r' \\')
    (EXT/'tables/comsol/partial_joint_paths.tex').write_text('\n'.join(lines)+'\n')
    j=next(j for j in protocol['requests'] if j['id']=='fine__left_only')
    config=dict(j['selection'],cases=['lof'])
    session=dict(schema=4,input_fingerprint=protocol['dataset_fingerprint'],selection=config,
                 display=dict(comparison='single',quantity='amplitude',scale='db',normalization='fixed',signal='lof',
                              db_range=40,db_top=0,reference_um_s=40,show_lof=False,weld=True,positions=False,match=True))
    write(CORE/'settings/partial_1-2_SSS_weld_sources_left_receivers.json',session)
    write(OUT/'saft_partial_weld_sources_left_receivers.json',session)
    print('Verified 78 signal-case metrics; collected 25 scientific configurations plus one replay control.')


if __name__=='__main__':main()
