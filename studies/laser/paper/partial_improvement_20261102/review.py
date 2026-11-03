"""Broader partial-weld review: native group statistics and controlled settings.

The production engine is unchanged. Combining disjoint receiver-group sums
before computing coherence exactly reproduces a joint-aperture reconstruction.
"""

# Locate portable scientific inputs without depending on the source checkout.
import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import data_root as _data_root, native_root as _native_root, executable as _executable, native_env as _native_env
import argparse
import csv
import itertools
import json
import os
import subprocess
from pathlib import Path
import numpy as np
from run import HERE, ROOT, STUDY, OUT, s

DEST=OUT/'review'
EXE=_executable('SAFTReview')
GROUPS={'left':[0],'weld':[1],'right':[2],'left_weld':[0,1],
        'sides':[0,2],'weld_right':[1,2],'all':[0,1,2]}
VIEWS={f'source{i}':([i],1) for i in range(5)}
VIEWS.update(weld_min1=([1,2,3],1),weld_min2=([1,2,3],2),all_min1=([0,1,2,3,4],1),all_min2=([0,1,2,3,4],2))


def base():
    c=json.loads((OUT/'admission_audit/saft_partial_weighted_1_5_angular15.json').read_text())['selection']
    c.update(sources=list(range(5)),receivers=list(range(331)),groups=[True]*3,
             spacing_mm=.2,minimum_sources=1,cases=['lof','reference','difference'])
    return c


def protocol():
    cat=json.loads((STUDY/'results/partial_geometry_paths_20261030/verification/native/path_catalogue.json').read_text())
    assert len(cat['paths'])==324
    jobs=[]
    for i,path in enumerate(cat['paths']):
        for g in range(3):
            jobs.append(dict(id=f'p{i:03d}g{g}',path_index=i,group=g,
                selection=dict(base(),path=path['id'],groups=[k==g for k in range(3)])))
    p=dict(purpose='Review partial-weld conclusions across the implemented path catalogue, physical receiver groups and fixed source subsets.',
        catalogue=cat,requests=jobs,receiver_group_combinations=GROUPS,source_views=VIEWS,
        fixed='Weighted 1-5 MHz, vy, 15% launch/receive gate, 5% wall gate, 1.5 mm source-receiver exclusion, direct LSR screening, 25% retained fraction, 8 real/effective receivers, CF exponent 1, RMS sources. 0.2 mm grid over x +/-7.5 mm and depth 16.65-25.25 mm.',
        recombination='Export independent native sum, power, weight, squared weight, accepted/eligible counts for each source and each of three disjoint receiver groups. Sum groups before applying count/effective gates, coherence and source RMS. No choice of the best group/source/path at an individual pixel.',
        cases=['lof','reference','difference'],scope='972 native configurations; 324 paths x 7 receiver-group combinations x 9 fixed source views. Includes 132 conventional, 48 bevel-reflection and 144 corner-timing hypotheses. Corner amplitudes/phases remain uncalibrated.',
        candidate_rule='Explore LOF contrast, support, pristine ratio and global peak location together. Finer and wider checks are required before promoting a new candidate. Report failures and low absolute amplitudes; this is retrospective, not blind validation.',
        dataset_fingerprint=s.dataset('partial')['input_fingerprint'],
        source_hashes={n:s.sha(_native_root()/n) for n in ('saft.cpp','saft.hpp','screening.cpp','dataset_geometry.cpp','groove_paths.cpp')})
    path=DEST/'protocol.json'
    from standalone import frozen_protocol
    return frozen_protocol(path,p,'source_hashes')


def execute(jobs,dest,exporter=True):
    dest.mkdir(parents=True,exist_ok=True)
    req=dest/'request.json';s.write(req,dict(requests=jobs))
    env=_native_env()
    with (dest/'compute.log').open('w') as log:
        subprocess.run([str(EXE if exporter else s.EXE),str(ROOT),str(_data_root()/'partial_weld'),str(req),str(dest)],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    for j in jobs:s.compact(dest/j['id'])


def locate(job_id):
    paths=list((DEST/'atlas').glob('b*/'+job_id+'/arrays.npz'))
    if len(paths)!=1:raise RuntimeError(f'Expected one result for {job_id}: {paths}')
    return paths[0]


def reconstruct(stats,sums,ids,minimum,c):
    accepted,eligible,w,w2,power,_,_=np.moveaxis(stats[ids],-1,0)
    total=sums[ids]
    eff=np.divide(w*w,w2,out=np.zeros_like(w),where=w2>0)
    valid=(accepted>=np.maximum(c['minimum_receivers'],np.ceil(c['receiver_fraction']*eligible)))&(eff>=c['minimum_effective_receivers']-1e-10)
    mean=np.divide(total,w,out=np.zeros_like(total),where=w>0)
    if c['coherence']:
        cf=np.clip(np.divide(abs(total)**2,w*power,out=np.zeros_like(w),where=w*power>0),0,1)
        mean*=cf**c['coherence_power']
    count=valid.sum(axis=0)
    if c['combination']=='rms':
        amplitude=np.sqrt(np.sum(np.where(valid,abs(mean)**2,0),axis=0)/np.maximum(count,1))
    else:
        amplitude=abs(np.sum(np.where(valid,mean,0),axis=0)/np.maximum(count,1))
    amplitude[count<minimum]=np.nan
    return amplitude,count,valid,mean


def stats_at(path_index):
    groups=[]
    for g in range(3):
        with np.load(locate(f'p{path_index:03d}g{g}')) as f:
            a=dict(f)
        groups.append(a)
    return groups


def verify():
    p=protocol();chosen=['1-2_SSS','1-2_SLL','3-2_SSLSS','1-1_LL']
    indices=[next(i for i,x in enumerate(p['catalogue']['paths']) if x['id']==name) for name in chosen]
    missing=[j for j in p['requests'] if j['path_index'] in indices and not list((DEST/'atlas').glob('b*/'+j['id']+'/arrays.npz'))]
    if missing:execute(missing,DEST/'atlas/bverify')
    jobs=[]
    for i in indices:
        for group,view in [('all','all_min2'),('sides','weld_min2'),('right','source3'),('left_weld','source1')]:
            ids,k=VIEWS[view]
            jobs.append(dict(id=f'p{i:03d}_{group}_{view}',path_index=i,group=group,view=view,
                selection=dict(base(),path=p['catalogue']['paths'][i]['id'],groups=[g in GROUPS[group] for g in range(3)],sources=ids,minimum_sources=k)))
    dest=DEST/'verification';pending=[j for j in jobs if not (dest/j['id']/'arrays.npz').exists()]
    if pending:execute(pending,dest,False)
    checks=[]
    for i in indices:
        groups=stats_at(i);n=len(groups[0]['points'])//2
        for job in [j for j in jobs if j['path_index']==i]:
            with np.load(dest/job['id']/'arrays.npz') as native:
                for case in p['cases']:
                    stats=sum(groups[g][case+'_2_stats'].reshape(5,n,7) for g in GROUPS[job['group']])
                    sums=sum(groups[g][case+'_2_sums'].reshape(5,n) for g in GROUPS[job['group']])
                    ids,k=VIEWS[job['view']];v,count,valid,means=reconstruct(stats,sums,ids,k,job['selection'])
                    np.testing.assert_allclose(v,native[case+'_2_amplitude'],rtol=2e-8,atol=2e-15,equal_nan=True)
                    np.testing.assert_array_equal(count,native[case+'_2_count'])
                    np.testing.assert_array_equal(valid.ravel(),native[case+'_2_valid'])
                    np.testing.assert_allclose(means.ravel(),native[case+'_2_means'],rtol=2e-8,atol=2e-15)
                    checks.append(dict(id=job['id'],case=case,max_error=float(np.nanmax(abs(v-native[case+'_2_amplitude']))) if np.isfinite(v).any() else 0))
    s.write(DEST/'verification.json',dict(status='passed',native_checks=len(jobs),signal_case_checks=len(checks),checks=checks,exporter_sha256=s.sha(EXE)))
    print('Group/source recombination verified against',len(jobs),'native reconstructions.',flush=True)


def atlas():
    p=protocol();assert json.loads((DEST/'verification.json').read_text())['status']=='passed'
    existing={f.parent.name for f in (DEST/'atlas').glob('b*/*/arrays.npz')}
    pending=[j for j in p['requests'] if j['id'] not in existing]
    offset=len(list((DEST/'atlas').glob('b*')))
    for start in range(0,len(pending),12):
        batch=pending[start:start+12]
        execute(batch,DEST/f'atlas/b{offset+start//12:03d}')
        print('atlas',len(existing)+start+len(batch),'/',len(p['requests']),flush=True)
    analyse()


def metric(v,points,distance):
    valid=np.isfinite(v);near=valid&(distance<=.5);back=valid&(distance>1)
    n=float(np.max(v[near])) if near.any() else 0
    b=float(np.max(v[back])) if back.any() else 0
    pk=int(np.argmax(np.where(valid,v,-np.inf))) if valid.any() else None
    return dict(cmax_db=float(20*np.log10(n/b)) if n>0 and b>0 else None,
        near_support_pct=float(valid[distance<=.5].mean()*100),support_pct=float(valid.mean()*100),
        near_um_s=n*1e6,background_um_s=b*1e6,peak_distance_mm=float(distance[pk]) if pk is not None else None,
        peak_x_mm=float(points[pk,0]*1000) if pk is not None else None,peak_depth_mm=float(points[pk,1]*1000) if pk is not None else None)


def analyse():
    p=protocol();c=base();rows=[];worst=0;validity_checks=0
    for i,path in enumerate(p['catalogue']['paths']):
        groups=stats_at(i);points=groups[0]['points'].reshape(-1,2);n=len(points);d=s.distances(points,'partial')
        for group,gs in GROUPS.items():
            values={}
            for case in p['cases']:
                stats=sum(groups[g][case+'_2_stats'].reshape(5,n,7) for g in gs)
                sums=sum(groups[g][case+'_2_sums'].reshape(5,n) for g in gs)
                for view,(ids,k) in VIEWS.items():
                    v,count,valid,means=reconstruct(stats,sums,ids,k,c)
                    values[case,view]=(v,metric(v,points,d))
                    if len(gs)==1 and view=='all_min1':
                        original=groups[gs[0]][case+'_2_amplitude']
                        np.testing.assert_allclose(v,original,rtol=2e-8,atol=2e-15,equal_nan=True)
                        if np.isfinite(v).any():worst=max(worst,float(np.nanmax(abs(v-original))))
                        validity_checks+=1
            for view in VIEWS:
                row=dict(path=path['id'],family=path['family'] or 'conventional',group=group,view=view)
                for case in p['cases']:
                    row.update({case+'_'+k:v for k,v in values[case,view][1].items()})
                a,b=values['lof',view][0],values['reference',view][0]
                common=(d<=.5)&np.isfinite(a)&np.isfinite(b)
                row['lof_over_pristine_near_db']=float(20*np.log10(a[common].max()/b[common].max())) if common.any() and b[common].max()>0 and a[common].max()>0 else None
                rows.append(row)
        if i%30==0:print('Analysed',i+1,'paths',flush=True)
    with (DEST/'atlas_metrics.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    ranked=[r for r in rows if r['lof_cmax_db'] is not None and r['lof_near_support_pct']>=20]
    ranked.sort(key=lambda r:r['lof_cmax_db'],reverse=True)
    s.write(DEST/'atlas_summary.json',dict(rows=len(rows),native_self_checks=validity_checks,max_self_error=worst,
        best=ranked[:150],positive_localised=[r for r in ranked if r['lof_cmax_db']>0 and r['lof_peak_distance_mm']<=.5]))
    for r in ranked[:30]:print(r['path'],r['group'],r['view'],round(r['lof_cmax_db'],2),round(r['lof_near_support_pct'],1),round(r['lof_peak_distance_mm'],3),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','verify','atlas','analyse']);cmd=ap.parse_args().command
    if cmd=='prepare':print('Prepared',len(protocol()['requests']),'native group reconstructions.')
    else:globals()[cmd]()
