"""Explore fixed weld-source/right-receiver apertures and other native paths.

Per-source complex means permit exact, fixed source-subset comparisons without
repeating propagation calculations. All paths, controls and failed results are
retained. No flaw coordinates enter reconstruction; metrics use the known flaw.
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
import itertools
import json
import csv
import numpy as np
from run import HERE, OUT, STUDY, s

DEST=OUT/'right_paths'
SUBSETS={'left_weld':([0],1),'centre_weld':([1],1),'right_weld':([2],1),
         'three_min1':([0,1,2],1),'three_min2':([0,1,2],2)}


def base():
    c=json.loads((OUT/'admission_audit/saft_partial_weighted_1_5_angular15.json').read_text())['selection']
    c.update(sources=[1,2,3],groups=[False,False,True],minimum_sources=1,
             spacing_mm=.2,cases=['lof','reference','difference'])
    return c


def paths():
    result=[]
    for before,after in [(1,1),(2,1),(1,2),(2,2),(1,3),(3,1),(2,3),(3,2)]:
        for modes in itertools.product('LS',repeat=before+after):
            result.append(f'{before}-{after}_'+''.join(modes))
    assert len(result)==132
    return result


def prepare():
    receivers=s.dataset('partial')['audit']['receivers']
    windows={name:[i for i,r in enumerate(receivers) if r['group']=='right' and lo<=r['x_m']*1000<hi]
             for name,lo,hi in [('all_right',0,100),('right_near',0,15),('right_middle',15,35),('right_far',35,100)]}
    jobs=[]
    for window,ids in windows.items():
        for path in paths():
            # Every conventional route for the full right-side aperture;
            # direct and single-reflection routes for each physical subaperture.
            if window!='all_right' and sum(map(int,path.split('_')[0].split('-')))>3:
                continue
            jobs.append(dict(id=window+'__'+path,window=window,
                             selection=dict(base(),path=path,receivers=ids)))
    p=dict(question='Do individual weld sources and right-side receiver apertures reveal useful alternative propagation paths?',
        scope='132 native conventional paths for all 157 right receivers, including mode conversions and up to five total ray legs; 20 direct/one-reflection routes in each of three fixed right-side windows. 192 native configurations.',
        source_views='Each weld source separately, plus the three-source RMS requiring one or two sources. The subset is fixed over the entire image, not selected per pixel.',
        fixed='Weighted 1-5 MHz, 15% joint angular gate, native geometry/reflection coefficients, 5% wall gate, CF exponent 1, 25% retained fraction, 8 accepted/effective receivers, LSR pulse screening, original complete filtering support, vy. Grid 0.2 mm over the same 15 x 8.6 mm field.',
        controls='LOF, pristine and waveform difference for every configuration. Preserve plane 1 (no pulse screening or CF) as a diagnostic, not as a like-for-like replacement.',
        selection='Known-flaw exploratory development. Rank fixed source/path/aperture choices; report coverage and off-flaw peak location with contrast. Confirm promising non-S|SS examples on finer and wider fields and pristine data before recommending.',
        limits='Conventional paths use the existing ray model; this is not a test of all physical wave mechanisms. No new COMSOL solve or acquisition.',
        windows=windows,requests=jobs,dataset_fingerprint=s.dataset('partial')['input_fingerprint'],
        native_sources={n:s.sha(_native_root()/n) for n in ('saft.cpp','saft.hpp','screening.cpp','dataset_geometry.cpp','groove_paths.cpp')})
    path=DEST/'protocol.json'
    from standalone import frozen_protocol
    return frozen_protocol(path,p,'native_sources')


def locate(stage,job_id):
    matches=list((DEST/stage).glob('batch_*/'+job_id+'/arrays.npz'))
    if len(matches)!=1:raise RuntimeError(f'Expected one result for {stage}/{job_id}: {matches}')
    return matches[0]


def image_for(a,case,subset,variant=2,source_count=3):
    ids,minimum=SUBSETS[subset]
    valid=a[f'{case}_{variant}_valid'].reshape(source_count,-1)[ids].astype(bool)
    means=a[f'{case}_{variant}_means'].reshape(source_count,-1)[ids]
    count=valid.sum(axis=0)
    value=np.sqrt(np.sum(np.where(valid,abs(means)**2,0),axis=0)/np.maximum(count,1))
    value[count<minimum]=np.nan
    return value,count


def summarise(jobs,stage):
    rows=[];checks=[]
    for job in jobs:
        with np.load(locate(stage,job['id'])) as archive:
            a=dict(archive)
            points=a['points'].reshape(-1,2)
            for variant in (1,2):
                for case in ('lof','reference','difference'):
                    expected,count=image_for(a,case,'three_min1',variant)
                    np.testing.assert_allclose(expected,a[f'{case}_{variant}_amplitude'],rtol=2e-14,atol=1e-18,equal_nan=True)
                    np.testing.assert_array_equal(count,a[f'{case}_{variant}_count'])
                    checks.append(float(np.nanmax(abs(expected-a[f'{case}_{variant}_amplitude']))) if np.isfinite(expected).any() else 0.)
                    for subset in SUBSETS:
                        values,count=image_for(a,case,subset,variant)
                        row=dict(id=job['id'],path=job['selection']['path'],window=job['window'],subset=subset,
                            case=case,variant=variant,**s.metrics(values,points,'partial'),
                            one_source_near_pct=float(np.mean((np.isfinite(values)&(count==1))[s.distances(points,'partial')<=.5])*100))
                        if case=='lof':
                            ref,_=image_for(a,'reference',subset,variant)
                            common=np.isfinite(values)&np.isfinite(ref)
                            near=common&(s.distances(points,'partial')<=.5)
                            row['lof_over_pristine_near_db']=float(20*np.log10(np.max(values[near])/np.max(ref[near]))) if near.any() and np.max(ref[near])>0 and np.max(values[near])>0 else None
                        rows.append(row)
    s.write(DEST/stage/'metrics.json',dict(rows=rows))
    s.write(DEST/stage/'verification.json',dict(status='passed',native_configurations=len(jobs),
        source_recombination_checks=len(checks),max_absolute_difference=max(checks)))
    flat=[]
    for r in rows:
        flat.append({k:r.get(k) for k in ('id','path','window','subset','case','variant','support_pct','near_support_pct','cmax_db','c99_db','one_source_near_pct','lof_over_pristine_near_db')}|
                    {name+'_'+key:r[name][key] if r[name] else None for name in ('peak','near','background') for key in ('um_s','distance_mm','x_mm','depth_mm')})
    with (DEST/stage/'metrics.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    shortlist=[r for r in rows if r['case']=='lof' and r['variant']==2 and r['cmax_db'] is not None]
    print('TOP SCREENED/COHERENCE RESULTS (ranked contrast, coverage retained)',flush=True)
    for r in sorted(shortlist,key=lambda x:x['cmax_db'],reverse=True)[:35]:
        print(r['id'],r['subset'],'C',round(r['cmax_db'],2),'support',round(r['near_support_pct'],1),
              'peak_dist',round(r['peak']['distance_mm'],3),'LOF/ref',round(r['lof_over_pristine_near_db'],2) if r['lof_over_pristine_near_db'] is not None else None,flush=True)


def run():
    p=prepare();folder=DEST/'survey'
    existing={q.parent.name for q in folder.glob('batch_*/*/arrays.npz')}
    pending=[j for j in p['requests'] if j['id'] not in existing]
    # Number new batches after existing ones so a resumed run never overwrites logs.
    offset=len(list(folder.glob('batch_*')))
    for start in range(0,len(pending),12):
        batch=pending[start:start+12];dest=folder/f'batch_{offset+start//12:03d}'
        s.execute('partial',dict(requests=[dict(id=j['id'],selection=j['selection']) for j in batch]),dest)
        for j in batch:s.compact(dest/j['id'])
        print('survey',len(existing)+start+len(batch),'/',len(p['requests']),'completed',flush=True)
    summarise(p['requests'],'survey')


def confirm(stage='confirmation'):
    p=json.loads((DEST/(stage+'_protocol.json')).read_text())
    folder=DEST/stage
    for job in p['requests']:
        target=folder/job['id']
        if not (target/'arrays.npz').exists():
            # Keep the temporary basename short for Windows MAX_PATH readers.
            s.execute('partial',dict(requests=[dict(id='image',selection=job['selection'])]),target/'execution')
            s.compact(target/'execution'/'image')
            # The nested location preserves the exact per-run benchmark and log.
            (target/'execution'/'image'/'arrays.npz').replace(target/'arrays.npz')
        print(stage,job['id'],'complete',flush=True)
    rows=[];checks=[]
    for job in p['requests']:
        with np.load(folder/job['id']/'arrays.npz') as archive:
            a=dict(archive)
            points=a['points'].reshape(-1,2)
            for case in ('lof','reference','difference'):
                v=a[case+'_2_amplitude']
                row=dict(id=job['id'],candidate=job['candidate'],check=job['check'],case=case,
                         selection=job['selection'],**s.metrics(v,points,'partial'))
                if case=='lof':
                    ref=a['reference_2_amplitude'];near=(s.distances(points,'partial')<=.5)&np.isfinite(v)&np.isfinite(ref)
                    row['lof_over_pristine_near_db']=float(20*np.log10(np.max(v[near])/np.max(ref[near]))) if near.any() and np.max(ref[near])>0 else None
                rows.append(row)
                if job['check'] in ('fine','replay'):
                    with np.load(locate('survey',job['survey_id'])) as original:
                        expected,_=image_for(original,case,job['subset'])
                        step=4 if job['check']=='fine' else 1
                        actual=v.reshape(len(a['z']),len(a['x']))[::step,::step].ravel()
                        np.testing.assert_allclose(actual,expected,rtol=2e-8,atol=2e-15,equal_nan=True)
                        checks.append(dict(id=job['id'],case=case,max_difference=float(np.nanmax(abs(actual-expected))) if np.isfinite(expected).any() else 0))
    s.write(folder/'metrics.json',dict(rows=rows))
    s.write(folder/'verification.json',dict(status='passed',source_subset_and_nested_grid_checks=checks))
    for r in rows:
        if r['case']=='lof':print(r['id'],'C',round(r['cmax_db'],2) if r['cmax_db'] is not None else None,
            'support',round(r['near_support_pct'],1),'peak',round(r['peak']['distance_mm'],3) if r['peak'] else None,
            'LOF/ref',round(r['lof_over_pristine_near_db'],2) if r['lof_over_pristine_near_db'] is not None else None,flush=True)


def verify_settings():
    folder=DEST/'preset_check';jobs=[];checks=[]
    for name in ('S_LL_right','SSL_SS_left','SSL_SS_two'):
        session=json.loads((DEST/f'settings/saft_partial_{name}.json').read_text())
        assert session['schema']==4 and session['input_fingerprint']==s.dataset('partial')['input_fingerprint']
        jobs.append(dict(id=name,selection=session['selection']))
    pending=[j for j in jobs if not (folder/j['id']/'arrays.npz').exists()]
    if pending:
        s.execute('partial',dict(requests=pending),folder)
        for j in pending:s.compact(folder/j['id'])
    for job in jobs:
        with np.load(folder/job['id']/'arrays.npz') as a,np.load(DEST/f'confirmation/{job["id"]}__fine/arrays.npz') as b:
            for field in ('amplitude','count','means','valid'):
                key='lof_2_'+field
                np.testing.assert_allclose(a[key],b[key],rtol=2e-8,atol=2e-15,equal_nan=True)
        checks.append(dict(id=job['id'],settings_sha256=s.sha(DEST/f'settings/saft_partial_{job["id"]}.json')))
    s.write(DEST/'settings_verification.json',dict(status='passed',checks=checks,
        method='Native LOF-only replay of each exported preset agrees with the three-case study; GUI interaction was not tested.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','run','summarise','confirm','reflection','verify_settings']);args=ap.parse_args()
    if args.command=='prepare':print('Frozen',len(prepare()['requests']),'native configurations.')
    elif args.command=='run':run()
    elif args.command=='confirm':confirm()
    elif args.command=='reflection':confirm('reflection')
    elif args.command=='verify_settings':verify_settings()
    else:summarise(prepare()['requests'],'survey')
