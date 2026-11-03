"""Focused follow-up: physics-based improvements and acquisition ablations."""

import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import executable as _executable, native_env as _native_env
from pathlib import Path
import json,copy,sys,subprocess,datetime,hashlib,csv,itertools
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
PREV=ROOT/'studies/laser/paper/weld_geometry_study_20261030'
sys.path.insert(0,str(PREV));import study as old
OUT=ROOT/'studies/laser/results/laser_weld_feasibility_20261030'
PATHS=old.CORE+['1-2_LLL','3-2_SSSSS']
def write(p,o):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(o,indent=2),encoding='utf-8')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def prepare():
    OUT.mkdir(parents=True,exist_ok=True)
    base=json.loads((PREV/'locked_protocol.json').read_text())['common_controls']
    factors=[('band1',{'low_MHz':1}),('weighted',{'processing':'weighted'}),('band1_weighted',{'low_MHz':1,'processing':'weighted'}),('plate_only',{'groups':[True,False,True]}),('band1_plate',{'low_MHz':1,'groups':[True,False,True]})]
    plan={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'basis':'Follow-up to known outcomes. Confirm 1-5 MHz as a candidate for suppressing low-frequency background; test response weighting and receiver groups. No parameter sweep over flaw coordinates, no flaw-centred crop, unchanged physical gates and support rules.',
       'candidate':'band1','candidate_basis':'Earlier one-factor evidence suggested a modest gain; validate at 0.1 mm, on reference controls, wider field and source/receiver group ablations. Exploratory confirmation on same known flaw, not held-out validation.',
       'paths':PATHS,'factors':factors,'base_protocol_sha256':sha(PREV/'locked_protocol.json')}
    lock=HERE/'followup_protocol.json'
    if not lock.exists():write(lock,plan)
    for kind in ['partial','full']:
        idx=[];jobs=[];d=old.dataset(kind)
        def add(path,name,changes,cases=None,spacing=.1):
            c=copy.deepcopy(base);c.update(path=path,receivers=list(range(len(d['audit']['receivers']))),spacing_mm=spacing,cases=cases or ['lof']);c.update(changes)
            idx.append({'index':len(jobs),'kind':kind,'path':path,'name':name,'selection':c});jobs.append({'selection':c})
        if kind=='partial':
            for p in PATHS:
                for name,changes in factors:add(p,name,changes,['lof','reference','difference'])
            for p in old.CORE:
                for name,groups in [('floor_only',[False,True,False]),('left_only',[True,False,False]),('right_only',[False,False,True])]:add(p,name,{'groups':groups})
            # Candidate validation: identical change in a finer and wider image.
            for name,changes,spacing in [('band1_fine',{'low_MHz':1},.05),('band1_wide',{'low_MHz':1,'x_min_mm':-45,'x_max_mm':45,'depth_min_mm':.25,'depth_max_mm':25.25},.25)]:add('1-2_SSS',name,changes,spacing=spacing)
            for name,groups in [('band1_floor',[False,True,False]),('band1_left',[True,False,False]),('band1_right',[False,False,True])]:add('1-2_SSS',name,{'low_MHz':1,'groups':groups})
        else:
            for p in old.CORE:
                for name,groups in [('plate_only',[True,False,True]),('floor_only',[False,True,False]),('left_only',[True,False,False]),('right_only',[False,False,True])]:add(p,name,{'groups':groups})
        write(OUT/kind/'index.json',idx);write(OUT/kind/'requests.json',{'verify_inputs':True,'requests':jobs});print(kind,len(jobs),'configs',sum(len(r['selection']['cases']) for r in idx),'case images')
def index(kind):return json.loads((OUT/kind/'index.json').read_text())
def load(kind,p,name):
    if name=='baseline':return old.load(kind,p)
    r=next(r for r in index(kind) if r['path']==p and r['name']==name)
    return np.load(OUT/kind/'numerical'/str(r['index'])/'images.npz')
def run(kind):
    exe=_executable('SAFTBatch')
    with (OUT/kind/'compute.log').open('w') as f:r=subprocess.run([str(exe),'--dataset',str(old.DATA[kind]),'--compute',str(OUT/kind/'requests.json'),'--output',str(OUT/kind/'numerical')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    assert r.returncode==0,(OUT/kind/'compute.log').read_text();print(kind,'completed',flush=True)
def source_subsets(a,kind):
    groups={'all':[0,1,2,3,4],'shoulders':[0,4],'central':[1,2,3]}
    groups.update({f'without_S{i}':[j for j in range(5) if j!=i] for i in range(5)})
    groups.update({f'pair_{i}_{j}':[i,j] for i,j in itertools.combinations(range(5),2)})
    vals=a['lof_means'][2];valid=a['lof_valid'][2].astype(bool);rows=[]
    for name,ids in groups.items():
        v=valid[ids];means=vals[ids];n=v.sum(axis=0);score=np.sqrt(np.sum(np.where(v,abs(means)**2,0),axis=0)/np.maximum(n,1));score[n<2]=np.nan
        if name=='all':assert np.allclose(score,a['lof_maps'][2],equal_nan=True,rtol=1e-13,atol=1e-17)
        rows.append({'source_group':name,'sources':ids,'metrics':old.metrics(score,a['points_m'],kind)})
    return rows
def summarise():
    rows=[];subs=[]
    for kind in ['partial','full']:
        for p in PATHS if kind=='partial' else old.CORE:
            a=old.load(kind,p)
            for case in ['lof','reference','difference']:rows.append({'kind':kind,'path':p,'name':'baseline','case':case,'metrics':old.metrics(a[case+'_maps'][2],a['points_m'],kind)})
        for r in index(kind):
            a=load(kind,r['path'],r['name'])
            for case in r['selection']['cases']:rows.append({'kind':kind,'path':r['path'],'name':r['name'],'case':case,'metrics':old.metrics(a[case+'_maps'][2],a['points_m'],kind)})
            if r['kind']=='partial' and r['path']=='1-2_SSS' and r['name']=='band1_fine':
                b=load(kind,r['path'],'band1');v=a['lof_maps'][2].reshape(len(a['depth_mm']),-1)[::2,::2].ravel();assert np.array_equal(v,b['lof_maps'][2],equal_nan=True)
        for p in old.CORE:
            for name in ['baseline']+(['band1'] if kind=='partial' and p=='1-2_SSS' else []):
                a=load(kind,p,name);subs.extend({'kind':kind,'path':p,'name':name,**r} for r in source_subsets(a,kind))
    write(OUT/'analysis.json',rows);write(OUT/'source_subsets.json',subs)
    flat=[]
    for r in rows:
        m=r['metrics'];flat.append({k:r[k] for k in ['kind','path','name','case']}|{k:m[k] for k in ['support_pct','near_support_pct','cmax_db','c99_db']}|{'near_um_s':m['near']['um_s'] if m['near'] else None,'peak_um_s':m['peak']['um_s'] if m['peak'] else None,'peak_distance_mm':m['peak']['distance_mm'] if m['peak'] else None})
    with (HERE/'results.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    print('Summarised',len(rows),'images and',len(subs),'source-subset assessments')
if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='run':run(sys.argv[2])
    elif sys.argv[1]=='summarise':summarise()
