"""Recorded second-stage confirmation, after seeing the first follow-up results."""

import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import executable as _executable, native_env as _native_env
from analyse import *
DEST=OUT/'confirmation'

def prepare():
    base=json.loads((PREV/'locked_protocol.json').read_text())['common_controls']
    plan={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'stage':'Confirmation after first follow-up; same known flaw, not independent validation.',
      'reason':'Response weighting improved partial S|SS Cmax at unchanged support; test finer/wider fields and acquisition groups. Preserve original followup_protocol.json, whose candidate was band1.',
      'primary_candidate':'weighted: 0.5-5 MHz with known receiver/wall response weighting',
      'secondary_candidate':'band1_weighted: 1-5 MHz and response weighting',
      'no_search':'No additional gate, margin, crop or pulse-screening search.',
      'first_stage_results_sha256':sha(OUT/'analysis.json')}
    if not (HERE/'confirmation_protocol.json').exists():write(HERE/'confirmation_protocol.json',plan)
    for kind in ['partial','full']:
        idx=[];jobs=[];d=old.dataset(kind)
        def add(p,name,changes):
            c=copy.deepcopy(base);c.update(path=p,receivers=list(range(len(d['audit']['receivers']))),spacing_mm=.1,cases=['lof']);c.update(changes)
            idx.append({'index':len(jobs),'kind':kind,'path':p,'name':name,'selection':c});jobs.append({'selection':c})
        if kind=='partial':
            p='1-2_SSS'
            for name,extra in [('weighted',{}),('band1_weighted',{'low_MHz':1})]:
                for suffix,changes in [('fine',{'spacing_mm':.05}),('wide',{'x_min_mm':-45,'x_max_mm':45,'depth_min_mm':.25,'depth_max_mm':25.25,'spacing_mm':.25})]:
                    add(p,name+'_'+suffix,{'processing':'weighted',**extra,**changes})
            for name,groups in [('plate',[True,False,True]),('floor',[False,True,False]),('left',[True,False,False]),('right',[False,False,True])]:add(p,'weighted_'+name,{'processing':'weighted','groups':groups})
            for name,src in [('shoulders',[0,4]),('central',[1,2,3]),('without_S3',[0,1,2,4]),('without_S4',[0,1,2,3])]:
                add(p,'weighted_'+name,{'processing':'weighted','sources':src})
            add(p,'weighted_central_plate',{'processing':'weighted','sources':[1,2,3],'groups':[True,False,True]})
            add(p,'baseline_central',{'sources':[1,2,3]})
        else:
            for p in old.CORE:
                add(p,'shoulders_plate',{'sources':[0,4],'groups':[True,False,True]})
            add('1-1_SS','baseline_shoulders',{'sources':[0,4]})
        write(DEST/kind/'index.json',idx);write(DEST/kind/'requests.json',{'verify_inputs':True,'requests':jobs})
        print(kind,len(jobs),'confirmation images',flush=True)

def index(kind):return json.loads((DEST/kind/'index.json').read_text())
def load(kind,p,name):
    r=next(r for r in index(kind) if r['path']==p and r['name']==name)
    return np.load(DEST/kind/'numerical'/str(r['index'])/'images.npz')
def run(kind):
    exe=_executable('SAFTBatch')
    with (DEST/kind/'compute.log').open('w') as f:r=subprocess.run([str(exe),'--dataset',str(old.DATA[kind]),'--compute',str(DEST/kind/'requests.json'),'--output',str(DEST/kind/'numerical')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    assert r.returncode==0,(DEST/kind/'compute.log').read_text();print(kind,'confirmation completed',flush=True)

def summarise():
    import analyse as first
    rows=[];checks=[]
    for kind in ['partial','full']:
        for r in index(kind):
            a=load(kind,r['path'],r['name']);m=old.metrics(a['lof_maps'][2],a['points_m'],kind)
            rows.append({k:r[k] for k in ['kind','path','name']}|{'case':'lof','metrics':m})
            if r['name'].endswith('_fine'):
                b=first.load(kind,r['path'],r['name'].removesuffix('_fine'));v=a['lof_maps'][2].reshape(len(a['depth_mm']),-1)[::2,::2].ravel()
                assert np.array_equal(v,b['lof_maps'][2],equal_nan=True);checks.append({'kind':kind,'name':r['name'],'nested_grid_exact':True})
            if r['name'] in ['baseline_central','baseline_shoulders','weighted_central','weighted_shoulders','weighted_without_S3','weighted_without_S4']:
                name,group=r['name'].split('_',1);b=first.load(kind,r['path'],name);ids=r['selection']['sources'];valid=b['lof_valid'][2,ids].astype(bool);n=valid.sum(axis=0)
                v=np.sqrt(np.sum(np.where(valid,abs(b['lof_means'][2,ids])**2,0),axis=0)/np.maximum(n,1));v[n<2]=np.nan
                assert np.allclose(v,a['lof_maps'][2],equal_nan=True,rtol=1e-13,atol=1e-17);checks.append({'kind':kind,'name':r['name'],'source_RMS_native_verified':True})
    write(DEST/'analysis.json',rows);write(DEST/'checks.json',checks)
    a=first.load('partial','1-2_SSS','weighted');write(DEST/'weighted_source_subsets.json',first.source_subsets(a,'partial'))
    for r in rows:
        m=r['metrics'];print(r['kind'],r['path'],r['name'],*[round(m[k],3) if m[k] is not None else None for k in ['near_support_pct','cmax_db','c99_db']],m['peak'])
if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='run':run(sys.argv[2])
    elif sys.argv[1]=='summarise':summarise()
