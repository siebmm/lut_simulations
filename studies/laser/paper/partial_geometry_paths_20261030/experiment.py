"""Fixed-setting test of finite-bevel and fixed-corner SAFT hypotheses."""
from pathlib import Path
import copy, datetime, hashlib, json, subprocess, sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
sys.path.insert(0, str(ROOT/'studies/laser/src'))
from standalone import executable, native_env
sys.path.insert(0,str(HERE.parent/'weld_geometry_study_20261030'))
import study as old
OUT=ROOT/'studies/laser/results/partial_geometry_paths_20261030'
DATA=old.DATA['partial']
PACKAGE=ROOT/'studies/laser/build/desktop_cpp/SAFTWorkbench'
BASE=json.loads((HERE.parent/'weld_geometry_study_20261030/locked_protocol.json').read_text())['common_controls']
def write(p,o):
    p.parent.mkdir(parents=True,exist_ok=True)
    p.write_text(json.dumps(o,indent=2),encoding='utf-8')
def exe():return executable('SAFTBatch')
def run_request(request,folder,program=None):
    folder.mkdir(parents=True,exist_ok=True)
    program=program or exe()
    write(folder/'executable_evidence.json',dict(executable=str(program),sha256=hashlib.sha256(program.read_bytes()).hexdigest()))
    with (folder/'compute.log').open('w') as log:
        p=subprocess.run([str(program),'--repository',str(ROOT),'--dataset',str(DATA),'--compute',str(request),'--output',str(folder)],cwd=ROOT,env=native_env(),stdout=log,stderr=subprocess.STDOUT)
    assert p.returncode==0,(folder/'compute.log').read_text()
def prepare():
    catalogue=json.loads((OUT/'verification/native/path_catalogue.json').read_text())
    paths=[p for p in catalogue['paths'] if p['family']]
    assert len(paths)==192
    controls=['1-1_SS','1-1_LL','1-2_SSS','2-1_SSS','1-2_LLL','3-2_SSSSS']
    jobs=[];index=[]
    for p in [{'id':p,'family':'conventional'} for p in controls]+paths:
        for processing,spacing in [('weighted',.1),('raw',.2 if p['family']!='conventional' else .1)]:
            c=copy.deepcopy(BASE)
            c.update(path=p['id'],processing=processing,spacing_mm=spacing,cases=['lof','reference','difference'],receivers=list(range(331)))
            index.append(dict(index=len(jobs),path=p['id'],family=p['family'],processing=processing,spacing_mm=spacing,selection=c))
            jobs.append(dict(selection=c))
    plan=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        basis='All six groove vertices and both sloping facets; both pure L and S; encounter before or after the image pixel; optional ID reflection before or after the geometry encounter; other branch direct or one ID reflection. No feature or setting chosen from LOF scores.',
        scope='Same-mode specular elastic reflection coefficients. Unit corner coefficient: timing-only diffraction hypothesis without a calibrated wedge amplitude/phase. Converted corner modes and Rayleigh-to-bulk diffraction are outside scope.',
        protocol='Common known-response weighting at 0.1 mm, plus a timing-mean sensitivity at 0.2 mm for every new path. Three matched signal cases for every job. Six conventional controls in both processing modes at 0.1 mm.',
        inference='Exploratory test on previously seen data, not blind validation. Report coverage and pristine controls with contrast. The strongest result across 192 new paths is descriptive, not a tuned recommendation.',
        base=BASE,catalogue=catalogue,requests=jobs,index=index)
    lock=HERE/'locked_protocol.json'
    if lock.exists():
        prev=json.loads(lock.read_text());assert prev['requests']==jobs;plan=prev
    else:write(lock,plan)
    write(OUT/'protocol_freeze.json',dict(sha256=hashlib.sha256(lock.read_bytes()).hexdigest(),protocol=plan))
    write(OUT/'requests.json',dict(requests=jobs));write(OUT/'index.json',index)
    print(len(jobs),'configurations;',len(jobs)*3,'signal-case reconstructions',flush=True)
def summarise():
    rows=[]
    for r in json.loads((OUT/'index.json').read_text()):
        with np.load(OUT/'numerical'/str(r['index'])/'images.npz') as a:
            for case in ['lof','reference','difference']:
                rows.append({k:r[k] for k in ['index','path','family','processing','spacing_mm']}|dict(case=case,metrics=old.metrics(a[case+'_maps'][2],a['points_m'],'partial')))
    write(OUT/'metrics.json',rows)
    for family in ['conventional','wall','corner']:
        subset=[r for r in rows if r['case']=='lof' and r['family']==family and r['processing']=='weighted']
        ranked=sorted([r for r in subset if r['metrics']['cmax_db'] is not None],key=lambda r:r['metrics']['cmax_db'],reverse=True)
        print(family,'supported',sum(r['metrics']['supported']>0 for r in subset),'of',len(subset))
        for r in ranked[:5]:print(r['path'],r['metrics']['cmax_db'],r['metrics']['near_support_pct'])
if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='run':run_request(OUT/'requests.json',OUT/'numerical');print('completed',flush=True)
    elif sys.argv[1]=='summarise':summarise()
