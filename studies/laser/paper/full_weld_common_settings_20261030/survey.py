"""Controlled test of every common direct-pulse-screen family selection."""

# Locate portable scientific inputs without depending on the source checkout.
import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import data_root as _data_root, native_root as _native_root, executable as _executable, native_env as _native_env
from pathlib import Path
import sys,json,copy,subprocess
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT/'studies/laser/paper/full_weld_four_paths_20261030'))
from analyse import read_blend,metrics,PATHS
OUT=ROOT/'studies/laser/results/full_weld_common_settings_20261030'

def prepare():
    m,_=read_blend(); base=m['layers'][0]['selection'];jobs=[];idx=[]
    for mode in ['off','L','S','R','LS','LR','SR','LSR']:
        for k,path in enumerate(PATHS):
            c=copy.deepcopy(base);c.update(path=path,pulse_screening=mode,spacing_mm=.2)
            idx.append({'index':len(jobs),'path':path,'screen':mode,'selection':c});jobs.append({'selection':c})
    (OUT/'survey_requests.json').write_text(json.dumps({'verify_inputs':True,'requests':jobs},indent=2))
    (OUT/'survey_index.json').write_text(json.dumps(idx,indent=2))
def run():
    exe=_executable('SAFTBatch')
    with (OUT/'survey_compute.log').open('w') as f:
        p=subprocess.run([str(exe),'--dataset',str(_data_root()/'full_weld'),'--compute',str(OUT/'survey_requests.json'),'--output',str(OUT/'survey')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError((OUT/'survey_compute.log').read_text())
def summarise():
    rows=[]
    for r in json.loads((OUT/'survey_index.json').read_text()):
        a=np.load(OUT/'survey'/str(r['index'])/'images.npz'); rows.append({k:r[k] for k in ['index','path','screen']}|metrics(a['lof_maps'][2],a['points_m']))
    (OUT/'survey_metrics.json').write_text(json.dumps(rows,indent=2))
    for mode in ['off','L','S','R','LS','LR','SR','LSR']:
        print(mode)
        for r in rows:
            if r['screen']!=mode:continue
            print(r['path'],*[round(r[k],3) if r[k] is not None else None for k in ['near_peak_um_s','contrast_peak_db','contrast_p99_db','support_pct','near_support_pct','global_peak_distance_mm']])
if __name__=='__main__':prepare();run();summarise()
