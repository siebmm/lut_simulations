"""Fixed joint source/receiver test using the current equal-mean SAFT protocol."""
from pathlib import Path
import copy
import json
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent/'filter_study_20261101'))
import study as s

OUT = ROOT/'studies/laser/results/partial_acquisition_20261101'
DATA = HERE.parent/'comsol_chapter_20261101/data'


def plan():
    records = json.loads((DATA/'images.json').read_text())
    base = copy.deepcopy(next(r['selection'] for r in records
                         if r['kind']=='partial' and r['id']=='baseline__1-2_SSS__fine'))
    assert base['processing']=='raw' and base['coherence'] and base['pulse_screening']=='LSR'
    jobs = []
    for path in s.PATHS:
        c = dict(base, path=path, spacing_mm=.2, sources=[1,2,3], groups=[True,False,True])
        jobs.append(dict(id='survey__'+path, selection=c, role='joint acquisition survey'))
    for name, groups in [('both_sides',[True,False,True]), ('left_only',[True,False,False]), ('right_only',[False,False,True])]:
        c = dict(base, sources=[1,2,3], groups=groups)
        jobs.append(dict(id='fine__'+name, selection=c, role='declared S|SS acquisition confirmation'))
    jobs.append(dict(id='control__all',selection=base,role='independent native reproduction of current paper'))
    return jobs


def main():
    if '--confirm-left' in sys.argv:
        return confirm_left()
    jobs = plan()
    protocol = dict(scope='18 existing paths with three weld-surface sources and both side-surface receiver groups at 0.2 mm; S|SS on both/left/right sides at 0.1 mm; one all-channel reproduction control.',
                    interpretation='Shoulder receivers mean the complete flat left/right side-surface groups, not narrow windows centred at source x=+/-12 mm.',
                    frozen='All geometry, filter, gates, coherence, screening and support rules fixed before examining this joint-acquisition test. No source/receiver IDs fitted to LOF contrast.',
                    cases=['lof','reference','difference'], requests=jobs,
                    scientific_configurations=21,reproduction_controls=1,
                    dataset_fingerprint=s.dataset('partial')['input_fingerprint'],executable_sha256=s.sha(s.EXE))
    p=HERE/'protocol.json'
    if p.exists(): assert json.loads(p.read_text())==protocol
    else:s.write(p,protocol)
    pending=[dict(id=j['id'],selection=j['selection']) for j in jobs if not (OUT/j['id']/'arrays.npz').exists()]
    if pending:
        print('Computing',len(pending),'configurations with the standard native filter.',flush=True)
        s.execute('partial',dict(requests=pending),OUT)
        for j in pending:s.compact(OUT/j['id'])
    rows=[]
    with np.load(OUT/'control__all/arrays.npz') as control:
        for job in jobs:
            with np.load(OUT/job['id']/'arrays.npz') as a:
                points=a['points'].reshape(-1,2)
                for case in ['lof','reference','difference']:
                    values=a[case+'_2_amplitude']
                    row=dict(id=job['id'],path=job['selection']['path'],spacing_mm=job['selection']['spacing_mm'],case=case,
                             **s.metrics(values,points,'partial'))
                    if job['id'].startswith('fine__'):
                        same=np.isfinite(values)&np.isfinite(control[case+'_2_amplitude'])
                        row['common_support_joint']=s.metrics(values,points,'partial',same)
                        row['common_support_all']=s.metrics(control[case+'_2_amplitude'],points,'partial',same)
                    rows.append(row)
                if job['id']=='control__all':
                    index=json.loads((DATA/'images.json').read_text())
                    old=next(r for r in index if r['kind']=='partial' and r['id']=='baseline__1-2_SSS__fine')
                    with np.load(DATA/(old['uid']+'.npz')) as previous:
                        for case in ['lof','reference','difference']:
                            np.testing.assert_allclose(a[case+'_2_amplitude'],previous[case+'_amplitude'].ravel(),rtol=2e-9,atol=1e-15,equal_nan=True)
                            np.testing.assert_array_equal(a[case+'_2_count'],previous[case+'_count'].ravel())
    s.write(OUT/'metrics.json',dict(rows=rows))
    s.write(OUT/'verification.json',dict(status='passed',configurations=len(jobs),scientific_configurations=21,
                                       native_reproduction_cases=3,protocol_sha256=s.sha(p),equal_receivers=True))
    for r in rows:
        if r['case']=='lof':print(r['id'],r['cmax_db'],r['near_support_pct'],r['peak'],flush=True)


def confirm_left():
    base=next(j['selection'] for j in plan() if j['id']=='fine__left_only')
    jobs=[]
    for name,groups in [('left_only',[True,False,False]),('both_sides',[True,False,True])]:
        for grid,changes in [('finer',dict(spacing_mm=.05)),('wide',dict(spacing_mm=.25,x_min_mm=-45,x_max_mm=45,depth_min_mm=.25,depth_max_mm=25.25))]:
            jobs.append(dict(id=name+'__'+grid,selection=dict(base,groups=groups,**changes)))
    protocol=dict(scope='Confirm left/both side S|SS with finer sampling and a wider field; fixed source, filter and admission rules.',
                  decision='Left-only is a candidate for this path/geometry. Check global peak localization, background contrast and coverage; do not infer whole-flaw inspection from a supported local maximum.',
                  requests=jobs,dataset_fingerprint=s.dataset('partial')['input_fingerprint'])
    p=HERE/'left_confirmation_protocol.json'
    if p.exists():assert json.loads(p.read_text())==protocol
    else:s.write(p,protocol)
    out=OUT/'left_confirmation'
    pending=[j for j in jobs if not (out/j['id']/'arrays.npz').exists()]
    if pending:
        print('Confirming left/both receiver sides on finer and wider fields.',flush=True)
        s.execute('partial',dict(requests=pending),out)
        for j in pending:s.compact(out/j['id'])
    rows=[]
    for j in jobs:
        with np.load(out/j['id']/'arrays.npz') as a:
            for case in ['lof','reference','difference']:
                rows.append(dict(id=j['id'],case=case,**s.metrics(a[case+'_2_amplitude'],a['points'].reshape(-1,2),'partial')))
    s.write(out/'metrics.json',dict(rows=rows))
    s.write(out/'verification.json',dict(status='passed',configurations=4,cases=12,protocol_sha256=s.sha(p)))
    for r in rows:
        if r['case']=='lof':print(r['id'],r['cmax_db'],r['near_support_pct'],r['peak'],flush=True)


if __name__=='__main__':main()
