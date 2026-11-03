"""Confirmation on finer/wider grids and physical acquisition groups.

Candidates are recorded before the confirmation output is evaluated.
Does not change the original survey protocol.
"""
import argparse, copy, json
import numpy as np
import study as s

def plan(kind):
    if kind=='full':
        choices=['baseline','bp_short_matched_guard','lp_short_matched_guard']
        paths=s.CORE
    else:
        choices=['baseline','band1_5','bp_short_matched_guard','lp_short_matched_guard','lp_short_initial_rest']
        paths=['1-2_SSS','1-1_SS','1-1_LL','3-2_SSSSS']
    jobs=[]
    for name in choices:
        for path in paths:
            c=s.base(kind,.1);c['path']=path
            jobs.append(dict(id=name+'__'+path+'__fine',variant=name,selection=c))
        # Acquisition and broad-field conclusions are anchored to the shared
        # band baseline, rather than choosing a separate filter per ablation.
        if name=='baseline':
            for path in (s.CORE if kind=='full' else ['1-2_SSS']):
                for label,changes in [('shoulders',dict(sources=[0,4])),('central',dict(sources=[1,2,3])),
                    ('plate',dict(groups=[True,False,True])),('floor',dict(groups=[False,True,False])),
                    ('shoulders_plate',dict(sources=[0,4],groups=[True,False,True]))]:
                    c=s.base(kind,.1);c.update(path=path,**changes)
                    jobs.append(dict(id=name+'__'+path+'__'+label,variant=name,selection=c))
                c=s.base(kind,.25);c.update(path=path,x_min_mm=-45,x_max_mm=45,depth_min_mm=.25,depth_max_mm=25.25)
                jobs.append(dict(id=name+'__'+path+'__wide',variant=name,selection=c))
        if name in ['baseline','lp_short_matched_guard']:
            path='1-1_SS' if kind=='full' else '1-2_SSS'
            c=s.base(kind,.05);c['path']=path
            jobs.append(dict(id=name+'__'+path+'__finer',variant=name,selection=c))
    return jobs

def run(kind):
    jobs=plan(kind);protocol=s.OUT/'confirmation'/kind/'protocol.json'
    if protocol.exists():assert json.loads(protocol.read_text())==jobs
    else:s.write(protocol,jobs)
    for name in dict.fromkeys(j['variant'] for j in jobs):
        group=[j for j in jobs if j['variant']==name]
        output=s.OUT/'confirmation'/kind
        if all((output/j['id']/'arrays.npz').exists() for j in group):continue
        spec=next(v for v in s.variants() if v['name']==name);folder=s.bank(kind,spec)
        req=[dict(id=j['id'],selection=j['selection'],signal_processing=dict(spec,bank=str(folder))) for j in group]
        print(kind,name,'confirmation',len(group),'jobs',flush=True)
        s.execute(kind,dict(requests=req),output)
        for j in group:s.compact(output/j['id'])
    summarise(kind)

def summarise(kind):
    output=s.OUT/'confirmation'/kind;rows=[]
    for job in json.loads((output/'protocol.json').read_text()):
        path=output/job['id']/'arrays.npz'
        if not path.exists():continue
        with np.load(path) as a:
            for case in ['lof','reference','difference']:
                rows.append(dict(id=job['id'],kind=kind,case=case,variant=job['variant'],path=job['selection']['path'],
                                 **s.metrics(a[case+'_2_amplitude'],a['points'].reshape(-1,2),kind)))
    s.write(output/'metrics.json',dict(rows=rows))
    for r in rows:
        if r['case']=='lof':print(r['id'],r['cmax_db'],r['near_support_pct'],flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('kind',choices=['full','partial']);args=parser.parse_args();run(args.kind)
