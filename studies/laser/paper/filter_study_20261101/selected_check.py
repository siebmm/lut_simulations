"""Independent checks of the full-weld shorter low-pass shear candidate."""
import copy,json
import numpy as np
import study as s

def main():
    kind='full';spec=next(v for v in s.variants() if v['name']=='lp_short_matched_guard')
    jobs=[]
    for label,changes in [('all',{}),('shoulders',dict(sources=[0,4])),('central',dict(sources=[1,2,3])),
        ('plate',dict(groups=[True,False,True])),('floor',dict(groups=[False,True,False])),
        ('shoulders_plate',dict(sources=[0,4],groups=[True,False,True])),
        ('gate15',dict(angular_gate=.15)),('gate25',dict(angular_gate=.25)),
        ('lowpass_launch_curve',dict(low_MHz=.0001)),
        ('wide',dict(spacing_mm=.25,x_min_mm=-45,x_max_mm=45,depth_min_mm=.25,depth_max_mm=25.25))]:
        c=s.base(kind,.1);c.update(path='1-1_SS',**changes);jobs.append(dict(id=label,selection=c))
    out=s.OUT/'selected_checks/full';s.write(out/'protocol.json',dict(
        selection_basis='Short low-pass candidate increased full-weld S|S contrast at 0.2, 0.1 and 0.05 mm without reducing near-flaw support; test physical groups, broad field and gate sensitivity before using it as a primary example.',
        curve_check='0.0001-5 MHz original long-filter launch curve approximates 0-5 MHz; receiving/wall physics unchanged. This checks coupling-band sensitivity, not a fitted new source.',
        spec=spec,jobs=jobs))
    folder=s.bank(kind,spec)
    req=[dict(**j,signal_processing=dict(spec,bank=str(folder))) for j in jobs]
    s.execute(kind,dict(requests=req),out)
    rows=[]
    for j in jobs:
        p=out/j['id'];s.compact(p)
        with np.load(p/'arrays.npz') as a:
            for case in ['lof','reference','difference']:
                rows.append(dict(id=j['id'],case=case,**s.metrics(a[case+'_2_amplitude'],a['points'].reshape(-1,2),kind)))
    s.write(out/'metrics.json',dict(rows=rows))
    for r in rows:
        if r['case']=='lof':print(r['id'],r['cmax_db'],r['near_support_pct'],flush=True)

if __name__=='__main__':main()
