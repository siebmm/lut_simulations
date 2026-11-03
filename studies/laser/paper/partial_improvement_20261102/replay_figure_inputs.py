"""Recompute only the grids used in the final review figures, then draw them."""
import json
import shutil
import numpy as np
import review as r
import review_confirm as c
from standalone import ROOT as INPUT_ROOT

def load(path):
    return json.loads(path.read_text())

def replay(jobs, destination):
    missing = [j for j in jobs if not (destination/j['id']/'arrays.npz').is_file()]
    if missing:
        print('Replaying figure inputs:', ', '.join(j['id'] for j in missing), flush=True)
        r.execute(missing, destination, False)

def main():
    verified=[]
    def check(folder, expected):
        with np.load(folder/'arrays.npz') as a:
            points=a['points'].reshape(-1,2)
            distances=r.s.distances(points,'partial')
            for case in ['lof','reference','difference']:
                actual=r.metric(a[case+'_2_amplitude'],points,distances)
                old=expected[case]
                for field in ['cmax_db','near_support_pct','support_pct']:
                    if old[field] is None:
                        assert actual[field] is None
                    else:
                        np.testing.assert_allclose(actual[field],old[field],rtol=2e-8,atol=1e-8)
                verified.append({'id':folder.name,'case':case,'passed':True})
    jobs = load(r.HERE/'protocol.json')['requests']
    replay([j for j in jobs if j['id']=='all__raw__min2'], r.OUT/'survey')
    jobs = load(r.DEST/'refinement/protocol.json')['requests']
    jobs = [j for j in jobs if j['id'] in ('muted_025','unmuted_025')]
    for job in jobs:
        saved = INPUT_ROOT/'output/main/partial'/job['id']/'arrays.npz'
        target = r.DEST/'refinement'/job['id']/'arrays.npz'
        if saved.exists() and not target.exists():
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(saved, target)
    replay(jobs, r.DEST/'refinement')
    rows=load(r.DEST/'refinement/metrics.json')['rows']
    for j in jobs:
        check(r.DEST/'refinement'/j['id'],{case:next(x for x in rows if x['id']==j['id'] and x['case']==case) for case in ['lof','reference','difference']})
    jobs = load(c.D/'protocol.json')['requests']
    jobs = [j for j in jobs if j['candidate'] in ('right_SLL_LL03','right_SSLSS_band15') and j['label'] in ('fine','wide')]
    replay([j for j in jobs if not list(c.D.glob('b*/'+j['id']+'/arrays.npz'))], c.D/'bfigure')
    rows=load(c.D/'metrics.json')['rows']
    for j in jobs:
        check(c.locate(j['id']).parent,next(x for x in rows if x['id']==j['id']))
    dest = r.DEST/'path_confirmation'
    jobs = load(dest/'protocol.json')['requests']
    jobs=[j for j in jobs if j['candidate'] in (0,9,8) and j['label']=='fine']
    replay(jobs, dest)
    rows=load(dest/'metrics.json')['rows']
    for j in jobs:
        check(dest/j['id'],next(x for x in rows if x['id']==j['id']))
    r.s.write(r.DEST/'figure_replay_check.json',{'status':'passed','checks':verified})
    import review_publication
    review_publication.figures()
    print('Five final review figures regenerated.', flush=True)

if __name__ == '__main__':
    main()
