"""Outcome-informed checks of candidates from the comprehensive settings audit."""
import argparse
import json
import numpy as np
import review as r
import settings_review as q

D=r.DEST/'confirmation'


def protocol():
    recipes=[
        ('broad_mute03','broad_SSS',dict(mutes=[True]*6,mute_halfwidth_us=.3)),
        ('broad_mute06','broad_SSS',dict(mutes=[True]*6,mute_halfwidth_us=.6)),
        ('broad_SS03','broad_SSS',dict(mutes=[False]*5+[True],mute_halfwidth_us=.3)),
        ('broad_CF2','broad_SSS',dict(coherence_power=2)),
        ('broad_left','broad_SSS',dict(aperture='left')),
        ('right_SLL_mute03','right_SLL',dict(mutes=[True]*6,mute_halfwidth_us=.3)),
        ('right_SLL_LL03','right_SLL',dict(mutes=[False]*3+[True,False,False],mute_halfwidth_us=.3)),
        ('right_SSLSS_band15','right_SSLSS',dict(low_MHz=1.5))]
    jobs=[]
    def add(candidate,label,c):
        jobs.append(dict(id=f'c{len(jobs):03d}',candidate=candidate,label=label,selection=c))
    for name,anchor,change in recipes:
        c=dict(q.anchors()[anchor],**change)
        for grid,fields in [('fine',dict(spacing_mm=.05)),('wide',dict(spacing_mm=.25,x_min_mm=-45,x_max_mm=45,depth_min_mm=.25,depth_max_mm=25.25))]:
            add(name,grid,dict(c,**fields))
    for anchor in ('broad_SSS','right_SLL'):
        for width in (.1,.2,.3,.4,.5,.6,.8):
            add(anchor,f'width{width:g}',dict(q.anchors()[anchor],spacing_mm=.1,mutes=[True]*6,mute_halfwidth_us=width))
        for aperture in ('even','odd'):
            add(anchor,aperture,dict(q.anchors()[anchor],spacing_mm=.1,mutes=[True]*6,mute_halfwidth_us=.3,aperture=aperture))
        add(anchor,'gate20_mute03',dict(q.anchors()[anchor],spacing_mm=.1,mutes=[True]*6,mute_halfwidth_us=.3,angular_gate=.2))
    p=dict(reason='The settings matrix reveals positive effects of fixed arrival mutes and higher coherence powers; verify field/grid dependence and sensitivity before interpreting increased contrast as an improvement.',
        definition='Fixed mutes are smooth tapers centred on specimen-only predicted arrivals; pulse-aware screens are a separate setting. Mask width is a hypothesis, not a measured packet duration. Coherence power changes nonlinear image weighting and does not add acoustic information.',
        baseline_metrics_sha256=r.s.sha(q.D/'native_metrics.json'),recipes=recipes,requests=jobs)
    path=D/'protocol.json'
    if path.exists():assert json.loads(path.read_text())==json.loads(json.dumps(p))
    else:r.s.write(path,p)
    return p


def run():
    p=protocol();done={f.parent.name for f in D.glob('b*/*/arrays.npz')}
    pending=[j for j in p['requests'] if j['id'] not in done];offset=len(list(D.glob('b*')))
    for start in range(0,len(pending),4):
        batch=pending[start:start+4];r.execute(batch,D/f'b{offset+start//4:03d}',False)
        print('confirm',len(done)+start+len(batch),'/',len(p['requests']),flush=True)
    analyse()


def locate(job_id):
    paths=list(D.glob('b*/'+job_id+'/arrays.npz'));assert len(paths)==1
    return paths[0]


def analyse():
    rows=[]
    for job in protocol()['requests']:
        with np.load(locate(job['id'])) as f:a=dict(f)
        points=a['points'].reshape(-1,2);d=r.s.distances(points,'partial')
        result={k:job[k] for k in ('id','candidate','label','selection')}
        for case in ('lof','reference','difference'):
            result[case]=r.metric(a[case+'_2_amplitude'],points,d)
        common=(d<=.5)&np.isfinite(a['lof_2_amplitude'])&np.isfinite(a['reference_2_amplitude'])
        l,p=a['lof_2_amplitude'][common],a['reference_2_amplitude'][common]
        result['lof_over_pristine_near_db']=float(20*np.log10(l.max()/p.max())) if common.any() and l.max()>0 and p.max()>0 else None
        rows.append(result)
    r.s.write(D/'metrics.json',dict(rows=rows))
    for x in rows:
        m=x['lof'];print(x['candidate'],x['label'],m['cmax_db'],m['near_support_pct'],m['peak_distance_mm'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','run','analyse']);cmd=ap.parse_args().command
    if cmd=='prepare':print('Prepared',len(protocol()['requests']),'follow-up checks.')
    else:globals()[cmd]()
