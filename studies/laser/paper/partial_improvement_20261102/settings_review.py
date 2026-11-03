"""Declared settings sensitivities for four partial-weld acquisition recipes."""
import argparse
import itertools
import json
import csv
import numpy as np
import review as r

D=r.DEST/'settings'


def anchors():
    b=r.base()
    return dict(
        broad_SSS=dict(b,path='1-2_SSS',minimum_sources=2),
        left_SSS=dict(b,path='1-2_SSS',sources=[1,2,3],groups=[True,False,False],minimum_sources=2),
        right_SLL=dict(b,path='1-2_SLL',sources=[3],groups=[False,False,True]),
        right_SSLSS=dict(b,path='3-2_SSLSS',sources=[1],groups=[False,False,True]))


def variations():
    v=[('baseline',{})]
    v += [(f'band{lo:g}_{hi:g}',dict(low_MHz=lo,high_MHz=hi)) for lo,hi in [(.5,3),(.5,5),(1,3),(1.5,5),(2,5)]]
    v += [(f'angular{g:g}',dict(angular_gate=g)) for g in (0,.05,.1,.2,.3)]
    v += [(f'wall{g:g}',dict(wall_gate=g)) for g in (0,.025,.1)]
    v += [(f'offset{g:g}',dict(min_offset_mm=g)) for g in (0,.75,3,5)]
    v += [('equal',dict(processing='raw')),('tangential',dict(component='vx'))]
    v += [('screen_'+x,dict(pulse_screening=x)) for x in ('off','L','S','R','LS','LR','SR')]
    v += [(f'reflect{bits:04b}',dict(pulse_reflections=[bool(bits&(1<<k)) for k in range(4)])) for bits in range(1,16)]
    v += [(f'tail{db}',dict(screen_tail_db=db)) for db in (20,40)]
    v += [(f'margin{t:g}',dict(screen_margin_us=t)) for t in (0,.05,.2)]
    v += [(f'mute{k}',dict(mutes=[i==k for i in range(6)])) for k in range(6)]
    v += [('mute_all',dict(mutes=[True]*6))]
    v += [(f'mute_all_width{w:g}',dict(mutes=[True]*6,mute_halfwidth_us=w)) for w in (.3,1.2)]
    v += [('aperture_'+x,dict(aperture=x)) for x in ('left','right','even','odd')]
    v += [('clean_case',dict(clean=True)),('clean_paired',dict(clean=True,repairs='paired'))]
    v += [('frozen_centre',dict(freeze=True,anchor_x_mm=0,anchor_depth_mm=20))]
    return v


def protocol():
    jobs=[]
    for anchor,c in anchors().items():
        for label,change in variations():
            jobs.append(dict(id=f'j{len(jobs):03d}',anchor=anchor,label=label,selection=dict(c,**change)))
    p=dict(purpose='Inspect every physical/processing control family across four established partial-weld recipes.',
        scope='A declared one-factor sensitivity matrix, plus all 16 reflected-screen combinations. Continuous ranges and arbitrary receiver masks have infinitely many combinations; this is not a global optimization or an exhaustive Cartesian search.',
        source_choices='All shots, three weld shots with left receivers, one right weld shot with right receivers, one left weld shot with right receivers.',
        diagnostics='vx is tangential motion, outside the assumed surface-normal optical receiver. Paired repair rejection uses matched-reference repair metadata. Frozen aperture uses the fixed weld-centre anchor (0,20) mm, not the flaw.',
        offline='Receiver fraction/count/effective floors, coherence exponent, minimum source count, fixed-source policy and coherent versus RMS source combination are recomputed from the exported sufficient statistics, then selected cases are checked natively.',
        source_hashes=r.protocol()['source_hashes'],requests=jobs)
    path=D/'protocol.json'
    if path.exists():assert json.loads(path.read_text())==p
    else:r.s.write(path,p)
    return p


def run():
    p=protocol();done={f.parent.name for f in D.glob('b*/*/arrays.npz')}
    missing=[j for j in p['requests'] if j['id'] not in done];offset=len(list(D.glob('b*')))
    for start in range(0,len(missing),10):
        r.execute(missing[start:start+10],D/f'b{offset+start//10:03d}')
        print('settings',len(done)+start+len(missing[start:start+10]),'/',len(p['requests']),flush=True)
    analyse()


def locate(job_id):
    matches=list(D.glob('b*/'+job_id+'/arrays.npz'))
    assert len(matches)==1,(job_id,matches)
    return matches[0]


def measure(a,c,label,anchor,job_id,changes=None):
    changes=changes or {};c=dict(c,**changes)
    points=a['points'].reshape(-1,2);d=r.s.distances(points,'partial');n=len(points);ns=len(c['sources'])
    k=len(c['sources']) if c['policy']=='fixed' else 1 if ns==1 else c['minimum_sources']
    row=dict(id=job_id,anchor=anchor,label=label);values={}
    for case in ('lof','reference','difference'):
        stats=a[case+'_2_stats'].reshape(ns,n,7);sums=a[case+'_2_sums'].reshape(ns,n)
        v,count,valid,means=r.reconstruct(stats,sums,list(range(ns)),k,c);values[case]=v
        row.update({case+'_'+key:value for key,value in r.metric(v,points,d).items()})
        if not changes:
            np.testing.assert_allclose(v,a[case+'_2_amplitude'],rtol=2e-8,atol=2e-15,equal_nan=True)
            np.testing.assert_array_equal(count,a[case+'_2_count'])
    common=(d<=.5)&np.isfinite(values['lof'])&np.isfinite(values['reference'])
    l,q=values['lof'][common],values['reference'][common]
    row['lof_over_pristine_near_db']=float(20*np.log10(l.max()/q.max())) if common.any() and l.max()>0 and q.max()>0 else None
    return row


def analyse():
    p=protocol();rows=[];offline=[]
    for job in p['requests']:
        with np.load(locate(job['id'])) as f:a=dict(f)
        rows.append(measure(a,job['selection'],job['label'],job['anchor'],job['id']))
        if job['label'] not in ('baseline','reflect1111'):continue
        changes=[]
        changes += [(f'fraction{q:g}',dict(receiver_fraction=q)) for q in (0,.1,.5,1)]
        changes += [(f'count{n}',dict(minimum_receivers=n)) for n in (1,4,16,32)]
        changes += [(f'effective{n}',dict(minimum_effective_receivers=n)) for n in (0,4,16,32)]
        changes += [('guards_relaxed',dict(receiver_fraction=0,minimum_receivers=1,minimum_effective_receivers=0))]
        changes += [(f'coherence{q:g}',dict(coherence_power=q)) for q in (0,.5,2)]
        changes += [('coherent_sources',dict(combination='coherent')),('all_sources_required',dict(policy='fixed'))]
        if len(job['selection']['sources'])>1:
            changes += [(f'minimum_sources{k}',dict(minimum_sources=k)) for k in range(1,len(job['selection']['sources'])+1) if k!=job['selection']['minimum_sources']]
        for label,change in changes:
            row=measure(a,job['selection'],job['label']+'__'+label,job['anchor'],job['id'],change)
            row['changes']=change;offline.append(row)
    for name,data in [('native',rows),('offline',offline)]:
        r.s.write(D/(name+'_metrics.json'),dict(rows=data))
        flat=[{k:v for k,v in row.items() if k!='changes'} for row in data]
        with (D/(name+'_metrics.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    r.s.write(D/'verification.json',dict(status='passed',native_configurations=len(p['requests']),source_recombination_cases=len(rows)*3,offline_views=len(offline)))
    for anchor in anchors():
        ranked=sorted([x for x in rows+offline if x['anchor']==anchor and x['lof_cmax_db'] is not None and x['lof_near_support_pct']>=40],key=lambda x:x['lof_cmax_db'],reverse=True)
        print(anchor,'top settings:',flush=True)
        for x in ranked[:10]:print(x['label'],round(x['lof_cmax_db'],2),round(x['lof_near_support_pct'],1),round(x['lof_peak_distance_mm'],3),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('command',choices=['prepare','run','analyse']);c=ap.parse_args().command
    if c=='prepare':print('Prepared',len(protocol()['requests']),'settings configurations.')
    else:globals()[c]()
