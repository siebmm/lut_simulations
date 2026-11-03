"""Step 6: explicit flat-plate multiview rays and complex free-wall response.

Coordinates are (x, OD depth). Global S polarization = (-dz, dx).
Modes before/after the image point are kept separate; no defect truth is used.
"""
from itertools import product
import json
from pathlib import Path
import numpy as np
from saft_apertures import reflection, ANGLES
from saft_imaging import sha
from study_paths import STUDY

CATALOG = STUDY.parent/'transducer_experiment/results/full_weld_simulations/sample1_dolphin_bevel/dolphin_paths/paths.json'
PAIRS = [(1,1),(2,1),(1,2),(2,2),(1,3),(3,1),(2,3),(3,2)]
CATALOG_PAIRS = {(1,1),(2,2),(1,3),(3,1),(2,3),(3,2)}
MIN_WALL_RESPONSE = .05


def catalogue():
    raw = json.loads(CATALOG.read_text())
    if {tuple(map(int,p['variant'].split('-'))) for p in raw['paths']} != CATALOG_PAIRS:
        raise ValueError('Dolphin catalogue topology changed; review mapping')
    views=[]
    for n,m in PAIRS:
        pair=f'{n}-{m}'
        labels=sorted({p['mode'] for p in raw['paths'] if p['variant']==pair})
        for letters in product('LS',repeat=n+m):
            mode=''.join(letters);before,after=mode[:n],mode[n:]
            views.append({'id':pair+'_'+mode,'pair':pair,'mode':mode,'before':before,'after':after,
                          'name':before+' | '+after,'reflections':n+m-2,'dolphin_labels':labels,
                          'catalogue_topology':bool(labels),
                          'route':route_labels(n,m)})
    return views,{'path':str(CATALOG.relative_to(STUDY.parent.parent)).replace('\\','/'),
                  'sha256':sha(CATALOG),'labels':sorted({p['mode'] for p in raw['paths']}),
                  'variants':len(raw['paths']),'leg_pairs':[f'{a}-{b}' for a,b in sorted(CATALOG_PAIRS)],
                  'scope':'Repository transcription of supplied geometric catalogue: illustrative offsets and old target coordinates are NOT imported. L/S assignments here are our 132 hypotheses, not verified Dolphin settings.'}


def walls(legs, h, outgoing=False):
    if legs not in [1,2,3]:raise ValueError('Expected 1, 2 or 3 legs on each side')
    return [] if legs==1 else [h] if legs==2 else ([0.,h] if outgoing else [h,0.])


def route_labels(n,m):
    label=lambda z:'OD' if z==0 else 'ID'
    return ['source']+[label(z) for z in walls(n,1.)]+['pixel']+[label(z) for z in walls(m,1.,True)]+['receiver']


def wall_coefficient(p, incident, outgoing, wall_depth, cl, cs):
    """Flux-normalized amplitude in a single consistent GLOBAL polarization basis.

    Reflection routine's local y points into material. At ID this coordinate
    reflection reverses the S basis, giving a minus sign for L/S conversion.
    Propagating modal flux is rho*c*cos(theta)*|amplitude|^2.
    """
    c={"L":cl,"S":cs};p=np.asarray(p)
    angle=np.arcsin(np.clip(p*c[incident],-1,1))
    reflected,_=reflection(angle,incident,cl,cs)
    ratio=np.sqrt(np.maximum(0,1-(p*c[outgoing])**2))/np.maximum(np.sqrt(np.maximum(0,1-(p*c[incident])**2)),1e-30)
    coefficient=reflected[...,0 if outgoing=='L' else 1]*np.sqrt(c[outgoing]/c[incident]*ratio)
    if wall_depth!=0 and incident!=outgoing:coefficient=-coefficient
    return coefficient


def branch(a,b,modes,h,cl,cs,outgoing=False,keep_points=False):
    """Stationary planar branch. All segment horizontal slownesses are equal."""
    if not modes or len(modes)>3 or any(m not in 'LS' for m in modes):raise ValueError('Invalid branch modes')
    a,b=np.broadcast_arrays(np.asarray(a,float),np.asarray(b,float))
    wall=walls(len(modes),h,outgoing)
    heights=[a[...,1]]+[np.full(a.shape[:-1],z) for z in wall]+[b[...,1]]
    hs=[np.abs(z1-z0) for z0,z1 in zip(heights[:-1],heights[1:])]
    if any(np.any(v<=0) for v in hs):raise ValueError('Nonpositive vertical leg')
    c={'L':cl,'S':cs};speeds=[c[m] for m in modes];dx=b[...,0]-a[...,0];distance=np.abs(dx)
    if len(set(modes))==1:
        p=dx/(np.hypot(dx,sum(hs))*speeds[0])
    else:
        lo=np.zeros_like(dx);hi=np.full_like(dx,(1-1e-13)/max(speeds))
        for _ in range(44):
            mid=(lo+hi)/2
            advance=sum(v*(mid*c0)/np.sqrt(1-(mid*c0)**2) for v,c0 in zip(hs,speeds))
            lo=np.where(advance<distance,mid,lo);hi=np.where(advance>=distance,mid,hi)
        p=np.sign(dx)*(lo+hi)/2
    offset=np.zeros_like(dx);time=np.zeros_like(dx);coefficient=np.ones_like(dx,dtype=complex)
    points=[a] if keep_points else None
    for j,(v,c0) in enumerate(zip(hs,speeds)):
        cosine=np.sqrt(np.maximum(0,1-(p*c0)**2))
        offset+=v*p*c0/cosine;time+=v/(c0*cosine)
        if j<len(wall):coefficient*=wall_coefficient(p,modes[j],modes[j+1],wall[j],cl,cs)
        if keep_points:points.append(b if j==len(modes)-1 else np.stack((a[...,0]+offset,heights[j+1]),axis=-1))
    return {'time':time,'p':p,'first_angle':np.arcsin(np.clip(p*speeds[0],-1,1)),
            'last_angle':np.arcsin(np.clip(p*speeds[-1],-1,1)),
            'wall_response':coefficient,'points':points}


def receiving_response(angle,mode,cl,cs):
    # Local OD y is depth, but saved COMSOL vy points upward.
    _,surface=reflection(angle,mode,cl,cs)
    _,normalizer=reflection(np.deg2rad(ANGLES),mode,cl,cs)
    return -surface[...,1]/np.max(abs(normalizer[:,1]))


def response_statistics(samples,mask,q):
    """Same admitted rays for timing-only and known-response weighted images.

    Bounded matched factor conj(q), normalized by sum|q|. Never inverse gain.
    q includes wall flux coefficients and normalized normal surface response.
    """
    weight=np.where(mask,abs(q),0)
    return {'count':mask.sum(axis=-1),'weight':weight.sum(axis=-1),'weight2':(weight**2).sum(axis=-1),
            'raw':np.sum(np.where(mask,samples,0),axis=-1),
            'weighted':np.sum(np.where(mask,samples*np.conj(q),0),axis=-1)}
