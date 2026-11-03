"""Independent Fermat minimisation, segment visibility and reciprocity checks."""
from experiment import *
from scipy.optimize import minimize
sys.path.insert(0,str(ROOT/'studies/laser/src'))
from verify_saft_datasets import PartialReference
from saft_apertures import reflection

def oracle(ref,a,b,via,mode='S',outgoing=False):
    # Millimetres keep the numerical optimisation well conditioned.
    a=np.array(a)*1000;b=np.array(b)*1000;pr=ref.profile*1000;h=ref.h*1000
    if not ref.steel(a/1000) or not ref.steel(b/1000):return None
    family,idx,*order=via.split(':');idx=int(idx);order=order[0] if order else ''
    back=np.array([[-75.,h],[75.,h]])
    faces=[];reflecting=[]
    if family=='wall':
        if order=='id-first':faces.append(back)
        faces.append(pr[idx:idx+2])
        if order=='id-last':faces.append(back)
        edges=np.array([f[1]-f[0] for f in faces]);origins=np.array([f[0] for f in faces])
        def objective(t):
            q=np.vstack([a,origins+t[:,None]*edges,b]);d=np.diff(q,axis=0);lens=np.linalg.norm(d,axis=1)
            u=d/np.maximum(lens[:,None],1e-30)
            return lens.sum(),np.sum((u[:-1]-u[1:])*edges,axis=1)
        opt=minimize(objective,np.full(len(faces),.5),jac=True,bounds=[(0,1)]*len(faces),method='L-BFGS-B',options=dict(ftol=1e-15,gtol=1e-10,maxiter=1000,maxls=40))
        if np.min(opt.x)<1e-7 or np.max(opt.x)>1-1e-7 or max(abs(objective(opt.x)[1]))>1e-5:return None
        q=np.vstack([a,origins+opt.x[:,None]*edges,b]);reflecting=list(zip(range(1,len(faces)+1),faces))
    else:
        corner=pr[idx]
        def bounce(u,v):
            t=(h-u[1])/(2*h-u[1]-v[1])
            return np.array([u[0]+t*(v[0]-u[0]),h])
        if order=='id-first':q=np.array([a,bounce(a,corner),corner,b]);reflecting=[(1,back)]
        elif order=='id-last':q=np.array([a,corner,bounce(corner,b),b]);reflecting=[(2,back)]
        else:q=np.array([a,corner,b])
    d=np.diff(q,axis=0);lens=np.linalg.norm(d,axis=1)
    if np.min(lens)<1e-7 or not all(ref.visible(u/1000,v/1000) for u,v in zip(q[:-1],q[1:])):return None
    if (not outgoing and d[0,1]<-1e-7) or (outgoing and d[-1,1]>1e-7):return None
    wall=1.+0j
    for i,face in reflecting:
        tangent=(face[1]-face[0])/np.linalg.norm(face[1]-face[0]);normal=np.array([-tangent[1],tangent[0]])
        if np.allclose(face[:,1],h):normal=-normal
        inc=d[i-1]/lens[i-1];out=d[i]/lens[i]
        if inc@normal>=-1e-9 or out@normal<=1e-9:return None
        if np.linalg.norm(out-inc+2*(inc@normal)*normal)>1e-6:return None
        wall*=reflection(np.array([np.arccos(np.clip(-inc@normal,0,1))]),mode,ref.cl,ref.cs)[0][0,int(mode=='S')]
    return dict(points=q/1000,time=lens.sum()/1000/(ref.cs if mode=='S' else ref.cl),wall=wall)

def main():
    ref=PartialReference(DATA);checks=[];expected=[]
    endpoints=list(ref.sources)+list(ref.xy[[0,78,156,157,165,173,330]])
    for x in [-.005,0,.005]:
      for z in [.017,.021,.025]:
       for endpoint in endpoints:
        for family,ids in [('wall',[2,4]),('corner',range(1,7))]:
         for idx in ids:
          for order in ['',':id-first',':id-last']:
           via=f'{family}:{idx}{order}';mode='S' if x else 'L'
           a=endpoint.tolist();b=[x,z];modes=mode*(3 if order else 2)
           value=oracle(ref,a,b,via,mode)
           checks.append(dict(a_m=a,b_m=b,modes=modes,via=via));expected.append(value)
           reverse=via.replace('id-first','id-last') if 'id-first' in via else via.replace('id-last','id-first')
           checks.append(dict(a_m=b,b_m=a,modes=modes,via=reverse,outgoing=True))
           expected.append(None if value is None else dict(points=value['points'][::-1],time=value['time'],wall=value['wall']))
    # A route through the open groove and a degenerate corner leg must fail.
    for a,b,via in [([0,0],[0,.021],'wall:2'),(ref.profile[3].tolist(),[0,.021],'corner:3')]:
        checks.append(dict(a_m=a,b_m=b,modes='SS',via=via));expected.append(None)
    cfg=copy.deepcopy(BASE);cfg.update(path='1-2_SSS',cases=['reference'],threads=1)
    request=OUT/'verification/request.json';write(request,dict(ray_checks=checks,requests=[dict(selection=cfg,points_m=[[0,.021]])]))
    run_request(request,OUT/'verification/native')
    actual=json.loads((OUT/'verification/native/geometry.json').read_text())['rays']
    errors=[];times=[];walls=[];mismatches=[]
    for i,(a,b) in enumerate(zip(actual,expected)):
        if a['valid']!=(b is not None):mismatches.append(dict(index=i,request=checks[i],native=a,expected=None if b is None else {k:(v.tolist() if isinstance(v,np.ndarray) else str(v)) for k,v in b.items()}));continue
        if b:
            errors.append(float(np.max(abs(np.array(a['points_m'])-b['points'])))*1000)
            times.append(abs(a['time_s']-b['time'])*1e6)
            walls.append(abs(complex(a['wall_real'],a['wall_imag'])-b['wall']))
    report=dict(rays=len(actual),valid=sum(a['valid'] for a in actual),max_vertex_error_mm=max(errors),max_time_error_us=max(times),max_complex_wall_error=max(walls),mismatches=mismatches)
    write(OUT/'verification/geometry_verification.json',report)
    print(json.dumps(report if not mismatches else dict(rays=len(actual),mismatches=mismatches[:5]),indent=2),flush=True)
    assert not mismatches and max(errors)<2e-5 and max(times)<1e-5 and max(walls)<2e-5
if __name__=='__main__':main()
