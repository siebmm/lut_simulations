"""Independent homogeneous ray tracing against the actual exposed surface."""

import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import executable as _executable, native_env as _native_env
from study import *

def profile(kind):return np.array(dataset(kind)['surface_profile_x_depth_m'])*1000
def depth_at(x,kind):
    p=profile(kind);same=np.isclose(p[:,0],x,atol=1e-9,rtol=0)
    if same.any():return p[same,1].min()
    return np.interp(x,p[:,0],p[:,1],left=np.inf,right=np.inf)
def visible(a,b,kind):
    h=dataset(kind)['material']['thickness_m']*1000
    if any(p[1]<depth_at(p[0],kind)-1e-7 or p[1]>h+1e-7 for p in [a,b]):return False
    if abs(b[0]-a[0])<1e-12:return True
    for x,z in profile(kind):
        t=(x-a[0])/(b[0]-a[0])
        if 1e-12<t<1-1e-12 and a[1]+t*(b[1]-a[1])<z-1e-7:return False
    return True
def traces(a,b,modes,kind,outgoing=False):
    a=np.array(a);b=np.array(b);h=dataset(kind)['material']['thickness_m']*1000
    assert len(set(modes))==1
    speed=dataset(kind)['material']['cs_m_s' if modes[0]=='S' else 'cl_m_s']/1000
    n=len(modes);pr=profile(kind);tops=sorted(set(pr[:-1,1][(np.diff(pr[:,0])>0)&(abs(np.diff(pr[:,1]))<1e-9)])) if n>=3 else [0]
    candidates=[]
    for top in tops:
        if n==3 and (a[1] if outgoing else b[1])<=top+1e-11:continue
        heights=[a[1]]+([] if n==1 else [h] if n==2 else ([top,h] if outgoing else [h,top]) if n==3 else [h,top,h])+[b[1]]
        dz=np.abs(np.diff(heights));total=dz.sum()
        if total==0:
            points=np.array([a,b]);length=abs(b[0]-a[0])
        else:
            points=np.column_stack([a[0]+(b[0]-a[0])*np.r_[0,np.cumsum(dz)]/total,heights]);length=np.hypot(b[0]-a[0],total)
        if any(not visible(u,v,kind) for u,v in zip(points[:-1],points[1:])):continue
        if n>=3:
            bounce=points[1 if n==3 and outgoing else 2]
            flat=any(q[0]>p[0] and abs(p[1]-top)<1e-9 and abs(q[1]-top)<1e-9 and p[0]-1e-7<=bounce[0]<=q[0]+1e-7 for p,q in zip(pr[:-1],pr[1:]))
            if not flat or abs(depth_at(bounce[0],kind)-top)>1e-7:continue
        candidates.append({'time_us':float(length/speed),'points_mm':points.tolist(),'top_mm':top})
    return sorted(candidates,key=lambda c:c['time_us'])

def check():
    folder=OUT/'geometry_checks';folder.mkdir(exist_ok=True)
    exe=_executable('SAFTBatch')
    summary={}
    for kind in DATA:
        d=dataset(kind);h=d['material']['thickness_m']*1000
        sources=[[s['x_mm'],h-s['xy_m'][1]*1000] for s in d['audit']['sources']]
        rids=[0,78,156,157,165,173,len(d['audit']['receivers'])-1]
        receivers=[[d['audit']['receivers'][r]['x_m']*1000,h-d['audit']['receivers'][r]['y_m']*1000] for r in rids]
        jobs=[];expected=[]
        for x in [-5,0,5]:
            for z in [17,21,25]:
                point=[x,z]
                for a,b,outgoing in [(s,point,False) for s in sources]+[(point,r,True) for r in receivers]:
                    for m in ['S','SS','SSS','L','LL','LLL']:
                        jobs.append({'a_m':(np.array(a)/1000).tolist(),'b_m':(np.array(b)/1000).tolist(),'modes':m,'outgoing':outgoing})
                        candidate=traces(a,b,m,kind,outgoing);expected.append(candidate[0] if candidate else None)
        config=copy.deepcopy(index(kind)[0]['selection']);config['cases']=['reference'];config['threads']=1
        request=folder/f'{kind}_request.json';write(request,{'ray_checks':jobs,'requests':[{'selection':config,'points_m':[[0,.021]]}]})
        p=subprocess.run([str(exe),'--dataset',str(DATA[kind]),'--compute',str(request),'--output',str(folder/kind)],cwd=ROOT,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        assert p.returncode==0,p.stdout
        actual=json.loads((folder/kind/'geometry.json').read_text())['rays'];errors=[];coordinates=[]
        for a,b in zip(actual,expected):
            assert a['valid']==bool(b),(a,b)
            if b:
                errors.append(abs(a['time_s']*1e6-b['time_us']));coordinates.append(np.max(abs(np.array(a['points_m'])*1000-b['points_mm'])))
        assert max(errors)<1e-8 and max(coordinates)<1e-7
        summary[kind]={'rays':len(actual),'valid_rays':sum(a['valid'] for a in actual),'max_time_error_us':max(errors),'max_vertex_error_mm':float(max(coordinates))}
    write(OUT/'geometry_verification.json',summary);print(json.dumps(summary,indent=2))
if __name__=='__main__':check()
