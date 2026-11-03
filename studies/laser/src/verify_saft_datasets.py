"""Independent partial-weld geometry, focusing, channel mapping and folder checks."""

# Locate portable scientific inputs without depending on the source checkout.
import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import data_root as _data_root, native_root as _native_root, executable as _executable, native_env as _native_env
from datetime import datetime, timezone
from itertools import product
import hashlib
import json
from pathlib import Path
import subprocess
import numpy as np
from scipy.optimize import brentq
from saft_imaging import kernel, filter_record
from saft_apertures import ANGLES, directivity, angular_gain, reflection
from study_paths import REPOSITORY


class PartialReference:
    def __init__(self, folder):
        self.folder = folder
        self.meta = json.loads((folder / 'dataset.json').read_text())
        self.profile = np.array(self.meta['surface_profile_x_depth_m'])
        self.audit = self.meta['audit']
        m = self.meta['material']
        self.h, self.cl, self.cs = m['thickness_m'], m['cl_m_s'], m['cs_m_s']
        self.sources = np.array([[s['xy_m'][0], self.h-s['xy_m'][1]] for s in self.audit['sources']])
        self.receivers = sorted(self.audit['receivers'], key=lambda r:int(r['receiver_id'][1:]))
        self.xy = np.array([[r['x_m'], self.h-r['y_m']] for r in self.receivers])
        self.signals = {}
        sp = self.audit['source_parameters']
        self.curves = directivity(('broad', .5, 5), self.cl, self.cs, sp['line_half_width_m'], sp['pulse_fwhm_s'])

    def depth(self, x):
        exact = self.profile[np.isclose(self.profile[:, 0], x, atol=1e-13, rtol=0)]
        if len(exact): return exact[:, 1].min()
        for a, b in zip(self.profile[:-1], self.profile[1:]):
            if a[0] < x < b[0]: return a[1]+(b[1]-a[1])*(x-a[0])/(b[0]-a[0])
        return np.inf

    def steel(self, p):
        return self.depth(p[0])-1e-10 <= p[1] <= self.h+1e-10

    def visible(self, a, b):
        if not self.steel(a) or not self.steel(b): return False
        if abs(b[0]-a[0]) < 1e-15: return True
        # Test inside each traversed surface interval, as well as its ends.
        cuts = sorted({0., 1., *[(x-a[0])/(b[0]-a[0]) for x in self.profile[:, 0]
                                 if 0 < (x-a[0])/(b[0]-a[0]) < 1]})
        probes = cuts + [(a+b)/2 for a,b in zip(cuts[:-1], cuts[1:])]
        return all(self.steel(a+t*(b-a)) for t in probes)

    def ray(self, a, b, modes, outgoing=False):
        a, b = np.array(a), np.array(b)
        if not self.steel(a) or not self.steel(b): return None
        candidates = []
        for top in ([0., .0164] if len(modes) == 3 else [0.]):
            if len(modes) == 3 and (a if outgoing else b)[1] <= top+1e-14: continue
            zs = [a[1]] + ([self.h] if len(modes)==2 else
                  ([top, self.h] if outgoing else [self.h, top]) if len(modes)==3 else []) + [b[1]]
            hs = abs(np.diff(zs)); speeds = np.array([self.cl if m=='L' else self.cs for m in modes])
            dx = b[0]-a[0]
            if len(modes)==1 and hs[0]<1e-14:
                if abs(dx)<1e-15 or not self.visible(a,b): continue
                return dict(time=abs(dx)/speeds[0], points=np.array([a,b]), first=np.sign(dx)*np.pi/2,
                            last=np.sign(dx)*np.pi/2, wall=1.+0j)
            if (hs<=0).any(): continue
            # Solve in dimensionless horizontal slowness; independent Brent solver.
            advance = lambda u: np.sum(hs*(u*speeds/speeds.max())/np.sqrt(1-(u*speeds/speeds.max())**2))
            u = brentq(lambda u:advance(u)-abs(dx), 0., 1-1e-14, xtol=2e-15)
            pc = np.sign(dx)*u*speeds/speeds.max(); co = np.sqrt(1-pc**2)
            xs = np.r_[a[0], a[0]+np.cumsum(hs*pc/co)]
            vertices = np.column_stack((xs, zs)); vertices[-1] = b
            if len(modes)==3:
                x = vertices[1 if outgoing else 2, 0]
                if not any(l[0]-1e-10<=x<=r[0]+1e-10 and r[0]>l[0] and
                           abs(l[1]-top)<1e-12 and abs(r[1]-top)<1e-12
                           for l,r in zip(self.profile[:-1],self.profile[1:])): continue
            if not all(self.visible(l,r) for l,r in zip(vertices[:-1],vertices[1:])): continue
            wall = 1.+0j
            for i in range(len(modes)-1):
                amp = reflection(np.array([np.arcsin(pc[i])]), modes[i], self.cl, self.cs)[0][0]
                q = amp[0 if modes[i+1]=='L' else 1]*np.sqrt(speeds[i+1]/speeds[i]*co[i+1]/co[i])
                if zs[i+1]==self.h and modes[i]!=modes[i+1]: q = -q
                wall *= q
            candidates.append(dict(time=np.sum(hs/(speeds*co)), points=vertices,
                                   first=np.arcsin(pc[0]), last=np.arcsin(pc[-1]), wall=wall))
        return min(candidates,key=lambda r:r['time']) if candidates else None

    def filtered(self, s, case, component):
        key = (s,case,component)
        if key not in self.signals:
            if case=='difference':
                self.signals[key] = self.filtered(s,'lof',component)-self.filtered(s,'reference',component)
            else:
                data = np.empty((331,5001))
                for group in ('left','weld','right'):
                    ids = [i for i,r in enumerate(self.receivers) if r['group']==group]
                    name = self.audit['sources'][s]['id']
                    table = np.loadtxt(self.folder/'sources'/name/case/'receivers'/f'{group}_{component}.csv', delimiter=',',skiprows=1)
                    data[ids] = table[:,1:].T
                self.signals[key] = filter_record(data,kernel(.5,5))
        return self.signals[key]

    def focus(self, path, points, component='vy'):
        n1,n2 = map(int,path.split('_')[0].split('-')); modes = path.split('_')[1]
        first,last = modes[:n1],modes[n1:]
        column = int(component=='vy')
        norm = max(abs(reflection(np.deg2rad(ANGLES),last[-1],self.cl,self.cs)[1][:,column]))
        shape = (5,len(points)); counts=np.zeros(shape); effective=np.zeros(shape)
        valid=np.zeros(shape,bool); means={case:np.zeros(shape,complex) for case in ('lof','reference','difference')}
        used=np.zeros((*shape,331),bool)
        for si,source in enumerate(self.sources):
            for pi,p in enumerate(points):
                inc = self.ray(source,p,first)
                if inc is None: continue
                gain = angular_gain(np.array([inc['first']]),self.curves[first[0]]['source'])[0]
                if gain<.2: continue
                rows=[]; eligible=0
                for ri,receiver in enumerate(self.xy):
                    rec=self.ray(p,receiver,last,True)
                    if rec is None: continue
                    surface=(-1 if column else 1)*reflection(np.array([rec['last']]),last[-1],self.cl,self.cs)[1][0,column]/norm
                    wall=inc['wall']*rec['wall']
                    if abs(wall)<.05 or abs(surface)<.2 or abs(np.rad2deg(rec['last']))>80 or np.linalg.norm(receiver-source)<.002-1e-14: continue
                    eligible+=1
                    time=self.audit['sources'][si]['source_trigger_s']+inc['time']+rec['time']
                    sample=time/1e-8; j=int(np.floor(sample)); alpha=sample-j
                    if j-328<0 or j+1+328>=5001: continue
                    rows.append((ri,j,alpha,wall*surface))
                if not rows: continue
                ids=np.array([r[0] for r in rows]); js=np.array([r[1] for r in rows]); alpha=np.array([r[2] for r in rows]); response=np.array([r[3] for r in rows])
                weights=abs(response); counts[si,pi]=len(rows); effective[si,pi]=weights.sum()**2/(weights**2).sum()
                valid[si,pi]=len(rows)>=max(8,np.ceil(.25*eligible)) and effective[si,pi]>=8-1e-10
                used[si,pi,ids]=True
                for case in means:
                    signal=self.filtered(si,case,component)
                    value=(1-alpha)*signal[ids,js]+alpha*signal[ids,js+1]
                    means[case][si,pi]=(value*np.conj(response)).sum()/weights.sum()
        source_count=valid.sum(axis=0); supported=source_count>=2
        maps={case:np.where(supported,np.sqrt(np.sum(np.where(valid,abs(m)**2,0),axis=0)/np.maximum(source_count,1)),np.nan) for case,m in means.items()}
        usage=(used & valid[:,:,None] & supported[None,:,None]).sum(axis=1)
        return dict(maps=maps,means=means,valid=valid,usage=usage,source_count=source_count)


def main():
    exe=_executable('SAFTBatch')
    if not exe.is_file():
        raise SystemExit('Build the native executables first; see README.md.')
    folder=_data_root()/'partial_weld'; ref=PartialReference(folder)
    out=REPOSITORY/'output/physics_check'; out.mkdir(parents=True,exist_ok=True)
    points=np.array([[2.5,20],[2.7,19],[0,17.5],[-2,21],[4,12],[-10,4],[10,23],[0,10],[3.1,8],[-3.2,9]])/1000
    modes=[''.join(m) for n in (1,2,3) for m in product('LS',repeat=n)]
    checks=[]; truth=[]
    for a in ref.sources:
        for p in points:
            for m in modes:
                for outgoing in (False,True):
                    aa,bb=(p,a) if outgoing else (a,p)
                    checks.append(dict(a_m=aa.tolist(),b_m=bb.tolist(),modes=m,outgoing=outgoing))
                    truth.append(ref.ray(aa,bb,m,outgoing))
    paths=['1-1_LL','1-1_SS','1-2_SSS','2-2_LSSL']
    jobs=[dict(selection=dict(path=p,cases=['lof','reference','difference']),points_m=points.tolist()) for p in paths]
    jobs.append(dict(selection=dict(path='1-1_LL',component='vx',cases=['lof','reference','difference']),points_m=points.tolist()))
    jobs.extend(dict(selection=dict(path=p,cases=['lof','difference'])) for p in ('1-1_LL','1-1_SS','2-2_LSSL'))
    jobs.append(dict(selection=dict(path='1-1_SS',pulse_screening='LSR',cases=['lof']),forbid_reference=True,clear_memory=True,
                     inspect=dict(source=2,receiver=314,has_point=True,x_mm=2.5,depth_mm=20)))
    jobs.append(dict(selection=dict(cases=['lof'],mutes=[False,False,False,False,True,False]),
                     inspect=dict(source=0,receiver=314,has_point=False)))
    request=out/'requests.json'; request.write_text(json.dumps(dict(requests=jobs,ray_checks=checks,verify_inputs=True)))
    subprocess.run([str(exe),'--repository',str(REPOSITORY),'--dataset',str(folder),'--compute',str(request),'--output',str(out)],env=_native_env(),check=True)
    actual=json.loads((out/'geometry.json').read_text())['rays']; errors=[]
    assert len(actual)==len(truth), 'Incomplete native ray export'
    for i,(got,want) in enumerate(zip(actual,truth)):
        assert got['valid']==(want is not None),(i,checks[i],got,want)
        if want is not None:
            errors.append(abs(got['time_s']-want['time']))
            np.testing.assert_allclose(got['time_s'],want['time'],rtol=1e-10,atol=1e-14)
            np.testing.assert_allclose(got['points_m'],want['points'],rtol=1e-9,atol=1e-11)
    amplitude_errors=[]
    for i,p in enumerate(paths+['1-1_LL']):
        print(f'Independent partial focusing: {p}, {"vx" if i==4 else "vy"}', flush=True)
        expected=ref.focus(p,points,'vx' if i==4 else 'vy')
        with np.load(out/str(i)/'images.npz') as data:
            for case in expected['maps']:
                for v in range(3):
                    np.testing.assert_allclose(data[case+'_maps'][v],expected['maps'][case],rtol=1e-8,atol=2e-14,equal_nan=True)
                    np.testing.assert_allclose(data[case+'_means'][v],expected['means'][case],rtol=1e-8,atol=2e-14)
                    np.testing.assert_array_equal(data[case+'_counts'][v],expected['source_count'])
                    np.testing.assert_array_equal(data[case+'_receiver_usage'][v],expected['usage'])
                    delta=abs(data[case+'_maps'][v]-expected['maps'][case]); amplitude_errors.extend(delta[np.isfinite(delta)].tolist())
    for i in range(5,9):
        with np.load(out/str(i)/'images.npz') as data:
            air=np.array([not ref.steel(p) for p in data['points_m']])
            assert np.isnan(data['lof_maps'][:,air]).all()
            assert data['lof_receiver_usage'].shape==(3,5,331)
    provenance=json.loads((out/'8/selection.json').read_text())
    previous=json.loads((out/'7/selection.json').read_text())
    refs=lambda p:{k:v for k,v in p['recordings_accessed'].items() if '/reference/' in k}
    assert refs(previous)==refs(provenance)
    trace=np.genfromtxt(out/'8/trace.csv',delimiter=',',names=True)
    raw=np.loadtxt(folder/'sources/bead_centre/lof/receivers/weld_vy.csv',delimiter=',',skiprows=1)[:,1]
    np.testing.assert_array_equal(trace['raw_m_s'],raw)
    # Unequal source/receiver depths break the LS/SL travel-time equality.
    nuisance=json.loads((out/'benchmark.json').read_text())['runs'][9]['background_arrival_us']
    trace=np.genfromtxt(out/'9/trace.csv',delimiter=',',names=True)
    centres=[nuisance[k][314] for k in (4,6)]
    assert abs(centres[0]-centres[1])>.6
    for centre in centres:
        assert trace['arrival_keep'][int(round(centre/.01))]==0
    # A separate empty application root proves that original run metadata,
    # models and other simulation outputs are not needed by a dataset folder.
    standalone=out/'isolated_root'; standalone.mkdir(exist_ok=True)
    one=out/'isolated.json'; one.write_text(json.dumps(dict(selection=dict(sources=[2]),points_m=[[.0025,.02]])))
    subprocess.run([str(exe),'--repository',str(standalone),'--dataset',str(folder),'--compute',str(one),'--output',str(out/'isolated')],env=_native_env(),check=True)
    invalid=out/'invalid.json'; invalid.write_text(json.dumps(dict(selection=dict(receivers=[331]))))
    rejected=subprocess.run([str(exe),'--repository',str(standalone),'--dataset',str(folder),'--compute',str(invalid),'--output',str(out/'invalid')],env=_native_env(),capture_output=True)
    assert rejected.returncode!=0 and b'R330' in rejected.stderr, rejected.stderr
    for name,count in [('full_weld',337),('partial_weld',331)]:
        d=_data_root()/name
        assert len(list(d.rglob('*.csv')))==60 and len(list(p for p in d.rglob('*') if p.is_file()))==61
        assert len(json.loads((d/'dataset.json').read_text())['audit']['receivers'])==count
    report=dict(status='passed',created_utc=datetime.now(timezone.utc).isoformat(),geometry_rays=len(checks),
        valid_rays=sum(t is not None for t in truth),max_time_error_s=max(errors),max_amplitude_error_m_s=max(amplitude_errors),
        focus_cases=5,full_grid_cases=5,receiver_count=331,checks=['Independent Brent ray times/vertices and groove visibility',
        'Independent NumPy focusing: LL, SS, S|SS, LS|SL and vx; LOF, reference and difference',
        'Exact admitted receiver usage counts','Full-grid air masks','LOF-only pulse screening forbids reference reads',
        'Recessed receiver raw-channel mapping','Unequal-depth LS and SL backwall muting',
        'Dataset runs with empty application root; receiver R331 rejected',
        'Two minimal 61-file datasets; original inputs preserved'],
        implementation_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sorted(_native_root().iterdir()) if p.suffix in ('.cpp','.hpp')},
        executable_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
    (out/'check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report|dict(implementation_sha256='recorded'),indent=2),flush=True)


if __name__=='__main__': main()
