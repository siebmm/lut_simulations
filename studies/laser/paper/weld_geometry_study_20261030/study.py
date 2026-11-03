"""Predeclared geometry-led full/partial-weld SAFT study. No LOF-based tuning."""

# Locate portable scientific inputs without depending on the source checkout.
import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import data_root as _data_root, native_root as _native_root, executable as _executable, native_env as _native_env
from pathlib import Path
import copy,hashlib,json,math,subprocess,sys,datetime,struct
import numpy as np
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[3]
OUT=ROOT/'studies/laser/results/weld_geometry_study_20261030'
DATA={kind:_data_root()/f'{kind}_weld' for kind in ['full','partial']}
PATHS=['1-1_SS','1-1_LL','1-1_LS','1-1_SL','1-2_SSS','2-1_SSS','1-2_LLL','2-1_LLL',
       '1-3_SSSS','2-2_SSSS','3-1_SSSS','1-3_LLLL','2-2_LLLL','3-1_LLLL',
       '2-3_SSSSS','3-2_SSSSS','2-3_LLLLL','3-2_LLLLL']
CORE=['1-1_SS','1-1_LL','1-2_SSS','2-1_SSS']
WIDE=[p for p in PATHS if set(p.split('_')[1])=={'S'} or p=='1-1_LL']
def label(path):
    split=int(path.split('-')[0]);modes=path.split('_')[1]
    return modes[:split]+' | '+modes[split:]
def dataset(kind):return json.loads((DATA[kind]/'dataset.json').read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,a):p.write_text(json.dumps(a,indent=2),encoding='utf-8')

def prepare():
    OUT.mkdir(parents=True,exist_ok=True)
    d=dataset('partial')
    # Design uses acquisition/surface metadata, never defect coordinates or image scores.
    floor=max(s['xy_m'][1] for s in d['audit']['sources'] if s['surface']=='weld')
    h=d['material']['thickness_m']*1000;floor_depth=h-floor*1000
    exclusion=math.ceil(max(s['footprint_half_width_m']*1000 for s in d['audit']['sources'])/.25)*.25
    lateral=math.ceil(2*max(abs(p[0]*1000) for p in d['surface_profile_x_depth_m'][1:-1])/.5)*.5
    zmin=round(floor_depth+.25,8);zmax=round(zmin+math.floor((h-.1-zmin)/.1+1e-8)*.1,8)
    base=dict(anchor_x_mm=0,anchor_depth_mm=0,angular_gate=.2,aperture='all',cases=['lof'],clean=False,
      coherence=True,coherence_power=1,combination='rms',component='vy',depth_min_mm=zmin,depth_max_mm=zmax,
      freeze=False,groups=[True]*3,high_MHz=5,low_MHz=.5,min_offset_mm=exclusion,minimum_effective_receivers=8,
      minimum_receivers=8,minimum_sources=2,mute_halfwidth_us=.6,mutes=[False]*6,policy='minimum',processing='raw',
      pulse_reflections=[False]*4,pulse_screening='LSR',receiver_fraction=.25,repairs='case',screen_margin_us=.1,
      screen_tail_db=30,sources=list(range(5)),spacing_mm=.1,threads=8,wall_gate=.05,x_min_mm=-lateral,x_max_mm=lateral)
    protocol={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'common_controls':base,'paths':PATHS,
      'selection_basis':'Frozen before this study reconstructs LOF images. Geometry/acquisition and established processing physics only. Historical flaw location is known; not a blind experiment.',
      'source_exclusion_basis':'Maximum prescribed source footprint half-width 1.425 mm, rounded up to the 0.25 mm acquisition increment: 1.5 mm.',
      'extent_basis':'Twice the groove half-opening rounded out to 0.5 mm; first centre 0.25 mm below the fill surface; final centre at least 0.1 mm above backwall.',
      'path_basis':'Every pure-L and pure-S split available with 2 to 5 total legs and at most 3 legs per branch; add both direct mode-converted scattering hypotheses.',
      'threshold_status':'20% angular, 5% wall, 8 receivers, 25% retained and 2 sources are predeclared engineering guards, not unique values derived from elasticity or optimized against the LOF.',
      'screening_basis':'All valid direct L/S/R families with common 30 dB tail, 0.1 us margin; no baseline reflected guard because it may overlap desired echoes. Reflected screens tested as diagnostics only.',
      'no_adaptive_selection':'Every declared path and sensitivity result is retained, including unsupported and off-flaw maxima. No metric selects replacement settings.',
      'fingerprints':{kind:dataset(kind)['input_fingerprint'] for kind in DATA}}
    protocol_path=HERE/'locked_protocol.json'
    if protocol_path.exists():
        old=json.loads(protocol_path.read_text());assert old['common_controls']==base and old['paths']==PATHS
        protocol=old
    else:write(protocol_path,protocol)
    factors=[('coarse',{}),('no_coherence',{'coherence':False}),('screen_off',{'pulse_screening':'off'}),
      ('screen_ssss',{'pulse_reflections':[False,False,False,True]}),('screen_reflected',{'pulse_reflections':[True]*4}),
      ('gate10',{'angular_gate':.1}),('gate30',{'angular_gate':.3}),('offset1',{'min_offset_mm':1}),('offset2',{'min_offset_mm':2}),
      ('receivers4',{'minimum_receivers':4,'minimum_effective_receivers':4}),('receivers16',{'minimum_receivers':16,'minimum_effective_receivers':16}),
      ('one_source',{'minimum_sources':1}),('three_sources',{'minimum_sources':3}),('band15',{'low_MHz':1}),
      ('tail20',{'screen_tail_db':20}),('tail40',{'screen_tail_db':40}),('shoulders',{'sources':[0,4]}),('floor_sources',{'sources':[1,2,3]})]
    factors.append(('exclude_repairs',{'clean':True}))
    for kind in ['partial','full']:
        jobs=[];index=[];d=dataset(kind)
        def add(path,name,changes=None,spacing=None,cases=None,inspect=None):
            c=copy.deepcopy(base);c.update(path=path,receivers=list(range(len(d['audit']['receivers']))));c.update(changes or {})
            if spacing:c['spacing_mm']=spacing
            if cases:c['cases']=cases
            j={'selection':c}
            if inspect:j['inspect']=inspect
            index.append({'index':len(jobs),'kind':kind,'path':path,'name':name,'selection':c});jobs.append(j)
        for path in PATHS:
            inspection={'source':2,'receiver':165 if kind=='partial' else 168,'signal':'lof','competitor':'2-3_SSSSS','has_point':False} if path==PATHS[0] else None
            add(path,'baseline',cases=['lof','reference','difference'],inspect=inspection)
        for path in PATHS:
            for name,changes in factors:
                if kind=='full' and name not in ['coarse','no_coherence','screen_reflected','exclude_repairs']:continue
                if kind=='partial' and name=='exclude_repairs':continue
                add(path,name,changes,spacing=.2)
            if kind=='partial' or path in CORE:add(path,'fine',spacing=.05)
            if path in WIDE:add(path,'overview',{'x_min_mm':-45,'x_max_mm':45,'depth_min_mm':.25,'depth_max_mm':25.25},spacing=.25)
        folder=OUT/kind;folder.mkdir(exist_ok=True)
        write(folder/'requests.json',{'verify_inputs':True,'requests':jobs});write(folder/'index.json',index)
        print(kind,len(jobs),'configurations;',sum(len(j['selection']['cases']) for j in jobs),'signal-case images',flush=True)
    write(OUT/'protocol_freeze.json',{'protocol_sha256':sha(protocol_path),'protocol':protocol,
          'raw_data_manifest_sha256':{k:sha(DATA[k]/'dataset.json') for k in DATA}})

def run(kind):
    exe=_executable('SAFTBatch')
    with (OUT/kind/'compute.log').open('w') as f:
        p=subprocess.run([str(exe),'--dataset',str(DATA[kind]),'--compute',str(OUT/kind/'requests.json'),
                          '--output',str(OUT/kind/'numerical')],cwd=ROOT,env=_native_env(),stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError((OUT/kind/'compute.log').read_text())
    print(kind,'completed',flush=True)

def index(kind):return json.loads((OUT/kind/'index.json').read_text())
def row(kind,path,name='baseline'):return next(r for r in index(kind) if r['path']==path and r['name']==name)
def load(kind,path,name='baseline'):return np.load(OUT/kind/'numerical'/str(row(kind,path,name)['index'])/'images.npz')
def tips(kind):
    d=dataset(kind);t=np.array(d['audit']['defect_benchmark_only']['tips_xy_mm']);t[:,1]=d['material']['thickness_m']*1000-t[:,1]
    return t[np.argsort(t[:,1])]
def distances(points,kind):
    a,b=tips(kind);v=b-a;p=points*1000;u=np.clip((p-a)@v/(v@v),0,1)
    return np.linalg.norm(p-a-u[:,None]*v,axis=1)
def metrics(values,points,kind):
    d=distances(points,kind);valid=np.isfinite(values);near=valid&(d<=.5);bg=valid&(d>1)
    def maximum(mask):
        if not mask.any():return None
        i=int(np.argmax(np.where(mask,values,-np.inf)));return {'um_s':float(values[i]*1e6),'x_mm':float(points[i,0]*1000),'depth_mm':float(points[i,1]*1000),'distance_mm':float(d[i])}
    pk=maximum(near);back=maximum(bg);b99=float(np.percentile(values[bg],99)*1e6) if bg.any() else None
    def db(a,b):return float(20*np.log10(a/b)) if a is not None and b is not None and a>0 and b>0 else None
    return {'supported':int(valid.sum()),'support_pct':float(100*valid.mean()),'near_support_pct':float(100*valid[d<=.5].mean()),
      'peak':maximum(valid),'near':pk,'background':back,'background_p99_um_s':b99,
      'cmax_db':db(pk['um_s'] if pk else None,back['um_s'] if back else None),'c99_db':db(pk['um_s'] if pk else None,b99)}

if __name__=='__main__':
    if sys.argv[1]=='prepare':prepare()
    elif sys.argv[1]=='run':run(sys.argv[2])
