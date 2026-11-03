"""Reproduce the four supplied full-weld recipes and controlled sensitivity runs.

Run from any directory with the study's Python environment. No COMSOL solves.
The known flaw is used only in post-reconstruction benchmark metrics/overlays.
"""

# Locate portable scientific inputs without depending on the source checkout.
import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import data_root as _data_root, native_root as _native_root, executable as _executable, native_env as _native_env
from pathlib import Path
import copy, hashlib, json, struct, subprocess, sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OUT = ROOT/'studies/laser/results/full_weld_four_paths_20261030'
DATA = _data_root()/'full_weld'
PATHS = ['1-1_SS','1-1_LL','1-2_SSS','2-1_SSS']
LABELS = ['S | S','L | L','S | SS','SS | S']

def read_blend():
    raw=(ROOT/'studies/laser/config/saft/full_weld/saft_paths_full_weld_all.saftblend').read_bytes()
    magic,n=struct.unpack_from('>II',raw); assert magic==0x53414642
    meta=json.loads(raw[8:8+n]); p=8+n; arrays=[]
    for _ in range(2+len(meta['layers'])):
        n,=struct.unpack_from('>I',raw,p); p+=4
        arrays.append(np.frombuffer(raw,dtype='>f8',count=n,offset=p).copy()); p+=8*n
    assert p==len(raw)
    return meta,arrays

def prepare():
    OUT.mkdir(parents=True,exist_ok=True)
    meta,arrays=read_blend(); jobs=[]; index=[]
    def add(k,name,changes=None,spacing=None,cases=None,inspect=False):
        c=copy.deepcopy(meta['layers'][k]['selection'])
        c.update(changes or {})
        if spacing: c['spacing_mm']=spacing
        if cases: c['cases']=cases
        j={'selection':c}
        if inspect:
            j['inspect']={'source':0,'receiver':94,'signal':'lof','competitor':'1-2_SSS',
                          'has_point':True,'x_mm':2.766365013,'depth_mm':19.71747506}
        index.append({'index':len(jobs),'path':PATHS[k],'name':name,'selection':c})
        jobs.append(j)
    for k in range(4):
        add(k,'baseline',cases=['lof','reference','difference'],inspect=k==0)
    standalone=json.loads((ROOT/'studies/laser/config/saft/full_weld/saft_settings_full_weld_all.json').read_text())['selection']
    add(3,'standalone',standalone,cases=['lof','reference','difference'])
    # Changes are one factor at a time relative to each coarse-grid baseline.
    variants=[('coarse',{}),('no_coherence',{'coherence':False}),
              ('no_screen',{'pulse_screening':'off'}),
              ('screen_LSR',{'pulse_screening':'LSR'}),
              ('exclude_repairs',{'clean':True}),
              ('gate10',{'angular_gate':.1}),('gate30',{'angular_gate':.3}),
              ('offset0',{'min_offset_mm':0}),('offset4',{'min_offset_mm':4}),
              ('wall0',{'wall_gate':0}),('wall10',{'wall_gate':.1}),
              ('band03',{'high_MHz':3}),('band15',{'low_MHz':1}),
              ('sources3',{'minimum_sources':3}),
              ('receivers16',{'minimum_receivers':16,'minimum_effective_receivers':16}),
              ('fraction50',{'receiver_fraction':.5}),
              ('tail20',{'screen_tail_db':20}),('tail40',{'screen_tail_db':40}),
              ('margin0',{'screen_margin_us':0}),('margin02',{'screen_margin_us':.2}),
              ('weighted',{'processing':'weighted'}),
              ('shoulders',{'sources':[0,4]})]
    for k in range(4):
        for name,changes in variants:
            if k>=2 and name in ['no_screen','tail20','tail40','margin0','margin02']: continue
            if k<2 and name in ['wall0','wall10']: continue
            if k==1 and name=='screen_LSR': continue
            add(k,name,changes,spacing=.2)
    for k in range(4):
        add(k,'fine',spacing=.05)
        add(k,'clean_fine',{'clean':True})
        add(k,'overview',{'x_min_mm':-45,'x_max_mm':45,'depth_min_mm':.25,'depth_max_mm':25.25},spacing=.25)
    (OUT/'requests.json').write_text(json.dumps({'verify_inputs':True,'requests':jobs},indent=2))
    (OUT/'index.json').write_text(json.dumps(index,indent=2))
    (HERE/'original_recipes.json').write_text(json.dumps({'blend_settings':meta['settings'],
        'layers':[{k:v for k,v in l.items() if k!='provenance'} for l in meta['layers']],
        'standalone_selection':standalone},indent=2))
    print(f'Prepared {len(jobs)} reconstructions',flush=True)

def run():
    exe=_executable('SAFTBatch')
    with (OUT/'compute.log').open('w') as log:
        result=subprocess.run([str(exe),'--dataset',str(DATA),'--compute',str(OUT/'requests.json'),
                               '--output',str(OUT/'numerical')],cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode: raise RuntimeError((OUT/'compute.log').read_text())

def load(row):
    return np.load(OUT/'numerical'/str(row['index'])/'images.npz')

def distances(x,z):
    p=np.column_stack([x,z]); a=np.array([2.83187473735,18.46747505998]); b=np.array([2.70085528914,20.96747505998])
    v=b-a; u=np.clip((p-a)@v/(v@v),0,1)
    return np.linalg.norm(p-a-u[:,None]*v,axis=1)

def metrics(a,points):
    x,z=points.T*1000; d=distances(x,z); valid=np.isfinite(a)
    near=valid&(d<=.5); bg=valid&(d>1.)
    im=int(np.nanargmax(a)); ip=int(np.nanargmax(np.where(near,a,np.nan))) if near.any() else im
    pk=float(np.max(a[near])) if near.any() else 0
    bgpk=float(np.max(a[bg])) if bg.any() else 0
    bg99=float(np.percentile(a[bg],99)) if bg.any() else 0
    def db(v,w):return float(20*np.log10(v/w)) if v>0 and w>0 else None
    return {'supported':int(valid.sum()),'support_pct':float(100*valid.mean()),'near_support_pct':float(100*valid[d<=.5].mean()),
            'global_peak_um_s':float(a[im]*1e6),'global_peak_x_mm':float(x[im]),'global_peak_depth_mm':float(z[im]),
            'global_peak_distance_mm':float(d[im]),'near_peak_um_s':pk*1e6,'near_peak_x_mm':float(x[ip]),'near_peak_depth_mm':float(z[ip]),
            'background_peak_um_s':bgpk*1e6,'background_p99_um_s':bg99*1e6,'contrast_peak_db':db(pk,bgpk),'contrast_p99_db':db(pk,bg99)}

def summarise():
    index=json.loads((OUT/'index.json').read_text()); results=[]; meta,arrays=read_blend()
    base={}; exact=[]
    for row in index:
        a=load(row); m=dict(index=row['index'],path=row['path'],name=row['name'],signals={})
        for case in row['selection']['cases']:
            m['signals'][case]=metrics(a[case+'_maps'][2],a['points_m'])
        if row['name']=='baseline':
            k=PATHS.index(row['path']); image=a['lof_maps'][2]
            exact.append({'path':row['path'],'equal_mask':bool(np.array_equal(np.isfinite(image),np.isfinite(arrays[2+k]))),
                          'max_abs_m_s':float(np.nanmax(np.abs(image-arrays[2+k])))})
            # Cross-case images use different nonlinear CF/record-specific masks: do not subtract them.
            common=np.isfinite(a['lof_maps'][2])&np.isfinite(a['reference_maps'][2])&np.isfinite(a['difference_maps'][2])
            m['common_supported']=int(common.sum())
            for case in row['selection']['cases']:
                m[case+'_common']=metrics(np.where(common,a[case+'_maps'][2],np.nan),a['points_m'])
            d=distances(*(a['points_m'].T*1000)); region=(d<=.5)&common
            ratio=np.linalg.norm(a['difference_maps'][2][region])/np.linalg.norm(a['lof_maps'][2][region])
            m['difference_to_lof_near_l2']=float(ratio)
            m['sources_near_hist']={str(int(v)):int(n) for v,n in zip(*np.unique(a['lof_counts'][2][d<=.5],return_counts=True))}
            m['receiver_use_per_source']=(a['lof_receiver_usage'][2]>0).sum(axis=1).tolist()
        if row['name']=='coarse':base[row['path']]=a['lof_maps'][2].copy()
        results.append(m)
    for r in results:
        row=index[r['index']]; a=load(row)['lof_maps'][2]; b=base[r['path']]
        if a.shape==b.shape:
            common=np.isfinite(a)&np.isfinite(b); both0=np.nan_to_num(a); base0=np.nan_to_num(b)
            r['baseline_support_retained_pct']=float(100*common.sum()/np.isfinite(b).sum())
            r['common_relative_l2_pct']=float(100*np.linalg.norm(a[common]-b[common])/np.linalg.norm(b[common])) if common.any() else None
    payload={'blend_reproduction':exact,'metrics_definition':{'near_mm':.5,'background_guard_mm':1.,'background_statistic':'peak and 99th percentile on supported saved grid'},'runs':results}
    (OUT/'analysis.json').write_text(json.dumps(payload,indent=2))
    print(json.dumps({'reproduction':exact,'baseline':[r for r in results if r['name'] in ('baseline','standalone')]},indent=2))

if __name__=='__main__':
    mode=sys.argv[1] if len(sys.argv)>1 else 'all'
    if mode in ('prepare','all'):prepare()
    if mode in ('run','all'):run()
    if mode in ('summarise','all'):summarise()
