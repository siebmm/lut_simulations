"""Full-weld study with one common reconstruction protocol for four paths."""

# Locate portable scientific inputs without depending on the source checkout.
import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import data_root as _data_root, native_root as _native_root, executable as _executable, native_env as _native_env
from pathlib import Path
import copy,hashlib,json,struct,subprocess,sys
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
OUT=ROOT/'studies/laser/results/full_weld_common_settings_20261030'
DATA=_data_root()/'full_weld'
PATHS=['1-1_SS','1-1_LL','1-2_SSS','2-1_SSS']
LABELS=['S | S','L | L','S | SS','SS | S']

def prepare():
    OUT.mkdir(exist_ok=True)
    protocol=HERE/'common_protocol.json'
    base=json.loads(protocol.read_text())['common_selection']
    jobs=[];idx=[]
    def add(k,name,changes=None,spacing=None,cases=None,inspect=False):
        c=copy.deepcopy(base);c['path']=PATHS[k];c.update(changes or {})
        if spacing:c['spacing_mm']=spacing
        if cases:c['cases']=cases
        job={'selection':c}
        if inspect:job['inspect']={'source':0,'receiver':94,'signal':'lof','competitor':'1-2_SSS','has_point':True,'x_mm':2.766365013,'depth_mm':19.71747506}
        idx.append({'index':len(jobs),'path':PATHS[k],'name':name,'selection':c});jobs.append(job)
    for k in range(4):add(k,'baseline',cases=['lof','reference','difference'],inspect=k==0)
    factors=[('coarse',{}),('no_coherence',{'coherence':False}),('exclude_repairs',{'clean':True}),
      ('gate10',{'angular_gate':.1}),('gate30',{'angular_gate':.3}),('offset0',{'min_offset_mm':0}),('offset4',{'min_offset_mm':4}),
      ('wall0',{'wall_gate':0}),('wall10',{'wall_gate':.1}),('band03',{'high_MHz':3}),('band15',{'low_MHz':1}),
      ('sources3',{'minimum_sources':3}),('receivers16',{'minimum_receivers':16,'minimum_effective_receivers':16}),
      ('fraction50',{'receiver_fraction':.5}),('tail20',{'screen_tail_db':20}),('tail40',{'screen_tail_db':40}),
      ('margin0',{'screen_margin_us':0}),('margin02',{'screen_margin_us':.2}),('weighted',{'processing':'weighted'}),
      ('shoulders',{'sources':[0,4]})]
    for k in range(4):
        for name,c in factors:add(k,name,c,spacing=.2)
        for screen in ['off','L','S','R','LS','LR','SR']:
            add(k,'no_screen' if screen=='off' else 'screen_'+screen,{'pulse_screening':screen},spacing=.2)
        add(k,'fine',spacing=.05)
        add(k,'clean_fine',{'clean':True})
        add(k,'overview',{'x_min_mm':-45,'x_max_mm':45,'depth_min_mm':.25,'depth_max_mm':25.25},spacing=.25)
    # Baseline recipes may differ only by path. No per-path control exceptions.
    normalized=[{k:v for k,v in j['selection'].items() if k!='path'} for j in jobs[:4]]
    assert all(v==normalized[0] for v in normalized)
    (OUT/'requests.json').write_text(json.dumps({'verify_inputs':True,'requests':jobs},indent=2))
    (OUT/'index.json').write_text(json.dumps(idx,indent=2))
    (HERE/'common_protocol.json').write_text(json.dumps({'common_selection':{k:v for k,v in base.items() if k!='path'},'paths':PATHS},indent=2))
    print('Prepared',len(jobs),'reconstructions',flush=True)

def run():
    exe=_executable('SAFTBatch')
    with (OUT/'compute.log').open('w') as f:
        p=subprocess.run([str(exe),'--dataset',str(DATA),'--compute',str(OUT/'requests.json'),'--output',str(OUT/'numerical')],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError((OUT/'compute.log').read_text())

def load(row):return np.load(OUT/'numerical'/str(row['index'])/'images.npz')
def distances(x,z):
    p=np.column_stack([x,z]);a=np.array([2.83187473735,18.46747505998]);b=np.array([2.70085528914,20.96747505998]);v=b-a
    u=np.clip((p-a)@v/(v@v),0,1);return np.linalg.norm(p-a-u[:,None]*v,axis=1)

def metrics(a,points):
    x,z=points.T*1000;d=distances(x,z);valid=np.isfinite(a);near=valid&(d<=.5);bg=valid&(d>1.)
    im=int(np.nanargmax(a)) if valid.any() else 0
    ip=int(np.nanargmax(np.where(near,a,np.nan))) if near.any() else im
    pk=float(np.max(a[near])) if near.any() else 0;bp=float(np.max(a[bg])) if bg.any() else 0;b99=float(np.percentile(a[bg],99)) if bg.any() else 0
    def db(v,w):return float(20*np.log10(v/w)) if v>0 and w>0 else None
    return {'supported':int(valid.sum()),'support_pct':float(100*valid.mean()),'near_support_pct':float(100*valid[d<=.5].mean()),
      'global_peak_um_s':float(a[im]*1e6),'global_peak_x_mm':float(x[im]),'global_peak_depth_mm':float(z[im]),'global_peak_distance_mm':float(d[im]),
      'near_peak_um_s':pk*1e6,'near_peak_x_mm':float(x[ip]),'near_peak_depth_mm':float(z[ip]),'background_peak_um_s':bp*1e6,'background_p99_um_s':b99*1e6,
      'contrast_peak_db':db(pk,bp),'contrast_p99_db':db(pk,b99)}

def summarise():
    idx=json.loads((OUT/'index.json').read_text());results=[];base={}
    for row in idx:
        a=load(row);r={'index':row['index'],'path':row['path'],'name':row['name'],'signals':{}}
        for case in row['selection']['cases']:r['signals'][case]=metrics(a[case+'_maps'][2],a['points_m'])
        if row['name']=='baseline':
            common=np.logical_and.reduce([np.isfinite(a[c+'_maps'][2]) for c in row['selection']['cases']]);r['common_supported']=int(common.sum())
            d=distances(*(a['points_m'].T*1000));near=d<=.5
            r['sources_near_hist']={str(int(k)):int(v) for k,v in zip(*np.unique(a['lof_counts'][2,near],return_counts=True))}
            r['receiver_use_per_source']=(a['lof_receiver_usage'][2]>0).sum(axis=1).tolist()
        if row['name']=='coarse':base[row['path']]=a['lof_maps'][2].copy()
        results.append(r)
    for r in results:
        a=load(idx[r['index']])['lof_maps'][2];b=base[r['path']]
        if a.shape==b.shape:
            common=np.isfinite(a)&np.isfinite(b)
            r['baseline_support_retained_pct']=float(100*common.sum()/np.isfinite(b).sum())
            r['common_relative_l2_pct']=float(100*np.linalg.norm(a[common]-b[common])/np.linalg.norm(b[common])) if common.any() else None
    # Verify source/receiver pulse windows are identical across all four routes.
    windows=[json.loads((OUT/'numerical'/str(k)/'pulse_windows.json').read_text()) for k in range(4)]
    shared_windows=all(w==windows[0] for w in windows);assert shared_windows
    checks=[]
    for k,path in enumerate(PATHS):
        a=load(idx[k]);fine=load(next(r for r in idx if r['path']==path and r['name']=='fine'))
        f=fine['lof_maps'][2].reshape(len(fine['depth_mm']),len(fine['x_mm']))[::2,::2].ravel()
        assert np.array_equal(np.isfinite(a['lof_maps'][2]),np.isfinite(f))
        err=float(np.nanmax(abs(a['lof_maps'][2]-f)));assert err<1e-17
        # Independent NumPy reassembly of the source-RMS amplitude from exported means.
        valid=a['lof_valid'][2].astype(bool);means=a['lof_means'][2];count=valid.sum(axis=0)
        rms=np.sqrt(np.sum(np.where(valid,abs(means)**2,0),axis=0)/np.maximum(count,1));rms[count<2]=np.nan
        rms_err=float(np.nanmax(abs(rms-a['lof_maps'][2])));assert rms_err<1e-17
        checks.append({'path':path,'shared_grid_error_m_s':err,'independent_rms_error_m_s':rms_err})
    common_mask=np.logical_and.reduce([np.isfinite(load(idx[k])['lof_maps'][2]) for k in range(4)])
    common_metrics=[metrics(np.where(common_mask,load(idx[k])['lof_maps'][2],np.nan),load(idx[k])['points_m']) for k in range(4)]
    for path in PATHS:
        ls=load(next(r for r in idx if r['path']==path and r['name']=='screen_LS'))['lof_maps'][2]
        lsr=load(next(r for r in idx if r['path']==path and r['name']=='coarse'))['lof_maps'][2]
        assert np.array_equal(ls,lsr,equal_nan=True)
        if path in PATHS[2:]:
            for name in ['no_screen','screen_L','screen_S','screen_R','screen_LS','screen_LR','screen_SR']:
                values=load(next(r for r in idx if r['path']==path and r['name']==name))['lof_maps'][2]
                assert np.array_equal(values,lsr,equal_nan=True)
    data={'common_protocol_verified':True,'identical_pulse_windows_verified':shared_windows,'numerical_checks':checks,
      'shared_all_path_support_metrics':common_metrics,'LS_LSR_primary_coarse_arrays_identical':True,
      'metrics_definition':{'near_mm':.5,'background_guard_mm':1.,'background_statistic':'peak and 99th percentile on supported grid'},'runs':results}
    (OUT/'analysis.json').write_text(json.dumps(data,indent=2))
    print(json.dumps({'checks':checks,'baseline':[r for r in results if r['name']=='baseline'],'shared_support':common_metrics},indent=2))

if __name__=='__main__':
    mode=sys.argv[1] if len(sys.argv)>1 else 'all'
    if mode in ['prepare','all']:prepare()
    if mode in ['run','all']:run()
    if mode in ['summarise','all']:summarise()
