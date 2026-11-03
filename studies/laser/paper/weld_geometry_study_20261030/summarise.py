"""Evaluate frozen reconstructions; ground truth enters here, after design."""

import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import executable as _executable, native_env as _native_env
from study import *
import csv,collections

def summarise(kind):
    idx=index(kind);runs=[];checks=[]
    for r in idx:
        a=load(kind,r['path'],r['name']);scores={c:metrics(a[c+'_maps'][2],a['points_m'],kind) for c in r['selection']['cases']}
        v=a['lof_valid'][2].astype(bool);m=a['lof_means'][2];count=v.sum(axis=0)
        rms=np.sqrt(np.sum(np.where(v,abs(m)**2,0),axis=0)/np.maximum(count,1));rms[count<r['selection']['minimum_sources']]=np.nan
        actual=a['lof_maps'][2];assert np.array_equal(np.isfinite(actual),np.isfinite(rms))
        err=float(np.nanmax(abs(actual-rms))) if np.isfinite(actual).any() else 0
        assert err<1e-17
        entry={**{k:r[k] for k in ['index','path','name']},'signals':scores,'rms_error_m_s':err}
        if r['name']=='baseline':
            near=distances(a['points_m'],kind)<=.5;valid=np.isfinite(actual)
            entry['near_source_count_histogram']={str(int(x)):int(n) for x,n in zip(*np.unique(a['lof_counts'][2,near],return_counts=True))}
            entry['source_qualified_fraction_pct']=[float(100*value[valid].mean()) if valid.any() else None for value in v]
            entry['unique_receivers_used_per_source']=(a['lof_receiver_usage'][2]>0).sum(axis=1).tolist()
            rejection=a['lof_pulse_rejection_pct'][2]
            entry['median_screened_pct']=float(np.nanmedian(rejection)) if np.isfinite(rejection).any() else None
            # Paired control on the same pixels; this includes differences in reference-derived windows.
            shared=np.logical_and.reduce([np.isfinite(a[c+'_maps'][2]) for c in ['lof','reference','difference']])
            entry['paired_case_common_mask']={c:metrics(np.where(shared,a[c+'_maps'][2],np.nan),a['points_m'],kind) for c in ['lof','reference','difference']}
        if r['name']=='fine':
            b=load(kind,r['path']);f=actual.reshape(len(a['depth_mm']),-1)[::2,::2].ravel();bval=b['lof_maps'][2]
            assert f.shape==bval.shape and np.array_equal(np.isfinite(f),np.isfinite(bval))
            error=float(np.nanmax(abs(f-bval))) if np.isfinite(f).any() else 0;assert error<1e-17
            checks.append({'path':r['path'],'shared_pixel_error_m_s':error})
        runs.append(entry)
    for r in runs:
        source=idx[r['index']]
        if source['selection']['spacing_mm']!=.2:continue
        a=load(kind,r['path'],r['name'])['lof_maps'][2];b=load(kind,r['path'],'coarse')['lof_maps'][2]
        mask=np.isfinite(a)&np.isfinite(b)
        r['coarse_common_l2_change_pct']=float(100*np.linalg.norm(a[mask]-b[mask])/np.linalg.norm(b[mask])) if mask.any() and np.linalg.norm(b[mask])>0 else None
    controls=[{k:v for k,v in r['selection'].items() if k!='path'} for r in idx[:18]]
    assert all(c==controls[0] for c in controls)
    windows=[json.loads((OUT/kind/'numerical'/str(i)/'pulse_windows.json').read_text()) for i in range(18)]
    assert all(w==windows[0] for w in windows)
    summary={case:{'count':len(w['windows']),'limited':sum(v['search_or_extent_limited'] for v in w['windows']),
                   'families':dict(collections.Counter(v['family'] for v in w['windows']))} for case,w in windows[0].items()}
    result={'kind':kind,'configurations':len(idx),'signal_case_images':sum(len(r['signals']) for r in runs),
      'identical_primary_controls_except_path':True,'identical_primary_windows':True,'pulse_windows':summary,
      'max_rms_error_m_s':max(r['rms_error_m_s'] for r in runs),'fine_grid_checks':checks,'runs':runs}
    write(OUT/kind/'analysis.json',result)
    print(kind,len(runs),'configs; RMS max',result['max_rms_error_m_s'],'fine grids',len(checks),flush=True)
    return result

def combine():
    results={kind:json.loads((OUT/kind/'analysis.json').read_text()) for kind in DATA}
    rows=[];matched=[]
    for kind,a in results.items():
        for r in a['runs']:
            for case,m in r['signals'].items():
                rows.append({'dataset':kind,'path':r['path'],'name':r['name'],'signal':case,
                  **{k:m[k] for k in ['support_pct','near_support_pct','cmax_db','c99_db']},
                  **{'peak_'+k:v for k,v in (m['peak'] or dict(um_s=None,x_mm=None,depth_mm=None,distance_mm=None)).items()},
                  'near_peak_um_s':m['near']['um_s'] if m['near'] else None,
                  'background_peak_um_s':m['background']['um_s'] if m['background'] else None,'background_p99_um_s':m['background_p99_um_s']})
    for path in PATHS:
        a=load('full',path);b=load('partial',path);assert np.array_equal(a['points_m'],b['points_m'])
        mask=np.isfinite(a['lof_maps'][2])&np.isfinite(b['lof_maps'][2])
        matched.append({'path':path,'full':metrics(np.where(mask,a['lof_maps'][2],np.nan),a['points_m'],'full'),
                        'partial':metrics(np.where(mask,b['lof_maps'][2],np.nan),b['points_m'],'partial')})
    with (HERE/'quantitative_results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    write(OUT/'matched_comparison.json',matched)
    exe=_executable('SAFTBatch')
    frozen=json.loads((OUT/'protocol_freeze.json').read_text());assert sha(HERE/'locked_protocol.json')==frozen['protocol_sha256']
    for kind in DATA:
        assert sha(DATA[kind]/'dataset.json')==frozen['raw_data_manifest_sha256'][kind]
        assert json.loads((OUT/kind/'numerical/input_verification.json').read_text())['status']=='passed'
    provenance={'protocol_sha256':sha(HERE/'locked_protocol.json'),'native_executable':str(exe.relative_to(ROOT)),
       'native_executable_sha256':sha(exe),'datasets':{k:{'manifest_sha256':sha(DATA[k]/'dataset.json'),'fingerprint':dataset(k)['input_fingerprint']} for k in DATA},
       'counts':{k:{q:results[k][q] for q in ['configurations','signal_case_images']} for k in DATA}}
    write(HERE/'reproducibility.json',provenance)
    print('Combined evidence:',len(rows),'case metrics',flush=True)

if __name__=='__main__':
    if sys.argv[1] in DATA:summarise(sys.argv[1])
    else:combine()
