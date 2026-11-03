"""Independent portable checks of reported metrics, source RMS and image sampling."""
from pathlib import Path
import json,hashlib
import numpy as np
HERE=Path(__file__).resolve().parent;DATA=HERE/'data'
def close(a,b):
    if a is None or b is None:assert a is None and b is None,(a,b)
    else:assert np.isclose(a,b,rtol=1e-11,atol=1e-11),(a,b)
def verify():
    rows=json.loads((DATA/'index.json').read_text());cases=0;grids=0;rms=0
    for r in rows:
        p=DATA/(r['uid']+'.npz');assert hashlib.sha256(p.read_bytes()).hexdigest()==r['sha256'];a=np.load(p)
        d=json.loads((DATA/(r['kind']+'_dataset.json')).read_text());tips=np.array(d['audit']['defect_benchmark_only']['tips_xy_mm']);tips[:,1]=d['material']['thickness_m']*1000-tips[:,1]
        points=a['points_m']*1000;u=tips[1]-tips[0];frac=np.minimum(1,np.maximum(0,np.einsum('ij,j->i',points-tips[0],u)/np.dot(u,u)));dist=np.linalg.norm(points-tips[0]-frac[:,None]*u,axis=1)
        for case,m in r['metrics'].items():
            v=a[case];ok=np.isfinite(v);near=ok&(dist<=.5);bg=ok&(dist>1);cases+=1
            assert int(ok.sum())==m['supported'];close(100*ok.mean(),m['support_pct']);close(100*ok[dist<=.5].mean(),m['near_support_pct'])
            for key,mask in [('peak',ok),('near',near),('background',bg)]:
                if not mask.any():assert m[key] is None;continue
                k=np.flatnonzero(mask)[np.argmax(v[mask])];expected=m[key]
                for field,val in [('um_s',v[k]*1e6),('x_mm',points[k,0]),('depth_mm',points[k,1]),('distance_mm',dist[k])]:close(val,expected[field])
            bp=np.quantile(v[bg],.99)*1e6 if bg.any() else None;close(bp,m['background_p99_um_s'])
            n=v[near].max() if near.any() else None
            for key,b in [('cmax_db',v[bg].max() if bg.any() else None),('c99_db',bp/1e6 if bp is not None else None)]:close(20*np.log10(n/b) if n is not None and b is not None and n>0 and b>0 else None,m[key])
        if 'lof_means' in a:
            good=a['lof_valid'].astype(bool);n=good.sum(axis=0);score=np.sqrt(np.sum(np.where(good,np.abs(a['lof_means'])**2,0),axis=0)/np.maximum(n,1));score[n<r['selection']['minimum_sources']]=np.nan
            assert np.allclose(score,a['lof'],equal_nan=True,rtol=1e-12,atol=1e-18);rms+=1
        if r['name']=='fine' or r['name'].endswith('_fine'):
            name='baseline' if r['name']=='fine' else r['name'].removesuffix('_fine')
            br=next(b for b in rows if b['kind']==r['kind'] and b['path']==r['path'] and b['name']==name);b=np.load(DATA/(br['uid']+'.npz'))
            assert np.array_equal(a['lof'].reshape(len(a['depth_mm']),-1)[::2,::2].ravel(),b['lof'],equal_nan=True);grids+=1
    new=[r for r in rows if r['stage'] in ['followup','confirmation']];report={'archived_configurations':len(rows),'independently_verified_case_metrics':cases,'independently_verified_source_RMS':rms,'exact_nested_grid_checks':grids,'new_configurations':len(new),'new_case_images':sum(len(r['metrics']) for r in new)}
    assert report['new_configurations']==82 and report['new_case_images']==142
    target=HERE/'verification/evidence_check.json';target.parent.mkdir(exist_ok=True);target.write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
if __name__=='__main__':verify()
