"""Independent array-to-score and coordinate checks for every archived image."""
from pathlib import Path
import json
import numpy as np
HERE=Path(__file__).resolve().parent;DATA=HERE/'data'

def main():
    index=json.loads((DATA/'images.json').read_text());assert len(index)==384
    survey=json.loads((DATA/'survey_metrics.json').read_text())['rows']
    others={('confirmation',k):json.loads((DATA/(k+'_confirmation_metrics.json')).read_text())['rows'] for k in ['full','partial']}
    others['selected_checks','full']=json.loads((DATA/'full_selected_metrics.json').read_text())['rows']
    checked=0;worst=0.
    for row in index:
        assert row['processing']=='raw' and row['selection']['processing']=='raw'
        kind=row['kind'];d=json.loads((DATA/(kind+'_dataset.json')).read_text())
        tips=np.array(d['audit']['defect_benchmark_only']['tips_xy_mm']);tips[:,1]=25.4-tips[:,1]
        with np.load(DATA/(row['uid']+'.npz')) as a:
            xx,zz=np.meshgrid(a['x_mm'],a['depth_mm']);p=np.column_stack([xx.ravel(),zz.ravel()])
            assert np.all(np.diff(a['x_mm'])>0) and np.all(np.diff(a['depth_mm'])>0)
            assert p[:,0].min()>=row['selection']['x_min_mm']-1e-9
            assert p[:,0].max()<=row['selection']['x_max_mm']+1e-9
            assert p[:,1].min()>=row['selection']['depth_min_mm']-1e-9
            assert p[:,1].max()<=row['selection']['depth_max_mm']+1e-9
            edge=tips[1]-tips[0];u=np.clip((p-tips[0])@edge/np.dot(edge,edge),0,1)
            distance=np.sqrt(((p-tips[0]-u[:,None]*edge)**2).sum(axis=1))
            for case in ['lof','reference','difference']:
                if row['stage']=='numerical':
                    expected=next(r for r in survey if r['kind']==kind and r['variant']==row['variant'] and r['path']==row['path'] and r['case']==case)
                else:expected=next(r for r in others[row['stage'],kind] if r['id']==row['id'] and r['case']==case)
                values=a[case+'_amplitude'].ravel();valid=np.isfinite(values)
                assert np.all(values[valid]>=0) and np.all(a[case+'_count'].ravel()[valid]>=2)
                near=valid&(distance<=.5);bg=valid&(distance>1)
                support=valid[distance<=.5].mean()*100
                np.testing.assert_allclose(support,expected['near_support_pct'],rtol=0,atol=1e-9)
                if near.any() and bg.any() and max(values[near])>0 and max(values[bg])>0:
                    contrast=20*np.log10(max(values[near])/max(values[bg]));worst=max(worst,abs(contrast-expected['cmax_db']))
                    np.testing.assert_allclose(contrast,expected['cmax_db'],rtol=0,atol=1e-8)
                else:assert expected['cmax_db'] is None
                if valid.any():
                    j=np.nanargmax(values);np.testing.assert_allclose([*p[j],distance[j]],
                        [expected['peak']['x_mm'],expected['peak']['depth_mm'],expected['peak']['distance_mm']],rtol=0,atol=1e-9)
                checked+=1
    result=dict(status='passed',configurations=len(index),signal_cases_checked=checked,max_contrast_difference_db=worst,
                coordinate_units='mm; bounds and native point correspondence verified',support='all finite image pixels have at least two qualified sources')
    (HERE/'verification').mkdir(parents=True, exist_ok=True)
    (HERE/'verification/data_check.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(result)

if __name__=='__main__':main()
