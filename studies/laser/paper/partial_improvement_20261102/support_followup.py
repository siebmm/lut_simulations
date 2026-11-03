"""Locate missing support and test one-source admission for the selected recipe."""
import copy
import json
import numpy as np
from run import HERE, OUT, s

DEST=OUT/'support_followup'
session=json.loads((OUT/'settings/saft_partial_weighted_1_5.json').read_text())
strict=copy.deepcopy(session['selection'])
strict['cases']=['lof','reference','difference']
relaxed=dict(strict,minimum_sources=1)
tips=np.array(s.dataset('partial')['audit']['defect_benchmark_only']['tips_xy_mm'])
tips[:,1]=25.4-tips[:,1]
line=tips[0]+np.linspace(0,1,101)[:,None]*(tips[1]-tips[0])
jobs=[dict(id='minimum_one_source',selection=relaxed),
      dict(id='exact_segment_support',selection=dict(strict,cases=['lof']),points_m=(line/1000).tolist())]
protocol=dict(question='Can the selected weighted 1-5 MHz recipe support the whole near-defect region, and is the actual defect line already supported?',
    controls='Only minimum_sources changes in the paired images. Exact segment points are used for an admission diagnostic only, not for tuning or modifying the full image.',
    requests=jobs)
p=DEST/'protocol.json'
if p.exists():assert json.loads(p.read_text())==protocol
else:s.write(p,protocol)
pending=[j for j in jobs if not (DEST/j['id']/'arrays.npz').exists()]
if pending:
    s.execute('partial',dict(requests=pending),DEST)
    for job in pending:s.compact(DEST/job['id'])
rows=[]
with np.load(OUT/'frequency/all__weighted__min2__band1_5/arrays.npz') as original,np.load(DEST/'minimum_one_source/arrays.npz') as changed:
    pts=original['points'].reshape(-1,2)
    np.testing.assert_array_equal(pts,changed['points'].reshape(-1,2))
    near=s.distances(pts,'partial')<=.5
    for case in ('lof','reference','difference'):
        key=case+'_2_amplitude';a,b=original[key],changed[key]
        keep=np.isfinite(a);added=np.isfinite(b)&~keep
        np.testing.assert_array_equal(a[keep],b[keep])
        assert np.all(changed[case+'_2_count'][added]==1)
        rows.append(dict(case=case,before=s.metrics(a,pts,'partial'),after=s.metrics(b,pts,'partial'),
            shared_amplitudes_exactly_unchanged=True,added_image_pixels=int(added.sum()),
            near_pixels=int(near.sum()),added_near_pixels=int((near&added).sum())))
    missing=near&~np.isfinite(original['lof_2_amplitude'])
    missing_data=dict(near_pixels=int(near.sum()),missing_pixels=int(missing.sum()),
        coordinates_mm=(pts[missing]*1000).tolist(),
        source_counts=original['lof_2_count'][missing].tolist(),
        qualifying_sources=[np.flatnonzero(v).tolist() for v in original['lof_2_valid'].reshape(5,-1)[:,missing].T])
with np.load(DEST/'exact_segment_support/arrays.npz') as a:
    line_data=dict(points=101,supported_points=int(np.isfinite(a['lof_2_amplitude']).sum()),
        minimum_qualifying_sources=int(a['lof_2_count'].min()),
        coordinates_mm=line.tolist(),source_counts=a['lof_2_count'].tolist())
session['selection']['minimum_sources']=1
s.write(DEST/'saft_partial_weighted_1_5_minimum_one_source.json',session)
s.write(DEST/'summary.json',dict(missing_neighbourhood=missing_data,exact_segment=line_data,comparisons=rows,
    interpretation='Near-defect support is a pixel-admission fraction in a 0.5 mm neighbourhood, not defect detection probability or percentage of defect length recovered. One-source admission retains existing supported values and adds less redundant measurements.'))
print(json.dumps(dict(missing_pixels=missing_data['missing_pixels'],near_pixels=missing_data['near_pixels'],
    exact_segment_supported=line_data['supported_points'],exact_segment_points=line_data['points'],
    minimum_segment_sources=line_data['minimum_qualifying_sources'],
    lof=rows[0]),indent=2),flush=True)
