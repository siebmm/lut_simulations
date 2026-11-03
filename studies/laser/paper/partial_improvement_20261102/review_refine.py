"""Final grid and fixed-source-subset checks for the backwall-SS mute candidate."""
import itertools,json
import numpy as np
import review as r
import review_confirm as c

D=r.DEST/'refinement'
p=c.protocol();candidate=next(j['selection'] for j in p['requests'] if j['candidate']=='broad_SS03' and j['label']=='fine')
jobs=[]
for name,change in [('muted',{}),('unmuted',dict(mutes=[False]*6))]:
    for grid,fields in [('025',dict(spacing_mm=.025)),('shifted',dict(spacing_mm=.05,x_min_mm=-7.475,x_max_mm=7.475,depth_min_mm=16.675,depth_max_mm=25.225))]:
        jobs.append(dict(id=name+'_'+grid,selection=dict(candidate,**change,**fields)))
protocol=dict(reason='Check convergence beyond 0.05 mm and sensitivity to a half-pixel grid shift, with a matched unmuted control. Keep the same physical acquisition, filter, gates and source threshold.',requests=jobs)
path=D/'protocol.json'
if path.exists():assert json.loads(path.read_text())==protocol
else:r.s.write(path,protocol)
for job in jobs:
    if not (D/job['id']/'arrays.npz').exists():r.execute([job],D,False)
    print('refinement',job['id'],'complete',flush=True)
rows=[]
for job in jobs:
    with np.load(D/job['id']/'arrays.npz') as a:
        points=a['points'].reshape(-1,2);d=r.s.distances(points,'partial')
        for case in ('lof','reference','difference'):
            rows.append(dict(id=job['id'],case=case,**r.metric(a[case+'_2_amplitude'],points,d)))
r.s.write(D/'metrics.json',dict(rows=rows))
selected=next(j for j in p['requests'] if j['candidate']=='broad_SS03' and j['label']=='fine')
with np.load(c.locate(selected['id'])) as f:a=dict(f)
points=a['points'].reshape(-1,2);d=r.s.distances(points,'partial');n=len(points)
subsets=[]
for size in range(1,6):
    for ids in itertools.combinations(range(5),size):
        for minimum in range(1,size+1):
            row=dict(sources=list(ids),minimum_sources=minimum)
            for case in ('lof','reference','difference'):
                valid=a[case+'_2_valid'].reshape(5,n)[list(ids)].astype(bool)
                means=a[case+'_2_means'].reshape(5,n)[list(ids)]
                count=valid.sum(axis=0)
                v=np.sqrt(np.sum(np.where(valid,abs(means)**2,0),axis=0)/np.maximum(count,1));v[count<minimum]=np.nan
                row[case]=r.metric(v,points,d)
            subsets.append(row)
r.s.write(D/'source_subsets.json',dict(scope='All 31 fixed nonempty source subsets and every feasible minimum source count: 80 views, three signal cases, native per-source means on the 0.05 mm grid.',rows=subsets))
for x in rows:
    if x['case']=='lof':print(x['id'],x['cmax_db'],x['near_support_pct'],x['peak_distance_mm'],flush=True)
print('Best broad-coverage source subsets:',flush=True)
for x in sorted([x for x in subsets if x['lof']['near_support_pct']>=90 and x['lof']['cmax_db'] is not None],key=lambda x:x['lof']['cmax_db'],reverse=True)[:8]:
    print(x['sources'],x['minimum_sources'],x['lof']['cmax_db'],x['lof']['near_support_pct'],flush=True)
