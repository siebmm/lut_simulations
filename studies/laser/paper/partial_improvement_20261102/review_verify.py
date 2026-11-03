"""Verify final native/offline and nested-grid results and export documented settings."""
import json
import numpy as np
import review as r
import review_confirm as c
import settings_review as q
D=r.DEST/'final_checks'
checks=[];metrics={}
for job in json.loads((D/'protocol.json').read_text())['requests']:
 with np.load(D/job['id']/'arrays.npz') as a:
  points=a['points'].reshape(-1,2);d=r.s.distances(points,'partial')
  metrics[job['id']]={case:r.metric(a[case+'_2_amplitude'],points,d) for case in ('lof','reference','difference')}
  if job['id']=='wide01':continue
  anchor=next(j for j in q.protocol()['requests'] if j['anchor']=='broad_SSS' and j['label']=='baseline')
  with np.load(q.locate(anchor['id'])) as b:
   n=len(points)
   for case in ('lof','reference','difference'):
    v,count,valid,means=r.reconstruct(b[case+'_2_stats'].reshape(5,n,7),b[case+'_2_sums'].reshape(5,n),list(range(5)),2,job['selection'])
    np.testing.assert_allclose(v,a[case+'_2_amplitude'],rtol=2e-8,atol=2e-15,equal_nan=True)
    np.testing.assert_array_equal(count,a[case+'_2_count'])
    checks.append(job['id']+'/'+case)
r.s.write(D/'metrics.json',metrics)
for muted in (True,False):
 prefix='muted' if muted else 'unmuted'
 finejob=next(j for j in c.protocol()['requests'] if j['candidate']=='broad_SS03' and j['label']=='fine')
 finepath=c.locate(finejob['id']) if muted else r.OUT/'admission_audit/angular_confirmation/fine/arrays.npz'
 with np.load(r.DEST/f'refinement/{prefix}_025/arrays.npz') as dense,np.load(finepath) as fine,np.load(r.DEST/f'refinement/{prefix}_shifted/arrays.npz') as shifted:
  for case in ('lof','reference','difference'):
   value=dense[case+'_2_amplitude'].reshape(len(dense['z']),len(dense['x']))
   for view,arr,slice_ in [('aligned',fine,(slice(None,None,2),slice(None,None,2))),('shifted',shifted,(slice(1,None,2),slice(1,None,2)))]:
    np.testing.assert_allclose(value[slice_].ravel(),arr[case+'_2_amplitude'],rtol=3e-7,atol=2e-12,equal_nan=True)
    checks.append(prefix+'/'+view+'/'+case)
template=json.loads((r.OUT/'settings/saft_partial_weighted_1_5.json').read_text())
settings={}
for name,job,spacing in [('saft_partial_weighted_SS_mute_025',finejob,.025),('saft_partial_weighted_SS_mute_05',finejob,.05),('saft_partial_right_S_LL_mute_05',next(j for j in c.protocol()['requests'] if j['candidate']=='right_SLL_LL03' and j['label']=='fine'),.05)]:
 result=dict(template,selection=dict(job['selection'],spacing_mm=spacing,cases=['lof']))
 path=r.DEST/'settings'/f'{name}.json';r.s.write(path,result);settings[path.name]=r.s.sha(path)
r.s.write(D/'verification.json',dict(status='passed',checks=checks,settings_sha256=settings,description='Three changed offline controls checked against independent native reconstruction; both aligned and half-pixel-shifted grids agree with their nested 0.025 mm points. Presets derive from the native-confirmed selection.'))
print('Passed',len(checks),'comparisons. Wide 0.1 mm:',metrics['wide01']['lof'])
