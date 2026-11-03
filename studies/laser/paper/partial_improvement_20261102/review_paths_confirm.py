"""Finer/wider checks of high-ranking atlas views, including pristine lookalikes."""
import csv,json
import numpy as np
import review as r
D=r.DEST/'path_confirmation'
CHOICES=[
 ('3-1_SLSL','left','source1'),('3-1_SLSS','left','source1'),
 ('3-2_LSSLL','left_weld','source1'),('3-1_LSSL','left_weld','source1'),
 ('3-1_SSSL','left','weld_min2'),('2-2_LSSL','sides','source4'),
 ('3-1_LSSS','weld','source1'),('3-1_SLLS','weld','source1'),
 ('3-1_SLSS','left_weld','source1'),('1-2_SSS','left','source3'),
 ('2-2_SLSS','sides','source4'),('2-1_SSS','left_weld','source1'),
 ('2-1_SLS','left_weld','source4'),
 ('exp_corner_5_L_rx_pixel-id_2','right','source4'),
 ('exp_corner_4_S_tx_id_1','left_weld','source0')]
jobs=[]
for candidate,(path,group,view) in enumerate(CHOICES):
 ids,k=r.VIEWS[view]
 c=dict(r.base(),path=path,groups=[g in r.GROUPS[group] for g in range(3)],sources=ids,minimum_sources=k)
 for label,fields in [('fine',dict(spacing_mm=.05)),('wide',dict(spacing_mm=.25,x_min_mm=-45,x_max_mm=45,depth_min_mm=.25,depth_max_mm=25.25))]:
  jobs.append(dict(id=f'a{len(jobs):03d}',candidate=candidate,label=label,path=path,group=group,view=view,selection=dict(c,**fields)))
p=dict(selection_reason='High local contrast, broad support, defect/pristine separation, or high-ranking pristine lookalikes; no changes to atlas settings except grid/field. Include two corner hypotheses with uncalibrated amplitudes/phases as diagnostics.',requests=jobs)
D.mkdir(exist_ok=True);r.s.write(D/'protocol.json',p)
for job in jobs:
 if not (D/job['id']/'arrays.npz').exists():r.execute([job],D,False)
 print(job['id'],job['path'],job['label'],'complete',flush=True)
rows=[]
for job in jobs:
 with np.load(D/job['id']/'arrays.npz') as a:
  points=a['points'].reshape(-1,2);d=r.s.distances(points,'partial')
  row={k:v for k,v in job.items() if k!='selection'}
  for case in ('lof','reference','difference'):row[case]=r.metric(a[case+'_2_amplitude'],points,d)
  common=(d<=.5)&np.isfinite(a['lof_2_amplitude'])&np.isfinite(a['reference_2_amplitude'])
  l,b=a['lof_2_amplitude'][common],a['reference_2_amplitude'][common]
  row['lof_over_pristine_near_db']=float(20*np.log10(l.max()/b.max())) if common.any() and l.max()>0 and b.max()>0 else None
  rows.append(row)
r.s.write(D/'metrics.json',dict(rows=rows))
for x in rows:print(x['path'],x['group'],x['view'],x['label'],x['lof']['cmax_db'],x['lof']['near_support_pct'],x['lof']['peak_distance_mm'],x['lof_over_pristine_near_db'],flush=True)
