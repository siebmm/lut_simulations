"""Controlled one-versus-two-source partial-weld audit; no filter/aperture tuning."""
from pathlib import Path
import copy,json,sys
import numpy as np
HERE=Path(__file__).resolve().parent
BASE=HERE.parent/'comsol_chapter_20261101'
ROOT=HERE.parents[3]
OUT=ROOT/'studies/laser/results/partial_source_requirement_20261101'
sys.path.insert(0,str(HERE.parent/'filter_study_20261101'))
import study as s

def main():
    existing=json.loads((BASE/'data/joint_acquisition/index.json').read_text())
    names={'all':'control__all','both':'fine__both_sides','left':'fine__left_only'}
    jobs=[]
    for name,oldid in names.items():
        cfg=copy.deepcopy(next(q['selection'] for q in existing if q['id']==oldid))
        assert cfg['minimum_sources']==2 and cfg['processing']=='raw'
        jobs.append(dict(id=name+'__one',selection=dict(cfg,minimum_sources=1),baseline=oldid))
        if name in ('all','left'):
            jobs.append(dict(id=name+'__replay_two',selection=cfg,baseline=oldid))
    for name,oldid in [('left_finer','left_only__finer'),('left_wide','left_only__wide'),
                       ('both_wide','both_sides__wide')]:
        cfg=copy.deepcopy(next(q['selection'] for q in existing if q['id']==oldid))
        jobs.append(dict(id=name+'__one',selection=dict(cfg,minimum_sources=1),baseline=oldid))
    protocol={'question':'What changes when one qualifying source is sufficient for the partial-weld S|SS image?',
      'frozen':'Change only minimum_sources from two to one. Keep source lists, physical receiver groups, filter, CF, LSR screening, receiver qualification and grid fixed within each pair.',
      'interpretation':'One-source support is a pixel-acceptance rule; it does not select a single laser or remove the other available source shots. Compare both the entire newly admitted field and the unchanged common support.',
      'cases':['lof','reference','difference'],'requests':jobs,
      'dataset_fingerprint':s.dataset('partial')['input_fingerprint'],'executable_sha256':s.sha(s.EXE)}
    OUT.mkdir(parents=True,exist_ok=True)
    p=OUT/'protocol.json'
    if p.exists():
        archived=json.loads(p.read_text())
        # A rebuild changes executable bytes, not the frozen scientific request.
        assert {k:v for k,v in archived.items() if k!='executable_sha256'} == {k:v for k,v in protocol.items() if k!='executable_sha256'}
        s.write(OUT/'replay_executable.json',{'sha256':s.sha(s.EXE)})
    else:s.write(p,protocol)
    pending=[dict(id=j['id'],selection=j['selection']) for j in jobs if not (OUT/j['id']/'arrays.npz').exists()]
    if pending:
        print('Computing',len(pending),'fixed source-threshold configurations.',flush=True)
        # A private request path avoids overwriting another study's request file.
        s.LOCAL=OUT/'native_request';s.LOCAL.mkdir(exist_ok=True)
        s.execute('partial',dict(requests=pending),OUT)
        for j in pending:s.compact(OUT/j['id'])
    rows=[];verification=[]
    for job in jobs:
        with np.load(OUT/job['id']/'arrays.npz') as a, np.load(BASE/f'data/joint_acquisition/{job["baseline"]}.npz') as b:
            points=a['points'].reshape(-1,2)
            xx,zz=np.meshgrid(b['x_mm'],b['depth_mm'])
            np.testing.assert_allclose(points*1000,np.c_[xx.ravel(),zz.ravel()],atol=1e-10)
            for case in protocol['cases']:
                v=a[case+'_2_amplitude'];old=b[case+'_amplitude'].ravel()
                count=a[case+'_2_count'];oldmask=np.isfinite(old);valid=np.isfinite(v)
                assert np.all(valid[oldmask])
                np.testing.assert_allclose(v[oldmask],old[oldmask],rtol=2e-9,atol=1e-15)
                np.testing.assert_array_equal(count[oldmask],b[case+'_count'].ravel()[oldmask])
                added=valid & ~oldmask
                if job['selection']['minimum_sources']==1:
                    assert np.all(count[added]==1)
                else:assert not added.any()
                row=dict(id=job['id'],case=case,baseline=job['baseline'],minimum_sources=job['selection']['minimum_sources'],
                    metrics=s.metrics(v,points,'partial'),baseline_metrics=s.metrics(old,points,'partial'),
                    common_metrics=s.metrics(v,points,'partial',oldmask),
                    newly_supported_metrics=s.metrics(v,points,'partial',added),
                    added_pixels=int(added.sum()),common_max_absolute_difference=float(np.max(np.abs(v[oldmask]-old[oldmask]))))
                rows.append(row)
                if case=='lof':
                    m=row['metrics'];print(job['id'],'C=',m['cmax_db'],'near=',m['near_support_pct'],'peak=',m['peak'],'added=',row['added_pixels'],flush=True)
                verification.append(dict(id=job['id'],case=case,common_support_unchanged=True,new_pixels_have_one_source=True))
    s.write(OUT/'metrics.json',{'rows':rows})
    s.write(OUT/'verification.json',{'status':'passed','signal_cases':len(rows),
        'new_scientific_configurations':6,'replay_controls':2,'checks':verification,
        'protocol_sha256':s.sha(p)})
    # Export a loadable one-source alternative; preserve the current two-source settings.
    session=json.loads((ROOT/'studies/laser/results/partial_acquisition_20261101/saft_partial_weld_sources_left_receivers.json').read_text())
    session['selection']['minimum_sources']=1
    s.write(OUT/'saft_partial_left_minimum_one_source.json',session)

if __name__=='__main__':main()
