"""Audit receiver/source qualification using unchanged native per-source sums.

Thresholds affect qualification, not waveforms or per-source means. Export their
sufficient statistics once, replay a declared matrix, and check native controls.
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
import argparse
import copy
import csv
import itertools
import json
import os
import subprocess
import numpy as np
from run import HERE, OUT, ROOT, STUDY, s

DEST=OUT/'admission_audit'
EXE=_executable('SAFTAdmissionAudit')
FRACTIONS=[0.,.1,.25,.5,1.]
MINIMUMS=[1,4,8,16,32]
EFFECTIVE=[0,4,8,16]
SOURCES=[1,2,3]


def anchors():
    result=[]
    for label,name,archive in [
        ('current_equal','baseline','survey/all__raw__min2'),
        ('weighted_1_5','weighted_1_5','frequency/all__weighted__min2__band1_5')]:
        c=json.loads((OUT/f'settings/saft_partial_{name}.json').read_text())['selection']
        c['cases']=['lof','reference','difference']
        result.append(dict(id=label,selection=c,archive=archive))
    return result


def protocol():
    p=dict(question='Does the 25% eligible-receiver retention rule, the count of 8, the effective count of 8, or the minimum source count limit partial-weld imaging?',
        interpretation='These are engineering admission guards, not established physical optima. This is retrospective sensitivity analysis on one known simulated defect.',
        definition='A source qualifies when accepted >= max(minimum_receivers, ceil(receiver_fraction * eligible)) and effective >= minimum_effective_receivers. A pixel requires minimum_sources qualifying sources.',
        fixed='Existing geometry/visibility, angular/wall/offset gates, finite filtering support, phase, pulse screening, coherence exponent, source and receiver lists, and image field remain fixed within each anchor.',
        replay='Per-source means and statistics are independent of these four qualification thresholds in the production engine. Apply the exact native comparisons and source RMS offline; verify selected outcomes with fresh native reconstruction.',
        fractions=FRACTIONS,minimum_receivers=MINIMUMS,minimum_effective_receivers=EFFECTIVE,
        minimum_sources=SOURCES,combinations_per_anchor=300,anchors=anchors(),
        source_hashes={n:s.sha(_native_root()/n) for n in ('saft.cpp','saft.hpp','screening.cpp','dataset_geometry.cpp','groove_paths.cpp')})
    path=DEST/'protocol.json'
    from standalone import frozen_protocol
    return frozen_protocol(path,p,'source_hashes')


def dump():
    p=protocol();pending=[dict(id=j['id'],selection=j['selection']) for j in p['anchors']
                         if not (DEST/'statistics'/j['id']/'arrays.npz').exists()]
    if pending:
        request=DEST/'native_request.json';s.write(request,dict(requests=pending))
        env=_native_env()
        with (DEST/'native.log').open('w') as log:
            subprocess.run([str(EXE),str(ROOT),str(_data_root()/'partial_weld'),str(request),str(DEST/'statistics')],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        for j in pending:s.compact(DEST/'statistics'/j['id'])
    print('Native sufficient statistics exported for both anchors.',flush=True)


def replay(a,case,q,n,e,minimum_sources):
    accepted=a[case+'_accepted'].reshape(5,-1)
    eligible=a[case+'_eligible'].reshape(5,-1)
    effective=a[case+'_effective'].reshape(5,-1)
    valid=(accepted>=np.maximum(n,np.ceil(q*eligible)))&(effective>=e-1e-10)
    means=a[case+'_means'].reshape(5,-1)
    count=valid.sum(axis=0)
    amplitude=np.sqrt(np.sum(np.where(valid,abs(means)**2,0),axis=0)/np.maximum(count,1))
    amplitude[count<minimum_sources]=np.nan
    return amplitude,count,valid


def score(v,points,distance):
    ok=np.isfinite(v);near=ok&(distance<=.5);back=ok&(distance>1)
    peak=int(np.argmax(np.where(ok,v,-np.inf))) if ok.any() else None
    a=float(np.max(v[near])) if near.any() else 0.
    b=float(np.max(v[back])) if back.any() else 0.
    return dict(cmax_db=float(20*np.log10(a/b)) if a>0 and b>0 else None,
        near_support_pct=float(100*np.mean(ok[distance<=.5])),support_pct=float(100*np.mean(ok)),
        near_peak_um_s=a*1e6,background_peak_um_s=b*1e6,
        global_peak_distance_mm=float(distance[peak]) if peak is not None else None,
        global_peak_x_mm=float(points[peak,0]*1000) if peak is not None else None,
        global_peak_depth_mm=float(points[peak,1]*1000) if peak is not None else None)


def analyse():
    p=protocol();rows=[];diagnostics=[];parity=[]
    for anchor in p['anchors']:
        with np.load(DEST/'statistics'/anchor['id']/'arrays.npz') as f:a=dict(f)
        points=a['points'].reshape(-1,2);distance=s.distances(points,'partial');near=distance<=.5
        with np.load(OUT/anchor['archive']/'arrays.npz') as f:old=dict(f)
        for case in ('lof','reference','difference'):
            base,basecount,basevalid=replay(a,case,.25,8,8,2)
            np.testing.assert_allclose(base,a[case+'_amplitude'],rtol=2e-14,atol=1e-18,equal_nan=True)
            np.testing.assert_allclose(base,old[case+'_2_amplitude'],rtol=2e-8,atol=2e-15,equal_nan=True)
            np.testing.assert_array_equal(basecount,old[case+'_2_count'])
            np.testing.assert_array_equal(basevalid.ravel(),a[case+'_valid'].astype(bool))
            parity.append(dict(anchor=anchor['id'],case=case,maximum_difference=float(np.nanmax(abs(base-old[case+'_2_amplitude'])))))
            expected=s.metrics(base,points,'partial')
            assert abs(score(base,points,distance)['cmax_db']-expected['cmax_db'])<1e-10
            accepted=a[case+'_accepted'].reshape(5,-1)
            eligible=a[case+'_eligible'].reshape(5,-1)
            effective=a[case+'_effective'].reshape(5,-1)
            missing=near&~np.isfinite(base)
            diagnostics.append(dict(anchor=anchor['id'],case=case,near_pixels=int(near.sum()),
                missing_near_pixels=int(missing.sum()),missing_points_mm=(points[missing]*1000).tolist(),
                accepted_at_missing=accepted[:,missing].T.tolist(),eligible_at_missing=eligible[:,missing].T.tolist(),
                effective_at_missing=effective[:,missing].T.tolist(),
                near_source_opportunities=int(5*near.sum()),
                near_zero_eligible_source_pixels=int(np.sum((eligible==0)[:,near])),
                near_fraction_failures_with_eligible=int(np.sum(((eligible>0)&(accepted<np.ceil(.25*eligible)))[:,near])),
                near_minimum_count_failures_with_eligible=int(np.sum(((eligible>0)&(accepted<8))[:,near])),
                near_effective_count_failures_with_eligible=int(np.sum(((eligible>0)&(effective<8-1e-10))[:,near]))))
            for q,n,e,k in itertools.product(FRACTIONS,MINIMUMS,EFFECTIVE,SOURCES):
                v,count,valid=replay(a,case,q,n,e,k)
                common=np.isfinite(base)&np.isfinite(v)
                row=dict(anchor=anchor['id'],case=case,receiver_fraction=q,minimum_receivers=n,
                    minimum_effective_receivers=e,minimum_sources=k,**score(v,points,distance),
                    gained_near_pixels=int(np.sum(near&np.isfinite(v)&~np.isfinite(base))),
                    lost_near_pixels=int(np.sum(near&~np.isfinite(v)&np.isfinite(base))),
                    changed_qualifying_source_pixels=int(np.sum(valid!=basevalid)),
                    maximum_common_amplitude_change_um_s=float(np.max(abs(v[common]-base[common]))*1e6) if common.any() else None,
                    near_one_source_pct=float(100*np.mean((np.isfinite(v)&(count==1))[near])))
                rows.append(row)
        print(anchor['id'],'analysed 300 combinations in three signal cases.',flush=True)
    s.write(DEST/'all_results.json',dict(rows=rows))
    with (DEST/'all_results.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    s.write(DEST/'qualification_diagnostics.json',dict(rows=diagnostics))
    s.write(DEST/'baseline_verification.json',dict(status='passed',checks=parity,
        native_executable_sha256=s.sha(EXE),exporter_sha256=s.sha(HERE/'admission_cli.cpp')))
    primary=[]
    for anchor in p['anchors']:
        for q,n,e,k,label in [
            (0.,8,8,2,'fraction off'),(.1,8,8,2,'10%'),(.25,8,8,2,'baseline 25%'),(.5,8,8,2,'50%'),(1.,8,8,2,'100%'),
            (.25,1,8,2,'minimum count 1'),(.25,4,4,2,'both count floors 4'),(.25,16,16,2,'both count floors 16'),
            (.25,8,0,2,'effective floor off'),(0.,1,0,2,'all receiver count guards relaxed'),
            (.25,8,8,1,'minimum one source'),(.25,8,8,3,'minimum three sources')]:
            row=next(r for r in rows if r['anchor']==anchor['id'] and r['case']=='lof' and
                     (r['receiver_fraction'],r['minimum_receivers'],r['minimum_effective_receivers'],r['minimum_sources'])==(q,n,e,k))
            primary.append(dict(label=label,**row))
            print(anchor['id'],label,'contrast',round(row['cmax_db'],3) if row['cmax_db'] is not None else None,
                  'near',round(row['near_support_pct'],2),'field',round(row['support_pct'],2),
                  'changed source-pixels',row['changed_qualifying_source_pixels'],flush=True)
    s.write(DEST/'primary_comparisons.json',dict(rows=primary))


def confirm():
    jobs=[]
    for a in anchors():
        for label,changes in [('fraction_off',dict(receiver_fraction=0)),
            ('fraction_one',dict(receiver_fraction=1)),
            ('all_relaxed',dict(receiver_fraction=0,minimum_receivers=1,minimum_effective_receivers=0)),
            ('counts16',dict(minimum_receivers=16,minimum_effective_receivers=16))]:
            jobs.append(dict(id=a['id']+'__'+label,selection=dict(a['selection'],**changes),anchor=a['id']))
    p=DEST/'confirmation_protocol.json'
    if p.exists():assert json.loads(p.read_text())==dict(requests=jobs)
    else:s.write(p,dict(requests=jobs))
    pending=[dict(id=j['id'],selection=j['selection']) for j in jobs if not (DEST/'confirmation'/j['id']/'arrays.npz').exists()]
    if pending:
        s.execute('partial',dict(requests=pending),DEST/'confirmation')
        for j in pending:s.compact(DEST/'confirmation'/j['id'])
    checks=[]
    for job in jobs:
        with np.load(DEST/'statistics'/job['anchor']/'arrays.npz') as f:a=dict(f)
        with np.load(DEST/'confirmation'/job['id']/'arrays.npz') as b:
            c=job['selection']
            for case in ('lof','reference','difference'):
                v,count,valid=replay(a,case,c['receiver_fraction'],c['minimum_receivers'],c['minimum_effective_receivers'],c['minimum_sources'])
                np.testing.assert_allclose(v,b[case+'_2_amplitude'],rtol=2e-8,atol=2e-15,equal_nan=True)
                np.testing.assert_array_equal(count,b[case+'_2_count'])
                np.testing.assert_array_equal(valid.ravel(),b[case+'_2_valid'].astype(bool))
                checks.append(dict(id=job['id'],case=case,max_absolute_difference=float(np.nanmax(abs(v-b[case+'_2_amplitude'])))))
    s.write(DEST/'confirmation_verification.json',dict(status='passed',native_configurations=len(jobs),signal_cases=len(checks),checks=checks))
    print('All threshold replays agree with independent native configurations.',flush=True)


def angular():
    import sys
    sys.path.insert(0,str(STUDY/'src'))
    from verify_saft_datasets import PartialReference
    from saft_apertures import directivity,angular_gain
    ref=PartialReference(_data_root()/'partial_weld');sp=ref.audit['source_parameters']
    curves=directivity(('audit',1,5),ref.cl,ref.cs,sp['line_half_width_m'],sp['pulse_fwhm_s'])
    diagnostics=json.loads((DEST/'qualification_diagnostics.json').read_text())['rows']
    missing=next(r for r in diagnostics if r['anchor']=='weighted_1_5' and r['case']=='lof')['missing_points_mm']
    rows=[]
    for point in missing:
        for si,source in enumerate(ref.sources):
            ray=ref.ray(source,np.array(point)/1000,'S')
            gain=float(angular_gain(np.array([ray['first']]),curves['S']['source'])[0]) if ray else None
            reason='blocked incoming geometry' if ray is None else 'below 20% launch gate' if gain<.2 else 'launch admitted'
            rows.append(dict(point_mm=point,source=si,source_id=ref.audit['sources'][si]['id'],launch_gain=gain,reason=reason))
    s.write(DEST/'missing_pixel_launch_audit.json',dict(rows=rows,
        method='Independent Python geometry/directivity implementation; incoming S leg, actual 1-5 MHz FIR response.'))
    base=next(a['selection'] for a in anchors() if a['id']=='weighted_1_5')
    jobs=[dict(id=f'gate{g:g}',selection=dict(base,angular_gate=g)) for g in (.15,.1,.05)]
    jobs.append(dict(id='gate0.1_fraction_off',selection=dict(base,angular_gate=.1,receiver_fraction=0)))
    p=dict(reason='Receiver count guards are inactive; all three weld sources at the 14 missing near pixels fail the 20% incoming-source angular gate. Test the existing GUI angular threshold in the recommended 1-5 MHz band.',
        controls='The production angular threshold applies jointly to launch and receive response. Keep two required sources, other count floors, screening, filter and field. The final control tests whether the 25% rule becomes active after angular admission changes.',requests=jobs)
    path=DEST/'angular/protocol.json'
    if path.exists():assert json.loads(path.read_text())==p
    else:s.write(path,p)
    pending=[j for j in jobs if not (DEST/'angular'/j['id']/'arrays.npz').exists()]
    if pending:
        s.execute('partial',dict(requests=pending),DEST/'angular')
        for job in pending:s.compact(DEST/'angular'/job['id'])
    results=[]
    for j in jobs:
        with np.load(DEST/'angular'/j['id']/'arrays.npz') as a:
            points=a['points'].reshape(-1,2);d=s.distances(points,'partial')
            for case in ('lof','reference','difference'):
                results.append(dict(id=j['id'],case=case,selection=j['selection'],**score(a[case+'_2_amplitude'],points,d)))
    s.write(DEST/'angular/metrics.json',dict(rows=results))
    for r in results:
        if r['case']=='lof':print(r['id'],r['cmax_db'],r['near_support_pct'],flush=True)


def angular_confirm():
    base=next(a['selection'] for a in anchors() if a['id']=='weighted_1_5')
    jobs=[]
    for label,changes in [('fine',dict(spacing_mm=.05)),
                          ('wide',dict(spacing_mm=.25,x_min_mm=-45,x_max_mm=45,depth_min_mm=.25,depth_max_mm=25.25))]:
        jobs.append(dict(id=label,selection=dict(base,angular_gate=.15,**changes)))
    path=DEST/'angular_confirmation/protocol.json'
    p=dict(reason='The 15% joint angular gate has slightly better 0.1 mm contrast and support than 20%; check finer/wider fields before suggesting it as an alternative.',requests=jobs)
    if path.exists():assert json.loads(path.read_text())==p
    else:s.write(path,p)
    pending=[j for j in jobs if not (DEST/'angular_confirmation'/j['id']/'arrays.npz').exists()]
    if pending:
        s.execute('partial',dict(requests=pending),DEST/'angular_confirmation')
        for j in pending:s.compact(DEST/'angular_confirmation'/j['id'])
    rows=[]
    for j in jobs:
        with np.load(DEST/'angular_confirmation'/j['id']/'arrays.npz') as a,np.load(OUT/'confirmation'/('weighted_1_5__'+j['id'])/'arrays.npz') as b:
            points=a['points'].reshape(-1,2);d=s.distances(points,'partial')
            np.testing.assert_array_equal(a['points'],b['points'])
            for case in ('lof','reference','difference'):
                rows.append(dict(id=j['id'],case=case,before=score(b[case+'_2_amplitude'],points,d),
                                 after=score(a[case+'_2_amplitude'],points,d)))
            if j['id']=='fine':
                with np.load(DEST/'angular/gate0.15/arrays.npz') as coarse:
                    for case in ('lof','reference','difference'):
                        v=a[case+'_2_amplitude'].reshape(len(a['z']),len(a['x']))[::2,::2].ravel()
                        np.testing.assert_allclose(v,coarse[case+'_2_amplitude'],rtol=2e-8,atol=2e-15,equal_nan=True)
    s.write(DEST/'angular_confirmation/metrics.json',dict(rows=rows,nested_grid='passed'))
    session=json.loads((OUT/'settings/saft_partial_weighted_1_5.json').read_text())
    session['selection']['angular_gate']=.15
    s.write(DEST/'saft_partial_weighted_1_5_angular15.json',session)
    for r in rows:
        if r['case']=='lof':print(r['id'],'before',r['before'],'after',r['after'],flush=True)


def report():
    p=protocol()
    rows=json.loads((DEST/'all_results.json').read_text())['rows']
    assert len(rows)==1800
    assert all(r['changed_qualifying_source_pixels']==0 for r in rows)
    # Source-count changes only mask additional pixels; receiver guards do not
    # change any qualifying source or any amplitude within a fixed source count.
    for anchor,case,k in itertools.product(('current_equal','weighted_1_5'),('lof','reference','difference'),SOURCES):
        group=[r for r in rows if (r['anchor'],r['case'],r['minimum_sources'])==(anchor,case,k)]
        assert len(group)==100
        assert all((r['cmax_db'],r['near_support_pct'],r['support_pct'])==
                   (group[0]['cmax_db'],group[0]['near_support_pct'],group[0]['support_pct']) for r in group)
    counts=[]
    for anchor in p['anchors']:
        with np.load(DEST/'statistics'/anchor['id']/'arrays.npz') as a:
            for case in ('lof','reference','difference'):
                accepted=a[case+'_accepted'];eligible=a[case+'_eligible'];effective=a[case+'_effective']
                populated=eligible>0
                np.testing.assert_array_equal(accepted,eligible)
                counts.append(dict(anchor=anchor['id'],case=case,
                    minimum_accepted_when_eligible=float(accepted[populated].min()),
                    minimum_effective_when_eligible=float(effective[populated].min()),
                    minimum_retained_fraction=float((accepted[populated]/eligible[populated]).min())))
    base=next(r for r in rows if (r['anchor'],r['case'],r['receiver_fraction'],r['minimum_receivers'],
              r['minimum_effective_receivers'],r['minimum_sources'])==('weighted_1_5','lof',.25,8,8,2))
    metrics=json.loads((DEST/'angular/metrics.json').read_text())['rows']
    angular=[dict(angular_gate=.2,near_support_pct=base['near_support_pct'],cmax_db=base['cmax_db'])]
    angular.extend(dict(angular_gate=r['selection']['angular_gate'],near_support_pct=r['near_support_pct'],cmax_db=r['cmax_db'])
                   for r in metrics if r['case']=='lof' and r['id']!='gate0.1_fraction_off')
    confirmations=json.loads((DEST/'angular_confirmation/metrics.json').read_text())
    checks=json.loads((DEST/'confirmation_verification.json').read_text())
    assert checks['status']=='passed' and confirmations['nested_grid']=='passed'
    summary=dict(status='passed',threshold_combinations=600,signal_case_evaluations=1800,
        receiver_guard_invariance='passed for all combinations at fixed source count',
        native_confirmation_configurations=checks['native_configurations'],counts=counts,angular_comparison=angular,
        finer_and_wider=[r for r in confirmations['rows'] if r['case']=='lof'])
    s.write(DEST/'summary.json',summary)
    table='\n'.join(f"  {100*r['angular_gate']:4.0f}%         {r['near_support_pct']:6.2f}%            {r['cmax_db']:5.2f} dB" for r in angular)
    body=f'''PARTIAL-WELD SAFT ADMISSION AUDIT - 2 NOVEMBER 2026

The 25% receiver-retention rule and the receiver-count floors do not restrict
the current S | SS reconstruction. Lowering them does not recover any pixels.
The angular response threshold does affect the missing near-defect coverage.

WHAT THE 25% RULE MEANS
A source contributes at a pixel when:
  accepted >= max(minimum_receivers, ceil(fraction * eligible))
  effective >= minimum_effective_receivers (native tolerance 1e-10)
The image then requires at least two contributing sources in the main comparison.
Eligible receivers have already passed geometry, angular response, wall,
source-receiver offset and aperture gates. Accepted receivers additionally have
complete filter/interpolation support and pass pulse screening. Thus 25% refers
to this eligible set for one source/pixel, not 25% of all 331 receivers.
Effective count is (sum weights)^2 / sum(weights squared); it measures weight
concentration, not a count of statistically independent observations.

TESTED RANGE
Fractions: 0, 10, 25, 50, 100%; minimum receivers: 1, 4, 8, 16, 32;
minimum effective receivers: 0, 4, 8, 16; minimum sources: 1, 2, 3.
All combinations were replayed for the equal-mean 0.5-5 MHz baseline and the
weighted 1-5 MHz candidate: 600 combinations, each in LOF, pristine and
waveform-difference cases. These are statistic replays, not new COMSOL solves.
Eight separate native reconstructions checked the replay (24 signal cases).
No qualifying source changed across receiver-guard choices. Where eligible
receivers exist, all survive; the minimum accepted count is 224. The weighted
effective count is at least 189.09. These conclusions apply to these fields
and recipes; they do not establish that the guards are always unnecessary.

THE ACTIVE LIMIT
On the 0.1 mm grid, 14 of 330 pixels within 0.5 mm of the defect lack two-source
support. The left shoulder's incoming ray is blocked by geometry, while all
three weld sources fall below the 20% launch-response threshold. Only the
right shoulder is admitted. The known defect line itself already had two-source
support at all 101 sampled locations in the preceding support audit.

Fresh native angular-threshold comparisons, weighted 1-5 MHz, two sources:
  Angular gate  Near support    Peak-to-background contrast
{table}
The Workbench angular threshold applies jointly to launch and receive response;
these tests change both. The 10% run with the fraction guard disabled is also
unchanged. Relaxing to 15% adds three near pixels on the 0.1 mm grid. On the
0.05 mm grid, coverage changes 95.82 -> 96.81%, contrast 9.08 -> 9.10 dB.
On the wider 0.25 mm grid, coverage changes 94.34 -> 96.23%, contrast
9.20 -> 9.12 dB. This is a small coverage gain with essentially unchanged
contrast, not an established contrast improvement. The 10% and 5% gates obtain
100% neighbourhood support but introduce substantially stronger background.

METRICS AND LIMITS
Near support is the fraction of sampled pixels within 0.5 mm of the prescribed
defect segment with a valid image value. It is not detection probability.
Contrast is 20 log10(near peak / distant background peak), with background
more than 1 mm from the segment. Unsupported pixels are excluded. Coverage
changes can expose new background, so contrast is not evaluated on a fixed mask.
The 15% choice is retrospective on one known simulated flaw; it is an optional
sensitivity setting, not a validated optimum or improved defect-tip sizing.
Geometry checks remain active. No COMSOL model, raw data or production engine
was changed. These are reconstruction changes using existing recordings.

FILES (relative to this README)
  saft_partial_weighted_1_5_angular15.json: optional Workbench preset.
  summary.json: compact verified comparison and receiver-count ranges.
  all_results.csv / all_results.json: all threshold replay outcomes.
  protocol.json: declared matrix, anchor settings and engine source hashes.
  statistics/: native per-source arrays and provenance for both anchors.
  qualification_diagnostics.json: source counts at missing near pixels.
  missing_pixel_launch_audit.json: independent ray/launch-gain diagnosis.
  baseline_verification.json / confirmation_verification.json: native parity.
  angular/: requests, arrays and metrics for 15%, 10%, 5% and fraction-off control.
  angular_confirmation/: finer/wider requests, arrays, metrics and grid check.

REPRODUCE FROM THE REPOSITORY ROOT
Prerequisites: existing partial-weld recordings, previous improvement study
settings/arrays, pinned Qt/MinGW toolchain and SAFTFilterStudy executable.
Build the exporter, which links the unchanged production engine:
  powershell -NoProfile -ExecutionPolicy Bypass -File studies/laser/paper/partial_improvement_20261102/build_admission.ps1
Run each command in this order with:
  python studies/laser/paper/partial_improvement_20261102/admission_audit.py COMMAND
Commands: prepare, dump, analyse, confirm, angular, angular_confirm, report.
Completed array exports are reused; frozen protocols must agree before reuse.
'''
    (DEST/'README.txt').write_text(body,encoding='utf-8')
    print('Saved verified admission summary, README and optional 15% preset.',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['prepare','dump','analyse','confirm','angular','angular_confirm','report'])
    cmd=parser.parse_args().command
    if cmd=='prepare':print('Frozen',protocol()['combinations_per_anchor'],'threshold combinations per anchor.',flush=True)
    else:globals()[cmd]()
