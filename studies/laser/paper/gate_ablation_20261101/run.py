"""Paired ablation of wall-response and source/receiver exclusion gates.

Reuses the checked production C++ study engine and published equal-mean
configurations. Never modifies those results or the shared signal-bank format.
Run --phase primary first, then --phase complete for all declared confirmations.
"""
from pathlib import Path
import argparse
import copy
import json
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STUDY = ROOT / 'studies/laser'
PAPER = STUDY / 'paper/comsol_chapter_20261101'
sys.path.insert(0, str(HERE.parent / 'filter_study_20261101'))
import study as s

OUT = STUDY / 'results/gate_ablation_20261101'
s.OUT = OUT
CONTROLS = {
    'wall_off': {'wall_gate': 0.0},
    'separation_off': {'min_offset_mm': 0.0},
    'both_off': {'wall_gate': 0.0, 'min_offset_mm': 0.0},
}
PRINCIPAL = {
    'full': 'selected_checks_full_all',
    'partial': 'confirmation_partial_baseline__1-2_SSS__fine',
}


def baseline_folder(row):
    origin = ('filter_study_20261101' if row['kind'] == 'full'
              else 'equal_receiver_study_20261101')
    path = STUDY / 'results' / origin / row['stage'] / row['kind']
    if row['stage'] == 'numerical':
        path /= row['variant']
    return path / row['id']


def is_primary(row):
    return (row['stage'] == 'numerical' and row['variant'] == 'baseline'
            or row['uid'] in PRINCIPAL.values())


def protocol():
    images = json.loads((PAPER / 'data/images.json').read_text())
    images = [r for r in images if r['stage'] != 'numerical'
              or r['variant'] == 'baseline']
    assert len(images) == 98
    for row in images:
        c = row['selection']
        assert c['processing'] == 'raw' and c['wall_gate'] == .05 and c['min_offset_mm'] == 1.5
    fixed = dict(
        purpose='Measure the independent and joint effects of the 5% wall gate and 1.5 mm source/receiver exclusion; all other settings fixed.',
        controls=CONTROLS, configurations=images,
        principal=PRINCIPAL,
        scope='22 baseline path-survey configurations plus all 76 finer/wider/filter/acquisition confirmations; three recorded signal cases each.',
        assessment='Compare complete untruncated image arrays and support, accepted-channel summaries and per-source coherent means. Known-flaw contrast/coverage/localisation are secondary diagnostics, not a tuning objective.',
        negligible_checks=dict(
            common_support_relative_l2_pct=1.0,
            common_support_max_change_pct_of_baseline_peak=1.0,
            support_changed_percentage_points=1.0,
            contrast_change_db=.1,
            note='Declared engineering indicators, not statistical confidence bounds. A clean primary image alone is insufficient to remove a gate globally.'),
        executable_sha256=s.sha(s.EXE),
        dataset_fingerprints={k:s.dataset(k)['input_fingerprint'] for k in ('full','partial')},
        baseline_sha256={r['uid']:s.sha(baseline_folder(r)/'arrays.npz') for r in images},
    )
    path = HERE / 'protocol.json'
    if path.exists():
        assert json.loads(path.read_text()) == fixed
    else:
        s.write(path, fixed)
    return fixed


def array_comparison(a, b):
    va, vb = np.isfinite(a), np.isfinite(b)
    common = va & vb
    d = b[common]-a[common]
    scale = float(np.max(np.abs(a[va]))) if va.any() else 0.0
    denom = float(np.linalg.norm(a[common]))
    return dict(
        exactly_equal=bool(np.array_equal(a,b,equal_nan=True)),
        baseline_support_pct=float(va.mean()*100),
        ablated_support_pct=float(vb.mean()*100),
        support_changed_pp=float(np.mean(va != vb)*100),
        gained_pixels=int(np.sum(vb & ~va)), lost_pixels=int(np.sum(va & ~vb)),
        common_pixels=int(common.sum()),
        common_relative_l2_pct=float(np.linalg.norm(d)/denom*100) if denom else None,
        common_max_change_pct_peak=float(np.max(np.abs(d))/scale*100) if d.size and scale else None,
    )


def summarise(p):
    rows=[]
    for cfg in p['configurations']:
        with np.load(baseline_folder(cfg)/'arrays.npz') as baseline:
            for label in CONTROLS:
                folder=OUT/cfg['kind']/(cfg['uid']+'__'+label)
                if not (folder/'arrays.npz').exists():
                    continue
                request=json.loads((folder/'request.json').read_text())
                expected=dict(cfg['selection'],**CONTROLS[label])
                assert request['selection']==expected
                with np.load(folder/'arrays.npz') as ablated:
                    np.testing.assert_array_equal(baseline['points'],ablated['points'])
                    for case in ('lof','reference','difference'):
                        stem=case+'_2'
                        a,b=baseline[stem+'_amplitude'],ablated[stem+'_amplitude']
                        before=s.metrics(a,baseline['points'].reshape(-1,2),cfg['kind'])
                        after=s.metrics(b,ablated['points'].reshape(-1,2),cfg['kind'])
                        deltas={name:(after[name]-before[name] if after[name] is not None and before[name] is not None else None)
                                for name in ('cmax_db','c99_db','near_support_pct','support_pct')}
                        rows.append(dict(uid=cfg['uid'],kind=cfg['kind'],path=cfg['path'],
                            variant=cfg['variant'],stage=cfg['stage'],id=cfg['id'],
                            control=label,case=case,principal=cfg['uid'] in PRINCIPAL.values(),
                            **array_comparison(a,b),before=before,after=after,delta=deltas,
                            source_counts_equal=bool(np.array_equal(baseline[stem+'_count'],ablated[stem+'_count'],equal_nan=True)),
                            max_source_count_change=float(np.max(np.abs(ablated[stem+'_count']-baseline[stem+'_count']))),
                            receiver_usage_equal=bool(np.array_equal(baseline[stem+'_receiver_usage'],ablated[stem+'_receiver_usage'])),
                            baseline_receiver_pixel_contributions=int(baseline[stem+'_receiver_usage'].sum()),
                            ablated_receiver_pixel_contributions=int(ablated[stem+'_receiver_usage'].sum()),
                            source_means_equal=bool(np.array_equal(baseline[stem+'_means'],ablated[stem+'_means'],equal_nan=True))))
    s.write(OUT/'comparison.json',dict(rows=rows,configurations_completed=len(rows)//9))
    principal=[r for r in rows if r['principal'] and r['case']=='lof']
    for r in principal:
        print(r['kind'],r['control'],'exact',r['exactly_equal'],
              'L2%',r['common_relative_l2_pct'],'mask pp',r['support_changed_pp'],
              'contrast',r['before']['cmax_db'],r['after']['cmax_db'],flush=True)
    print('Compared',len(rows),'paired signal-case images.',flush=True)


def run(phase):
    p=protocol()
    selected=[r for r in p['configurations'] if phase=='complete' or is_primary(r)]
    for kind in ('full','partial'):
        for variant in dict.fromkeys(r['variant'] for r in selected if r['kind']==kind):
            configs=[r for r in selected if r['kind']==kind and r['variant']==variant]
            jobs=[]
            for row in configs:
                for label,changes in CONTROLS.items():
                    job_id=row['uid']+'__'+label
                    if not (OUT/kind/job_id/'arrays.npz').exists():
                        jobs.append(dict(id=job_id,selection=dict(row['selection'],**changes)))
            if not jobs:
                continue
            spec=configs[0]['signal_processing']
            assert all(r['signal_processing']==spec for r in configs)
            bank=s.bank(kind,spec)
            for job in jobs:
                job['signal_processing']=dict(spec,bank=str(bank))
            # Repeat the principal baseline in the same process to check the
            # archived comparator and external-bank settings independently.
            for row in configs:
                if row['uid'] in PRINCIPAL.values():
                    jobs.append(dict(id=row['uid']+'__repeat',selection=row['selection'],
                                     signal_processing=dict(spec,bank=str(bank))))
            print(kind,variant,len(jobs),'jobs start',flush=True)
            s.execute(kind,dict(requests=jobs),OUT/kind)
            for job in jobs:
                s.compact(OUT/kind/job['id'])
            for row in configs:
                if row['uid'] in PRINCIPAL.values():
                    with np.load(baseline_folder(row)/'arrays.npz') as a, np.load(OUT/kind/(row['uid']+'__repeat')/'arrays.npz') as b:
                        for key in a.files:
                            np.testing.assert_allclose(a[key],b[key],rtol=1e-12,atol=1e-18,equal_nan=True)
            print(kind,variant,'done',flush=True)
            summarise(p)
    for row in p['configurations']:
        assert s.sha(baseline_folder(row)/'arrays.npz') == p['baseline_sha256'][row['uid']]
    summarise(p)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--phase',choices=['primary','complete','summarise'],default='primary')
    args=parser.parse_args()
    if args.phase=='summarise':
        summarise(protocol())
    else:
        run(args.phase)
