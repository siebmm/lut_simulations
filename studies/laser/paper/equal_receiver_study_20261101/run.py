"""Repeat the declared partial-weld survey with equal receiver averaging.

The path/filter/grid/acquisition matrix is fixed before inspecting outputs.
Original weighted study products are preserved in their existing directory.
"""
from pathlib import Path
import copy
import json
import sys
import time
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent / 'filter_study_20261101'))
import study as s
import confirm

ORIGINAL = s.OUT
s.OUT = ROOT / 'studies/laser/results/equal_receiver_study_20261101'
original_base = s.base


def base(kind, spacing=.2):
    config = original_base(kind, spacing)
    config['processing'] = 'raw'
    return config


s.base = base


def main():
    s.OUT.mkdir(parents=True, exist_ok=True)
    protocol = json.loads((s.HERE / 'protocol.json').read_text())
    protocol['base']['partial'] = base('partial')
    protocol['created_utc'] = time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())
    protocol['receiver_combination'] = 'Equal receiver mean after travel-time alignment for both welds; no modeled amplitude weights or receiving/wall phase corrections. Geometry, physical gates, pulse screening, coherence exponent and all other controls retained.'
    protocol['scope'] = 'Recompute all 252 partial-weld survey and 28 confirmation configurations; retain the 104 existing full-weld configurations.'
    protocol['selection'] = 'User-selected simpler common receiver combination, fixed before the repeat outputs. All 18 paths and 14 filters retained. Baseline S|SS stays the declared principal partial-weld example, subject to its supported evidence.'
    path = HERE / 'protocol.json'
    if path.exists():
        previous = json.loads(path.read_text())
        assert previous['base'] == protocol['base'] and previous['variants'] == protocol['variants']
        protocol = previous
    else:
        s.write(path, protocol)
    jobs = confirm.plan('partial')
    target = s.OUT / 'confirmation/partial/protocol.json'
    if target.exists():
        assert json.loads(target.read_text()) == jobs
    else:
        s.write(target, jobs)
    s.write(s.OUT / 'provenance.json', dict(
        protocol_sha256=s.sha(path), executable_sha256=s.sha(s.EXE),
        source_study_sha256=s.sha(s.HERE / 'study.py'),
        confirmation_driver_sha256=s.sha(s.HERE / 'confirm.py'),
        dataset_fingerprint=s.dataset('partial')['input_fingerprint']))
    for spec in protocol['variants']:
        name = spec['name']
        output = s.OUT / 'numerical/partial' / name
        survey = [dict(id=p, selection=dict(base('partial'), path=p)) for p in s.PATHS]
        group = [dict(id=j['id'], selection=j['selection']) for j in jobs if j['variant'] == name]
        pending = [(output, survey), (s.OUT / 'confirmation/partial', group)]
        if all(all((out/j['id']/'arrays.npz').exists() for j in requests) for out, requests in pending):
            continue
        bank = s.bank('partial', spec)
        for out, requests in pending:
            requests = [dict(**j, signal_processing=dict(spec, bank=str(bank))) for j in requests
                        if not (out/j['id']/'arrays.npz').exists()]
            if not requests:
                continue
            print(name, out.parent.name, len(requests), 'jobs', flush=True)
            s.execute('partial', dict(requests=requests), out)
            for job in requests:
                s.compact(out/job['id'])
            print(name, out.parent.name, 'complete', flush=True)
    confirm.summarise('partial')
    rows = []
    for spec in protocol['variants']:
        for p in s.PATHS:
            folder = s.OUT / 'numerical/partial' / spec['name'] / p
            with np.load(folder/'arrays.npz') as a, np.load(s.OUT/'numerical/partial/baseline'/p/'arrays.npz') as b:
                points = a['points'].reshape(-1, 2)
                mask = np.isfinite(b['lof_2_amplitude'])
                for case in ['lof', 'reference', 'difference']:
                    values = a[case+'_2_amplitude']
                    rows.append(dict(kind='partial', variant=spec['name'], path=p, case=case,
                        **s.metrics(values, points, 'partial'),
                        common_baseline_support=s.metrics(values, points, 'partial', mask)))
    s.write(s.OUT/'metrics.json', dict(rows=rows))
    # Independent prior native run of the same equal-mean conventional route.
    import csv
    control = list(csv.DictReader((ROOT/'studies/laser/paper/comsol_chapter_20261101/data/geometry_extension/all_path_metrics.csv').open()))
    actual = json.loads((s.OUT/'confirmation/partial/metrics.json').read_text())['rows']
    for case in ['lof', 'reference', 'difference']:
        old = next(r for r in control if r['path']=='1-2_SSS' and r['processing']=='raw' and r['case']==case)
        new = next(r for r in actual if r['id']=='baseline__1-2_SSS__fine' and r['case']==case)
        np.testing.assert_allclose(new['cmax_db'], float(old['Cmax_dB']), rtol=0, atol=1e-8)
        np.testing.assert_allclose(new['near_support_pct'], float(old['near_support_pct']), rtol=0, atol=1e-9)
    s.write(s.OUT/'verification/reproduction.json', dict(status='passed',
        configurations=280, signal_cases=840, all_processing='raw',
        independent_baseline_cases=3, contrast_tolerance_db=1e-8))
    print('Equal-receiver survey and confirmations complete; independent baseline reproduced.', flush=True)


if __name__ == '__main__':
    if '--verify' in sys.argv:
        s.verify()
    else:
        main()
