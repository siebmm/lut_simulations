"""Recompute the published full/partial images from checked receiver histories."""
from pathlib import Path
import argparse
import json
import os
import sys
import numpy as np

from replay_workspace import ROOT, prepare_workspace


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', choices=['both', 'full', 'partial'], default='both')
    args = parser.parse_args()
    work, env = prepare_workspace()
    STUDY = work/'studies/laser'
    os.environ['SAFT_REPOSITORY'] = env['SAFT_REPOSITORY']
    sys.path[:0] = [str(STUDY/'paper/filter_study_20261101'), str(STUDY/'src')]
    import study as s
    from standalone import executable
    dest = ROOT / 'output/main'
    s.OUT = dest
    s.LOCAL = ROOT / 'build/filter_work'
    s.EXE = executable('SAFTFilterStudy')
    if not s.EXE.exists():
        raise SystemExit('Build the native executables first; see README.md.')
    s.LOCAL.mkdir(parents=True, exist_ok=True)
    # Each process finds its Qt/C++ runtime beside the executable.
    checks = []
    for kind in ['full', 'partial']:
        if args.case not in ('both', kind):
            continue
        if kind == 'full':
            saved = STUDY / 'results/filter_study_20261101/selected_checks/full'
            protocol = json.loads((saved / 'protocol.json').read_text())
            job = next(j for j in protocol['jobs'] if j['id'] == 'all')
            raw = s.LOCAL / 'raw/full'
            s.execute('full', {'dump_raw': True}, raw)
            print('Full-weld receiver hashes verified; filtering histories.', flush=True)
            bank = s.bank('full', protocol['spec'])
            jobs = [dict(job, signal_processing=dict(protocol['spec'], bank=str(bank)))]
            expected = [r for r in json.loads((saved / 'metrics.json').read_text())['rows'] if r['id'] == 'all']
        else:
            saved = STUDY / 'results/partial_improvement_20261102/review/refinement'
            protocol = json.loads((saved / 'protocol.json').read_text())
            jobs = [j for j in protocol['requests'] if j['id'] in ('muted_025', 'unmuted_025')]
            expected = [r for r in json.loads((saved / 'metrics.json').read_text())['rows'] if r['id'] in ('muted_025', 'unmuted_025')]
        out = dest / kind
        print(f'Reconstructing {kind}: ' + ', '.join(j['id'] for j in jobs), flush=True)
        s.execute(kind, {'requests': jobs}, out)
        for job in jobs:
            folder = out / job['id']
            s.compact(folder)
            with np.load(folder / 'arrays.npz') as arrays:
                points = arrays['points'].reshape(-1, 2)
                for case in ['lof', 'reference', 'difference']:
                    actual = s.metrics(arrays[case + '_2_amplitude'], points, kind)
                    old = next(r for r in expected if r['id'] == job['id'] and r['case'] == case)
                    paired = {
                        'cmax_db': (actual['cmax_db'], old['cmax_db']),
                        'near_support_pct': (actual['near_support_pct'], old['near_support_pct']),
                        'support_pct': (actual['support_pct'], old['support_pct']),
                        'near_um_s': (actual['near']['um_s'], old['near']['um_s'] if kind == 'full' else old['near_um_s']),
                        'background_um_s': (actual['background']['um_s'], old['background']['um_s'] if kind == 'full' else old['background_um_s']),
                        'peak_distance_mm': (actual['peak']['distance_mm'], old['peak']['distance_mm'] if kind == 'full' else old['peak_distance_mm']),
                    }
                    for name, (a, b) in paired.items():
                        if not np.isclose(a, b, rtol=2e-8, atol=1e-8):
                            raise AssertionError((kind, job['id'], case, name, a, b))
                    checks.append(dict(kind=kind, id=job['id'], case=case, passed=True,
                                       metrics=actual, comparison={k: dict(actual=a, archived=b) for k, (a, b) in paired.items()}))
                    if case == 'lof':
                        print(f"{kind}/{job['id']}: {actual['cmax_db']:.8f} dB, {actual['near_support_pct']:.8f}% near support: PASS", flush=True)
    previous = ROOT / 'output/reproduction_check.json'
    if args.case != 'both' and previous.exists():
        checks = [c for c in json.loads(previous.read_text())['checks'] if c['kind'] != args.case] + checks
    report = dict(status='passed', scope='Fresh reconstruction of published grids, plus matched unmuted partial control; no COMSOL solve.',
                  executable=str(s.EXE.relative_to(ROOT)) if s.EXE.is_relative_to(ROOT) else str(s.EXE), checks=checks)
    s.write(ROOT / 'output/reproduction_check.json', report)
    print('Saved output/reproduction_check.json', flush=True)


if __name__ == '__main__':
    main()
