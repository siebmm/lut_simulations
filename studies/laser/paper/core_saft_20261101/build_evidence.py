"""Select and independently verify the common-protocol LOF-only evidence.

This is a presentation subset of existing reconstructions, not a new FE study.
No parameters are optimized, and the existing papers/results are not modified.
"""
from pathlib import Path
import hashlib
import json
import shutil
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCE = HERE.parent / 'comsol_chapter_20261101'
DATA = HERE / 'data'
PATHS = ['1-1_LL', '1-1_SS', '1-2_SSS', '2-1_SSS']


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def metrics(a, dataset):
    xx, zz = np.meshgrid(a['x_mm'], a['depth_mm'])
    p = np.column_stack([xx.ravel(), zz.ravel()])
    tips = np.array(dataset['audit']['defect_benchmark_only']['tips_xy_mm'])
    tips[:, 1] = 25.4-tips[:, 1]
    edge = tips[1]-tips[0]
    u = np.clip((p-tips[0])@edge/(edge@edge), 0, 1)
    distance = np.linalg.norm(p-tips[0]-u[:, None]*edge, axis=1)
    v = a['amplitude'].ravel()*1e6
    valid = np.isfinite(v)
    near, back = valid & (distance <= .5), valid & (distance > 1.)
    def peak(mask):
        if not mask.any():
            return None
        j = np.argmax(np.where(mask, v, -np.inf))
        return dict(um_s=float(v[j]), x_mm=float(p[j, 0]), depth_mm=float(p[j, 1]), distance_mm=float(distance[j]))
    n, b = peak(near), peak(back)
    return dict(cmax_db=float(20*np.log10(n['um_s']/b['um_s'])) if n and b and min(n['um_s'], b['um_s']) > 0 else None,
                near_support_pct=float(valid[distance <= .5].mean()*100),
                support_pct=float(valid.mean()*100), peak=peak(valid), near=n, background=b)


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    index = json.loads((SOURCE/'data/images.json').read_text())
    wanted = {('full', f'baseline__{path}__fine') for path in PATHS}
    for kind, paths in [('full', ['1-1_LL', '1-1_SS']), ('partial', ['1-2_SSS'])]:
        for path in paths:
            for group in ['fine', 'shoulders', 'central', 'plate', 'floor', 'shoulders_plate']:
                wanted.add((kind, f'baseline__{path}__{group}'))
    wanted |= {('full', 'baseline__1-1_LL__wide'), ('partial', 'baseline__1-2_SSS__wide')}
    rows = [r for r in index if r['stage'] == 'confirmation' and (r['kind'], r['id']) in wanted]
    assert len(rows) == len(wanted) == 22
    records = []
    for kind in ['full', 'partial']:
        shutil.copy2(SOURCE/f'data/{kind}_dataset.json', DATA/f'{kind}_dataset.json')
    for r in rows:
        c = r['selection']
        assert r['variant'] == 'baseline' and c['processing'] == 'raw'
        assert c['coherence'] and c['coherence_power'] == 1 and c['pulse_screening'] == 'LSR'
        assert not c['clean'] and not any(c['mutes']) and not any(c['pulse_reflections']) and not c['freeze']
        source = SOURCE/'data'/(r['uid']+'.npz')
        with np.load(source) as a:
            fields = {name: a[name] for name in ['x_mm', 'depth_mm']}
            fields.update({name: a['lof_'+name] for name in ['amplitude', 'count', 'cf', 'retained', 'screened']})
        assert np.all(fields['count'][np.isfinite(fields['amplitude'])] >= 2)
        np.savez_compressed(DATA/(r['uid']+'.npz'), **fields)
        d = json.loads((DATA/(r['kind']+'_dataset.json')).read_text())
        measured = metrics(fields, d)
        expected = next(v for v in json.loads((SOURCE/'data'/(r['kind']+'_confirmation_metrics.json')).read_text())['rows']
                        if v['id'] == r['id'] and v['case'] == 'lof')
        for key in ['cmax_db', 'near_support_pct', 'support_pct']:
            if expected[key] is None:
                assert measured[key] is None
            else:
                np.testing.assert_allclose(measured[key], expected[key], atol=1e-9, rtol=0)
        if measured['peak']:
            for key, value in measured['peak'].items():
                np.testing.assert_allclose(value, expected['peak'][key], atol=1e-9, rtol=0)
        records.append(dict(r, metrics=measured, source_array_sha256=sha(source)))
    main = [r for r in records if r['id'].endswith('__fine')]
    assert len(main) == 5
    ignored = {'path', 'receivers', 'x_min_mm', 'x_max_mm', 'depth_min_mm', 'depth_max_mm'}
    common = {k:v for k,v in main[0]['selection'].items() if k not in ignored}
    assert all({k:v for k,v in r['selection'].items() if k not in ignored} == common for r in main)
    assert all(r['signal_processing'] == main[0]['signal_processing'] for r in main)
    write(DATA/'images.json', records)
    write(DATA/'protocol.json', dict(
        definition='LOF-only basic SAFT: common bandpass, equal receiver mean, receiver coherence, direct-arrival pulse screening and source RMS.',
        common_controls=dict(common, cases=['lof']), filter=main[0]['signal_processing'],
        principal_paths={'full': PATHS, 'partial': ['1-2_SSS']},
        presentation='One common filter and admission protocol; no response weighting, image fusion or reference-subtracted image is presented.',
        scope='22 baseline LOF-only reconstructions; the separately indexed data/joint_acquisition supplement supplies the left-side result and its controls.',
        provenance='Reuses verified existing COMSOL exports and native SAFT reconstructions; no new COMSOL solves.',
        dataset_fingerprints={k:json.loads((DATA/f'{k}_dataset.json').read_text())['input_fingerprint'] for k in ['full','partial']}))
    settings = HERE/'settings'; settings.mkdir(exist_ok=True)
    for r in main:
        dataset = json.loads((DATA/(r['kind']+'_dataset.json')).read_text())
        write(settings/(r['kind']+'_'+r['path']+'.json'), dict(
            schema=4, input_fingerprint=dataset['input_fingerprint'],
            selection=dict(r['selection'], cases=['lof']),
            display=dict(comparison='single', quantity='amplitude', scale='db',
                         normalization='fixed', signal='lof', db_range=40, db_top=0,
                         reference_um_s=10 if r['kind']=='full' else 40,
                         show_lof=False, weld=True, positions=False, match=True)))
    figs = HERE/'figures/core_saft'; figs.mkdir(parents=True, exist_ok=True)
    for name in ['geometry', 'source_calibration']:
        for ext in ['pdf', 'png']:
            shutil.copy2(SOURCE/f'figures/comsol/{name}.{ext}', figs/f'{name}.{ext}')
        target = HERE/'verification/figure_layout'; target.mkdir(parents=True, exist_ok=True)
        shutil.copy2(SOURCE/f'verification/figure_layout/{name}.json', target/f'{name}.json')
    if not (HERE/'core_heading_style.tex').exists():
        shutil.copy2(SOURCE/'comsol_heading_style.tex', HERE/'core_heading_style.tex')
    shutil.copy2(SOURCE/'comsol_references.bib', HERE/'core_references.bib')
    shutil.copy2(SOURCE/'figure_layout.py', HERE/'figure_layout.py')
    write(HERE/'verification/data_check.json', dict(status='passed', images=22, principal_images=5,
        signal='lof', equal_receiver_averaging=True, coherence=True, pulse_screening='LSR', common_filter=True,
        metrics_independently_recomputed=True, no_new_finite_element_solves=True))
    write(DATA/'sha256.json', {p.name:sha(p) for p in sorted(DATA.iterdir()) if p.is_file() and p.name != 'sha256.json'})
    print('Verified 22 LOF-only images; five principal paths have identical processing controls.')


if __name__ == '__main__':
    main()
