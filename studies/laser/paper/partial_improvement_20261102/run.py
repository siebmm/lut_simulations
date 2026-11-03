"""Explore partial-weld SAFT settings with fixed fields and pristine controls.

Uses the unchanged production native engine. All search results are retained;
the known flaw is used only for retrospective assessment, never for focusing.
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
import json
import io
import sys
import time
import zipfile
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
STUDY = ROOT / 'studies/laser'
OUT = STUDY / 'results/partial_improvement_20261102'
sys.path.insert(0, str(HERE.parent / 'filter_study_20261101'))
import study as s

s.OUT = OUT
s.LOCAL = STUDY / 'build/partial_improvement_20261102'


def base():
    c = s.base('partial', .1)
    c.update(path='1-2_SSS', processing='raw')
    return c


def plan():
    apertures = {
        'all': dict(sources=[0, 1, 2, 3, 4], groups=[True, True, True]),
        'weld_sides': dict(sources=[1, 2, 3], groups=[True, False, True]),
        'weld_left': dict(sources=[1, 2, 3], groups=[True, False, False]),
    }
    jobs = []
    def add(name, changes, role):
        jobs.append(dict(id=name, selection=dict(base(), **changes), role=role))
    for aperture, changes in apertures.items():
        for processing in ('raw', 'weighted'):
            for minimum in (2, 1):
                add(f'{aperture}__{processing}__min{minimum}',
                    dict(changes, processing=processing, minimum_sources=minimum),
                    'receiver response / physical aperture / source-support comparison')
    for aperture, changes in apertures.items():
        for gate in (.1, .05):
            add(f'{aperture}__weighted__gate{gate:g}',
                dict(changes, processing='weighted', angular_gate=gate),
                'relax angular response admission while retaining response weights and two-source requirement')
    for aperture, minimum in [('all', 2), ('weld_sides', 2), ('weld_left', 1)]:
        for name, screening in [('tail20', dict(screen_tail_db=20)),
                                ('tail40', dict(screen_tail_db=40)),
                                ('margin005', dict(screen_margin_us=.05))]:
            add(f'{aperture}__weighted__min{minimum}__{name}',
                dict(apertures[aperture], processing='weighted', minimum_sources=minimum, **screening),
                'direct-packet screen sensitivity; same L/S/R predictions and filter')
    return jobs


def protocol():
    p = dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
        purpose='Improve the partial-weld result without hiding coverage loss or stronger off-flaw indications.',
        status='Exploratory development on a known simulated defect; not blind validation.',
        fixed='S|SS, vy, nominal 0.5-5 MHz 657-tap filter, identical 15 mm by 8.6 mm field, 0.1 mm grid, RMS sources, CF exponent 1, all declared physical gates except explicit ablations.',
        no_flaw_input='No defect coordinates, flaw-specific delays, pixel-specific parameter selection, or reference subtraction in LOF reconstructions.',
        evaluation='Near-flaw peak / strongest background more than 1 mm from the segment; near support within 0.5 mm. Also pristine/difference controls, full-field maximum location, common support and along-segment extent.',
        candidate_rule='Prefer positive improvements with at least 90% near-flaw support and a global maximum within 0.5 mm of the flaw. Report one-source-only coverage separately. Confirm selected candidates with wider/finer fields and physical source/receiver removal.',
        requests=plan(), dataset_fingerprint=s.dataset('partial')['input_fingerprint'],
        native_sources={name:s.sha(_native_root()/name) for name in ('saft.cpp','saft.hpp','screening.cpp','dataset_geometry.cpp','groove_paths.cpp')})
    target = HERE / 'protocol.json'
    if target.exists():
        old = json.loads(target.read_text())
        assert old['requests'] == p['requests']
        return old
    s.write(target, p)
    return p


def execute(jobs, stage):
    dest = OUT / stage
    pending = [dict(id=j['id'], selection=j['selection']) for j in jobs
               if not (dest/j['id']/'arrays.npz').exists()]
    if pending:
        print(f'{stage}: {len(pending)} configurations, three signal cases each', flush=True)
        s.execute('partial', dict(requests=pending), dest)
        for job in pending:
            s.compact(dest/job['id'])
    summarise(jobs, stage)


def summarise(jobs, stage):
    rows = []
    baseline = OUT/'survey/all__raw__min2/arrays.npz'
    with np.load(baseline) as b:
        for job in jobs:
            with np.load(OUT/stage/job['id']/'arrays.npz') as a:
                points = a['points'].reshape(-1, 2)
                matching = points.shape == b['points'].reshape(-1,2).shape and np.allclose(points,b['points'].reshape(-1,2),rtol=0,atol=1e-12)
                for case in ('lof','reference','difference'):
                    values = a[case+'_2_amplitude']
                    metrics = s.metrics(values, points, 'partial')
                    row = dict(id=job['id'], stage=stage, case=case, selection=job['selection'], **metrics)
                    if matching:
                        mask = np.isfinite(b[case+'_2_amplitude']) & np.isfinite(values)
                        row['common_candidate'] = s.metrics(values,points,'partial',mask)
                        row['common_baseline'] = s.metrics(b[case+'_2_amplitude'],points,'partial',mask)
                    near = s.distances(points,'partial') <= .5
                    row['near_one_source_pct'] = float(100*np.mean(np.isfinite(values[near]) & (a[case+'_2_count'][near]==1)))
                    rows.append(row)
    s.write(OUT/stage/'metrics.json',dict(rows=rows))
    fields=['id','case','cmax_db','c99_db','near_support_pct','near_one_source_pct','support_pct']
    with (OUT/stage/'metrics.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields+['peak_distance_mm','near_um_s','background_um_s'])
        w.writeheader()
        for row in rows:
            w.writerow({k:row[k] for k in fields} | dict(
                peak_distance_mm=row['peak']['distance_mm'] if row['peak'] else None,
                near_um_s=row['near']['um_s'] if row['near'] else None,
                background_um_s=row['background']['um_s'] if row['background'] else None))
    for r in sorted((r for r in rows if r['case']=='lof'),key=lambda r:r['cmax_db'] if r['cmax_db'] is not None else -999,reverse=True):
        print(r['id'], 'Cmax',round(r['cmax_db'],3) if r['cmax_db'] is not None else None,
              'support',round(r['near_support_pct'],1),'one source',round(r['near_one_source_pct'],1),
              'peak distance',round(r['peak']['distance_mm'],3) if r['peak'] else None,flush=True)


def verify():
    checks=[]
    for name, archive, member in [
        ('all__raw__min2','comsol_chapter_evidence.zip','data/confirmation_partial_baseline__1-2_SSS__fine.npz'),
        ('all__weighted__min2','laser_weld_feasibility_sources.zip','data/partial__1-2_SSS__weighted.npz')]:
        folder = 'comsol_chapter_20261101' if archive == 'comsol_chapter_evidence.zip' else 'laser_weld_feasibility_20261030'
        old = HERE.parent / folder / member
        payload = old.read_bytes()
        with np.load(OUT/'survey'/name/'arrays.npz') as a,np.load(io.BytesIO(payload)) as b:
            worst=0
            xx,zz=np.meshgrid(b['x_mm'],b['depth_mm'])
            np.testing.assert_allclose(a['points'].reshape(-1,2)*1000,np.c_[xx.ravel(),zz.ravel()],rtol=0,atol=1e-9)
            for case in ('lof','reference','difference'):
                want=b[case+'_amplitude'].ravel() if case+'_amplitude' in b else b[case].ravel()
                actual=a[case+'_2_amplitude']
                np.testing.assert_allclose(actual,want,rtol=2e-8,atol=2e-15,equal_nan=True)
                worst=max(worst,float(np.nanmax(abs(actual-want))))
                if case+'_count' in b:
                    np.testing.assert_array_equal(a[case+'_2_count'],b[case+'_count'].ravel())
                if case+'_means' in b:
                    np.testing.assert_allclose(a[case+'_2_means'],b[case+'_means'].ravel(),rtol=2e-8,atol=2e-15)
                    np.testing.assert_array_equal(a[case+'_2_valid'],b[case+'_valid'].ravel())
            checks.append(dict(id=name,status='passed',maximum_absolute_difference=worst,archive=str(old.relative_to(ROOT)),member=member,archive_sha256=s.sha(old)))
    s.write(OUT/'verification.json',dict(status='passed',baseline_reproduction=checks,
        executable_sha256=s.sha(s.EXE),driver_sha256=s.sha(Path(__file__)),protocol_sha256=s.sha(HERE/'protocol.json')))
    print('Both archived baselines reproduced in all three signal cases.',flush=True)


def frequency_plan():
    anchors = [j for j in plan() if j['id'] in
               ('all__raw__min2','all__weighted__min2','weld_sides__weighted__min1')]
    jobs=[]
    for lo,hi in [(.5,3),(1,3),(1,5),(1.5,5),(2,5)]:
        for anchor in anchors:
            jobs.append(dict(id=f'{anchor["id"]}__band{lo:g}_{hi:g}',
                             selection=dict(anchor['selection'],low_MHz=lo,high_MHz=hi)))
    p=dict(reason='After the first settings comparison: separate lower-frequency structural clutter from higher-frequency flaw responses, within the existing 5 MHz design range.',
           status='Exploratory second stage after viewing the first 27 outcomes.',
           controls='Same 657-tap kernel, 3.28 us support, CF exponent 1, LSR screen and field. Native source directivity curves and packet windows are recomputed consistently for each declared band. This is a complete band-setting comparison, not a filter-only ablation.',
           assessment='Report absolute amplitudes as well as contrast; a small high-frequency signal is not an experimentally validated detection.',
           requests=jobs)
    target=HERE/'frequency_protocol.json'
    if target.exists():assert json.loads(target.read_text())==p
    else:s.write(target,p)
    return jobs


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('command',choices=['prepare','run','verify','confirm','frequency'])
    args=parser.parse_args()
    if args.command=='prepare':
        p=protocol();print('Frozen',len(p['requests']),'configurations',flush=True)
    elif args.command=='run':execute(protocol()['requests'],'survey')
    elif args.command=='verify':verify()
    elif args.command=='frequency':execute(frequency_plan(),'frequency')
    elif args.command=='confirm':
        p=json.loads((HERE/'confirmation_protocol.json').read_text())
        execute(p['requests'],'confirmation')
