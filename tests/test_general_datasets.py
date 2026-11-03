"""Contract and physics regressions for arbitrary supported 2D acquisitions."""
import contextlib
import copy
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import numpy as np
from scipy import optimize, signal

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
sys.path.insert(0, str(ROOT / 'studies/laser/src'))
from dataset_io import load_dataset, read_json, validate_geometry, visible
from make_example import generate
from standalone import executable, native_env


def write(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


@contextlib.contextmanager
def changed(path, value):
    before = path.read_bytes()
    try:
        write(path, value)
        yield
    finally:
        path.write_bytes(before)


class DatasetTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        (ROOT / 'build').mkdir(exist_ok=True)
        cls.temporary = tempfile.TemporaryDirectory(dir=ROOT / 'build')
        cls.work = Path(cls.temporary.name)
        cls.dataset = generate(cls.work / 'dataset')
        cls.acquisition = read_json(cls.dataset / 'acquisition.json')
        cls.geometry = read_json(cls.dataset / 'geometry.json')
        cls.counter = 0

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def selection(self):
        return dict(path='1-1_LL', component='vz', cases=['measurement', 'baseline', 'difference'],
                    low_MHz=.5, high_MHz=5, minimum_sources=1, minimum_receivers=1,
                    minimum_effective_receivers=1, receiver_fraction=0, min_offset_mm=0, threads=2)

    def native(self, request, dataset=None, tool='SAFTBatch', success=True):
        program = executable(tool)
        self.assertTrue(program.is_file(), 'Build all four native tools before running the integration suite')
        type(self).counter += 1
        folder = self.work / f'run{self.counter}'
        folder.mkdir()
        path = folder / 'request.json'
        write(path, request)
        dataset = dataset or self.dataset
        if tool == 'SAFTBatch':
            command = [str(program), '--repository', str(self.work), '--dataset', str(dataset),
                       '--compute', str(path), '--output', str(folder / 'result')]
        else:
            command = [str(program), str(self.work), str(dataset), str(path), str(folder / 'result')]
        result = subprocess.run(command, env=native_env(), capture_output=True, text=True, timeout=180)
        if success:
            self.assertEqual(result.returncode, 0, result.stderr)
            return folder / 'result'
        self.assertNotEqual(result.returncode, 0, 'Invalid input was accepted')
        return result.stderr

    def point_request(self, points=None):
        return dict(requests=[dict(id='check', selection=self.selection(),
                                  points_m=points if points is not None else [[.003, .012]])])

    def test_contract_and_peak(self):
        *_, report = load_dataset(self.dataset)
        self.assertEqual((report['sources'], report['receivers'], report['samples']), (7, 13, 2401))
        selection = self.selection() | dict(x_min_mm=2.5, x_max_mm=3.5, depth_min_mm=11.5,
                                           depth_max_mm=12.5, spacing_mm=.1)
        result = self.native(dict(verify_inputs=True, requests=[dict(selection=selection)]))
        with np.load(result / '0/images.npz') as arrays:
            image = arrays['measurement_maps'][2]
            self.assertTrue(np.isfinite(image).all())
            np.testing.assert_allclose(arrays['points_m'][image.argmax()], [.003, .012], atol=1e-12)
            self.assertEqual(arrays['measurement_receiver_usage'].shape, (3, 7, 13))

    def test_all_four_interfaces(self):
        selection = self.selection() | dict(x_min_mm=2.9, x_max_mm=3.1, depth_min_mm=11.9,
                                           depth_max_mm=12.1, spacing_mm=.1)
        request = dict(requests=[dict(id='check', selection=selection)])
        folders = {name: self.native(request, tool=name) for name in
                   ('SAFTFilterStudy', 'SAFTReview', 'SAFTAdmissionAudit', 'SAFTBatch')}
        with np.load(folders['SAFTBatch'] / '0/images.npz') as batch:
            for case in selection['cases']:
                expected = batch[case+'_maps'][2]
                self.assertTrue(np.isfinite(expected).all())
                for name in ('SAFTFilterStudy', 'SAFTReview', 'SAFTAdmissionAudit'):
                    suffix = '_amplitude.f64' if name == 'SAFTAdmissionAudit' else '_2_amplitude.f64'
                    actual = np.fromfile(folders[name] / 'check' / (case+suffix), dtype='<f8')
                    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=1e-15)

    def test_independent_filter_interpolation_and_focus(self):
        points = np.array([[.003, .012], [.005, .015], [-.010, .012]])
        output = self.native(self.point_request(points.tolist()))
        a = self.acquisition
        dt, start = a['sampling']['dt_s'], a['sampling']['start_s']
        half = a['filter_half_samples']
        n = np.arange(-half, half+1)
        width, centre = 4.5e6, 2.75e6
        base = 2*width*dt*np.sinc(width*n*dt)*np.kaiser(2*half+1, .1102*(55-8.7))
        kernel = base*np.exp(2j*np.pi*centre*n*dt)*2/base.sum()
        polygons, _ = validate_geometry(self.geometry)
        means = np.zeros((7, len(points)), dtype=complex)
        valid = np.zeros(means.shape, dtype=bool)
        for si, source in enumerate(a['sources']):
            values = np.loadtxt(self.dataset / f'histories/{source["id"]}_measurement_vz.csv', delimiter=',', skiprows=1)[:, 1:].T
            filtered = signal.fftconvolve(values, kernel[None, :], axes=-1, mode='full')[:, half:half+2401]
            source_point = np.array(source['position_m'])
            for pi, p in enumerate(points):
                samples = []
                for ri, receiver in enumerate(a['receivers']):
                    r = np.array(receiver['position_m'])
                    if not visible(source_point, p, polygons) or not visible(p, r, polygons):
                        continue
                    t = source['trigger_s']+(np.linalg.norm(p-source_point)+np.linalg.norm(p-r))/6000
                    f = (t-start)/dt
                    j = int(np.floor(f))
                    if j-half < 0 or j+1+half >= 2401:
                        continue
                    samples.append(filtered[ri, j]*(1-f+j)+filtered[ri, j+1]*(f-j))
                if samples:
                    means[si, pi] = np.mean(samples)
                    valid[si, pi] = True
        with np.load(output / '0/images.npz') as arrays:
            np.testing.assert_allclose(arrays['measurement_means'][2], means, rtol=1e-9, atol=1e-15)
            np.testing.assert_array_equal(arrays['measurement_valid'][2], valid)
            expected = np.sqrt(np.sum(abs(means)**2, axis=0)/np.maximum(valid.sum(axis=0), 1))
            expected[valid.sum(axis=0) == 0] = np.nan
            np.testing.assert_allclose(arrays['measurement_maps'][2], expected, rtol=1e-9, atol=1e-15)
            self.assertTrue(np.isnan(arrays['measurement_maps'][2, -1]), 'Void was reconstructed as material')

    def test_reflected_and_converted_travel_times(self):
        a, b = np.array([.010, .002]), np.array([.016, .009])
        left, right = np.array([.023, .027]), np.array([-.020, .022])
        edge = right-left
        checks = [dict(a_m=a.tolist(), b_m=b.tolist(), modes=m, via='bottom') for m in ('LL', 'LS', 'SL', 'SS')]
        request = self.point_request() | dict(ray_checks=checks)
        result = self.native(request)
        rays = read_json(result / 'geometry.json')['rays']
        for modes, ray in zip(('LL', 'LS', 'SL', 'SS'), rays):
            c1, c2 = [6000 if m == 'L' else 3200 for m in modes]
            def derivative(u):
                p = left+u*edge
                return np.dot(p-a, edge)/(np.linalg.norm(p-a)*c1)+np.dot(p-b, edge)/(np.linalg.norm(p-b)*c2)
            u = optimize.brentq(derivative, 0, 1, xtol=1e-14)
            p = left+u*edge
            expected = np.linalg.norm(p-a)/c1+np.linalg.norm(p-b)/c2
            self.assertTrue(ray['valid'])
            self.assertAlmostEqual(ray['time_s'], expected, delta=1e-14)
            np.testing.assert_allclose(ray['points_m'][1], p, atol=1e-11)

    def test_rotated_translated_geometry_and_weighted_components(self):
        acquisition = copy.deepcopy(self.acquisition)
        acquisition['response'] = dict(model='laser_free_surface', line_half_width_m=.0005, pulse_fwhm_s=6.5e-9)
        points = np.array([[.003, .012], [.006, .014]])
        request = self.point_request(points.tolist())
        request['requests'][0]['selection']['component'] = 'vx'
        with changed(self.dataset / 'acquisition.json', acquisition):
            original = self.native(request)
            rotated = self.work / 'rotated'
            shutil.copytree(self.dataset, rotated)
        rotation = np.array([[0., -1.], [1., 0.]])
        offset = np.array([.04, -.08])
        geometry = copy.deepcopy(self.geometry)
        for loop in geometry['loops']:
            loop['vertices_m'] = (np.array(loop['vertices_m']) @ rotation.T+offset).tolist()
        for key in ('sources', 'receivers'):
            for endpoint in acquisition[key]:
                endpoint['position_m'] = (rotation @ endpoint['position_m']+offset).tolist()
                endpoint['inward_normal'] = (rotation @ endpoint['inward_normal']).tolist()
        for source in acquisition['sources']:
            for case in acquisition['cases']:
                vx = rotated / f'histories/{source["id"]}_{case}_vx.csv'
                vz = rotated / f'histories/{source["id"]}_{case}_vz.csv'
                x, z = np.loadtxt(vx, delimiter=',', skiprows=1), np.loadtxt(vz, delimiter=',', skiprows=1)
                header = vx.read_text().splitlines()[0]
                x[:, 1:], z[:, 1:] = -z[:, 1:].copy(), x[:, 1:].copy()
                for path, array in ((vx, x), (vz, z)):
                    np.savetxt(path, array, delimiter=',', fmt='%.17g', header=header, comments='')
        write(rotated / 'geometry.json', geometry)
        write(rotated / 'acquisition.json', acquisition)
        load_dataset(rotated)
        request['requests'][0]['selection']['component'] = 'vz'
        request['requests'][0]['points_m'] = (points @ rotation.T+offset).tolist()
        actual = self.native(request, dataset=rotated)
        with np.load(original / '0/images.npz') as old, np.load(actual / '0/images.npz') as new:
            self.assertTrue(np.isfinite(old['measurement_maps'][2]).all())
            np.testing.assert_allclose(new['measurement_maps'], old['measurement_maps'], rtol=2e-10, atol=1e-14)

    def test_large_receiver_count_and_long_record(self):
        for name, sources, receivers, samples in (('many', 1, 350, 1001), ('long', 2, 3, 9001)):
            dataset = generate(self.work / name, sources, receivers, samples)
            result = self.native(self.point_request(), dataset=dataset)
            with np.load(result / '0/images.npz') as arrays:
                self.assertEqual(arrays['measurement_receiver_usage'].shape, (3, sources, receivers))
                self.assertTrue(np.isfinite(arrays['measurement_maps']).all())

    def test_sampling_in_trace_export_and_pulse_screening(self):
        request = self.point_request()
        job = request['requests'][0]
        job['selection']['pulse_screening'] = 'LS'
        job['inspect'] = dict(source=0, receiver=0, signal='baseline', has_point=True,
                              x_mm=3, depth_mm=12, competitor='1-1_SS')
        result = self.native(request)
        trace = np.genfromtxt(result / '0/trace.csv', delimiter=',', names=True)
        np.testing.assert_allclose(trace['time_us'], -2+np.arange(2401)*.02, atol=1e-12)
        windows = read_json(result / '0/pulse_windows.json')['baseline']['windows']
        self.assertEqual({w['source'] for w in windows}, set(range(7)))
        first = next(w for w in windows if w['source'] == 0 and w['receiver'] == 0 and w['family'] == 'L surface')
        self.assertAlmostEqual(first['peak_time_us'], .5, delta=.021)
        geometry = copy.deepcopy(self.geometry)
        geometry.pop('screening_boundaries')
        with changed(self.dataset / 'geometry.json', geometry):
            job['selection']['pulse_reflections'] = [True, False, False, False]
            self.assertIn('screening_boundaries', self.native(request, success=False))

    def test_external_processed_bank_uses_dataset_dimensions(self):
        selection = self.selection() | dict(x_min_mm=2.9, x_max_mm=3.1, depth_min_mm=11.9,
                                           depth_max_mm=12.1, spacing_mm=.1)
        request = dict(requests=[dict(id='check', selection=selection)])
        baseline = self.native(request)
        fingerprint = read_json(baseline / 'benchmark.json')['provenance']['input_fingerprint']
        bank = self.work / 'bank'
        bank.mkdir()
        write(bank / 'bank.json', dict(dataset_fingerprint=fingerprint))
        half, dt = 128, 20e-9
        n = np.arange(-half, half+1)
        base = 2*4.5e6*dt*np.sinc(4.5e6*n*dt)*np.kaiser(2*half+1, .1102*(55-8.7))
        kernel = base*np.exp(2j*np.pi*2.75e6*n*dt)*2/base.sum()
        for si in range(7):
            for case in ('measurement', 'baseline'):
                values = np.loadtxt(self.dataset / f'histories/shot{si}_{case}_vz.csv', delimiter=',', skiprows=1)[:, 1:].T
                processed = signal.fftconvolve(values, kernel[None, :], axes=-1, mode='full')[:, half:half+2401]
                processed.astype('<c16').tofile(bank / f'{case}_{si}.c128')
        request['requests'][0]['signal_processing'] = dict(bank=str(bank), support_half_samples=half)
        actual = self.native(request, tool='SAFTFilterStudy')
        with np.load(baseline / '0/images.npz') as arrays:
            for case in selection['cases']:
                values = np.fromfile(actual / f'check/{case}_2_amplitude.f64', dtype='<f8')
                np.testing.assert_allclose(values, arrays[case+'_maps'][2], rtol=1e-9, atol=1e-15)

    def test_bad_time_header_and_nonfinite_samples_rejected(self):
        path = self.dataset / self.acquisition['recordings'][0]['path']
        original = path.read_bytes()
        for kind in ('time', 'header', 'nan', 'extra'):
            try:
                lines = original.decode().splitlines()
                if kind == 'header':
                    lines[0] = lines[0].replace('sensor0', 'wrong')
                elif kind == 'extra':
                    lines.append(lines[-1])
                else:
                    row = lines[1].split(',')
                    row[0 if kind == 'time' else 1] = '0.001' if kind == 'time' else 'nan'
                    lines[1] = ','.join(row)
                path.write_text('\n'.join(lines)+'\n')
                with self.assertRaises(ValueError):
                    load_dataset(self.dataset)
                request = self.point_request()
                request['requests'][0]['selection']['component'] = 'vx'
                self.native(request, success=False)
            finally:
                path.write_bytes(original)

    def test_missing_record_and_nyquist_rejected(self):
        acquisition = copy.deepcopy(self.acquisition)
        acquisition['recordings'].pop()
        with changed(self.dataset / 'acquisition.json', acquisition):
            with self.assertRaises(ValueError):
                load_dataset(self.dataset)
            self.native(self.point_request(), success=False)
        request = self.point_request()
        request['requests'][0]['selection']['high_MHz'] = 30
        self.assertIn('Nyquist', self.native(request, success=False))

    def test_invalid_geometry_and_paths_rejected(self):
        for transform in ('winding', 'crossing', 'outside_hole', 'unknown_boundary'):
            geometry = copy.deepcopy(self.geometry)
            if transform == 'winding':
                geometry['loops'][0]['vertices_m'].reverse()
            elif transform == 'crossing':
                geometry['loops'][0]['vertices_m'][1:3] = list(reversed(geometry['loops'][0]['vertices_m'][1:3]))
            elif transform == 'outside_hole':
                geometry['loops'][1]['vertices_m'] = (np.array(geometry['loops'][1]['vertices_m'])+.2).tolist()
            else:
                geometry['paths'][0]['source']['boundaries'] = ['missing']
            with self.assertRaises(ValueError):
                validate_geometry(geometry)
            with changed(self.dataset / 'geometry.json', geometry):
                self.native(self.point_request(), success=False)

    def test_no_reference_required(self):
        acquisition = copy.deepcopy(self.acquisition)
        acquisition.pop('difference')
        acquisition['cases'] = ['measurement']
        acquisition['recordings'] = [r for r in acquisition['recordings'] if r['case'] == 'measurement']
        with changed(self.dataset / 'acquisition.json', acquisition):
            load_dataset(self.dataset)
            request = self.point_request()
            request['requests'][0]['selection']['cases'] = ['measurement']
            self.native(request)
            request['requests'][0]['selection']['cases'] = ['difference']
            self.native(request, success=False)


if __name__ == '__main__':
    unittest.main()
