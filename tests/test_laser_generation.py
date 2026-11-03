"""Check local-source analysis against the independent retained COMSOL solutions."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ANALYSIS = ROOT/'models/laser_generation/analyze.py'
spec = importlib.util.spec_from_file_location('laser_calibration_analysis', ANALYSIS)
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


class LaserGenerationTests(unittest.TestCase):
    def test_both_reference_fluences_reproduce_comsol_comparisons(self):
        reference = json.loads((ROOT/'studies/laser/paper/comsol_chapter_20261101/data/source_calibration/provenance.json').read_text())
        fields = ['least_squares_dipole_scale', 'fit_relative_L2_error',
                  'peak_temperature_rise_K', 'thermal_1D_analytic_relative_L2_error',
                  'displacement_derivative_vs_velocity_relative_L2_error',
                  'final_thermal_energy_relative_error']
        with tempfile.TemporaryDirectory() as folder:
            for fluence in ('10', '100'):
                with self.subTest(fluence=fluence):
                    output = Path(folder)/fluence
                    subprocess.run([sys.executable, str(ANALYSIS), '--stored-fluence', fluence,
                                    '--output', str(output)], check=True, capture_output=True)
                    actual = json.loads((output/'calibration_comparison.json').read_text())
                    expected = reference['runs'][fluence]['historical_statistics']
                    for field in fields:
                        np.testing.assert_allclose(actual[field], expected[field], rtol=1e-8, atol=1e-12, err_msg=field)
                    self.assertEqual(actual['output_samples'], 451)
                    self.assertTrue((output/'calibration_waveform_comparison.png').is_file())
                    self.assertTrue((output/'calibration_temperature.png').is_file())

    def test_comsol_export_time_and_component_order(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'trace.txt'
            path.write_text('% x y u @ t=0 v @ t=0 u @ t=2e-9 v @ t=2e-9\n'
                            '0.001 -0.002 10 20 11 21\n'
                            '0.002 0 30 40 31 41\n')
            time, xy, values = analysis.read_export(path, 2)
            np.testing.assert_array_equal(time, [0, 2e-9])
            np.testing.assert_array_equal(xy, [[0.001, -0.002], [0.002, 0]])
            np.testing.assert_array_equal(values, [[[10, 20], [11, 21]], [[30, 40], [31, 41]]])
            path.write_text(path.read_text().replace('v @ t=0', 'v @ t=1e-9'))
            with self.assertRaisesRegex(AssertionError, 'grouped by output time'):
                analysis.read_export(path, 2)

    def test_incomplete_simulation_cannot_be_reported_as_a_calibration(self):
        with tempfile.TemporaryDirectory() as folder:
            exports = Path(folder)/'exports'
            exports.mkdir()
            (exports/'calibration_parameters.json').write_text('{"complete":false}')
            result = subprocess.run([sys.executable, str(ANALYSIS), '--exports', str(exports),
                                     '--output', str(Path(folder)/'results')], capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn('Finish both source solves', result.stderr)
            self.assertFalse((Path(folder)/'results/calibration_comparison.json').exists())


if __name__ == '__main__':
    unittest.main()
