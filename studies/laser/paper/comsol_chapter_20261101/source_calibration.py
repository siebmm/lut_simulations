"""Independently verify the retained thermoelastic calibration exports."""
from pathlib import Path
import hashlib
import json
import re
import sys
import numpy as np
from scipy.integrate import quad
from scipy.signal import butter, sosfiltfilt

HERE = Path(__file__).resolve().parent
DATA = HERE / 'data/source_calibration'
PARAMS = dict(E_Pa=212442565550.3707, nu=0.2889775930054178,
              rho_kg_m3=7859.97504274123, k_W_m_K=52., cp_J_kg_K=480.,
              alpha_K_inv=11e-6, ambient_K=293.15, w0_m=0.0005,
              pulse_fwhm_s=6.5e-9, pulse_centre_s=30e-9,
              production_fluence_J_m2=100.,
              production_multiplier=1.0844342935972646)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()






def verify():
    evidence = json.loads((DATA / 'provenance.json').read_text())
    p = evidence['parameters']
    reports = {}
    sigma = p['pulse_fwhm_s'] / (2 * np.sqrt(2 * np.log(2)))
    for label, run in evidence['runs'].items():
        fluence = float(label)
        path = DATA / run['array_file']
        assert sha(path) == run['array_sha256']
        with np.load(path) as f:
            t, xy = f['time_s'], f['arc_xy_m']
            dt = np.diff(t)
            np.testing.assert_allclose(dt, 2e-9, rtol=1e-9, atol=1e-16)
            assert len(t) == 451 and xy.shape == (7, 2)
            assert f['surface_xy_m'].shape == (5, 2)
            sos = butter(4, 5e6, fs=1/dt[0], output='sos')
            thermal = sosfiltfilt(sos, f['thermal_arc'][:, :, 2:4], axis=1)
            dipole = sosfiltfilt(sos, f['dipole_arc'][:, :, 2:4], axis=1)
            theta = np.rad2deg(np.arctan2(xy[:, 0], -xy[:, 1]))
            fit = (theta >= 14) & (theta <= 76)
            window = (t >= 0.12e-6) & (t <= 0.82e-6)
            a, b = dipole[fit][:, window].ravel(), thermal[fit][:, window].ravel()
            gain = float(np.dot(a, b) / np.dot(a, a))
            residual = float(np.linalg.norm(gain*a-b) / np.linalg.norm(b))
            per_angle = [float(np.linalg.norm(gain*d[window]-v[window]) /
                               np.linalg.norm(v[window])) for d, v in zip(dipole, thermal)]
            derivative = sosfiltfilt(sos, np.gradient(f['thermal_arc'][:, :, :2], dt[0], axis=1), axis=1)
            derivative_error = float(np.linalg.norm(derivative[:, 10:-10]-thermal[:, 10:-10]) /
                                     np.linalg.norm(thermal[:, 10:-10]))
            def surface_temperature(tt):
                if tt <= 0:
                    return 0.
                integrand = lambda z: 2*fluence*np.exp(-0.5*((tt-z*z-p['pulse_centre_s'])/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
                integral = quad(integrand, 0, np.sqrt(tt),
                                points=[np.sqrt(max(tt-p['pulse_centre_s'], 0))], epsabs=1e-6)[0]
                return integral / np.sqrt(np.pi*p['k_W_m_K']*p['rho_kg_m3']*p['cp_J_kg_K'])
            analytic = np.array([surface_temperature(tt) for tt in t])
            centre = int(np.argmin(f['surface_xy_m'][:, 0]))
            thermal_error = float(np.linalg.norm(f['thermal_surface'][centre, :, 4]-p['ambient_K']-analytic) / np.linalg.norm(analytic))
            energy_expected = fluence*p['w0_m']*np.sqrt(np.pi/2)
            np.testing.assert_allclose(f['energy'][:, 0], t, rtol=0, atol=1e-15)
            energy_ratio = float(f['energy'][-1, 1]/energy_expected)
            np.testing.assert_allclose(energy_ratio, f['energy'][-1, 2], rtol=1e-12)
            report = dict(least_squares_dipole_scale=gain, fit_relative_L2_error=residual,
                per_angle_deg=theta.tolist(), per_angle_relative_L2_error=per_angle,
                peak_temperature_rise_K=float(f['thermal_surface'][:, :, 4].max()-p['ambient_K']),
                thermal_1D_analytic_relative_L2_error=thermal_error,
                displacement_derivative_vs_velocity_relative_L2_error=derivative_error,
                final_full_width_thermal_energy_J_per_m=float(f['energy'][-1, 1]),
                final_thermal_energy_relative_error=abs(energy_ratio-1))
            for key, actual in report.items():
                np.testing.assert_allclose(actual, run['historical_statistics'][key], rtol=1e-8, atol=1e-12, err_msg=key)
            reports[label] = report
    a0 = p['E_Pa']*p['alpha_K_inv']*p['production_fluence_J_m2'] / ((1-p['nu'])*p['rho_kg_m3']*p['cp_J_kg_K']*p['w0_m'])
    a = a0*p['production_multiplier']
    report = {'runs': reports, 'production_A0_Pa': a0, 'production_A_Pa': a,
              'production_peak_lobe_traction_Pa': 2*np.exp(-0.5)*a,
              'status': 'All fit, withheld-angle, temperature, energy and velocity checks reproduce the historical exports.'}
    (HERE / 'verification').mkdir(parents=True, exist_ok=True)
    (HERE / 'verification/source_calibration_check.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    verify()
