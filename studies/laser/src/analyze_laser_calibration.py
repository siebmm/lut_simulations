"""Validate COMSOL receiver exports and compare the 0-5 MHz outgoing field.

Run after RunLaserCalibration.class; uses NumPy, SciPy, and Matplotlib.
"""
from study_paths import STUDY
import json
import re
from pathlib import Path
import numpy as np
from scipy.signal import butter, sosfiltfilt
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_export(filename, nvar):
    path = Path(filename)
    lines = path.read_text(encoding='utf-8-sig').splitlines()
    header = next(line for line in lines if line.startswith('%') and '@ t=' in line)
    column_times = np.array([float(t) for t in re.findall(r'@\s*t=([+\-\d.eE]+)', header)])
    a = np.loadtxt(path, comments='%')
    a = np.atleast_2d(a)
    assert a.shape[1] == 2 + len(column_times), (filename, a.shape, len(column_times))
    assert len(column_times) % nvar == 0
    t = column_times.reshape(-1, nvar)
    assert np.all(t == t[:, :1]), 'Expected expressions grouped by output time'
    assert np.isfinite(a).all(), f'Nonfinite data in {filename}'
    assert len(np.unique(a[:, :2], axis=0)) == len(a), 'Duplicate receiver coordinates'
    return t[:, 0], a[:, :2], a[:, 2:].reshape(len(a), -1, nvar)


def main():
    parameters_path = STUDY / 'results/thermoelastic_source_calibration/exports/calibration_parameters.json'
    if parameters_path.exists():
        parameters = json.loads(parameters_path.read_text(encoding='utf-8'))
        if not parameters.get('complete'):
            raise ValueError('Calibration run is incomplete; finish both solves before comparing exports')
        fluence = float(parameters['absorbed_fluence_J_m2'])
    else:
        # Preserved exports predate run metadata and were solved at 10 J/m2.
        fluence = 10.0
    if not np.isfinite(fluence) or fluence <= 0:
        raise ValueError('Absorbed fluence must be finite and positive')
    t, xy, thermal = read_export(STUDY / 'results/thermoelastic_source_calibration/exports/calibration_thermal_arc.txt', 4)
    td, xyd, dipole = read_export(STUDY / 'results/thermoelastic_source_calibration/exports/calibration_dipole_arc.txt', 4)
    ts, xys, surface = read_export(STUDY / 'results/thermoelastic_source_calibration/exports/calibration_thermal_surface.txt', 5)
    tsd, xysd, surface_dipole = read_export(STUDY / 'results/thermoelastic_source_calibration/exports/calibration_dipole_surface.txt', 4)
    assert np.allclose(t, td, atol=1e-15) and np.allclose(xy, xyd, atol=1e-12)
    assert np.allclose(t, ts, atol=1e-15) and np.allclose(t, tsd, atol=1e-15)
    assert np.allclose(xys, xysd, atol=1e-12)
    assert len(xy) == 7 and len(xys) == 5
    dt = np.diff(t)
    assert np.allclose(dt, dt[0], rtol=1e-9, atol=1e-16)
    sos = butter(4, 5e6, fs=1/dt[0], output='sos')
    tv = sosfiltfilt(sos, thermal[:, :, 2:4], axis=1)
    dv = sosfiltfilt(sos, dipole[:, :, 2:4], axis=1)
    derivative_velocity = sosfiltfilt(sos, np.gradient(thermal[:, :, :2], dt[0], axis=1), axis=1)
    derivative_error = np.linalg.norm(derivative_velocity[:, 10:-10]-tv[:, 10:-10])/np.linalg.norm(tv[:, 10:-10])
    sv = sosfiltfilt(sos, surface[:, :, 3], axis=1)
    sdv = sosfiltfilt(sos, surface_dipole[:, :, 3], axis=1)
    # Exclude epicentre (dipole weak/null longitudinal response) and surface
    # from fitting, then report their errors independently.
    theta = np.rad2deg(np.arctan2(xy[:, 0], -xy[:, 1]))
    fit_receivers = (theta >= 14) & (theta <= 76)
    fit_times = (t >= 0.12e-6) & (t <= 0.82e-6)
    a = dv[fit_receivers][:, fit_times, :].ravel()
    b = tv[fit_receivers][:, fit_times, :].ravel()
    assert np.linalg.norm(a) > 0 and np.linalg.norm(b) > 0
    scale = float(np.dot(a, b)/np.dot(a, a))
    fit_error = float(np.linalg.norm(scale*a-b)/np.linalg.norm(b))
    errors = []
    for i in range(len(xy)):
        ref = tv[i, fit_times, :]
        err = np.linalg.norm(scale*dv[i, fit_times, :]-ref)/np.linalg.norm(ref)
        errors.append(float(err))
    temperature_max = float(np.max(surface[:, :, 4]))
    # Independent 1D semi-infinite heat-conduction check at the beam centre.
    # sqrt(t-s) substitution removes the Green-function endpoint singularity.
    sigma = 6.5e-9/(2*np.sqrt(2*np.log(2)))
    def expected_temperature(tt):
        if tt <= 0:
            return 0.
        f = lambda z: 2*fluence*np.exp(-.5*((tt-z*z-30e-9)/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
        split = np.sqrt(max(tt-30e-9, 0))
        val = quad(f, 0, np.sqrt(tt), points=[split], epsabs=1e-6)[0]
        return val/np.sqrt(np.pi*52*7859.97504274123*480)
    analytic_temperature = np.array([expected_temperature(tt) for tt in t])
    centre = int(np.argmin(xys[:, 0]))
    thermal_error = np.linalg.norm(surface[centre, :, 4]-293.15-analytic_temperature)/np.linalg.norm(analytic_temperature)
    stats = {
        'absorbed_fluence_J_m2': fluence,
        'output_samples': len(t), 'dt_s': float(dt[0]),
        'bulk_receivers': len(xy), 'surface_receivers': len(xys),
        'peak_surface_temperature_K': temperature_max,
        'peak_temperature_rise_K': temperature_max-293.15,
        'thermal_1D_analytic_relative_L2_error': float(thermal_error),
        'displacement_derivative_vs_velocity_relative_L2_error': float(derivative_error),
        'peak_thermal_arc_displacement_m': float(np.max(np.abs(thermal[:, :, :2]))),
        'peak_thermal_surface_vertical_displacement_m': float(np.max(np.abs(surface[:, :, 1]))),
        'filter': 'Zero-phase fourth-order Butterworth low-pass, 5 MHz cutoff',
        'fit_receivers_deg_from_downward_vertical': theta[fit_receivers].tolist(),
        'fit_window_s': [0.12e-6, 0.82e-6],
        'least_squares_dipole_scale': scale,
        'fit_relative_L2_error': fit_error,
        'per_angle_deg': theta.tolist(),
        'per_angle_relative_L2_error': errors,
        'status': 'Numerical source comparison; not experimental calibration or convergence certification',
    }
    assert temperature_max > 293.151, 'Laser heating not detected'
    assert stats['peak_thermal_arc_displacement_m'] > 1e-16, 'No outgoing wave detected'
    energy_path = Path(STUDY / 'results/thermoelastic_source_calibration/exports/calibration_energy_check.txt')
    if energy_path.exists():
        energy = np.loadtxt(energy_path, comments='%')
        assert energy.shape == (len(t), 3) and np.isfinite(energy).all()
        assert np.allclose(energy[:, 0], t, atol=1e-15)
        stats['final_full_width_thermal_energy_J_per_m'] = float(energy[-1, 1])
        stats['final_thermal_energy_over_absorbed_energy'] = float(energy[-1, 2])
        stats['final_thermal_energy_relative_error'] = float(abs(energy[-1, 2]-1))
    Path(STUDY / 'results/thermoelastic_source_calibration/calibration_comparison.json').write_text(json.dumps(stats, indent=2)+'\n')
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for ax, angle in zip(axes.flat[:5], [15, 30, 45, 60, 75]):
        i = int(np.argmin(abs(theta-angle)))
        ax.plot(t*1e6, tv[i, :, 1]*1e6, label='Thermoelastic')
        ax.plot(t*1e6, scale*dv[i, :, 1]*1e6, '--', label='Fitted dipole')
        ax.set(title=f'{angle} degrees, r = 2 mm', xlabel='Time (us)', ylabel='Vertical velocity (um/s)')
        ax.grid(alpha=.2)
    axes.flat[0].legend(fontsize=8)
    ax = axes.flat[5]
    i = int(np.argmax(xys[:, 0]))
    ax.plot(t*1e6, sv[i]*1e6, label='Thermoelastic')
    ax.plot(t*1e6, scale*sdv[i]*1e6, '--', label='Fitted dipole')
    ax.set(title='Surface receiver x = 2 mm (not fitted)', xlabel='Time (us)', ylabel='Vertical velocity (um/s)')
    ax.grid(alpha=.2)
    fig.suptitle(f'5 MHz low-pass comparison | fitted scale {scale:.4f} | bulk relative error {fit_error:.1%}')
    fig.savefig(STUDY / 'results/thermoelastic_source_calibration/calibration_waveform_comparison.png', dpi=170)
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.plot(t*1e9, surface[centre, :, 4]-293.15, label='COMSOL')
    ax.plot(t*1e9, analytic_temperature, '--', label='1D semi-infinite heat conduction')
    ax.set(xlim=(0, 120), xlabel='Time (ns)', ylabel='Surface temperature rise (K)',
           title=f'Laser spot centre: absorbed fluence {fluence:g} J/m2')
    ax.grid(alpha=.2)
    ax.legend()
    fig.savefig(STUDY / 'results/thermoelastic_source_calibration/calibration_temperature.png', dpi=170)
    print(json.dumps(stats, indent=2))


if __name__ == '__main__':
    main()
