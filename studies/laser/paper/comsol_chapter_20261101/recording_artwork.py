"""Recording figures with data-free legend rows.

Default: regenerate from compact archived plotting data. --collect refreshes
those data from the original acquisition CSVs in the NDT repository. Filtering,
normalization and arrival hypotheses match the original recording figures.
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
import argparse, hashlib, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib import patheffects as pe
from figure_layout import legend_strip
from make_artifacts import save

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DATA = HERE / 'data/recordings'
COLORS = ['#0072b2', '#d55e00', '#009e73', '#aa4499']
STYLES = ['--', ':', '-.', '-']
PATHS = ['S | S', 'L | L', 'S | SS', 'SS | S']


def collect():
    sys.path.insert(0, str(ROOT / 'studies/laser/src'))
    sys.path.insert(0, str(HERE.parent / 'weld_geometry_study_20261030'))
    from saft_imaging import kernel, filter_record
    from saft_path_physics import branch
    from geometry_checks import traces
    DATA.mkdir(parents=True, exist_ok=True)
    fingerprints = {}
    check = {}
    for kind, source, receiver in [('full', 'shoulder_left', 94), ('partial', 'bead_left', 330)]:
        directory = _data_root() / (kind + '_weld')
        d = json.loads((directory / 'dataset.json').read_text())
        archived = json.loads((HERE / 'data' / (kind + '_dataset.json')).read_text())
        assert d == archived
        receivers = d['audit']['receivers']
        t = np.arange(5001) * .01
        x = np.array([r['x_m'] * 1000 for r in receivers])
        groups = np.array([r['group'] for r in receivers])
        raw = {}
        for case in (['lof'] if kind == 'full' else ['lof', 'reference']):
            raw[case] = np.empty((len(receivers), len(t)))
            for group in ['left', 'weld', 'right']:
                p = directory / 'sources' / source / case / 'receivers' / (group + '_vy.csv')
                with p.open() as stream:
                    ids = [int(v[1:]) for v in stream.readline().strip().split(',')[1:]]
                raw[case][ids] = np.loadtxt(p, delimiter=',', skiprows=1)[:, 1:].T
                fingerprints[p.relative_to(ROOT).as_posix()] = hashlib.sha256(p.read_bytes()).hexdigest()
        h = kernel(.5, 5)
        filt = {case: filter_record(values, h) for case, values in raw.items()}
        ref = abs(filt['lof']).max()
        spectrum = abs(np.fft.rfft(raw['lof'][receiver] * np.hanning(len(t)), 8192))
        values = dict(time_us=t, receiver_x_mm=x, receiver_groups=groups,
                      envelope_db=20*np.log10(np.maximum(abs(filt['lof']) / ref, 1e-12)),
                      bscan_reference_um_s=ref*1e6, receiver=receiver,
                      frequency_MHz=np.fft.rfftfreq(8192, 1e-8)/1e6,
                      spectrum_db=20*np.log10(np.maximum(spectrum / spectrum.max(), 1e-8 if kind == 'full' else 1e-9)))
        if kind == 'partial':
            filt['difference'] = filt['lof'] - filt['reference']
        for case, array in filt.items():
            values[case + '_trace_um_s'] = array[receiver] * 1e6
        if kind == 'full':
            mat = d['material']; height = mat['thickness_m']
            sx = d['audit']['sources'][0]['x_mm']
            tips = np.array(d['audit']['defect_benchmark_only']['tips_xy_mm'])
            tips[:, 1] = height*1000 - tips[:, 1]
            point = tips.mean(axis=0) / 1000
            paths = []
            for before, after in [('S', 'S'), ('L', 'L'), ('S', 'SS'), ('SS', 'S')]:
                inc = branch([sx/1000, 0], point, before, height, mat['cl_m_s'], mat['cs_m_s'])
                rec = branch(point, np.column_stack([x/1000, np.zeros(len(x))]), after,
                             height, mat['cl_m_s'], mat['cs_m_s'], outgoing=True)
                paths.append(.03 + (inc['time'] + rec['time']) * 1e6)
            values['pixel_paths_us'] = np.array(paths)
            # Preserve the original trace's independently evaluated arrival values.
            old = json.loads((ROOT / 'studies/laser/results/full_weld_common_settings_20261030/figure_evidence.json').read_text())
            values['trace_paths_us'] = np.array(old['ray_times_us'])
            values['direct_paths_us'] = np.array([.03 + abs(x-sx)*1000/mat[key]
                                                  for key in ['cl_m_s', 'cs_m_s', 'cr_m_s']])
            trace = np.genfromtxt(ROOT / 'studies/laser/results/full_weld_common_settings_20261030/numerical/0/trace.csv', delimiter=',', names=True)
            native = trace['filtered_real_m_s'] + 1j*trace['filtered_imag_m_s']
            error = float(np.max(abs(filt['lof'][receiver] - native)))
            assert error < 1e-14
            check[kind] = {'filter_max_error_against_native_m_s': error}
        else:
            modes = ['LL', 'SS', 'SSSS']
            values['return_modes'] = np.array(modes)
            arrivals = [[.03 + q['time_us'] for q in traces([-1.5,16.4], [x[receiver],16.4], m, kind)] for m in modes]
            assert all(len(v) == 1 for v in arrivals)
            values['return_times_us'] = np.array(arrivals).ravel()
            old = json.loads((ROOT / 'studies/laser/results/weld_geometry_study_20261030/figure_evidence.json').read_text())['signals']
            for mode, times in zip(modes, arrivals):
                np.testing.assert_allclose(times, old['predicted_specimen_return_times_us'][mode], rtol=0, atol=1e-10)
            check[kind] = {'return_times_match_original': True}
        np.testing.assert_allclose(ref*1e6, old['bscan_reference_um_s'], rtol=1e-12)
        check[kind]['bscan_reference_um_s'] = float(ref*1e6)
        np.savez_compressed(DATA / (kind + '.npz'), **values)
    provenance = dict(description='Plotting arrays regenerated from the original recording CSVs with unchanged 0.5-5 MHz kernel, spectrum window and arrival hypotheses.',
                      inputs_sha256=fingerprints, checks=check,
                      arrays_sha256={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in DATA.glob('*.npz')})
    (DATA / 'provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    (HERE / 'verification/recording_figure_check.json').write_text(json.dumps({'status': 'passed', **provenance}, indent=2)+'\n')


def bscan(ax, values, title):
    x, time = values['receiver_x_mm'], values['time_us']
    for group in ['left', 'weld', 'right']:
        ids = np.flatnonzero(values['receiver_groups'] == group)
        im = ax.pcolormesh(x[ids], time, values['envelope_db'][ids].T,
                          cmap='magma', vmin=-60, vmax=0, shading='nearest', rasterized=True)
    ax.set(xlim=(x.min(), x.max()), ylim=(30,0), xlabel='Receiver $x$ (mm)',
           ylabel=r'Time ($\mu$s)', title=title)
    return im


def spectrum(ax, values, title):
    ax.plot(values['frequency_MHz'], values['spectrum_db'], c='#555555', lw=.7)
    ax.axvspan(.5, 5, color='#dddddd', alpha=.55, zorder=-1)
    ax.set(xlim=(0,8), ylim=(-100,0), xlabel='Frequency (MHz)',
           ylabel='dB below spectral peak', title=title)


def full():
    v = np.load(DATA / 'full.npz'); t = v['time_us']
    fig = plt.figure(figsize=(7,6.5), layout='constrained')
    gs = fig.add_gridspec(5,2, height_ratios=[1,.12,1.55,.12,.12], width_ratios=[1.6,1])
    ax = fig.add_subplot(gs[0,0])
    ax.plot(t, v['lof_trace_um_s'].real, c='#222222', lw=.65, label=r'Filtered signed $v_y$')
    ax.set(xlim=(0,30), xlabel=r'Time ($\mu$s)', ylabel=r'$v_y$ ($\mu$m/s)',
           title=r'(a) S0 $\rightarrow$ R94 ($x_r=-20.05$ mm)')
    for time, col, ls in zip(v['trace_paths_us'], COLORS, STYLES):
        ax.axvline(time, c=col, ls=ls, lw=.9)
    legend_strip(fig, gs[1,0], *ax.get_legend_handles_labels(), fontsize=7)
    spectrum(fig.add_subplot(gs[0,1]), v, '(b) Raw-record Hann spectrum')
    legend_strip(fig, gs[1,1], [Patch(fc='#dddddd', label='0.5-5 MHz band')], fontsize=7)
    ax = fig.add_subplot(gs[2,:])
    im = bscan(ax, v, '(c) Filtered B-scan; pixel paths focus at the LOF centre')
    path_handles = []
    for times, col, ls, label in zip(v['pixel_paths_us'], COLORS, STYLES, PATHS):
        line, = ax.plot(v['receiver_x_mm'], times, c=col, ls=ls, lw=1.2, label=label)
        path_handles.append(line)
    direct_handles = []
    for times, ls, mode in zip(v['direct_paths_us'], [':','--','-'], ['L','S','Rayleigh']):
        ax.plot(v['receiver_x_mm'], times, c='white', ls=ls, lw=.65, alpha=.8)
        # A thin dark halo keeps the white key legible on the paper background.
        direct_handles.append(Line2D([], [], c='white', ls=ls, lw=.9,
            path_effects=[pe.Stroke(linewidth=2.1, foreground='#777777'), pe.Normal()],
            label='Direct '+mode))
    legend_strip(fig, gs[3,:], path_handles, ncol=4, fontsize=7, title='Pixel paths in (a) and (c)', title_fontsize=7)
    legend_strip(fig, gs[4,:], direct_handles, ncol=3, fontsize=7, title='White direct-arrival guides in (c)', title_fontsize=7)
    ref = float(v['bscan_reference_um_s'])
    fig.colorbar(im, ax=ax, pad=.02, label=f'Envelope (dB)\n0 dB = {ref:.1f} '+r'$\mu$m/s', ticks=[-60,-40,-20,0])
    save(fig, 'full_recordings')


def partial():
    v = np.load(DATA / 'partial.npz'); t = v['time_us']
    fig = plt.figure(figsize=(7,7.2), layout='constrained')
    gs = fig.add_gridspec(5,2, height_ratios=[1,.23,1,.15,1.35])
    ax = fig.add_subplot(gs[0,:])
    cases = [('lof','#111111','-','LOF'), ('reference',COLORS[1],'--','Pristine'),
             ('difference',COLORS[0],':','Difference')]
    for case, col, ls, label in cases[:2]:
        ax.plot(t, v[case+'_trace_um_s'].real, c=col, ls=ls, lw=.7, label=label)
    ax.axvspan(0,3.28,color='#dddddd')
    ax.set(xlim=(0,18), xlabel=r'Time ($\mu$s)', ylabel=r'$v_y$ ($\mu$m/s)',
           title=r'(a) S1 at $x=-1.5$ mm $\rightarrow$ R330 at $x=1.94$ mm; both on weld surface')
    handles, labels = ax.get_legend_handles_labels()
    for mode, time, col in zip(v['return_modes'], v['return_times_us'], [COLORS[0],COLORS[2],COLORS[3]]):
        line = ax.axvline(time,c=col,ls=':',lw=.9)
        handles.append(line); labels.append(f'{mode}: {time:.3f} '+r'$\mu$s')
    legend_strip(fig, gs[1,:], handles, labels, ncol=5, fontsize=7)
    ax = fig.add_subplot(gs[2,0])
    for case, col, ls, label in cases:
        ax.plot(t, abs(v[case+'_trace_um_s']), c=col, ls=ls, lw=.8, label=label)
    ax.axvspan(0,3.28,color='#dddddd')
    ax.set(xlim=(0,18), xlabel=r'Time ($\mu$s)', ylabel=r'Envelope ($\mu$m/s)',title='(b) Waveform controls')
    legend_strip(fig, gs[3,0], *ax.get_legend_handles_labels(), ncol=3, fontsize=7, columnspacing=1)
    spectrum(fig.add_subplot(gs[2,1]), v, '(c) Raw-record Hann spectrum')
    legend_strip(fig, gs[3,1], [Patch(fc='#dddddd', label='0.5-5 MHz band')], fontsize=7)
    ax = fig.add_subplot(gs[4,:])
    im = bscan(ax, v, '(d) S1 B-scan; weld surface is 16.4 mm below side surfaces')
    ax.set_xlim(-44,44)
    fig.colorbar(im, ax=ax, pad=.02, label='Envelope (dB)', ticks=[-60,-40,-20,0])
    save(fig, 'partial_recordings')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(); parser.add_argument('--collect', action='store_true')
    if parser.parse_args().collect: collect()
    full(); partial()
