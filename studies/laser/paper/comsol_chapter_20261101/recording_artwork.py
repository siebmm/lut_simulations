"""Recording figures with data-free legend rows.

Regenerate from the retained recording arrays. Filtering,
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
    full(); partial()
