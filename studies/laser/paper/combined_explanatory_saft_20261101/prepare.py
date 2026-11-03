"""Draw the paper's path illustration and collect its explanatory figure panels."""
from pathlib import Path
import hashlib, json, re, shutil, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / 'comsol_chapter_20261101'
EXPL = HERE.parent / 'explanatory_saft_20261101'
REPORTS = HERE.parents[1] / 'results/reports'
FIG = HERE / 'figures/combined_explanatory'
QA = HERE / 'verification'
FIG.mkdir(parents=True, exist_ok=True)
QA.mkdir(exist_ok=True)
sys.path.insert(0, str(BASE))
import geometry_artwork as g

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

sources = {BASE / 'figures/comsol/geometry.pdf': 'geometry.pdf'}
for n in ['source_calibration', 'coupling', 'full_common_paths', 'primary_full',
          'partial_joint', 'partial_source_threshold']:
    sources[EXPL / f'figures/explanatory/{n}.pdf'] = f'{n}.pdf'
for src, dest in sources.items():
    shutil.copyfile(src, FIG / dest)

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                     'axes.labelsize': 8, 'axes.titlesize': 9,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'pdf.fonttype': 42, 'axes.linewidth': .6})
fig = plt.figure(figsize=(7.1, 4.5), facecolor='white')
axes = [fig.add_axes([.08, .21, .40, .68]), fig.add_axes([.58, .21, .40, .68])]
colors = ['#d55e00', '#0072b2', '#a12c52']
for ax, kind in zip(axes, ['full', 'partial']):
    ax.add_patch(Rectangle((-75, 0), 150, 25.4, fc='#f3f3f3', ec='#555555', lw=.6))
    g.outline(ax, kind, flaw=False, heat=False)
    ax.set(xlim=(-24, 7), ylim=(27.5, -2), xlabel='$x$ (mm)', ylabel='Depth (mm)')
    ax.set_aspect('equal', adjustable='box')
    ax.set_xticks([-20, -10, 0])
    ax.set_yticks([0, 10, 20, 25.4])
axes[0].set_title('(a) Full weld: S | S', loc='left', pad=10)
axes[1].set_title('(b) Partial weld: S | SS', loc='left', pad=10)
full_s = np.array([-12., 0.]); full_p = np.array([0., 20.]); full_r = np.array([-20.0534, 0.])
part_s = np.array([0., 16.4]); part_p = np.array([2., 20.]); part_r = np.array([-12., 0.])
bounce = np.array([part_p[0] + (part_r[0]-part_p[0])*(25.4-part_p[1])/(50.8-part_p[1]), 25.4])
for ax, s, p, r, ret in [(axes[0], full_s, full_p, full_r, [full_p, full_r]),
                        (axes[1], part_s, part_p, part_r, [part_p, bounce, part_r])]:
    ax.plot(*np.array([s, p]).T, color=colors[0], lw=1.5, zorder=8)
    ax.plot(*np.array(ret).T, color=colors[1], lw=1.4, zorder=8)
    ax.scatter(*s, c=colors[0], marker='v', s=35, edgecolors='white', lw=.5, zorder=10)
    ax.scatter(*r, c=colors[1], marker='s', s=24, edgecolors='white', lw=.5, zorder=10)
    ax.scatter(*p, c='black', s=19, zorder=10)
    ax.annotate('Trial pixel', xy=p, xytext=(p[0]-12, p[1]+2), fontsize=7, zorder=12,
                arrowprops={'arrowstyle': '-', 'lw': .6, 'shrinkB': 4})
axes[0].text(-12, -1, 'S0', ha='center', fontsize=7, zorder=12)
axes[0].text(-20.0534, -1, 'R94', ha='center', fontsize=7, zorder=12)
axes[1].text(-12, -1, 'Receiver', ha='center', fontsize=7, zorder=12)
axes[1].text(1, 15.4, 'Source', fontsize=7, zorder=12)
axes[1].plot(*np.array([part_p, part_r]).T, color=colors[2], ls='--', lw=1, zorder=8)
axes[1].text(-1, 8, 'Air', fontsize=8, ha='center', color='#555555', zorder=12)
axes[1].annotate('Backwall bounce', xy=bounce, xytext=(-19, 25.0), fontsize=7, zorder=12,
                 arrowprops={'arrowstyle': '-', 'lw': .6, 'shrinkB': 4})
handles = [Line2D([], [], c=colors[0], lw=1.5, label='Incident route'),
           Line2D([], [], c=colors[1], lw=1.5, label='Receiving route'),
           Line2D([], [], c=colors[2], lw=1, ls='--', label='Blocked direct return')]
legend = fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.53, .055),
                    ncol=3, fontsize=8, frameon=False)
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
lb = legend.get_window_extent(renderer)
assert all(not lb.overlaps(ax.get_window_extent(renderer)) for ax in axes)
for ax in axes:
    for obj in [ax.title, ax.xaxis.label, ax.yaxis.label, *ax.texts,
                *ax.get_xticklabels(), *ax.get_yticklabels()]:
        if not obj.get_text(): continue
        b = obj.get_window_extent(renderer)
        assert b.x0 >= 0 and b.y0 >= 0 and b.x1 <= fig.bbox.width and b.y1 <= fig.bbox.height, obj.get_text()
    q = ax.transData.transform([[0, 20], [1, 20], [0, 21]])
    assert np.isclose(np.linalg.norm(q[1]-q[0]), np.linalg.norm(q[2]-q[0]))
fig.savefig(FIG/'path_examples.pdf', facecolor='white')
fig.savefig(QA/'path_examples.png', dpi=180, facecolor='white')
plt.close(fig)
cS = 3.237967635
full_delay = .030 + (np.linalg.norm(full_s-full_p)+np.linalg.norm(full_r-full_p))/cS
part_delay = .030 + sum(np.linalg.norm(b-a) for a,b in [(part_s,part_p),(part_p,bounce),(bounce,part_r)])/cS
assert abs(full_delay-15.98) < .01 and abs(part_delay-11.751) < .001

print('Regenerated combined path diagram and copied source artwork.')
