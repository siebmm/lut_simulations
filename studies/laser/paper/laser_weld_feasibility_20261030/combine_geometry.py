"""One specimen/geometry illustration, redrawn from the original model coordinates.

Run this script directly. It exports PNG (600 dpi), vector PDF and editable SVG
to results/reports/figures, without replacing the two source illustrations.
"""
from pathlib import Path
import hashlib
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np

from figures import DATA, dataset, grid_box, outline, specimen, tips

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
OUTPUT = ROOT / 'studies/laser/results/reports/figures'
NAME = 'laser_weld_combined_geometry'
BLUE, GREEN, ORANGE, RED = '#0072b2', '#009e73', '#d55e00', '#b2182b'


def acquisition(ax, kind, detail=False):
    data = dataset(kind)
    receivers = data['audit']['receivers']
    for ids, color in [(range(314), BLUE), (range(314, len(receivers)), GREEN)]:
        xy = np.array([[receivers[j]['x_m'] * 1000,
                        25.4 - receivers[j]['y_m'] * 1000] for j in ids])
        ax.scatter(*xy.T, s=7 if detail else 2.5, color=color, linewidths=0,
                   zorder=5, clip_on=True)
    sources = np.array([[v['x_mm'], 25.4 - v['xy_m'][1] * 1000]
                        for v in data['audit']['sources']])
    ax.scatter(*sources.T, s=29 if detail else 19, marker='v', color=ORANGE,
               edgecolors='white', linewidths=.45, zorder=9, clip_on=True)
    return sources


def panel(ax, label, title):
    ax.set_title(f'({label}) {title}', loc='left', pad=9, fontsize=9)
    ax.set_aspect('equal', adjustable='box')
    ax.tick_params(length=3, width=.6, pad=3)
    for spine in ax.spines.values():
        spine.set_color('#555555')
        spine.set_linewidth(.55)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                         'axes.labelsize': 8, 'axes.titlesize': 9,
                         'xtick.labelsize': 7, 'ytick.labelsize': 7,
                         'pdf.fonttype': 42, 'ps.fonttype': 42,
                         'svg.fonttype': 'none', 'hatch.linewidth': .55})
    fig = plt.figure(figsize=(7.2, 3.85), facecolor='white')
    # Shared left-column scales retain the specimen's actual proportions.
    full = fig.add_axes([.085, .635, .545, .245])
    partial = fig.add_axes([.085, .275, .545, .245])
    detail = fig.add_axes([.745, .265, .235, .66])
    boxes = {}
    for ax, kind, letter in [(full, 'full', 'a'), (partial, 'partial', 'b')]:
        specimen(ax, kind)
        boxes[kind] = grid_box(ax, kind)
        sources = acquisition(ax, kind)
        ax.set(xlim=(-78, 78), ylim=(29, -8), xlabel='$x$ (mm)',
               ylabel='Depth (mm)')
        ax.set_xticks([-75, -60, 0, 60, 75])
        ax.set_yticks([0, 25.4] if kind == 'full' else [0, 16.4, 25.4])
        panel(ax, letter, 'Full weld' if kind == 'full' else 'Partial weld (9 mm fill)')
        for index, x in [(0, -20), (4, 20)]:
            ax.annotate(f'S{index}', xy=sources[index], xytext=(x, -2.7),
                        fontsize=6.5, ha='center', va='bottom', zorder=10,
                        arrowprops={'arrowstyle': '-', 'lw': .45, 'color': '#555555',
                                    'shrinkA': 1, 'shrinkB': 4})
        if kind == 'full':
            ax.set_xlabel('')  # The two overview panels share one x-axis label.
            ax.text(0, -2.7, 'S1-S3', fontsize=6.5, ha='center', va='bottom', zorder=10)
        else:
            ax.annotate('S1-S3', xy=sources[2], xytext=(-21, 13.2),
                        fontsize=6.5, ha='center', zorder=10,
                        arrowprops={'arrowstyle': '-', 'lw': .5, 'shrinkB': 4})

    # One enlarged view provides the groove floor, continuous J-root and LOF.
    detail.set_facecolor('#f0f0f0')
    outline(detail, 'partial', heat=False)
    acquisition(detail, 'partial', detail=True)
    root = np.array(dataset('partial')['geometry']['right'][:2]) * 1000
    detail.plot(root[:, 0], 25.4 - root[:, 1], color='#777777', lw=.65, ls=':', zorder=3)
    detail.axhline(25.4, color='#333333', lw=.7)
    flaw = tips('partial')
    detail.annotate('LOF\n18.467 to\n20.967 mm', xy=flaw.mean(axis=0), xytext=(-1.7, 18.4),
                    fontsize=7, ha='left', va='top', zorder=10,
                    arrowprops={'arrowstyle': '->', 'lw': .7, 'color': '#333333',
                                'connectionstyle': 'arc3,rad=0.05', 'shrinkB': 3})
    detail.text(0, 15.6, 'Floor: 16.4 mm', fontsize=7, ha='center', zorder=10)
    detail.text(-4.0, 25.12, 'Backwall', fontsize=7, zorder=10)
    detail.set(xlim=(-4.5, 4.5), ylim=(26.0, 14.5), xlabel='$x$ (mm)', ylabel='Depth (mm)')
    detail.set_xticks([-4, 0, 4])
    detail.set_yticks([16.4, 20, 25.4])
    panel(detail, 'c', 'Partial-weld detail')

    handles = [Line2D([], [], ls='', marker='v', markersize=5.5, color=ORANGE,
                      label='Laser sources'),
               Line2D([], [], ls='', marker='o', markersize=3.5, color=BLUE,
                      label='Plate receivers'),
               Line2D([], [], ls='', marker='o', markersize=3.5, color=GREEN,
                      label='Central receivers'),
               Line2D([], [], color=RED, ls='--', lw=1.2, label='Known LOF'),
               Patch(fill=False, edgecolor='#222222', ls='--', lw=.75, label='SAFT grid'),
               Patch(facecolor='#e1e1e1', edgecolor='#aaaaaa', hatch='///',
                     label='15 mm absorbers')]
    # Matplotlib fills legend columns first: acquisition above geometry symbols.
    handles = [handles[j] for j in (0, 3, 1, 4, 2, 5)]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.52, .07),
               ncol=3, frameon=False, fontsize=7.2, columnspacing=2.1,
               handlelength=2.0, labelspacing=.7)

    fig.canvas.draw()
    scales = []
    for ax in (full, partial, detail):
        q = ax.transData.transform([[0, 20], [1, 20], [0, 21]])
        scale = np.linalg.norm(q[1] - q[0]) / np.linalg.norm(q[2] - q[0])
        assert np.isclose(scale, 1, rtol=1e-12)
        scales.append(float(scale))
        weld = ax.lines[0].get_xydata()
        assert np.array_equal(weld[0], weld[-1])
    assert np.allclose(boxes['full'].get_bbox().bounds, [-6, 5, 12, 20.25])
    assert np.allclose(boxes['partial'].get_bbox().bounds, [-7.5, 16.65, 15, 8.6])
    assert np.allclose(tips('full'), flaw)
    pixels = np.asarray(fig.canvas.buffer_rgba())
    for ax, point in [(full, (40, -1)), (partial, (40, -1)), (detail, (0, 15))]:
        x, y = ax.transData.transform(point)
        assert np.all(pixels[pixels.shape[0] - 1 - int(round(y)), int(round(x)), :3] == 255)
    for suffix in ('png', 'pdf', 'svg'):
        fig.savefig(OUTPUT / f'{NAME}.{suffix}', dpi=600, facecolor='white',
                    bbox_inches='tight', pad_inches=.05)
    plt.close(fig)
    caption = (
        'Specimen geometry and laser-ultrasonic acquisition layout for (a) the full weld and '
        '(b) the partial weld, with (c) an enlarged view of the partial-weld floor, J-bevel '
        'and lack of fusion (LOF). The specimens are 150 mm wide and 25.4 mm thick; hatched '
        'regions denote 15 mm absorbing strips. The partial weld has 9 mm fill, placing its '
        'exposed floor 16.4 mm below the original outer surface (OD), from which all depths '
        'are measured. Orange triangles mark the shoulder sources S0/S4 and central sources '
        'S1-S3. Blue and green markers denote plate and central receivers, respectively '
        '(337 total for the full weld; 331 for the partial weld). Black dashed rectangles '
        'show the SAFT grids: x = -6 to 6 mm, depth = 5-25.25 mm for the full weld; '
        'x = -7.5 to 7.5 mm, depth = 16.65-25.25 mm for the partial weld. '
        'The red dashed line and open circles mark the prescribed LOF and its tips at '
        '18.467-20.967 mm depth. Grey internal boundaries are bonded except along the LOF; '
        'the dotted root centreline denotes a bonded domain partition.\n'
    )
    (OUTPUT / f'{NAME}_caption.txt').write_text(caption, encoding='utf-8')
    verification = {'status': 'passed', 'x_to_depth_scale_ratios': scales,
                    'closed_J_bevel': True, 'air_is_white': True,
                    'full_grid_mm': list(boxes['full'].get_bbox().bounds),
                    'partial_grid_mm': list(boxes['partial'].get_bbox().bounds),
                    'same_LOF_coordinates_in_both_models': True,
                    'input_sha256': {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                     for p in [DATA/'full_dataset.json', DATA/'partial_dataset.json',
                                               DATA/'geometry_display_extents.json']}}
    (HERE/'verification/combined_geometry_figure_check.json').write_text(
        json.dumps(verification, indent=2)+'\n', encoding='utf-8')
    print(OUTPUT / f'{NAME}.png')


if __name__ == '__main__':
    main()
