from pathlib import Path
import hashlib,json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch,Rectangle
from figure_layout import legend_strip,verify_layout
HERE=Path(__file__).resolve().parent;DATA=HERE/'data';FIG=HERE/'figures/comsol'
OUTPUT=FIG;NAME='geometry'
BLUE,GREEN,ORANGE,RED='#0072b2','#009e73','#d55e00','#b2182b'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.6})
def save(fig,name):
    name='physics' if name=='F02_physics' else name
    verify_layout(fig,name)
    FIG.mkdir(parents=True,exist_ok=True)
    for suffix in ['pdf','svg','png']:fig.savefig(FIG/(name+'.'+suffix),bbox_inches='tight',pad_inches=.05,dpi=600)
    plt.close(fig)

def dataset(kind):return json.loads((DATA/f'{kind}_dataset.json').read_text())

def tips(kind):
    d=dataset(kind);t=np.array(d['audit']['defect_benchmark_only']['tips_xy_mm']);t[:,1]=25.4-t[:,1];return t[np.argsort(t[:,1])]

def profile(kind):return np.array(dataset(kind)['surface_profile_x_depth_m'])*1000

def outline(ax,kind,flaw=True,heat=True):
    d=dataset(kind);w=np.array(d['geometry']['weld'])*1000
    # Dataset vertices describe a polygon: its closing root segment is implicit.
    # A plain open polyline omits the final left-arc-to-centre segment.
    w=np.vstack([w,w[0]])
    ax.plot(w[:,0],25.4-w[:,1],c='#9a9a9a' if heat else '#777777',lw=.6,zorder=3)
    if kind=='partial':
        p=profile(kind);ax.fill_between(p[:,0],-3,p[:,1],color='white',lw=0,zorder=4);ax.plot(*p.T,c='#555555',lw=.6,zorder=4)
    if flaw:
        t=tips(kind);c='#82f5ff' if heat else '#b2182b';line,=ax.plot(*t.T,c=c,lw=.9,ls=(0,(3,1.5)),zorder=6)
        if heat:line.set_path_effects([pe.Stroke(linewidth=1.8,foreground='#111111'),pe.Normal()])
        ax.scatter(*t.T,s=12,fc='none',ec=c,lw=.7,zorder=7)

def grid_box(ax,kind):
    b=json.loads((DATA/'geometry_display_extents.json').read_text())[kind]['bounds']
    rect=Rectangle((b['x_min_mm'],b['depth_min_mm']),b['x_max_mm']-b['x_min_mm'],b['depth_max_mm']-b['depth_min_mm'],fill=False,ec='#222222',lw=.7,ls='--',zorder=8)
    rect.set_gid('imaging-extent');ax.add_patch(rect);return rect

def specimen(ax,kind):
    ax.set_facecolor('white')
    ax.add_patch(Rectangle((-75,0),150,25.4,fc='#f0f0f0',ec='#333333',lw=.7,zorder=0))
    for x in [-75,60]:ax.add_patch(Rectangle((x,0),15,25.4,fc='#e1e1e1',ec='#aaaaaa',lw=.5,hatch='///',zorder=1))
    outline(ax,kind,heat=False)

def physics():
    # Curves are archived, not refitted to any reconstructed indication.
    curves=np.load(DATA/'directivity.npz');angles=curves['angle_deg']
    fig=plt.figure(figsize=(7.1,2.85),layout='constrained')
    gs=fig.add_gridspec(2,3,width_ratios=[1,1,1.1],height_ratios=[1,.16])
    axs=[fig.add_subplot(gs[0,i]) for i in range(3)]
    for i,key in enumerate(['source','vy']):
        for mode,col,ls in [('L','#0072b2','-'),('S','#d55e00','--')]:axs[i].plot(angles,curves[mode+'_'+key],color=col,ls=ls,label=mode,lw=1.2)
        axs[i].axhline(.2,c='#555555',ls=':',lw=.8,label='20% gate');axs[i].set(xlim=(0,80),ylim=(0,1.04),xlabel='Angle from normal (deg)',title=['(a) Dipole launch','(b) Normal reception'][i]);axs[i].grid(alpha=.15)
    legend_strip(fig,gs[1,:2],*axs[0].get_legend_handles_labels(),ncol=3,fontsize=7)
    axs[0].set_ylabel('Separately normalized magnitude')
    ax=axs[2];ax.set_facecolor('#f3f3f3');outline(ax,'partial',False,False)
    s=np.array([1.5,16.4]);p=np.array([0.,21.]);r=np.array([-12.,0.]);dz=np.array([25.4-p[1],25.4]);bounce=np.array([p[0]+(r[0]-p[0])*dz[0]/dz.sum(),25.4])
    ax.plot(*np.array([s,p]).T,c='#d55e00',lw=1.4,label='Incident S');ax.plot(*np.array([p,bounce,r]).T,c='#0072b2',lw=1.3,ls='--',label='Reflected SS')
    ax.scatter(*s,s=32,marker='v',c='#d55e00',zorder=8);ax.scatter(*r,s=24,marker='s',c='#0072b2',zorder=8);ax.scatter(*p,s=22,c='black',zorder=8)
    ax.annotate('Pixel',xy=p,xytext=(5,21),fontsize=7,arrowprops={'arrowstyle':'-','lw':.6});ax.axhline(25.4,c='#333333',lw=.7)
    ax.set(xlim=(-16,10),ylim=(27,-2),xlabel='$x$ (mm)',ylabel='Depth (mm)',title='(c) S | SS example')
    legend_strip(fig,gs[1,2],*ax.get_legend_handles_labels(),ncol=1,fontsize=7)
    save(fig,'F02_physics')

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

def geometry():
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
    detail.text(0, 15.45, 'Weld surface\n16.4 mm', fontsize=7, ha='center', va='center', zorder=10)
    detail.text(-4.0, 25.12, 'Backwall', fontsize=7, zorder=10)
    detail.set(xlim=(-4.5, 4.5), ylim=(26.0, 14.5), xlabel='$x$ (mm)', ylabel='Depth (mm)')
    detail.set_xticks([-4, 0, 4])
    detail.set_yticks([16.4, 20, 25.4])
    panel(detail, 'c', 'Partial-weld detail')

    handles = [Line2D([], [], ls='', marker='v', markersize=5.5, color=ORANGE,
                      label='Laser sources'),
               Line2D([], [], ls='', marker='o', markersize=3.5, color=BLUE,
                      label='Side-surface receivers'),
               Line2D([], [], ls='', marker='o', markersize=3.5, color=GREEN,
                      label='Weld-surface receivers'),
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
    verify_layout(fig,NAME)
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
        'are measured. Orange triangles mark the shoulder sources S0/S4 and weld-surface sources '
        'S1-S3. Blue and green markers denote side-surface and weld-surface receivers, respectively '
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

if __name__=='__main__':
    geometry();physics()
