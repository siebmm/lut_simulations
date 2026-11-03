"""Explanatory diagrams and coverage maps; never alter archived image amplitudes."""
from pathlib import Path
import hashlib, json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.colors import ListedColormap, BoundaryNorm

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / 'comsol_chapter_20261101'
sys.path.insert(0, str(BASE))
import geometry_artwork as g
OUT = HERE / 'figures/explanatory'
QA = HERE / 'verification'
OUT.mkdir(parents=True, exist_ok=True)
QA.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8,
                     'axes.labelsize': 8, 'axes.titlesize': 9,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7,
                     'legend.fontsize': 7, 'pdf.fonttype': 42,
                     'ps.fonttype': 42, 'axes.linewidth': .6})
INPUTS = {}
CHECKS = {}

def record(path):
    INPUTS[str(path.relative_to(BASE))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path

def save(fig, name):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    for ax in fig.axes:
        for artist in [ax.title, ax.xaxis.label, ax.yaxis.label,
                       *ax.get_xticklabels(), *ax.get_yticklabels(), *ax.texts]:
            if not artist.get_visible() or not artist.get_text():
                continue
            bb = artist.get_window_extent(renderer)
            if bb.x0 < 0 or bb.y0 < 0 or bb.x1 > fig.bbox.x1 or bb.y1 > fig.bbox.y1:
                outside.append(artist.get_text())
    for legend in fig.legends:
        bb = legend.get_window_extent(renderer)
        assert fig.bbox.contains(bb.x0, bb.y0) and fig.bbox.contains(bb.x1, bb.y1)
        assert not any(bb.overlaps(ax.get_window_extent(renderer)) for ax in fig.axes)
    assert not outside, (name, outside)
    CHECKS[name] = {'text_inside_canvas': True, 'legends_outside_data': True}
    for suffix in ('pdf', 'png'):
        fig.savefig(OUT / f'{name}.{suffix}', dpi=600, facecolor='white')
    plt.close(fig)

def geometry(kind):
    record(BASE / f'data/{kind}_dataset.json')
    record(BASE / 'data/geometry_display_extents.json')
    fig = plt.figure(figsize=(7.2, 5.2))
    wide = fig.add_axes([.09, .67, .86, .25])
    acquisition = fig.add_axes([.09, .21, .53, .35])
    detail = fig.add_axes([.72, .21, .23, .35])
    for ax in (wide, acquisition, detail):
        g.specimen(ax, kind)
        g.acquisition(ax, kind, detail=ax is not wide)
        ax.set_aspect('equal', adjustable='box')
        ax.set_ylabel('Depth (mm)')
        ax.set_xlabel('$x$ (mm)')
    g.grid_box(wide, kind)
    g.grid_box(acquisition, kind)
    wide.set(xlim=(-78, 78), ylim=(29, -5), title='(a) Complete specimen and acquisition')
    wide.set_xticks([-75, -60, -40, 0, 40, 60, 75])
    wide.set_yticks([0, 25.4])
    acquisition.set(xlim=(-16, 16), ylim=(27, -3), title='(b) Source positions and imaging field')
    acquisition.set_xticks([-12, 0, 12])
    acquisition.set_yticks([0, 16.4, 25.4] if kind == 'partial' else [0, 10, 20, 25.4])
    src = g.acquisition(acquisition, kind, detail=True)
    for i in (0, 4):
        acquisition.annotate(f'S{i}', src[i], xytext=(0, 6), textcoords='offset points',
                             ha='center', fontsize=7,zorder=12)
    if kind == 'partial':
        acquisition.annotate('S1-S3', src[2], xytext=(-8, 12), fontsize=7,
                             arrowprops={'arrowstyle': '-', 'lw': .6},zorder=12)
    else:
        acquisition.text(0, -1.1, 'S1-S3', ha='center', fontsize=7,zorder=12)
    detail.set(xlim=(-4.5, 4.5), ylim=(26, 14.5), title='(c) Embedded LOF')
    detail.set_xticks([-4, 0, 4]); detail.set_yticks([16.4, 20, 25.4])
    flaw = g.tips(kind)
    detail.annotate('LOF', flaw.mean(axis=0), xytext=(-3.4, 19), fontsize=8,
                    arrowprops={'arrowstyle': '->', 'lw': .65},zorder=12)
    if kind == 'partial':
        detail.text(0, 15.5, 'Weld surface', ha='center', fontsize=7,zorder=12)
    detail.text(-4, 25.1, 'Backwall', fontsize=7)
    handles = [Line2D([], [], ls='', marker='v', color=g.ORANGE, label='Sources'),
               Line2D([], [], ls='', marker='o', ms=3, color=g.BLUE, label='Side receivers'),
               Line2D([], [], ls='', marker='o', ms=3, color=g.GREEN, label='Weld receivers'),
               Line2D([], [], color=g.RED, ls='--', label='Prescribed LOF'),
               Patch(fill=False, ec='black', ls='--', label='SAFT field'),
               Patch(fc='#ddd', ec='#aaa', hatch='///', label='Absorbing strips')]
    fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(.52, .035),
               ncol=3, frameon=False, columnspacing=2.5, labelspacing=.8)
    save(fig, f'geometry_{kind}')

def coupling():
    a = np.load(record(BASE / 'data/directivity.npz'))
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.75))
    fig.subplots_adjust(left=.09, right=.98, bottom=.27, top=.86, wspace=.25)
    for i, key in enumerate(('source', 'vy')):
        ax = axs[i]
        for mode, col, ls in [('L', g.BLUE, '-'), ('S', g.ORANGE, '--')]:
            ax.plot(a['angle_deg'], a[f'{mode}_{key}'], c=col, ls=ls, lw=1.2, label=mode)
        ax.axhline(.2, c='#555', ls=':', lw=.9, label='20% admission level')
        ax.set(xlim=(0, 80), ylim=(0, 1.04), xlabel='Angle from normal (degrees)',
               title=['(a) Source launch', '(b) Normal reception'][i])
        ax.grid(alpha=.15)
    axs[0].set_ylabel('Separately normalized amplitude')
    fig.legend(*axs[0].get_legend_handles_labels(), loc='lower center',
               bbox_to_anchor=(.5, .015), ncol=3, frameon=False)
    save(fig, 'coupling')

def timing(kind):
    fig = plt.figure(figsize=(7.2, 3.65))
    ax = fig.add_axes([.09, .25, .40, .60])
    bx = fig.add_axes([.60, .29, .37, .53])
    g.specimen(ax, kind)
    # No actual LOF in this illustrative calculation: remove its red artists.
    for line in list(ax.lines):
        if line.get_color() == g.RED: line.remove()
    for coll in list(ax.collections): coll.remove()
    cs, cl = 3.237967635, 5.943603925
    if kind == 'full':
        s, p, r = np.array([-12., 0]), np.array([0., 20]), np.array([-20.0533975316, 0])
        path = np.array([s, p, r])
        ax.plot(*path.T, c=g.BLUE, lw=1.3)
        ax.set(xlim=(-24, 7), ylim=(27, -3))
        receivers = np.linspace(-43.5534, 43.5534, 500)
        lengths = np.linalg.norm(p-s) + np.sqrt((receivers-p[0])**2+p[1]**2)
        for speed, col, ls, label in [(cs,g.BLUE,'-','S | S'),(cl,g.ORANGE,'--','L | L')]:
            bx.plot(receivers,.03+lengths/speed,c=col,ls=ls,label=label)
        ts = .03 + (np.linalg.norm(p-s)+np.linalg.norm(r-p))/cs
        tl = .03 + (np.linalg.norm(p-s)+np.linalg.norm(r-p))/cl
        bx.scatter(r[0],ts,c=g.BLUE,s=17,zorder=3)
        bx.annotate(f'R94: {ts:.2f} us', (r[0],ts), xytext=(-38,21.5), fontsize=7,
                    arrowprops={'arrowstyle':'-','lw':.6})
        bx.set(xlabel='Receiver $x$ (mm)',ylabel='Predicted arrival (us)',
               xlim=(-44,44),ylim=(6,24),title='(b) One pixel, many receiver delays')
        bx.set_xticks([-40,-20,0,20,40])
        CHECKS['full_example'] = dict(S_delay_us=ts,L_delay_us=tl,
                                     incident_mm=float(np.linalg.norm(p-s)),return_mm=float(np.linalg.norm(r-p)))
    else:
        s, p, r = np.array([0.,16.4]), np.array([2.,20.]), np.array([-12.,0.])
        bounce=np.array([p[0]+(r[0]-p[0])*(25.4-p[1])/(50.8-p[1]-r[1]),25.4])
        ax.plot(*np.array([s,p]).T,c=g.ORANGE,lw=1.4,zorder=8)
        ax.plot(*np.array([p,bounce,r]).T,c=g.BLUE,lw=1.3,zorder=8)
        ax.plot(*np.array([p,r]).T,c='#777',ls=':',lw=1.1,zorder=8)
        ax.set(xlim=(-15,8),ylim=(27,-3))
        incident=np.linalg.norm(p-s)
        direct=np.linalg.norm(r-p)
        reflected=np.linalg.norm(bounce-p)+np.linalg.norm(r-bounce)
        for y,length in [(1,reflected),(0,direct)]:
            bx.barh(y,incident/cs,left=.03,height=.32,color=g.ORANGE)
            bx.barh(y,length/cs,left=.03+incident/cs,height=.32,
                    color=g.BLUE if y else '#aaa',hatch=None if y else '//')
            bx.text(.03+(incident+length)/cs+.15,y,f'{.03+(incident+length)/cs:.2f} us',va='center',fontsize=7)
        bx.set(xlim=(0,14),ylim=(-.6,1.65),xlabel='Predicted arrival (us)',
               yticks=[0,1],yticklabels=['S | S','S | SS'],title='(b) Timing and geometric admission')
        bx.text(0,.29,'Rejected: direct return crosses air',fontsize=7,color='#555')
        bx.text(0,1.3,'Valid reflected route',fontsize=7)
        bx.grid(axis='x',alpha=.15)
        CHECKS['partial_example'] = dict(incident_mm=float(incident),return_mm=float(reflected),
            reflected_delay_us=float(.03+(incident+reflected)/cs),direct_delay_us=float(.03+(incident+direct)/cs),
            bounce_x_mm=float(bounce[0]),direct_floor_crossing_x_mm=-.52,
            reflected_floor_crossing_x_mm=float(bounce[0]+(r[0]-bounce[0])*9/25.4))
    for coord,marker,color,label,offset in [(s,'v',g.ORANGE,'Source',(4,8)),
                                           (r,'s',g.BLUE,'Receiver',(-3,-11)),
                                           (p,'o','black','Trial pixel',(7,0))]:
        ax.scatter(*coord,marker=marker,c=color,s=28,zorder=10)
        ax.annotate(label,coord,xytext=offset,textcoords='offset points',fontsize=7,va='center',zorder=12,
                    bbox={'facecolor':'white','edgecolor':'none','alpha':.8,'pad':.5})
    ax.set_aspect('equal',adjustable='box')
    ax.set(xlabel='$x$ (mm)',ylabel='Depth (mm)',title='(a) Illustrative ray geometry')
    if kind == 'full':
        fig.legend(*bx.get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.5,.07),
                   frameon=False,ncol=2)
    else:
        fig.legend(handles=[Line2D([],[],c=g.ORANGE,label='Incident S'),
                            Line2D([],[],c=g.BLUE,label='Backwall-assisted SS return'),
                            Line2D([],[],c='#777',ls=':',label='Blocked direct return')],
                   loc='lower center',bbox_to_anchor=(.5,.07),frameon=False,ncol=3)
    save(fig,f'timing_{kind}')

def coverage():
    fig, axs=plt.subplots(1,2,figsize=(7.2,3.5))
    fig.subplots_adjust(left=.08,right=.87,bottom=.20,top=.84,wspace=.25)
    cmap=ListedColormap(plt.get_cmap('cividis')([.10,.38,.66,.94]))
    cmap.set_bad('#d4d4d4');norm=BoundaryNorm(np.arange(1.5,6),4)
    t=g.tips('partial');counts=[]
    for ax, stem, title in zip(axs,('control__all','fine__left_only'),
            ('(a) All sources and receivers','(b) Weld sources, left receivers')):
        d=np.load(record(BASE/f'data/joint_acquisition/{stem}.npz'))
        x,z=d['x_mm'],d['depth_mm'];a=d['lof_amplitude'];n=d['lof_count']
        assert a.shape==n.shape==(len(z),len(x))
        valid=np.isfinite(a)
        assert np.all(n[valid]>=2)
        xx,zz=np.meshgrid(x,z);points=np.stack([xx,zz],axis=-1)
        v=t[1]-t[0];u=np.clip(np.sum((points-t[0])*v,axis=-1)/np.dot(v,v),0,1)
        near=np.linalg.norm(points-(t[0]+u[...,None]*v),axis=-1)<=.5
        fraction=100*np.mean(valid[near]);counts.append(float(fraction))
        im=ax.pcolormesh(x,z,np.ma.array(n,mask=~valid),cmap=cmap,norm=norm,shading='nearest',rasterized=True)
        g.outline(ax,'partial',flaw=False,heat=False)
        ax.set(xlim=(-7.5,7.5),ylim=(25.3,16.4),xlabel='$x$ (mm)',title=title)
        ax.set_aspect('equal',adjustable='box');ax.set_ylabel('Depth (mm)')
        ax.plot([7.85,7.85],[t[0,1],t[1,1]],color=g.RED,lw=1,clip_on=False)
        for tip in t:
            ax.plot([7.55,7.85],[tip[1],tip[1]],color=g.RED,lw=1,clip_on=False)
        ax.text(.5,-.30,f'Near-flaw support: {fraction:.1f}%',transform=ax.transAxes,ha='center',fontsize=8)
    cax=fig.add_axes([.91,.30,.018,.38]);fig.colorbar(im,cax=cax,ticks=[2,3,4,5])
    cax.set_ylabel('Qualified sources',labelpad=6)
    assert np.allclose(counts,[95.7575757576,48.4848484848])
    CHECKS['coverage_data']={'near_support_percent':counts,'count_is_qualifying_source_count':True}
    save(fig,'partial_coverage')

if __name__=='__main__':
    geometry('full');geometry('partial');coupling();timing('full');timing('partial');coverage()
    (QA/'new_figures.json').write_text(json.dumps({'checks':CHECKS,'input_sha256':INPUTS},indent=2)+'\n')
    print('Created six vector figures with 600 dpi previews; geometry, timing and coverage checks passed.')
