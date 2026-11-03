"""Publication figures and tables from the portable LOF-only evidence."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from figure_layout import verify_layout, legend_strip

HERE = Path(__file__).resolve().parent
DATA = HERE/'data'; FIG = HERE/'figures/core_saft'; TABLE = HERE/'tables/core_saft'
INDEX = json.loads((DATA/'images.json').read_text())
DS = {k:json.loads((DATA/f'{k}_dataset.json').read_text()) for k in ['full','partial']}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':8.5,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.6})
CM = plt.get_cmap('magma').copy(); CM.set_bad('#d2d2d2')


def row(kind, path, group='fine'):
    return next(r for r in INDEX if r['kind']==kind and r['id']==f'baseline__{path}__{group}')


def label(path, tex=False):
    n=int(path.split('-')[0]); modes=path.split('_')[1]
    return ('$'+modes[:n]+r'\mid '+modes[n:]+'$') if tex else modes[:n]+' | '+modes[n:]


def outline(ax, kind, marker=True):
    d=DS[kind]; w=np.array(d['geometry']['weld'])*1000; w=np.vstack([w,w[0]])
    ax.plot(w[:,0],25.4-w[:,1],c='#999999',lw=.65,zorder=3)
    if kind=='partial':
        p=np.array(d['surface_profile_x_depth_m'])*1000
        ax.fill_between(p[:,0],-10,p[:,1],color='white',lw=0,zorder=4)
        ax.plot(*p.T,c='#555555',lw=.6,zorder=4)
    if marker:
        tips=np.array(d['audit']['defect_benchmark_only']['tips_xy_mm']); tips[:,1]=25.4-tips[:,1]; tips=tips[np.argsort(tips[:,1])]
        ax.plot([1.015,1.04,1.04,1.015],[tips[0,1],tips[0,1],tips[1,1],tips[1,1]],
                c='#777777',lw=.65,transform=ax.get_yaxis_transform(),clip_on=False)
        x=float(tips[:,0].mean())
        ax.plot([x,x],[1.005,1.035],c='#777777',lw=.7,transform=ax.get_xaxis_transform(),clip_on=False)


def save(fig,name):
    verify_layout(fig,name)
    for ext in ['pdf','png']:
        fig.savefig(FIG/f'{name}.{ext}',bbox_inches='tight',pad_inches=.05,dpi=600)
    plt.close(fig)


def heat(ax,r,ref,title,count=False):
    with np.load(DATA/(r['uid']+'.npz')) as a:
        x,z=a['x_mm'],a['depth_mm']; dx=np.diff(x).mean(); dz=np.diff(z).mean()
        value=np.where(np.isfinite(a['amplitude']),a['count'],np.nan) if count else 20*np.log10(np.maximum(a['amplitude']*1e6/ref,1e-15))
        cmap=plt.get_cmap('cividis').copy() if count else CM
        if count:cmap.set_bad('#d2d2d2')
        im=ax.imshow(value,origin='upper',extent=[x[0]-dx/2,x[-1]+dx/2,z[-1]+dz/2,z[0]-dz/2],
                     vmin=2 if count else -40,vmax=5 if count else 0,cmap=cmap,interpolation='nearest',aspect='equal')
        ax.set(xlim=(x[0],x[-1]),ylim=(z[-1],z[0]),xlabel='$x$ (mm)')
    outline(ax,r['kind']); ax.set_title(title,loc='left',pad=10)
    ax.set_xticks([-6,0,6] if r['kind']=='full' else [-7.5,0,7.5]); ax.set_yticks([5,10,15,20,25] if r['kind']=='full' else [17,21,25])
    return im


def full_paths():
    paths=['1-1_LL','1-1_SS','1-2_SSS','2-1_SSS']
    fig,axes=plt.subplots(1,4,figsize=(7.1,3.9),layout='constrained')
    for i,(ax,path) in enumerate(zip(axes,paths)):
        im=heat(ax,row('full',path),10,f'({chr(97+i)}) '+label(path))
        if i==0:ax.set_ylabel('Depth below outer surface (mm)')
        else:ax.set_yticklabels([])
    fig.colorbar(im,ax=axes,orientation='horizontal',pad=.04,shrink=.9,aspect=45,ticks=[-40,-30,-20,-10,0],
                 label=r'Image amplitude (dB); 0 dB = 10 $\mu$m/s')
    save(fig,'full_paths')


def partial_result():
    fig,axes=plt.subplots(2,2,figsize=(7.1,5.0),layout='constrained')
    index=json.loads((DATA/'joint_acquisition/index.json').read_text())
    names=[('control__all','All sources / all receivers'),
           ('fine__both_sides','Weld sources / both side receivers'),
           ('fine__left_only','Weld sources / left side receivers'),
           ('fine__right_only','Weld sources / right side receivers')]
    for i,(ax,(identifier,title)) in enumerate(zip(axes.flat,names)):
        with np.load(DATA/'joint_acquisition'/(identifier+'.npz')) as a:
            x,z=a['x_mm'],a['depth_mm']; dx=np.diff(x).mean(); dz=np.diff(z).mean()
            value=20*np.log10(np.maximum(a['lof_amplitude']*1e6/40,1e-15))
            im=ax.imshow(value,origin='upper',extent=[x[0]-dx/2,x[-1]+dx/2,z[-1]+dz/2,z[0]-dz/2],
                         vmin=-40,vmax=0,cmap=CM,interpolation='nearest',aspect='equal')
            ax.set(xlim=(x[0],x[-1]),ylim=(z[-1],z[0]),xlabel='$x$ (mm)',xticks=[-7.5,0,7.5],yticks=[17,21,25])
        m=next(r for r in index if r['id']==identifier)['metrics'][0]
        ax.set_title(f'({chr(97+i)}) {title}\n'+f"$C_{{\\max}}$ = {m['cmax_db']:.2f} dB; near support {m['near_support_pct']:.1f}%",loc='left',pad=10,fontsize=8)
        outline(ax,'partial')
        if i%2==0:ax.set_ylabel('Depth below outer surface (mm)')
    fig.colorbar(im,ax=axes,orientation='horizontal',pad=.04,shrink=.85,aspect=45,
                 ticks=[-40,-30,-20,-10,0],label=r'Amplitude (dB); 0 dB = 40 $\mu$m/s')
    save(fig,'partial_result')


def segment(ax,a,b,color):
    ax.plot([a[0],b[0]],[a[1],b[1]],c=color,lw=1.1,zorder=5)
    a,b=np.array(a),np.array(b)
    ax.annotate('',xy=.4*a+.6*b,xytext=.6*a+.4*b,arrowprops=dict(arrowstyle='->',color=color,lw=1))


def paths():
    fig=plt.figure(figsize=(7.1,3.05),layout='constrained'); gs=fig.add_gridspec(2,2,height_ratios=[1,.17])
    p=np.array([2.,20.]); colors=['#d55e00','#0072b2']
    for i,kind in enumerate(['full','partial']):
        ax=fig.add_subplot(gs[0,i]); ax.add_patch(Rectangle((-25,0),50,25.4,fc='#f0f0f0',ec='#777777',lw=.6))
        outline(ax,kind,marker=False)
        src=np.array([-12.,0.]) if i==0 else np.array([0.,16.4]); rec=np.array([-20.,0.])
        segment(ax,src,p,colors[0])
        if i:
            b=np.array([p[0]+(rec[0]-p[0])*(25.4-p[1])/(50.8-p[1]),25.4])
            segment(ax,p,b,colors[1]); segment(ax,b,rec,colors[1])
        else:segment(ax,p,rec,colors[1])
        ax.scatter(*src,marker='v',c=colors[0],s=30,zorder=6); ax.scatter(*rec,marker='s',c=colors[1],s=24,zorder=6)
        ax.scatter(*p,c='black',s=17,zorder=6)
        ax.annotate('Trial pixel',xy=p,xytext=(-18,22.5),fontsize=7,arrowprops=dict(arrowstyle='-',lw=.6))
        ax.set(xlim=(-25,25),ylim=(27,-2),xlabel='$x$ (mm)',xticks=[-20,0,20],yticks=[0,10,20,25.4],
               title=['(a) Full weld: L | L','(b) Partial weld: S | SS'][i]);ax.set_aspect('equal')
        if i==0:ax.set_ylabel('Depth below outer surface (mm)')
    handles=[Line2D([],[],c=colors[0],marker='v',label='Source to pixel'),Line2D([],[],c=colors[1],marker='s',label='Pixel to receiver')]
    legend_strip(fig,gs[1,:],handles,ncol=2,fontsize=8)
    save(fig,'paths')


def tables():
    TABLE.mkdir(parents=True,exist_ok=True)
    lines=[]
    for kind,paths in [('full',['1-1_LL','1-1_SS','1-2_SSS','2-1_SSS']),('partial',['1-2_SSS'])]:
        for path in paths:
            r=row(kind,path)['metrics']
            if kind=='partial':
                joint=json.loads((DATA/'joint_acquisition/index.json').read_text())
                r=next(q for q in joint if q['id']=='fine__left_only')['metrics'][0]
            pk=r['peak']
            lines.append(f"{kind.capitalize()} & {label(path,True)} & {pk['um_s']:.2f} & {r['cmax_db']:.2f} & {r['near_support_pct']:.1f} & {pk['distance_mm']:.3f}"+r' \\')
    (TABLE/'results.tex').write_text('\n'.join(lines)+'\n')
    lines=[]
    names={'fine':'All sources / all receivers','shoulders':'Shoulder sources / all receivers','central':'Weld-surface sources / all receivers',
           'plate':'All sources / side-surface receivers','floor':'All sources / weld-surface receivers','shoulders_plate':'Shoulder sources / side-surface receivers'}
    for group,name in names.items():
        values=[]
        for kind,path in [('full','1-1_LL'),('full','1-1_SS'),('partial','1-2_SSS')]:
            m=row(kind,path,group)['metrics']; c='--' if m['cmax_db'] is None else f"{m['cmax_db']:.2f}"
            values += [c, f"{m['near_support_pct']:.1f}"]
        lines.append(name+' & '+' & '.join(values)+r' \\')
    (TABLE/'acquisition.tex').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':
    paths(); full_paths(); partial_result(); tables()
    print('Generated three new figures and two numerical tables from verified LOF-only arrays.')
