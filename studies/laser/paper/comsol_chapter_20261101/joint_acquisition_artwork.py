"""Regenerate the joint-acquisition figure from the portable evidence."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from figure_layout import verify_layout
HERE=Path(__file__).resolve().parent
DATA=HERE/'data'; FIG=HERE/'figures/comsol'
DS={k:json.loads((DATA/f'{k}_dataset.json').read_text()) for k in ['full','partial']}
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':8.5,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.6})
CM=plt.get_cmap('magma').copy(); CM.set_bad('#d2d2d2')


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
    save(fig,'partial_joint')

if __name__=='__main__':partial_result()
