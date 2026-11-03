"""Regenerate the five partial-weld figures from checked reconstructions."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import review as r

F=r.DEST/'figures';F.mkdir(parents=True,exist_ok=True)

def load(path):return json.loads(path.read_text(encoding='utf-8'))
def arr(path):
    with np.load(path) as arrays:return dict(arrays)

confirmation=load(r.CONFIRMATION/'metrics.json')['rows']
def confirmed(name,label):return next(x for x in confirmation if x['candidate']==name and x['label']==label)
def ca(name,label='fine'):return arr(r.locate_confirmation(confirmed(name,label)['id']))
tips=np.array(r.s.dataset('partial')['audit']['defect_benchmark_only']['tips_xy_mm'])
tips[:,1]=25.4-tips[:,1];tips=tips[np.argsort(tips[:,1])]
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':10,'pdf.fonttype':42,'savefig.dpi':200})
cmap=plt.get_cmap('magma').copy();cmap.set_bad('#cccccc')

def panel(ax,a,case,title,wide=False):
    v=a[case+'_2_amplitude'];db=20*np.log10(np.maximum(v.reshape(len(a['z']),len(a['x']))/60e-6,1e-9))
    im=ax.pcolormesh(a['x'],a['z'],np.ma.masked_invalid(db),shading='nearest',rasterized=True,cmap=cmap,vmin=-50,vmax=0)
    ax.set_aspect('equal');ax.set(xlabel='x (mm)',ylabel='Depth (mm)',title=title);ax.invert_yaxis()
    xlo,xhi=float(a['x'][0]),float(a['x'][-1]);zlo,zhi=float(a['z'][0]),float(a['z'][-1])
    if wide:
        profile=np.array(r.s.dataset('partial')['surface_profile_x_depth_m'])*1000
        ax.fill_between(profile[:,0],0,profile[:,1],color='white',zorder=2)
        ax.plot(profile[:,0],profile[:,1],color='#888888',lw=.6)
        ax.set(xlim=(-45,45),ylim=(25.4,0))
    for x,z in tips:
        ax.plot([xlo-.02*(xhi-xlo),xlo],[z,z],color='#0077bb',lw=1.4,clip_on=False)
        ax.plot([x,x],[zhi,zhi+.025*(zhi-zlo)],color='#0077bb',lw=1.4,clip_on=False)
    if not wide:ax.set(xlim=(xlo,xhi),ylim=(zhi,zlo))
    return im

def figures():
    selections=[('Equal mean, 0.5–5 MHz\nOriginal baseline',r.OUT/'survey/all__raw__min2/arrays.npz'),
      ('Weighted, 1–5 MHz\nNo fixed mute',r.DEST/'refinement/unmuted_025/arrays.npz'),
      ('Weighted, 1–5 MHz\nNarrow SS mute',r.DEST/'refinement/muted_025/arrays.npz')]
    fig,axes=plt.subplots(2,3,figsize=(13.7,7.4),layout='constrained')
    for col,(name,p) in enumerate(selections):
        a=arr(p);points=a['points'].reshape(-1,2);d=r.s.distances(points,'partial')
        for row,case in enumerate(('lof','reference')):
            m=r.metric(a[case+'_2_amplitude'],points,d)
            title=(name if row==0 else 'Pristine, identical settings')+f'\nC = {m["cmax_db"]:.2f} dB; near peak = {m["near_um_s"]:.2f} µm/s'
            im=panel(axes[row,col],a,case,title)
    fig.colorbar(im,ax=axes,pad=.012,shrink=.85,label='dB relative to 60 µm/s (shared scale)')
    fig.savefig(F/'partial_reassessment.pdf');fig.savefig(F/'partial_reassessment.png');plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(11,3.6),layout='constrained')
    for name,label,color in [('broad_SSS','All-source S | SS','#0077bb'),('right_SLL','S3/right S | LL','#bb5511')]:
        rows=[x for x in confirmation if x['candidate']==name and x['label'].startswith('width')]
        x=[x['selection']['mute_halfwidth_us'] for x in rows]
        axes[0].plot(x,[x['lof']['cmax_db'] for x in rows],'o-',label=label,color=color)
        axes[1].plot(x,[x['lof']['near_support_pct'] for x in rows],'o-',label=label,color=color)
    for ax in axes:ax.set_xlabel('Mute half-width (µs)');ax.grid(alpha=.2);ax.legend(fontsize=9)
    axes[0].set_ylabel('Maximum-background contrast (dB)');axes[1].set(ylabel='Near-flaw support (%)',ylim=(0,103))
    fig.savefig(F/'mute_sensitivity.pdf');fig.savefig(F/'mute_sensitivity.png');plt.close(fig)
    choices=[('right_SLL_LL03','S3 / right: S | LL'),('right_SSLSS_band15','S1 / right: SSL | SS')]
    fig,axes=plt.subplots(2,2,figsize=(10,7.3),layout='constrained')
    for col,(name,label) in enumerate(choices):
        a=ca(name)
        for row,case in enumerate(('lof','reference')):
            m=confirmed(name,'fine')[case]
            im=panel(axes[row,col],a,case,(label if row==0 else 'Pristine control')+f'\nC = {m["cmax_db"]:.2f} dB; near peak = {m["near_um_s"]:.3g} µm/s')
    fig.colorbar(im,ax=axes,shrink=.85,pad=.012,label='dB relative to 60 µm/s')
    fig.savefig(F/'right_paths_controls.pdf');fig.savefig(F/'right_paths_controls.png');plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(11,6.5),layout='constrained')
    for ax,(name,label) in zip(axes,choices):
        m=confirmed(name,'wide')['lof'];im=panel(ax,ca(name,'wide'),'lof',label+f' — C = {m["cmax_db"]:.2f} dB; peak distance = {m["peak_distance_mm"]:.2f} mm',True)
    fig.colorbar(im,ax=axes,shrink=.85,pad=.015,label='dB relative to 60 µm/s')
    fig.savefig(F/'right_paths_wide.pdf');fig.savefig(F/'right_paths_wide.png');plt.close(fig)
    pathrows=load(r.DEST/'path_confirmation/metrics.json')['rows']
    fig,axes=plt.subplots(2,3,figsize=(13.7,7.4),layout='constrained')
    for col,(index,title) in enumerate([(0,'S1 / left: SLS | L'),(9,'S3 / left: S | SS'),(8,'S1 / left + weld: SLS | S')]):
        row=next(x for x in pathrows if x['candidate']==index and x['label']=='fine');a=arr(r.DEST/'path_confirmation'/row['id']/'arrays.npz')
        for ri,case in enumerate(('lof','reference')):
            m=row[case];im=panel(axes[ri,col],a,case,(title if ri==0 else 'Pristine control')+f'\nC = {m["cmax_db"]:.2f} dB; near peak = {m["near_um_s"]:.2f} µm/s')
    fig.colorbar(im,ax=axes,pad=.012,shrink=.85,label='dB relative to 60 µm/s (shared scale)')
    fig.savefig(F/'atlas_controls.pdf');fig.savefig(F/'atlas_controls.png');plt.close(fig)

if __name__ == '__main__':
    figures()
