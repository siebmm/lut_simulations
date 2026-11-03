"""Turn the fixed minimum-source audit into a figure, tables and explanatory text."""
from pathlib import Path
import json,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent;BASE=HERE.parent/'comsol_chapter_20261101'
OUT=HERE.parents[1]/'results/partial_source_requirement_20261101'
sys.path.insert(0,str(BASE))
import geometry_artwork as g

def main():
    rows=json.loads((OUT/'metrics.json').read_text())['rows']
    def row(id,case='lof'):return next(r for r in rows if r['id']==id and r['case']==case)
    lines=[]
    for tag,label in [('all','All sources / all receivers'),('both','Weld sources / both sides'),('left','Weld sources / left side')]:
        r=row(tag+'__one')
        for threshold,key in [(2,'baseline_metrics'),(1,'metrics')]:
            m=r[key]
            lines.append(f"{label} & {threshold} & {m['cmax_db']:.2f} & {m['near_support_pct']:.1f} & {m['peak']['distance_mm']:.3f}"+r' \\')
    (HERE/'source_threshold_table.tex').write_text('\n'.join(lines)+'\n')
    lines=[]
    for id,label in [('all__one','All / standard'),('both__one','Weld, both / standard'),
        ('left__one','Weld, left / standard'),('left_finer__one','Weld, left / finer'),
        ('left_wide__one','Weld, left / wide'),('both_wide__one','Weld, both / wide')]:
        rs=[row(id,c) for c in ('lof','reference','difference')]
        lines.append(label+' & '+' & '.join(f"{r['metrics']['cmax_db']:.2f}" for r in rs)+
                     f" & {rs[0]['metrics']['near_support_pct']:.1f}"+r' \\')
    (HERE/'source_threshold_controls.tex').write_text('\n'.join(lines)+'\n')

    old=np.load(BASE/'data/joint_acquisition/fine__left_only.npz')
    new=np.load(OUT/'left__one/arrays.npz')
    x,z=old['x_mm'],old['depth_mm'];shape=(len(z),len(x))
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
       'axes.titlesize':9,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42})
    fig,axs=plt.subplots(2,3,figsize=(7.2,4.7))
    fig.subplots_adjust(left=.09,right=.90,bottom=.12,top=.91,hspace=.42,wspace=.29)
    cmap=plt.get_cmap('magma').copy();cmap.set_bad('#d4d4d4');t=g.tips('partial')
    for i,threshold in enumerate((2,1)):
        for j,case in enumerate(('lof','reference','difference')):
            ax=axs[i,j]
            values=old[case+'_amplitude'] if i==0 else new[case+'_2_amplitude'].reshape(shape)
            db=20*np.log10(np.maximum(values,1e-30)/(60e-6))
            im=ax.pcolormesh(x,z,db,cmap=cmap,vmin=-45,vmax=0,shading='nearest',rasterized=True)
            g.outline(ax,'partial',flaw=False,heat=False)
            ax.set(xlim=(-7.5,7.5),ylim=(25.3,16.4),xlabel='$x$ (mm)')
            ax.set_aspect('equal',adjustable='box');ax.set_xticks([-5,0,5]);ax.set_yticks([18,21,24])
            ax.set_title(f"({chr(97+3*i+j)}) {['LOF only','Pristine','Waveform difference'][j]}",fontsize=8)
            if j==0:ax.set_ylabel(f'Minimum {threshold} sources\nDepth (mm)')
            ax.plot([7.8,7.8],[t[0,1],t[1,1]],c=g.RED,lw=.9,clip_on=False)
            for tip in t:ax.plot([7.5,7.8],[tip[1],tip[1]],c=g.RED,lw=.9,clip_on=False)
    cax=fig.add_axes([.94,.24,.015,.48]);fig.colorbar(im,cax=cax,ticks=[0,-15,-30,-45])
    cax.set_ylabel('Amplitude (dB)',fontsize=8)
    fig.text(.5,.04,'Shared reference: 0 dB = 60 µm/s; grey = unsupported',ha='center',fontsize=8)
    for suffix in ('pdf','png'):
        fig.savefig(HERE/f'figures/explanatory/partial_source_threshold.{suffix}',dpi=600)
    plt.close(fig)
    print('Minimum-source comparison figure and two tables generated from verified arrays.')

if __name__=='__main__':main()
