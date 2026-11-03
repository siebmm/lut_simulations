"""Scientific comparison figures and compact tables from native exports."""
from pathlib import Path
import json
import numpy as np
HERE=Path(__file__).resolve().parent
if (HERE/'data/dataset.json').exists():
    OUT=HERE/'data'
    class PortableDataset:
        @staticmethod
        def dataset(kind):return json.loads((OUT/'dataset.json').read_text())
        @staticmethod
        def tips(kind):
            d=PortableDataset.dataset(kind);t=np.array(d['audit']['defect_benchmark_only']['tips_xy_mm']);t[:,1]=25.4-t[:,1]
            return t[np.argsort(t[:,1])]
    old=PortableDataset
    def write(p,o):
        p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(o,indent=2),encoding='utf-8')
else:
    from experiment import *
import csv, math, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.patches import Rectangle

FIG=HERE/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,
 'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
 'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.6})
CM=plt.get_cmap('magma').copy();CM.set_bad('#c9c9c9')
def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',pad_inches=.05)
    fig.savefig(FIG/(name+'.png'),bbox_inches='tight',pad_inches=.05,dpi=400);plt.close(fig)
def outline(ax,flaw=True):
    d=old.dataset('partial');w=np.array(d['geometry']['weld'])*1000;w=np.vstack([w,w[0]])
    ax.plot(w[:,0],25.4-w[:,1],c='#999999',lw=.65,zorder=3)
    if flaw:
        tips=old.tips('partial');line,=ax.plot(*tips.T,c='#92f7ef',lw=.9,ls=(0,(3,2)),zorder=5)
        line.set_path_effects([pe.Stroke(linewidth=1.7,foreground='#151515'),pe.Normal()])
        ax.scatter(*tips.T,s=11,fc='none',ec='#92f7ef',lw=.7,zorder=6)
def heat(ax,index,case,reference=80,db=60):
    with np.load(OUT/'numerical'/str(index)/'images.npz') as a:
        x,z=a['x_mm'],a['depth_mm'];dx=x[1]-x[0];dz=z[1]-z[0]
        v=20*np.log10(np.maximum(a[case+'_maps'][2].reshape(len(z),len(x))*1e6/reference,1e-15))
    im=ax.imshow(v,extent=[x[0]-dx/2,x[-1]+dx/2,z[-1]+dz/2,z[0]-dz/2],cmap=CM,vmin=-db,vmax=0,interpolation='nearest',aspect='equal')
    outline(ax);ax.set(xlim=(-7.5,7.5),ylim=(25.25,16.65),xticks=[-7.5,0,7.5],yticks=[17,21,25],xlabel='$x$ (mm)');return im
def geometry():
    req=json.loads((OUT/'verification/request.json').read_text())['ray_checks']
    rays=json.loads((OUT/'verification/native/geometry.json').read_text())['rays']
    examples=[]
    for via in ['wall:4:id-last','corner:4:id-first']:
        candidates=[(q,r) for q,r in zip(req,rays) if q['via']==via and not q.get('outgoing') and r['valid'] and q['b_m']==[.005,.021]]
        target=[.012,0] if via.startswith('wall') else [0,.0164]
        examples.append(min(candidates,key=lambda qr:np.linalg.norm(np.array(qr[0]['a_m'])-target)))
    pr=np.array(old.dataset('partial')['surface_profile_x_depth_m'])*1000
    fig,axs=plt.subplots(1,2,figsize=(7.1,2.65),layout='constrained')
    for ax,(q,r),title in zip(axs,examples,['(a) Finite-bevel receiving branch','(b) Fixed-corner receiving branch']):
        ax.set_facecolor('white');ax.fill_between(pr[:,0],pr[:,1],25.4,color='#eeeeee',zorder=2);ax.plot(*pr.T,c='#555555',lw=.8,zorder=3)
        ax.axhline(25.4,c='#444444',lw=.8);p=np.array(r['points_m'])*1000
        ax.plot(*p.T,c='#be5c21',lw=1.2,zorder=4);ax.scatter(*p[0],marker='o',fc='white',ec='#be5c21',s=30,zorder=5)
        ax.text(p[0,0]+1,p[0,1]+1.2,'receiver',fontsize=7)
        ax.scatter(*p[-1],marker='x',c='black',s=25,zorder=5);ax.text(p[-1,0]+1,p[-1,1]+.2,'p',fontsize=9)
        if 'corner' in q['via']:
            c=pr[4];ax.scatter(*c,s=22,fc='white',ec='#be5c21',zorder=6);ax.annotate('$C_R$',xy=c,xytext=(7,13),arrowprops={'arrowstyle':'-','lw':.6})
        else:
            c=p[1];ax.scatter(*c,s=18,c='#be5c21',zorder=6);ax.annotate('reflection',xy=c,xytext=(-18,7),arrowprops={'arrowstyle':'-','lw':.6})
        ax.set(xlim=(-20,24),ylim=(26,-1),xticks=[-20,0,20],yticks=[0,16.4,25.4],xlabel='$x$ (mm)',title=title);ax.set_aspect('equal')
    axs[0].set_ylabel('Depth below OD (mm)');save(fig,'geometry')

def main():
    rows=json.loads((OUT/'metrics.json').read_text());index=json.loads((OUT/'index.json').read_text())
    def get(path,processing='weighted',case='lof'):return next(r for r in rows if r['path']==path and r['processing']==processing and r['case']==case)
    def best(family):
        candidates=[r for r in rows if r['family']==family and r['processing']=='weighted' and r['case']=='lof' and r['metrics']['cmax_db'] is not None]
        return max(candidates,key=lambda r:r['metrics']['cmax_db']) if candidates else None
    corner_best=best('corner');wall=best('wall');baseline=get('1-2_SSS');summaries=[]
    # An explicitly labelled coverage-qualified illustration, not a replacement protocol.
    corner=max([r for r in rows if r['family']=='corner' and r['processing']=='weighted' and r['case']=='lof' and r['metrics']['near_support_pct']>=80 and r['metrics']['cmax_db'] is not None],key=lambda r:r['metrics']['cmax_db'])
    for family in ['conventional','corner','wall']:
        rr=[r for r in rows if r['family']==family and r['processing']=='weighted' and r['case']=='lof']
        summaries.append(dict(family=family,paths=len(rr),any_support=sum(r['metrics']['supported']>0 for r in rr),near_support=sum(r['metrics']['near_support_pct']>0 for r in rr),
          peak_near=sum(r['metrics']['peak'] is not None and r['metrics']['peak']['distance_mm']<=.5 for r in rr),positive_contrast=sum((r['metrics']['cmax_db'] or -1)>0 for r in rr)))
    selected=[('Conventional, timing mean',get('1-2_SSS','raw')),('Conventional, response weighted',baseline),('C1: best corner contrast',corner_best),('C2: corner, broad coverage',corner)]
    if wall:selected.append(('B1: best bevel contrast',wall))
    def fmt(x):return '--' if x is None else f'{x:.2f}'
    table=[]
    for name,r in selected:
        m=r['metrics'];table.append(name+' & '+' & '.join([fmt(m['cmax_db']),fmt(m['c99_db']),f"{m['near_support_pct']:.1f}",f"{m['support_pct']:.1f}",fmt(m['peak']['distance_mm'] if m['peak'] else None)])+r' \\')
    (HERE/'table_scores.tex').write_text('\n'.join(table)+'\n')
    (HERE/'table_families.tex').write_text('\n'.join(r['family'].capitalize()+' & '+' & '.join(str(r[k]) for k in ['paths','any_support','near_support','peak_near','positive_contrast'])+r' \\' for r in summaries)+'\n')
    summary=dict(baseline=baseline,corner_best=corner_best,corner=corner,wall=wall,families=summaries,corner_reference=get(corner['path'],case='reference'),corner_difference=get(corner['path'],case='difference'),corner_raw=get(corner['path'],'raw'),
       catalogue=json.loads((OUT/'verification/native/path_catalogue.json').read_text()))
    write(OUT/'summary.json',summary)
    with (HERE/'all_path_metrics.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['path','family','processing','spacing_mm','case','support_pct','near_support_pct','Cmax_dB','C99_dB','peak_distance_mm'])
        for r in rows:
            m=r['metrics'];w.writerow([r[k] for k in ['path','family','processing','spacing_mm','case']]+[m['support_pct'],m['near_support_pct'],m['cmax_db'],m['c99_db'],m['peak']['distance_mm'] if m['peak'] else None])
    geometry()
    fig,axs=plt.subplots(2,3,figsize=(7.1,4.2),layout='constrained')
    for j,(name,r) in enumerate([('$S\\mid SS$',baseline),('C2: corner timing',corner)]):
        for i,(case,title) in enumerate([('lof','LOF only'),('reference','Pristine reference'),('difference','LOF - reference')]):
            im=heat(axs[j,i],r['index'],case);axs[j,i].set_title(f'({chr(97+j*3+i)}) {title}')
            if i==0:axs[j,i].set_ylabel(name+'\nDepth (mm)')
    fig.colorbar(im,ax=axs.ravel().tolist(),orientation='horizontal',pad=.04,shrink=.85,aspect=45,label=r'Common scale: dB relative to 80 $\mu$m/s',ticks=[-60,-40,-20,0]);save(fig,'common_scale')
    peaks=[get(corner['path'],case=case)['metrics']['peak']['um_s'] for case in ['lof','reference','difference']]
    ref=10**math.ceil(math.log10(max(peaks)))
    fig,axs=plt.subplots(1,3,figsize=(7.1,2.5),layout='constrained')
    for ax,case,title in zip(axs,['lof','reference','difference'],['(a) LOF only','(b) Pristine reference','(c) LOF - reference']):
        im=heat(ax,corner['index'],case,reference=ref,db=40);ax.set_title(title)
    axs[0].set_ylabel('Depth (mm)');fig.colorbar(im,ax=axs.ravel().tolist(),orientation='horizontal',pad=.05,shrink=.85,aspect=45,label=f'Expanded corner scale: dB relative to {ref:g} '+r'$\mu$m/s',ticks=[-40,-20,0]);save(fig,'corner_expanded')
    fig,axs=plt.subplots(1,2,figsize=(7.1,2.8),layout='constrained')
    for family,color,marker in [('corner','#9a5e97','o'),('wall','#d09145','s')]:
        rr=[r for r in rows if r['family']==family and r['processing']=='weighted' and r['case']=='lof' and r['metrics']['cmax_db'] is not None]
        axs[0].scatter([r['metrics']['near_support_pct'] for r in rr],[r['metrics']['cmax_db'] for r in rr],s=16,c=color,marker=marker,alpha=.7,label=family.capitalize())
    bm=baseline['metrics'];axs[0].scatter([bm['near_support_pct']],[bm['cmax_db']],c='black',s=55,marker='*',label='$S\\mid SS$')
    axs[0].axhline(0,c='#555555',lw=.6);axs[0].set(xlabel='Near-flaw supported pixels (%)',ylabel='$C_{\\max}$ (dB)',title='(a) Every evaluable weighted route',xlim=(-3,103));axs[0].legend(loc='lower left')
    labels=['Corner','Bevel'];rr=summaries[1:]
    axs[1].bar(labels,[r['paths'] for r in rr],color='#dddddd',label='Total hypotheses')
    axs[1].bar(labels,[r['any_support'] for r in rr],color='#9b9b9b',label='Any image support')
    axs[1].bar(labels,[r['near_support'] for r in rr],color='#494949',label='Any near-flaw support')
    axs[1].set(ylabel='Number of paths',title='(b) Coverage before contrast');axs[1].legend(loc='upper right');save(fig,'all_paths')
    for name,r in [('conventional_S_SS',baseline),('corner_C2_exploratory',corner)]:
        dest=OUT/'review_settings'/f'{name}.json';dest.parent.mkdir(exist_ok=True)
        settings=json.loads((OUT/'numerical'/str(r['index'])/'settings.json').read_text())
        # Retain native display defaults; this file is an image recipe, not a fitted amplitude scale.
        settings['display']={}
        write(dest,settings)
    print(json.dumps({k:v for k,v in summary.items() if k!='catalogue'},indent=2))
if __name__=='__main__':main()
