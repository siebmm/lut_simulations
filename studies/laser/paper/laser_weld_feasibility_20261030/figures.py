"""Portable publication artwork, exclusively from archived numerical arrays."""
from pathlib import Path
import json, numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.patches import Rectangle,Patch
from matplotlib.lines import Line2D
HERE=Path(__file__).resolve().parent;DATA=HERE/'data';FIG=HERE/'figures';FIG.mkdir(exist_ok=True)
INDEX=json.loads((DATA/'index.json').read_text());CORE=['1-1_SS','1-1_LL','1-2_SSS','2-1_SSS']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.6})
CM=plt.get_cmap('magma').copy();CM.set_bad('#c9c9c9');REF=80e-6
def row(kind,path,name='baseline'):return next(r for r in INDEX if r['kind']==kind and r['path']==path and r['name']==name)
def load(kind,path,name='baseline'):return np.load(DATA/(row(kind,path,name)['uid']+'.npz'))
def metric(kind,path,name='baseline',case='lof'):return row(kind,path,name)['metrics'][case]
def label(p):n=int(p.split('-')[0]);m=p.split('_')[1];return m[:n]+' | '+m[n:]
def dataset(kind):return json.loads((DATA/f'{kind}_dataset.json').read_text())
def tips(kind):
    d=dataset(kind);t=np.array(d['audit']['defect_benchmark_only']['tips_xy_mm']);t[:,1]=25.4-t[:,1];return t[np.argsort(t[:,1])]
def profile(kind):return np.array(dataset(kind)['surface_profile_x_depth_m'])*1000
def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',pad_inches=.05)
    fig.savefig(FIG/(name+'.png'),bbox_inches='tight',pad_inches=.05,dpi=600);plt.close(fig)
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
def heat(ax,kind,path,name='baseline',case='lof',title=None,wide=False,db=60):
    a=load(kind,path,name);x,z=a['x_mm'],a['depth_mm'];dx=np.diff(x).mean();dz=np.diff(z).mean()
    v=20*np.log10(np.maximum(a[case].reshape(len(z),len(x))/REF,1e-15))
    im=ax.imshow(v,extent=[x[0]-dx/2,x[-1]+dx/2,z[-1]+dz/2,z[0]-dz/2],cmap=CM,vmin=-db,vmax=0,aspect='equal',interpolation='nearest')
    outline(ax,kind);ax.set(xlim=(x[0],x[-1]),ylim=(z[-1],z[0]),xlabel='$x$ (mm)',title=title or label(path))
    ax.set_xticks([-40,-20,0,20,40] if wide else [-7.5,0,7.5]);ax.set_yticks([0,10,20,25.4] if wide else [17,21,25]);return im
def colorbar(fig,axes,im):fig.colorbar(im,ax=axes,orientation='horizontal',pad=.045,shrink=.85,aspect=45,label=r'Amplitude score (dB); common 0 dB = 80 $\mu$m/s',ticks=np.linspace(im.norm.vmin,0,4 if im.norm.vmin==-60 else 3))

def grid_box(ax,kind):
    b=json.loads((DATA/'geometry_display_extents.json').read_text())[kind]['bounds']
    rect=Rectangle((b['x_min_mm'],b['depth_min_mm']),b['x_max_mm']-b['x_min_mm'],b['depth_max_mm']-b['depth_min_mm'],fill=False,ec='#222222',lw=.7,ls='--',zorder=8)
    rect.set_gid('imaging-extent');ax.add_patch(rect);return rect

def specimen(ax,kind):
    ax.set_facecolor('white')
    ax.add_patch(Rectangle((-75,0),150,25.4,fc='#f0f0f0',ec='#333333',lw=.7,zorder=0))
    for x in [-75,60]:ax.add_patch(Rectangle((x,0),15,25.4,fc='#e1e1e1',ec='#aaaaaa',lw=.5,hatch='///',zorder=1))
    outline(ax,kind,heat=False)

def verify_air_shading(fig,axs):
    """Check exported pixel colours above the OD, away from labels/markers."""
    fig.canvas.draw();rgba=np.asarray(fig.canvas.buffer_rgba());checks=[]
    for ax in axs:
        x,y=ax.transData.transform([40,-1]);rgb=rgba[rgba.shape[0]-1-int(round(y)),int(round(x)),:3]
        assert np.all(rgb==255),rgb.tolist()
        checks.append({'sample_x_depth_mm':[40,-1],'air_rgb':rgb.tolist()})
    return checks

def geometry():
    fig,axs=plt.subplots(1,3,figsize=(7.1,2.9),layout='constrained',gridspec_kw={'width_ratios':[1,1,.85]})
    for i,kind in enumerate(['full','partial']):
        ax=axs[i];specimen(ax,kind)
        d=dataset(kind);r=d['audit']['receivers'];s=d['audit']['sources']
        for group,ids,col in [('plate',range(314),'#0072b2'),('central',range(314,len(r)),'#009e73')]:
            xy=np.array([[r[j]['x_m']*1000,25.4-r[j]['y_m']*1000] for j in ids]);ax.scatter(*xy.T,s=3,c=col,zorder=5)
        xy=np.array([[v['x_mm'],25.4-v['xy_m'][1]*1000] for v in s]);ax.scatter(*xy.T,s=25,marker='v',c='#d55e00',ec='white',lw=.3,zorder=8)
        for j in [0,4]:ax.text(xy[j,0],xy[j,1]-1.1,f'S{j}',ha='center',fontsize=7,zorder=10)
        ax.annotate('S1-S3',xy=xy[2],xytext=(25,9 if kind=='full' else 13),arrowprops={'arrowstyle':'-','lw':.6},fontsize=7,ha='center')
        grid_box(ax,kind)
        ax.set(xlim=(-78,78),ylim=(28,-3),xlabel='$x$ (mm)',title=f'({chr(97+i)}) {kind.capitalize()} weld');ax.set_xticks([-75,0,75]);ax.set_yticks([0,16.4,25.4])
        if i==0:ax.set_ylabel('Depth below OD (mm)')
    ax=axs[2];ax.set_facecolor('#f0f0f0');outline(ax,'partial',heat=False);t=tips('partial')
    # The actual model has bonded steel below the meeting point of the J arcs.
    # Show its centreline as a domain partition, not an open root gap.
    root=np.array(dataset('partial')['geometry']['right'][:2])*1000
    ax.plot(root[:,0],25.4-root[:,1],color='#777777',lw=.6,ls=':',zorder=3)
    ax.axhline(16.4,c='#555555',ls=':',lw=.7);ax.axhline(25.4,c='#333333',lw=.7)
    ax.annotate('LOF\n18.467-20.967 mm',xy=t.mean(axis=0),xytext=(-4,22.6),arrowprops={'arrowstyle':'->','lw':.7},fontsize=7)
    ax.text(-4.1,16.1,'Floor: 16.4 mm',fontsize=7,zorder=10);ax.text(-4.1,25.05,'Backwall',fontsize=7,zorder=10)
    ax.set(xlim=(-4.5,4.5),ylim=(26,14.5),xlabel='$x$ (mm)',title='(c) Partial-weld detail');ax.set_xticks([-4,0,4]);ax.set_yticks([16.4,20,25.4])
    ax.set_aspect('equal',adjustable='box')
    # Verify the rendered drawing itself, including its physical aspect ratio.
    air_checks=verify_air_shading(fig,axs[:2]);q=ax.transData.transform([[0,20],[1,20],[0,21]])
    scale=np.linalg.norm(q[1]-q[0])/np.linalg.norm(q[2]-q[0])
    assert np.isclose(scale,1,rtol=1e-12)
    boxes=[[p for p in panel.patches if p.get_gid()=='imaging-extent'] for panel in axs[:2]]
    assert [len(p) for p in boxes]==[1,1]
    assert boxes[0][0].get_bbox().bounds==(-6,5,12,20.25)
    assert np.allclose(boxes[1][0].get_bbox().bounds,(-7.5,16.65,15,8.6))
    for panel in axs:
        xy=panel.lines[0].get_xydata();assert np.array_equal(xy[0],xy[-1])
    audit={'closed_weld_polygons':True,'full_weld_box_x_z_width_height_mm':list(boxes[0][0].get_bbox().bounds),'partial_weld_box_x_z_width_height_mm':list(boxes[1][0].get_bbox().bounds),
           'whole_specimen_width_mm':150,'absorber_regions_x_mm':[[-75,-60],[60,75]],
           'above_surface_air_checks':air_checks,
           'detail_x_to_depth_scale_ratio':float(scale),'root_closure_depth_mm':float(25.4-root[1,1]),
           'root_below_closure':'Bonded centreline partition, not an open crack.',
           'source_geometry':'Unmodified archived dataset polygon vertices.'}
    dest=HERE/'verification';dest.mkdir(exist_ok=True)
    (dest/'geometry_figure_check.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    save(fig,'F01_geometry')

def whole_specimen():
    fig,axs=plt.subplots(2,1,figsize=(7.1,3.9),layout='constrained')
    for i,kind in enumerate(['full','partial']):
        ax=axs[i];specimen(ax,kind);grid_box(ax,kind)
        d=dataset(kind);r=d['audit']['receivers'];s=d['audit']['sources']
        for ids,col in [(range(314),'#0072b2'),(range(314,len(r)),'#009e73')]:
            xy=np.array([[r[j]['x_m']*1000,25.4-r[j]['y_m']*1000] for j in ids]);ax.scatter(*xy.T,s=3,c=col,zorder=5)
        xy=np.array([[v['x_mm'],25.4-v['xy_m'][1]*1000] for v in s]);ax.scatter(*xy.T,s=23,marker='v',c='#d55e00',ec='white',lw=.3,zorder=8)
        b=json.loads((DATA/'geometry_display_extents.json').read_text())[kind]['bounds']
        txt=f'({chr(97+i)}) {kind.capitalize()} weld: '+('saved full-weld grid' if kind=='full' else 'partial-weld study grid')
        ax.set(xlim=(-78,78),ylim=(29,-4),xlabel='$x$ (mm)',ylabel='Depth (mm)',title=txt)
        ax.set_aspect('equal',adjustable='box');ax.set_xticks([-75,-60,0,60,75]);ax.set_yticks([0,16.4,25.4] if kind=='partial' else [0,25.4])
    handles=[Line2D([],[],ls='',marker='v',color='#d55e00',label='Sources'),Line2D([],[],ls='',marker='.',color='#0072b2',label='Plate receivers'),Line2D([],[],ls='',marker='.',color='#009e73',label='Central receivers'),Line2D([],[],c='#b2182b',ls='--',label='Known LOF'),Patch(fc='#e1e1e1',ec='#aaaaaa',hatch='///',label='15 mm absorbers')]
    fig.legend(handles=handles,loc='outside lower center',ncol=5,frameon=False,fontsize=7)
    air_checks=verify_air_shading(fig,axs)
    (HERE/'verification/whole_specimen_figure_check.json').write_text(json.dumps({'above_surface_air_checks':air_checks,'equal_coordinate_scales':all(ax.get_aspect()==1 for ax in axs)},indent=2),encoding='utf-8')
    save(fig,'F07_whole_specimen')

def physics():
    # Curves are archived, not refitted to any reconstructed indication.
    curves=np.load(DATA/'directivity.npz');angles=curves['angle_deg']
    fig,axs=plt.subplots(1,3,figsize=(7.1,2.8),layout='constrained',gridspec_kw={'width_ratios':[1,1,1.1]})
    for i,key in enumerate(['source','vy']):
        for mode,col,ls in [('L','#0072b2','-'),('S','#d55e00','--')]:axs[i].plot(angles,curves[mode+'_'+key],color=col,ls=ls,label=mode,lw=1.2)
        axs[i].axhline(.2,c='#555555',ls=':',lw=.8);axs[i].set(xlim=(0,80),ylim=(0,1.04),xlabel='Angle from normal (deg)',title=['(a) Dipole launch','(b) Normal reception'][i]);axs[i].grid(alpha=.15);axs[i].legend(frameon=False,loc='upper left')
    axs[0].set_ylabel('Separately normalized magnitude')
    ax=axs[2];ax.set_facecolor('#f3f3f3');outline(ax,'partial',False,False)
    s=np.array([1.5,16.4]);p=np.array([0.,21.]);r=np.array([-12.,0.]);dz=np.array([25.4-p[1],25.4]);bounce=np.array([p[0]+(r[0]-p[0])*dz[0]/dz.sum(),25.4])
    ax.plot(*np.array([s,p]).T,c='#d55e00',lw=1.4,label='Incident S');ax.plot(*np.array([p,bounce,r]).T,c='#0072b2',lw=1.3,ls='--',label='Reflected SS')
    ax.scatter(*s,s=32,marker='v',c='#d55e00',zorder=8);ax.scatter(*r,s=24,marker='s',c='#0072b2',zorder=8);ax.scatter(*p,s=22,c='black',zorder=8)
    ax.annotate('Pixel',xy=p,xytext=(5,21),fontsize=7,arrowprops={'arrowstyle':'-','lw':.6});ax.axhline(25.4,c='#333333',lw=.7)
    ax.set(xlim=(-16,10),ylim=(27,-2),xlabel='$x$ (mm)',ylabel='Depth (mm)',title='(c) S | SS example');ax.legend(frameon=False,loc='upper right',bbox_to_anchor=(1,.77),fontsize=6.5)
    save(fig,'F02_physics')

def atlas():
    fig,axs=plt.subplots(2,4,figsize=(7.1,3.5),layout='constrained')
    for i,kind in enumerate(['full','partial']):
        for j,p in enumerate(CORE):
            ax=axs[i,j];im=heat(ax,kind,p,title=f'({chr(97+i*4+j)}) '+label(p));m=metric(kind,p)
            ax.text(.025,.04,f'$C_{{max}}$ = {m["cmax_db"]:.1f} dB',transform=ax.transAxes,color='white',fontsize=6.3,bbox={'fc':'black','ec':'none','alpha':.7,'pad':1.5})
            if j==0:ax.set_ylabel(kind.capitalize()+' weld\nDepth (mm)')
            else:ax.set_yticklabels([])
            if i==0:ax.set_xlabel('');ax.set_xticklabels([])
    colorbar(fig,axs,im);save(fig,'F03_core_images')

def improvement():
    names=['baseline','weighted','band1_weighted'];titles=['Equal mean, 0.5-5 MHz','Response-weighted, 0.5-5 MHz','Response-weighted, 1-5 MHz']
    fig,axs=plt.subplots(2,3,figsize=(7.1,4.1),layout='constrained')
    for i,case in enumerate(['lof','reference']):
        for j,(name,title) in enumerate(zip(names,titles)):
            ax=axs[i,j];im=heat(ax,'partial','1-2_SSS',name,case,title=f'({chr(97+i*3+j)}) '+title if i==0 else f'({chr(97+i*3+j)}) Pristine control',db=40);m=metric('partial','1-2_SSS',name,case)
            ax.text(.025,.04,f'$I_N$ = {m["near"]["um_s"]:.2f} $\\mu$m/s',transform=ax.transAxes,color='white',fontsize=7,bbox={'fc':'black','ec':'none','alpha':.7,'pad':1.5})
            if j==0:ax.set_ylabel(('LOF present' if i==0 else 'No LOF')+'\nDepth (mm)')
            else:ax.set_yticklabels([])
            if i==0:ax.set_xlabel('');ax.set_xticklabels([])
    colorbar(fig,axs,im);save(fig,'F04_improvement')

def wide():
    fig,ax=plt.subplots(figsize=(7.1,2.7),layout='constrained')
    im=heat(ax,'partial','1-2_SSS','weighted_wide',title='Partial weld: response-weighted S | SS, 0.5-5 MHz',wide=True,db=40)
    ax.set_ylabel('Depth (mm)');ax.add_patch(Rectangle((-7.5,16.65),15,8.6,fill=False,ec='#444444',lw=.8,ls='--',zorder=8))
    colorbar(fig,ax,im);save(fig,'F05_wide_field')

def acquisition():
    specs=[('full','1-1_SS','baseline','Full S | S'),('full','1-1_LL','baseline','Full L | L'),('full','1-2_SSS','baseline','Full S | SS'),('full','2-1_SSS','baseline','Full SS | S'),('partial','1-2_SSS','weighted','Partial S | SS\n(response-weighted)')]
    src=json.loads((DATA/'source_subsets.json').read_text());ws=json.loads((DATA/'weighted_source_subsets.json').read_text())
    groups=['All sources','Shoulders','Central three','All receivers','Both plates','Central surface']
    values=np.zeros((5,6,2));mask=np.zeros((5,6),bool)
    for i,(kind,p,name,_) in enumerate(specs):
        for j,group in enumerate(['all','shoulders','central']):
            rr=next(v for v in (ws if name=='weighted' else src) if v['source_group']==group and (name=='weighted' or (v['kind']==kind and v['path']==p and v['name']==name)))
            m=rr['metrics'];values[i,j]=[m['near_support_pct'],m['cmax_db'] if m['cmax_db'] is not None else np.nan]
        for j,n in enumerate([name,'weighted_plate' if name=='weighted' else 'plate_only','weighted_floor' if name=='weighted' else 'floor_only'],3):
            m=metric(kind,p,n);values[i,j]=[m['near_support_pct'],m['cmax_db'] if m['cmax_db'] is not None else np.nan]
    fig,axs=plt.subplots(1,2,figsize=(7.1,3.5),layout='constrained')
    for j,(title,cm,vmin,vmax) in enumerate([('(a) Near-flaw support (%)','cividis',0,100),('(b) Maximum-background contrast (dB)','RdBu',-15,15)]):
        cmap=plt.get_cmap(cm).copy();cmap.set_bad('#e5e5e5');im=axs[j].imshow(values[:,:,j],cmap=cmap,vmin=vmin,vmax=vmax,aspect='auto');axs[j].set_title(title,fontsize=8)
        axs[j].set_xticks(range(6),groups,rotation=42,ha='right');axs[j].set_yticks(range(5),[r[3] for r in specs] if j==0 else ['']*5);axs[j].axvline(2.5,c='white',lw=2)
        for k in range(5):
            for l in range(6):
                v=values[k,l,j];t='n/a' if not np.isfinite(v) else f'{v:.0f}' if j==0 else f'{v:.1f}';white=(j==0 and v<50) or(j==1 and abs(v)>8)
                axs[j].text(l,k,t,ha='center',va='center',color='white' if white else 'black',fontsize=7)
        fig.colorbar(im,ax=axs[j],orientation='horizontal',pad=.03,shrink=.85,ticks=[0,50,100] if j==0 else [-15,0,15])
    save(fig,'F06_acquisition')
if __name__=='__main__':
    for f in [geometry,physics,atlas,improvement,wide,acquisition,whole_specimen]:f();print(f.__name__,flush=True)
