"""Publication artwork from float64 arrays; no GUI screenshots or LOF tuning."""
from study import *
from geometry_checks import traces,profile
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.patches import Rectangle
from matplotlib.colors import BoundaryNorm
from scipy.ndimage import map_coordinates
sys.path.insert(0,str(ROOT/'studies/laser/src'))
from saft_imaging import kernel,filter_record
from saft_apertures import directivity,ANGLES

FIG=HERE/'figures';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,'axes.titlesize':9,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.6})
AN={k:json.loads((OUT/k/'analysis.json').read_text()) for k in DATA}
CM=plt.get_cmap('magma').copy();CM.set_bad('#d4d4d4')
COLORS=['#0072b2','#d55e00','#009e73','#aa4499','#a6761d','#cc79a7']
PRIMARY_REF=max(r['signals']['lof']['peak']['um_s'] for a in AN.values() for r in a['runs'] if r['name']=='baseline' and r['signals']['lof']['peak'])/1e6
CONTROL_REF=max(m['peak']['um_s'] for a in AN.values() for r in a['runs'] if r['name']=='baseline' for m in r['signals'].values() if m['peak'])/1e6
def result(kind,path,name='baseline'):return next(r for r in AN[kind]['runs'] if r['path']==path and r['name']==name)
def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',pad_inches=.04)
    fig.savefig(FIG/(name+'.png'),bbox_inches='tight',pad_inches=.04,dpi=600);plt.close(fig)
def outline(ax,kind,flaw=True):
    d=dataset(kind);w=np.array(d['geometry']['weld'])*1000;h=d['material']['thickness_m']*1000
    ax.plot(w[:,0],h-w[:,1],c='#949494',lw=.6,zorder=3)
    if kind=='partial':
        p=profile(kind);ax.fill_between(p[:,0],-.5,p[:,1],color='white',lw=0,zorder=4);ax.plot(*p.T,c='#666666',lw=.7,zorder=4)
    if flaw:
        t=tips(kind);line,=ax.plot(*t.T,c='#7ff9ff',lw=.9,ls=(0,(3,1.5)),zorder=5)
        line.set_path_effects([pe.Stroke(linewidth=1.8,foreground='#111111'),pe.Normal()])
        ax.scatter(*t.T,s=10,fc='none',ec='#7ff9ff',lw=.6,zorder=6)
def image(ax,kind,path,case='lof',name='baseline',reference=None,db=60,title=None,values=None):
    a=load(kind,path,name);x,z=a['x_mm'],a['depth_mm'];v=a[case+'_maps'][2] if values is None else values
    ref=reference or PRIMARY_REF;dx=np.diff(x).mean();dz=np.diff(z).mean()
    data=20*np.log10(np.maximum(v.reshape(len(z),len(x))/ref,1e-15))
    im=ax.imshow(data,extent=[x[0]-dx/2,x[-1]+dx/2,z[-1]+dz/2,z[0]-dz/2],cmap=CM,vmin=-db,vmax=0,aspect='equal',interpolation='nearest')
    outline(ax,kind);ax.set_xlim(x[0],x[-1]);ax.set_ylim(z[-1],z[0]);ax.set_xticks([-7.5,0,7.5]);ax.set_yticks([17,21,25]);ax.set_xlabel('$x$ (mm)')
    ax.set_title(title or label(path),pad=4)
    if not np.isfinite(v).any():ax.text(.5,.5,'No supported pixels',transform=ax.transAxes,ha='center',va='center',fontsize=7,color='#333333')
    return im

def geometry():
    fig=plt.figure(figsize=(7,5.2),layout='constrained');gs=fig.add_gridspec(2,2,width_ratios=[1.9,1])
    for i,kind in enumerate(['full','partial']):
        ax=fig.add_subplot(gs[i,0]);d=dataset(kind);h=d['material']['thickness_m']*1000
        ax.add_patch(Rectangle((-50,0),100,h,fc='#f4f4f4',ec='#555555',lw=.6));outline(ax,kind,False)
        r=np.array([[r['x_m']*1000,h-r['y_m']*1000] for r in d['audit']['receivers']]);s=np.array([[s['x_mm'],h-s['xy_m'][1]*1000] for s in d['audit']['sources']])
        ax.scatter(*r.T,s=4,color='#666666',zorder=5);ax.scatter(*s.T,s=28,marker='v',c='#d55e00',zorder=6)
        ax.add_patch(Rectangle((-7.5,16.65),15,8.6,fill=False,ec='#0072b2',lw=1,ls='--',zorder=6))
        ax.plot(*tips(kind).T,color='#ba1222',lw=2,zorder=7)
        ax.set(xlim=(-46,46),ylim=(28,-2),xlabel='$x$ (mm)',ylabel='Depth below OD (mm)',title=f'({chr(97+i)}) '+('Full weld: 337 receivers' if kind=='full' else 'Partial weld: 331 receivers'))
    ax=fig.add_subplot(gs[:,1]);kind='partial';ax.set_facecolor('#f4f4f4');outline(ax,kind,False)
    t=tips(kind);ax.plot(*t.T,c='#ba1222',lw=1.6);ax.scatter(*t.T,c='#ba1222',s=16)
    ax.axhline(16.4,c='#555555',ls=':',lw=.7);ax.text(-4.1,16.05,'Fill surface: 16.4 mm',fontsize=7,zorder=9)
    ax.annotate('Upper tip\n18.467 mm',xy=t[0],xytext=(-4.1,18.3),arrowprops={'arrowstyle':'->','lw':.6},fontsize=7)
    ax.annotate('Lower tip\n20.967 mm',xy=t[1],xytext=(-4.1,21.7),arrowprops={'arrowstyle':'->','lw':.6},fontsize=7)
    ax.axhline(25.4,c='black',lw=.8);ax.text(-4.1,25.,'Backwall: 25.4 mm',fontsize=7)
    ax.set(xlim=(-4.6,4.6),ylim=(26,14.5),xlabel='$x$ (mm)',ylabel='Depth below OD (mm)',title='(c) Geometry and known LOF')
    save(fig,'G01_geometry')

def rays():
    paths=[p for p in PATHS if set(p.split('_')[1])=={'S'}]
    fig,axs=plt.subplots(4,2,figsize=(7,7.5),layout='constrained');kind='partial';d=dataset(kind);h=d['material']['thickness_m']*1000
    a=np.array([-1.5,16.4]);r=d['audit']['receivers'][330];b=np.array([r['x_m']*1000,h-r['y_m']*1000]);p=np.array([0.,20.95]);ev=[]
    for ax,path in zip(axs.flat,paths):
        before,after=label(path).split(' | ');left=traces(a,p,before,kind);right=traces(p,b,after,kind,True)
        outline(ax,kind,False)
        for candidates,c,ls in [(left,'#d55e00','-'),(right,'#0072b2','--')]:
            if candidates:
                xy=np.array(candidates[0]['points_mm']);ax.plot(*xy.T,c=c,ls=ls,lw=1.2)
        ax.scatter([a[0],b[0]],[a[1],b[1]],s=20,marker='v',c=['#d55e00','#0072b2'],zorder=6);ax.scatter(*p,s=14,c='black',zorder=6)
        ax.axhline(h,c='#555555',lw=.6);time=.03+left[0]['time_us']+right[0]['time_us'] if left and right else None
        ev.append({'path':path,'time_us':time,'source_mm':a.tolist(),'receiver_mm':b.tolist(),'point_mm':p.tolist()})
        ax.set(xlim=(-5,5),ylim=(26,15.7),xlabel='$x$ (mm)',ylabel='Depth (mm)',title=label(path)+(f': {time:.2f} µs' if time else ': no route'))
    save(fig,'G02_rays');return ev

def physics():
    d=dataset('partial');mat=d['material'];curves=directivity(('broad',.5,5),mat['cl_m_s'],mat['cs_m_s'])
    fig,ax=plt.subplots(1,3,figsize=(7,2.6),layout='constrained')
    for j,field in enumerate(['source','vy']):
        for mode,c,ls in [('L',COLORS[0],'-'),('S',COLORS[1],'--')]:ax[j].plot(ANGLES,curves[mode][field],c=c,ls=ls,label=mode)
        ax[j].axhline(.2,c='#444444',ls=':',lw=.8);ax[j].set(xlim=(0,80),ylim=(0,1.05),xlabel='Angle from normal (°)',title=['(a) Dipole launch','(b) Normal reception'][j]);ax[j].grid(alpha=.15)
    ax[0].set_ylabel('Separately normalized magnitude');ax[0].legend(frameon=False)
    x=np.linspace(.25,9,200)
    for m,c in [('L',COLORS[0]),('S',COLORS[1])]:
        speed=mat['cl_m_s' if m=='L' else 'cs_m_s']/1000
        ax[2].plot(x,.03+2*x/speed,c=c,label=m)
    ax[2].axhspan(0,3.28,color='#dddddd');ax[2].axhline(3.28,c='#555555',ls=':',lw=.8)
    ax[2].set(xlim=(0,9),ylim=(0,6),xlabel='Distance below fill (mm)',ylabel='Round-trip time (µs)',title='(c) Vertical timing example')
    save(fig,'G03_physics')

def atlas(kind,case='lof'):
    ref=PRIMARY_REF if case=='lof' else CONTROL_REF
    for part in range(2):
        fig,axs=plt.subplots(3,3,figsize=(7,5.8),layout='constrained')
        for i,(ax,path) in enumerate(zip(axs.flat,PATHS[part*9:(part+1)*9])):
            im=image(ax,kind,path,case,reference=ref,title=f'{part*9+i+1}. '+label(path));m=result(kind,path)['signals'][case]
            pk=m['peak']['um_s'] if m['peak'] else None;txt=f'Peak {pk:.3g} µm/s' if pk is not None else ''
            ax.text(.025,.04,txt,transform=ax.transAxes,color='white',fontsize=6.5,bbox={'fc':'black','ec':'none','alpha':.65,'pad':1})
            if i%3==0:ax.set_ylabel('Depth (mm)')
            else:ax.set_yticklabels([])
            if i<6:ax.set_xlabel('');ax.set_xticklabels([])
        fig.colorbar(im,ax=axs,orientation='horizontal',shrink=.85,pad=.02,aspect=42,label=f'Amplitude score (dB); 0 dB = {ref*1e6:.3f} µm/s',ticks=[-60,-40,-20,0])
        save(fig,f'{kind}_{case}_atlas_{part+1}')

def comparison():
    fig,axs=plt.subplots(1,3,figsize=(7,6),layout='constrained',gridspec_kw={'width_ratios':[1,1,1]})
    for j,(field,title,vmin,vmax,cmap) in enumerate([('c99_db','(a) $C_{99}$ (dB)',-30,30,'RdBu'),('near_support_pct','(b) Near-flaw support (%)',0,100,'cividis'),('distance','(c) Global peak distance (mm)',0,10,'magma_r')]):
        values=np.full((18,2),np.nan)
        for i,path in enumerate(PATHS):
            for k,kind in enumerate(['full','partial']):
                m=result(kind,path)['signals']['lof'];v=m['peak']['distance_mm'] if field=='distance' and m['peak'] else m.get(field)
                if v is not None:values[i,k]=v
        cm=plt.get_cmap(cmap).copy();cm.set_bad('#eeeeee');im=axs[j].imshow(values,cmap=cm,vmin=vmin,vmax=vmax,aspect='auto')
        axs[j].set_xticks([0,1],['Full','Partial']);axs[j].set_yticks(range(18),[label(p) for p in PATHS] if j==0 else ['']*18);axs[j].set_title(title)
        for i in range(18):
            for k in range(2):
                v=values[i,k];txt='—' if np.isnan(v) else f'{v:.0f}' if field=='near_support_pct' else f'{v:.1f}'
                white=(j==0 and abs(v)>15) or(j==1 and v<50) or(j==2 and v>5)
                axs[j].text(k,i,txt,ha='center',va='center',fontsize=7,c='white' if white else 'black')
        fig.colorbar(im,ax=axs[j],orientation='horizontal',pad=.04,shrink=.85)
    save(fig,'G04_comparison')

def source_support():
    fig,axs=plt.subplots(1,2,figsize=(7,6),layout='constrained')
    for ax,kind in zip(axs,['full','partial']):
        values=np.array([[v if v is not None else np.nan for v in result(kind,p)['source_qualified_fraction_pct']] for p in PATHS])
        cm=plt.get_cmap('cividis').copy();cm.set_bad('#eeeeee');im=ax.imshow(values,cmap=cm,vmin=0,vmax=100,aspect='auto')
        ax.set_xticks(range(5),['S0','S1','S2','S3','S4']);ax.set_yticks(range(18),[label(p) for p in PATHS] if kind=='full' else ['']*18);ax.set_title(kind.capitalize()+' weld');ax.set_xlabel('Source')
        for i in range(18):
            for j in range(5):
                v=values[i,j];ax.text(j,i,'—' if np.isnan(v) else f'{v:.0f}',ha='center',va='center',fontsize=6.5,c='black' if np.isnan(v) or v>55 else 'white')
    fig.colorbar(im,ax=axs,orientation='horizontal',shrink=.7,label='Qualifies at supported image pixels (%)',pad=.03)
    save(fig,'G05_source_support')

def receiver_usage():
    paths=[p for p in PATHS if set(p.split('_')[1])=={'S'}];fig,axs=plt.subplots(4,2,figsize=(7,6),layout='constrained')
    for i,(ax,path) in enumerate(zip(axs.flat,paths)):
        a=load('partial',path);n=np.isfinite(a['lof_maps'][2]).sum();usage=a['lof_receiver_usage'][2]*100/max(n,1)
        im=ax.imshow(usage,cmap='cividis',vmin=0,vmax=100,extent=[-.5,330.5,4.5,-.5],aspect='auto',interpolation='nearest')
        for v in [156.5,313.5]:ax.axvline(v,c='white',lw=.7,ls='--')
        ax.set_title(label(path));ax.set_xticks([0,157,330]);ax.set_yticks(range(5),['S0','S1','S2','S3','S4']);ax.set_xlabel('Receiver ID (left / right / fill)')
        if i%2==0:ax.set_ylabel('Source')
    fig.colorbar(im,ax=axs,orientation='horizontal',shrink=.8,pad=.03,label='Used at supported image pixels (%)')
    save(fig,'G06_receiver_usage')

def sensitivity():
    choices=[('no_coherence','CF off'),('screen_off','Direct off'),('screen_ssss','+ SSSS'),('screen_reflected','+ All echoes'),('gate30','Gate 30%'),('three_sources','3 sources')]
    fig,axs=plt.subplots(1,2,figsize=(7,6.7),layout='constrained')
    for j,field in enumerate(['delta','coverage']):
        values=np.full((18,len(choices)),np.nan)
        for i,path in enumerate(PATHS):
            base=result('partial',path,'coarse')['signals']['lof']['c99_db']
            for k,(name,_) in enumerate(choices):
                m=result('partial',path,name)['signals']['lof']
                v=(m['c99_db']-base) if field=='delta' and m['c99_db'] is not None and base is not None else m['near_support_pct'] if field=='coverage' else None
                if v is not None:values[i,k]=v
        cm=plt.get_cmap('RdBu' if j==0 else 'cividis').copy();cm.set_bad('#eeeeee');im=axs[j].imshow(values,cmap=cm,vmin=-20 if j==0 else 0,vmax=20 if j==0 else 100,aspect='auto')
        axs[j].set_xticks(range(len(choices)),[q[1] for q in choices],rotation=60,ha='right');axs[j].set_yticks(range(18),[label(p) for p in PATHS] if j==0 else ['']*18);axs[j].set_title(['(a) Change in $C_{99}$ (dB)','(b) Near-flaw support (%)'][j])
        for i in range(18):
            for k in range(len(choices)):
                v=values[i,k];txt='—' if np.isnan(v) else f'{v:+.0f}' if j==0 else f'{v:.0f}'
                c='black' if np.isnan(v) or(j==0 and abs(v)<10) or(j==1 and v>55) else 'white'
                axs[j].text(k,i,txt,ha='center',va='center',fontsize=6,c=c)
        fig.colorbar(im,ax=axs[j],orientation='horizontal',shrink=.85,pad=.03)
    save(fig,'G07_sensitivity')

def controls():
    paths=['1-1_SS','1-1_LL','1-2_SSS','1-3_LLLL','2-2_SSSS','3-2_SSSSS']
    fig,axs=plt.subplots(6,3,figsize=(7,9),layout='constrained')
    for i,path in enumerate(paths):
        for j,case in enumerate(['lof','reference','difference']):
            im=image(axs[i,j],'partial',path,case,reference=CONTROL_REF,title=['LOF only','Reference','LOF − reference'][j] if i==0 else ' ')
            if j==0:axs[i,j].set_ylabel(label(path)+'\nDepth (mm)')
            else:axs[i,j].set_yticklabels([])
            if i<5:axs[i,j].set_xticklabels([]);axs[i,j].set_xlabel('')
    fig.colorbar(im,ax=axs,orientation='horizontal',shrink=.85,pad=.02,label=f'Amplitude score (dB); 0 dB = {CONTROL_REF*1e6:.3f} µm/s',ticks=[-60,-40,-20,0])
    save(fig,'G08_controls')

def wide():
    chosen=['1-1_SS','1-1_LL','1-2_SSS','2-1_SSS','2-2_SSSS','3-2_SSSSS']
    ref=max(result(kind,p,'overview')['signals']['lof']['peak']['um_s'] for kind in DATA for p in WIDE if result(kind,p,'overview')['signals']['lof']['peak'])/1e6
    for part in range(2):
        fig,axs=plt.subplots(3,2,figsize=(7,5.3),layout='constrained')
        for i,path in enumerate(chosen[part*3:part*3+3]):
            for j,kind in enumerate(['full','partial']):
                ax=axs[i,j];im=image(ax,kind,path,name='overview',reference=ref,title=kind.capitalize()+': '+label(path))
                ax.set_xticks([-40,0,40]);ax.set_yticks([0,10,20]);ax.set_ylabel('Depth (mm)' if j==0 else '')
                ax.add_patch(Rectangle((-7.5,16.65),15,8.6,fill=False,ec='white',lw=.6,ls='--',zorder=7))
        fig.colorbar(im,ax=axs,orientation='horizontal',pad=.02,shrink=.85,label=f'Amplitude score (dB); 0 dB = {ref*1e6:.2f} µm/s',ticks=[-60,-40,-20,0])
        save(fig,f'G09_wide_{part+1}')
    return ref

def signals():
    kind='partial';d=dataset(kind);nr=len(d['audit']['receivers']);records={}
    for case in ['lof','reference']:
        raw=np.empty((nr,5001))
        for group in ['left','weld','right']:
            p=DATA[kind]/'sources/bead_left'/case/'receivers'/f'{group}_vy.csv'
            ids=[int(n[1:]) for n in p.open().readline().strip().split(',')[1:]];raw[ids]=np.loadtxt(p,delimiter=',',skiprows=1)[:,1:].T
        records[case]=raw
    filt={case:filter_record(raw,kernel(.5,5)) for case,raw in records.items()};filt['difference']=filt['lof']-filt['reference']
    r=330;x=d['audit']['receivers'][r]['x_m']*1000;t=np.arange(5001)*.01
    fig=plt.figure(figsize=(7,6.2),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1,1,1.35])
    ax=fig.add_subplot(gs[0,:]);
    for case,color,ls in [('lof','#111111','-'),('reference',COLORS[1],'--')]:ax.plot(t,filt[case][r].real*1e6,c=color,ls=ls,lw=.7,label=case.capitalize())
    ax.axvspan(0,3.28,color='#dddddd');ax.set(xlim=(0,18),ylabel='$v_y$ (µm/s)',title=f'(a) S1 at $x=-1.5$ mm → R330 at $x={x:.2f}$ mm; both on fill surface');ax.legend(frameon=False,ncol=2)
    a=[-1.5,16.4];b=[x,16.4];centres={}
    for modes,c in [('LL',COLORS[0]),('SS',COLORS[2]),('SSSS',COLORS[3])]:
        cands=traces(a,b,modes,kind);centres[modes]=[.03+q['time_us'] for q in cands]
        for q,tt in enumerate(centres[modes]):
            ax.axvline(tt,c=c,ls=':',lw=.9);ax.text(tt+.08,.97,modes,rotation=90,c=c,transform=ax.get_xaxis_transform(),va='top',fontsize=7)
    ax=fig.add_subplot(gs[1,0])
    for case,c,ls in [('lof','#222222','-'),('reference',COLORS[1],'--'),('difference',COLORS[0],':')]:ax.plot(t,abs(filt[case][r])*1e6,c=c,ls=ls,lw=.8,label=case)
    ax.axvspan(0,3.28,color='#dddddd');ax.set(xlim=(0,18),xlabel='Time (µs)',ylabel='Envelope (µm/s)',title='(b) Waveform controls');ax.legend(frameon=False,fontsize=6)
    ax=fig.add_subplot(gs[1,1]);f=np.fft.rfftfreq(8192,1e-8)/1e6;y=abs(np.fft.rfft(records['lof'][r]*np.hanning(5001),8192));y/=y.max()
    ax.plot(f,20*np.log10(np.maximum(y,1e-9)),c='#555555',lw=.7);ax.axvspan(.5,5,color='#dddddd',zorder=-1);ax.set(xlim=(0,8),ylim=(-100,0),xlabel='Frequency (MHz)',ylabel='dB below peak',title='(c) Raw-record Hann spectrum')
    ax=fig.add_subplot(gs[2,:]);rpos=np.array([v['x_m']*1000 for v in d['audit']['receivers']]);envelope=abs(filt['lof']);ref=envelope.max()
    for group in ['left','weld','right']:
        ids=[i for i,r0 in enumerate(d['audit']['receivers']) if r0['group']==group]
        im=ax.pcolormesh(rpos[ids],t,20*np.log10(np.maximum(envelope[ids].T/ref,1e-12)),cmap='magma',vmin=-60,vmax=0,shading='nearest',rasterized=True)
    ax.set(xlim=(-44,44),ylim=(30,0),xlabel='Receiver $x$ (mm)',ylabel='Time (µs)',title='(d) S1 B-scan; centre group is 16.4 mm below plate receivers')
    fig.colorbar(im,ax=ax,pad=.02,label='Envelope (dB)',ticks=[-60,-40,-20,0]);save(fig,'G10_signals')
    return {'source':1,'receiver':330,'predicted_specimen_return_times_us':centres,'bscan_reference_um_s':float(ref*1e6)}

def partial_wide():
    ref=max(result(kind,p,'overview')['signals']['lof']['peak']['um_s'] for kind in DATA for p in WIDE if result(kind,p,'overview')['signals']['lof']['peak'])/1e6
    fig,axs=plt.subplots(3,3,figsize=(7,4.6),layout='constrained')
    for i,(ax,path) in enumerate(zip(axs.flat,WIDE)):
        im=image(ax,'partial',path,name='overview',reference=ref)
        ax.set_xticks([-40,0,40]);ax.set_yticks([0,10,20]);ax.set_ylabel('Depth (mm)' if i%3==0 else '')
        ax.add_patch(Rectangle((-7.5,16.65),15,8.6,fill=False,ec='white',lw=.6,ls='--',zorder=7))
    fig.colorbar(im,ax=axs,orientation='horizontal',pad=.02,shrink=.85,label=f'Amplitude score (dB); 0 dB = {ref*1e6:.2f} µm/s',ticks=[-60,-40,-20,0]);save(fig,'G12_partial_wide')

def profiles_figure():
    paths=CORE+['2-2_SSSS','3-2_SSSSS'];fig,axs=plt.subplots(3,2,figsize=(7,6),layout='constrained')
    z=np.linspace(17.5,22,901);tt=tips('partial');x=tt[0,0]+(z-tt[0,1])*(tt[1,0]-tt[0,0])/(tt[1,1]-tt[0,1])
    for ax,path in zip(axs.flat,paths):
        for kind,c,ls in [('full',COLORS[0],'-'),('partial',COLORS[1],'--')]:
            a=load(kind,path);v=map_coordinates(a['lof_maps'][2].reshape(len(a['depth_mm']),-1),[(z-a['depth_mm'][0])/.1,(x-a['x_mm'][0])/.1],order=1,mode='nearest')
            ax.plot(z,v*1e6,c=c,ls=ls,label=kind.capitalize())
        ax.axvspan(tt[0,1],tt[1,1],color='#dddddd',zorder=-1,alpha=.5);ax.set(xlim=(17.5,22),xlabel='Depth along extended LOF line (mm)',ylabel='Score (µm/s)',title=label(path));ax.grid(alpha=.15)
    axs[0,0].legend(frameon=False);save(fig,'G11_profiles')

if __name__=='__main__':
    geometry();ray_evidence=rays();physics()
    for kind in DATA:
        for case in ['lof','reference','difference']:atlas(kind,case)
    comparison();source_support();receiver_usage();sensitivity();controls();wide_ref=wide();signal_evidence=signals();profiles_figure();partial_wide()
    write(OUT/'figure_evidence.json',{'primary_reference_um_s':PRIMARY_REF*1e6,'control_reference_um_s':CONTROL_REF*1e6,'wide_reference_um_s':wide_ref*1e6,'ray_examples':ray_evidence,'signals':signal_evidence})
    print('Figures generated',flush=True)
