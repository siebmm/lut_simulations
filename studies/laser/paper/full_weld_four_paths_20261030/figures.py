"""Publication figures from exported float64 reconstructions, not GUI screenshots."""
from analyse import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.patches import Rectangle
from matplotlib.colors import BoundaryNorm
from scipy.ndimage import map_coordinates
sys.path.insert(0,str(ROOT/'studies/laser/src'))
from saft_apertures import directivity,ANGLES
from saft_path_physics import branch
from saft_imaging import kernel,filter_record

FIG=HERE/'figures'; FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,
 'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'ps.fonttype':42,
 'axes.linewidth':.65,'lines.linewidth':1.1,'savefig.dpi':600})
index=json.loads((OUT/'index.json').read_text())
analysis=json.loads((OUT/'analysis.json').read_text())
dataset=json.loads((DATA/'dataset.json').read_text())
mat=dataset['material']; H=mat['thickness_m']*1000
S=np.array([s['x_mm'] for s in dataset['audit']['sources']]); R=np.array([r['x_m']*1000 for r in dataset['audit']['receivers']])
TIPS=np.array([[2.83187473735,18.46747505998],[2.70085528914,20.96747505998]])
P=TIPS.mean(axis=0)
palette=['#0072b2','#d55e00','#009e73','#aa4499']
cmap=plt.get_cmap('magma').copy();cmap.set_bad('#dddddd')
def row(k,name='baseline'):return next(r for r in index if r['path']==PATHS[k] and r['name']==name)
def arrays(k,name='baseline'):return load(row(k,name))
def save(fig,name):
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',pad_inches=.04)
    fig.savefig(FIG/(name+'.png'),bbox_inches='tight',pad_inches=.04,dpi=600)
    plt.close(fig)
def outline(ax,color='#a0a0a0',lw=.65):
    w=np.array(dataset['geometry']['weld'])*1000
    ax.plot(w[:,0],H-w[:,1],color=color,lw=lw,zorder=3)
def flaw(ax,markers=True):
    line,=ax.plot(*TIPS.T,color='#7ff9ff',lw=1.1,ls=(0,(3,1.5)),zorder=5)
    line.set_path_effects([pe.Stroke(linewidth=2.,foreground='#111111'),pe.Normal()])
    if markers:ax.scatter(*TIPS.T,s=11,facecolors='none',edgecolors='#7ff9ff',linewidths=.7,zorder=6)
def img(ax,a,values,reference,db_range=40,zoom=False,title=None):
    x,z=a['x_mm'],a['depth_mm']; dx=np.diff(x).mean(); dz=np.diff(z).mean()
    values=values.reshape(len(z),len(x)); v=20*np.log10(np.maximum(values/reference,1e-15))
    im=ax.imshow(v,extent=[x[0]-dx/2,x[-1]+dx/2,z[-1]+dz/2,z[0]-dz/2],cmap=cmap,vmin=-db_range,vmax=0,interpolation='nearest',aspect='equal')
    outline(ax);flaw(ax)
    if zoom:ax.set_xlim(1.1,4.5);ax.set_ylim(22,17.4);ax.set_xticks([1.5,2.5,3.5,4.5]);ax.set_yticks([18,19,20,21,22])
    else:ax.set_xlim(x[0],x[-1]);ax.set_ylim(z[-1],z[0]);ax.set_xticks([-6,0,6]);ax.set_yticks([5,10,15,20,25])
    ax.set_xlabel('$x$ (mm)')
    if title:ax.set_title(title,pad=5)
    return im

def geometry():
    fig,ax=plt.subplots(1,2,figsize=(7,3.15),gridspec_kw={'width_ratios':[2.15,1]},layout='constrained')
    a=ax[0];a.add_patch(Rectangle((-75,0),150,H,fc='#f5f5f5',ec='black',lw=.8))
    for side in [-75,60]:a.add_patch(Rectangle((side,0),15,H,fc='#dddddd',hatch='////',ec='#b5b5b5',lw=.4))
    outline(a,'#666666');a.scatter(R,np.zeros(len(R)),s=5,color='#777777',label='337 receivers',zorder=4)
    a.scatter(S,np.zeros(5),s=30,marker='v',color='#d55e00',zorder=5,label='Five sources')
    a.add_patch(Rectangle((-6,5),12,20.2,fill=False,ec='#0072b2',lw=1,ls='--'))
    a.annotate('Saved SAFT grid',xy=(-6,9),xytext=(-50,9),arrowprops={'arrowstyle':'->','lw':.7},fontsize=8)
    a.plot(*TIPS.T,color='#c00000',lw=2)
    a.set(xlim=(-77,77),ylim=(31,-6),xlabel='$x$ (mm)',ylabel='Depth $z$ (mm)',title='(a) Full-weld model and acquisition')
    a.legend(loc='lower center',ncol=2,frameon=False,bbox_to_anchor=(.5,-.35))
    a.text(0,29,'Sources: $x=-12,-1.5,0,1.5,12$ mm',ha='center',fontsize=7)
    a=ax[1];outline(a,'#777777'); a.plot(*TIPS.T,color='#c00000',lw=1.8);a.scatter(*TIPS.T,color='#c00000',s=15)
    a.annotate('Upper tip\n(2.832, 18.467)',xy=TIPS[0],xytext=(.05,17.6),arrowprops={'arrowstyle':'->','lw':.6},fontsize=7)
    a.annotate('Lower tip\n(2.701, 20.967)',xy=TIPS[1],xytext=(.05,22.2),arrowprops={'arrowstyle':'->','lw':.6},fontsize=7)
    a.text(3.12,19.6,'Known LOF\n2.50 mm height',fontsize=7,va='center')
    a.set(xlim=(0,5.8),ylim=(24,16.65),xlabel='$x$ (mm)',ylabel='Depth $z$ (mm)',title='(b) Right fusion boundary')
    save(fig,'Fig01_geometry')

def rays():
    fig,axs=plt.subplots(2,2,figsize=(7,5),layout='constrained'); times=[]
    for k,ax in enumerate(axs.flat):
        before,after=[('S','S'),('L','L'),('S','SS'),('SS','S')][k]
        inc=branch([S[0]/1000,0],P/1000,before,H/1000,mat['cl_m_s'],mat['cs_m_s'],keep_points=True)
        rec=branch(P/1000,[R[94]/1000,0],after,H/1000,mat['cl_m_s'],mat['cs_m_s'],outgoing=True,keep_points=True)
        for b,c,style,lab in [(inc,'#d55e00','-','Source to pixel'),(rec,'#0072b2','--','Pixel to receiver')]:
            pts=np.array(b['points'])*1000;ax.plot(*pts.T,color=c,ls=style,lw=1.4,label=lab)
            for q in range(len(pts)-1):
                mid=(pts[q]+pts[q+1])/2; vec=(pts[q+1]-pts[q])*.12
                ax.annotate('',xy=mid+vec,xytext=mid-vec,arrowprops={'arrowstyle':'->','color':c,'lw':1.})
        t=.03+(inc['time']+rec['time'])*1e6;times.append(float(t))
        ax.axhline(0,c='#555555',lw=.8);ax.axhline(H,c='#555555',lw=.8);outline(ax,'#cccccc');ax.plot(*TIPS.T,c='black',lw=1.4)
        ax.scatter([S[0],R[94]],[0,0],marker='v',s=25,c=['#d55e00','#0072b2'],zorder=5)
        ax.scatter(*P,s=18,c='black',zorder=5)
        ax.set(xlim=(-25,7),ylim=(27,-2),xlabel='$x$ (mm)',ylabel='Depth (mm)',title=f'({chr(97+k)}) {LABELS[k]}: {t:.2f} µs')
        ax.set_aspect('equal');ax.set_yticks([0,10,20,25.4])
    fig.legend(*axs[0,0].get_legend_handles_labels(),ncol=2,loc='outside lower center',frameon=False)
    save(fig,'Fig02_paths');return times

def physics():
    curves=directivity(('broad',.5,5),mat['cl_m_s'],mat['cs_m_s'])
    fig,axs=plt.subplots(1,3,figsize=(7,2.65),layout='constrained')
    for i,field in enumerate(['source','vy']):
        for mode,c,ls in [('L',palette[0],'-'),('S',palette[1],'--')]:axs[i].plot(ANGLES,curves[mode][field],c=c,ls=ls,label=mode)
        axs[i].axhline(.2,c='#555555',ls=':',label='20% gate');axs[i].set(xlim=(0,80),ylim=(0,1.05),xlabel='Angle from normal (°)',title=['(a) Dipole launch','(b) Normal-velocity reception'][i])
        axs[i].set_xticks([0,20,40,60,80]);axs[i].grid(alpha=.15)
    axs[0].set_ylabel('Separately normalized magnitude');axs[0].legend(frameon=False)
    n=32768;h=kernel(.5,5);f=np.fft.fftfreq(n,1e-8)/1e6;mag=abs(np.fft.fft(h,n))/2
    m=(f>=0)&(f<=7);axs[2].plot(f[m],20*np.log10(np.maximum(mag[m],1e-8)),c='#222222')
    axs[2].axvspan(.5,5,color='#dddddd',alpha=.5,zorder=-1);axs[2].set(xlim=(0,7),ylim=(-80,3),xlabel='Frequency (MHz)',ylabel='Filter gain (dB)',title='(c) Quadrature FIR')
    axs[2].grid(alpha=.15); save(fig,'Fig03_response')
    return curves

def mainimages():
    fig,axs=plt.subplots(1,4,figsize=(7,4.1),layout='constrained'); ref=max(np.nanmax(arrays(k)['lof_maps'][2]) for k in range(4))
    for k,ax in enumerate(axs):
        a=arrays(k);im=img(ax,a,a['lof_maps'][2],ref,title=f'({chr(97+k)}) {LABELS[k]}')
        ax.text(.04,.025,f'Peak {np.nanmax(a["lof_maps"][2])*1e6:.2f}',transform=ax.transAxes,c='white',fontsize=7,bbox={'fc':'black','ec':'none','alpha':.6,'pad':1.5})
        if k==0:ax.set_ylabel('Depth (mm)')
        else:ax.set_yticklabels([])
    fig.colorbar(im,ax=axs,shrink=.75,pad=.02,label=f'Amplitude score (dB); 0 dB = {ref*1e6:.3f} µm/s',ticks=[-40,-30,-20,-10,0])
    save(fig,'Fig04_lof_grid')
    # Deliberately common physical reference across all twelve panels.
    ref=max(np.nanmax(arrays(k)[c+'_maps'][2]) for k in range(4) for c in ['lof','reference','difference'])
    fig,axs=plt.subplots(3,4,figsize=(7,7.6),layout='constrained')
    for j,(case,label) in enumerate([('lof','LOF only'),('reference','Reference'),('difference','LOF − reference')]):
        for k in range(4):
            a=arrays(k);im=img(axs[j,k],a,a[case+'_maps'][2],ref,zoom=True,title=LABELS[k] if j==0 else None)
            if k==0:axs[j,k].set_ylabel(label+'\nDepth (mm)')
            else:axs[j,k].set_yticklabels([])
            if j<2:axs[j,k].set_xlabel('');axs[j,k].set_xticklabels([])
    fig.colorbar(im,ax=axs,shrink=.7,pad=.025,label=f'Amplitude score (dB); 0 dB = {ref*1e6:.3f} µm/s',ticks=[-40,-30,-20,-10,0])
    save(fig,'Fig05_reference_control')
    fig,axs=plt.subplots(4,1,figsize=(7,8.2),layout='constrained')
    ref=max(np.nanmax(arrays(k,'overview')['lof_maps'][2]) for k in range(4))
    for k,ax in enumerate(axs.flat):
        a=arrays(k,'overview');im=img(ax,a,a['lof_maps'][2],ref,db_range=60,title=f'({chr(97+k)}) {LABELS[k]}')
        ax.set_xticks([-40,-20,0,20,40]);ax.set_yticks([0,10,20]);ax.add_patch(Rectangle((-6,5),12,20.2,fill=False,ec='white',ls='--',lw=.65,zorder=4))
        ax.text(.02,.06,f'Peak {np.nanmax(a["lof_maps"][2])*1e6:.1f} µm/s',transform=ax.transAxes,c='white',fontsize=7,bbox={'fc':'black','ec':'none','alpha':.6,'pad':1.5})
        ax.set_ylabel('Depth (mm)')
    fig.colorbar(im,ax=axs,shrink=.8,label=f'Amplitude score (dB); 0 dB = {ref*1e6:.2f} µm/s',ticks=[-60,-40,-20,0])
    save(fig,'Fig06_extended_field')

def ablation():
    fig,axs=plt.subplots(2,3,figsize=(7,7),layout='constrained')
    for k in range(2):
        a=arrays(k,'no_screen');b=arrays(k,'coarse');data=[a['lof_maps'][1],a['lof_maps'][2],b['lof_maps'][2]]
        ref=np.nanmax(data[0])
        for j in range(3):
            im=img(axs[k,j],a,data[j],ref,db_range=60,title=['Filtered mean','+ Coherence','+ Pulse screen'][j] if k==0 else None)
            if j==0:axs[k,j].set_ylabel(LABELS[k]+'\nDepth (mm)')
            else:axs[k,j].set_yticklabels([])
        fig.colorbar(im,ax=axs[k,:],shrink=.8,pad=.03,label=f'dB; 0 dB = {ref*1e6:.2f} µm/s',ticks=[-60,-40,-20,0])
    save(fig,'Fig07_background_controls')

def coverage():
    fig,axs=plt.subplots(2,4,figsize=(7,5.4),layout='constrained',gridspec_kw={'height_ratios':[3,1.1]})
    cm=plt.get_cmap('viridis',4).copy();cm.set_bad('#dddddd')
    for k in range(4):
        a=arrays(k);counts=np.where(np.isfinite(a['lof_maps'][2]),a['lof_counts'][2],np.nan).reshape(len(a['depth_mm']),-1)
        im=axs[0,k].imshow(counts,cmap=cm,norm=BoundaryNorm([1.5,2.5,3.5,4.5,5.5],4),extent=[-6.05,6.05,25.25,4.95],aspect='equal',interpolation='nearest')
        outline(axs[0,k]);flaw(axs[0,k]);axs[0,k].set_title(LABELS[k]);axs[0,k].set_xticks([-6,0,6]);axs[0,k].set_xlabel('$x$ (mm)');axs[0,k].set_xlim(-6,6);axs[0,k].set_ylim(25.2,5)
        usage=a['lof_receiver_usage'][2]*100/np.isfinite(a['lof_maps'][2]).sum()
        iu=axs[1,k].imshow(usage,cmap='cividis',vmin=0,vmax=100,extent=[-.5,336.5,4.5,-.5],aspect='auto',interpolation='nearest')
        axs[1,k].set_xticks([0,168,336]);axs[1,k].set_xlabel('Receiver ID');axs[1,k].set_yticks(range(5));axs[1,k].set_yticklabels(['S0','S1','S2','S3','S4'])
        if k==0:axs[0,k].set_ylabel('Depth (mm)');axs[1,k].set_ylabel('Source')
        else:axs[0,k].set_yticklabels([]);axs[1,k].set_yticklabels([])
    fig.colorbar(im,ax=axs[0],shrink=.8,label='Qualifying sources',ticks=[2,3,4,5])
    fig.colorbar(iu,ax=axs[1],label='Used pixels (%)',ticks=[0,50,100])
    save(fig,'Fig08_coverage')

def sensitivity():
    names=['no_coherence','no_screen','screen_LSR','exclude_repairs','gate10','gate30','offset0','offset4','wall0','wall10','band03','band15','sources3','receivers16','fraction50','tail20','tail40','margin0','margin02','weighted','shoulders']
    labels=['Coherence off','Pulse screen off','Screen L + S + R','Exclude repairs','Angular gate 10%','Angular gate 30%','Exclusion 0 mm','Exclusion 4 mm','Wall gate 0%','Wall gate 10%','Band 0.5–3 MHz','Band 1–5 MHz','Minimum 3 sources','Minimum 16 receivers','Eligible fraction 50%','Pulse tail 20 dB','Pulse tail 40 dB','Timing margin 0 µs','Timing margin 0.2 µs','Response weighting','Shoulder sources only']
    delta=np.full((len(names),4),np.nan);support=delta.copy();textvalues={}
    for k,path in enumerate(PATHS):
        base=next(r for r in analysis['runs'] if r['path']==path and r['name']=='coarse')['signals']['lof']
        for j,name in enumerate(names):
            r=next((r for r in analysis['runs'] if r['path']==path and r['name']==name),None)
            if r:
                m=r['signals']['lof'];support[j,k]=m['support_pct']
                if m['contrast_p99_db'] is not None:delta[j,k]=m['contrast_p99_db']-base['contrast_p99_db']
    fig,axs=plt.subplots(1,2,figsize=(7,6.8),layout='constrained',gridspec_kw={'width_ratios':[1,1]})
    for q,(v,title,cmap0,vmin,vmax) in enumerate([(delta,'(a) Change in $C_{99}$ (dB)','RdBu',-15,15),(support,'(b) Supported grid (%)','cividis',0,100)]):
        cm=plt.get_cmap(cmap0).copy();cm.set_bad('#eeeeee');im=axs[q].imshow(v,cmap=cm,vmin=vmin,vmax=vmax,aspect='auto')
        axs[q].set_xticks(range(4),LABELS,rotation=0,fontsize=7);axs[q].set_yticks(range(len(names)),labels if q==0 else ['']*len(names));axs[q].set_title(title)
        for j in range(len(names)):
            for k in range(4):
                val=v[j,k];label='—' if np.isnan(val) else f'{val:+.1f}' if q==0 else f'{val:.0f}'
                color='black' if np.isnan(val) or (q==0 and abs(val)<8) or(q==1 and val>55) else 'white'
                axs[q].text(k,j,label,ha='center',va='center',fontsize=7,color=color)
        fig.colorbar(im,ax=axs[q],orientation='horizontal',shrink=.8,pad=.04)
    save(fig,'Fig09_sensitivity')

def signals(times):
    raw=np.empty((337,5001))
    for group in ['left','weld','right']:
        path=DATA/'sources/shoulder_left/lof/receivers'/f'{group}_vy.csv'
        header=path.open().readline().strip().split(','); ids=[int(s[1:]) for s in header[1:]]
        raw[ids]=np.loadtxt(path,delimiter=',',skiprows=1)[:,1:].T
    filtered=filter_record(raw,kernel(.5,5)); envelope=abs(filtered); ref=envelope.max()
    # Cross-check the independent Python convolution against the C++ trace.
    trace=np.genfromtxt(OUT/'numerical/0/trace.csv',delimiter=',',names=True)
    err=np.max(abs(filtered[94]-(trace['filtered_real_m_s']+1j*trace['filtered_imag_m_s'])))
    fig=plt.figure(figsize=(7,5.7),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.55],width_ratios=[1.6,1])
    ax=fig.add_subplot(gs[0,0]);t=np.arange(5001)*.01
    ax.plot(t,filtered[94].real*1e6,c='#222222',lw=.65,label='Filtered signed $v_y$')
    ax.set(xlim=(0,30),xlabel='Time (µs)',ylabel='$v_y$ (µm/s)',title='(a) S0 → R94 ($x_r=-20.05$ mm)')
    for k,tt in enumerate(times):ax.axvline(tt,c=palette[k],ls=['--',':','-.','-'][k],lw=.9)
    ax.legend(frameon=False,loc='upper right')
    ax=fig.add_subplot(gs[0,1]);w=np.hanning(5001);f=np.fft.rfftfreq(8192,1e-8)/1e6
    y=abs(np.fft.rfft(raw[94]*w,8192)); y=y/y.max()
    ax.plot(f,20*np.log10(np.maximum(y,1e-8)),c='#555555',lw=.7);ax.axvspan(.5,5,color='#dddddd',alpha=.55,zorder=-1)
    ax.set(xlim=(0,8),ylim=(-100,0),xlabel='Frequency (MHz)',ylabel='dB below spectral peak',title='(b) Raw-record Hann spectrum')
    ax=fig.add_subplot(gs[1,:]);vals=20*np.log10(np.maximum(envelope/ref,1e-12))
    for group in ['left','weld','right']:
        ids=np.array([i for i,r in enumerate(dataset['audit']['receivers']) if r['group']==group])
        im=ax.pcolormesh(R[ids],t,vals[ids].T,cmap='magma',vmin=-60,vmax=0,shading='nearest',rasterized=True)
    for k,(before,after) in enumerate([('S','S'),('L','L'),('S','SS'),('SS','S')]):
        inc=branch([S[0]/1000,0],P/1000,before,H/1000,mat['cl_m_s'],mat['cs_m_s'])
        rec=branch(P/1000,np.column_stack([R/1000,np.zeros(337)]),after,H/1000,mat['cl_m_s'],mat['cs_m_s'],outgoing=True)
        ax.plot(R,.03+(inc['time']+rec['time'])*1e6,c=palette[k],ls=['--',':','-.','-'][k],lw=1.2,label=LABELS[k])
    for speed,ls in [(mat['cl_m_s'],':'),(mat['cs_m_s'],'--'),(mat['cr_m_s'],'-')]:
        ax.plot(R,.03+abs(R-S[0])*1000/speed,c='white',ls=ls,lw=.65,alpha=.8)
    ax.set(xlim=(R.min(),R.max()),ylim=(30,0),xlabel='Receiver $x$ (mm)',ylabel='Time (µs)',title='(c) Filtered B-scan; coloured curves focus at the LOF centre')
    ax.legend(loc='lower left',ncol=4,facecolor='white',framealpha=.9)
    fig.colorbar(im,ax=ax,pad=.02,label=f'Envelope (dB)\n0 dB = {ref*1e6:.1f} µm/s',ticks=[-60,-40,-20,0])
    save(fig,'Fig10_signals')
    return {'independent_filter_max_error_m_s':float(err),'trace_source':0,'trace_receiver':94,'ray_times_us':times,'bscan_reference_um_s':float(ref*1e6)}

def profiles():
    fig,axs=plt.subplots(1,2,figsize=(7,3.1),layout='constrained')
    d=np.linspace(17,22.5,1101);xx=TIPS[0,0]+(d-TIPS[0,1])*(TIPS[1,0]-TIPS[0,0])/(TIPS[1,1]-TIPS[0,1])
    for q,case in enumerate(['lof','difference']):
        for k in range(4):
            a=arrays(k);im=a[case+'_maps'][2].reshape(len(a['depth_mm']),-1)
            vals=map_coordinates(im,[(d-5)/.1,(xx+6)/.1],order=1,mode='nearest')
            axs[q].plot(d,vals*1e6,c=palette[k],ls=['-','--','-.',':'][k],label=LABELS[k])
        for tip in TIPS[:,1]:axs[q].axvline(tip,c='#555555',ls=':',lw=.8)
        axs[q].axvspan(TIPS[0,1],TIPS[1,1],color='#dddddd',alpha=.35,zorder=-1)
        axs[q].set(xlim=(17,22.5),xlabel='Depth on extended LOF line (mm)',ylabel='Amplitude score (µm/s)',title=['(a) LOF only','(b) LOF − reference'][q]);axs[q].grid(alpha=.15)
    axs[0].legend(frameon=False,ncol=2);save(fig,'Fig11_profiles')

if __name__=='__main__':
    geometry();times=rays();curves=physics();mainimages();ablation();coverage();sensitivity();evidence=signals(times);profiles()
    (OUT/'figure_evidence.json').write_text(json.dumps(evidence,indent=2));print(json.dumps(evidence,indent=2))
