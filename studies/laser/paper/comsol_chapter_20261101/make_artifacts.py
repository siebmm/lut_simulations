"""Collect numerical evidence (--collect), then regenerate the chapter artwork.

Default plotting uses only the portable data/ folder. No simulations or fitting.
"""

# Locate portable scientific inputs without depending on the source checkout.
import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import data_root as _data_root, native_root as _native_root, executable as _executable, native_env as _native_env
from pathlib import Path
import argparse,copy,hashlib,json,shutil,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects as pe
from matplotlib.lines import Line2D
from figure_layout import legend_strip,verify_layout

HERE=Path(__file__).resolve().parent; DATA=HERE/'data'; FIG=HERE/'figures/comsol'
ROOT=HERE.parents[3]; STUDY=ROOT/'studies/laser'; RESULTS=STUDY/'results/filter_study_20261101'
PARTIAL_RESULTS=STUDY/'results/equal_receiver_study_20261101'
PARTIAL_REFERENCE=40.
CORE=['1-1_SS','1-1_LL','1-2_SSS','2-1_SSS']
CM=plt.get_cmap('magma').copy();CM.set_bad('#cfcfcf')
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.labelsize':8,
 'axes.titlesize':8.5,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,
 'pdf.fonttype':42,'ps.fonttype':42,'axes.linewidth':.6})

def write(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
def label(p):n=int(p.split('-')[0]);m=p.split('_')[1];return m[:n]+' | '+m[n:]

def collect():
    DATA.mkdir(exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
    index=[]
    for stage in ['numerical','confirmation','selected_checks']:
        for kind in ['full','partial']:
            source=PARTIAL_RESULTS if kind=='partial' else RESULTS
            for p in sorted((source/stage/kind).rglob('arrays.npz')):
                job=json.loads((p.parent/'request.json').read_text());config=job['selection']
                spec=job.get('signal_processing',{});uid=f'{stage}_{kind}_'+p.parent.relative_to(source/stage/kind).as_posix().replace('/','__')
                assert config['processing']=='raw', 'Publication images must use equal receiver averaging'
                with np.load(p) as a:
                    shape=(len(a['z']),len(a['x']));fields={'x_mm':a['x'],'depth_mm':a['z']}
                    xx,zz=np.meshgrid(a['x'],a['z'])
                    np.testing.assert_allclose(np.column_stack([xx.ravel(),zz.ravel()])/1000,a['points'].reshape(-1,2),rtol=0,atol=1e-14)
                    for case in ['lof','reference','difference']:
                        for attr in ['amplitude','count','cf','retained','screened']:
                            fields[case+'_'+attr]=a[case+'_2_'+attr].reshape(shape)
                    np.savez_compressed(DATA/(uid+'.npz'),**fields)
                index.append(dict(uid=uid,stage=stage,kind=kind,id=job['id'],path=config['path'],
                                  processing=config['processing'],variant=spec.get('name'),selection=config,signal_processing={k:v for k,v in spec.items() if k!='bank'}))
    write(DATA/'images.json',index)
    for kind in ['full','partial']:shutil.copy2(_data_root()/f'{kind}_weld/dataset.json',DATA/f'{kind}_dataset.json')
    assert len(index)==384 and sum(r['kind']=='partial' for r in index)==280
    for relative,name in [('verification/filter_audit.json','filter_audit.json'),
         ('confirmation/full/metrics.json','full_confirmation_metrics.json'),
         ('selected_checks/full/metrics.json','full_selected_metrics.json')]:
        shutil.copy2(RESULTS/relative,DATA/name)
    survey=json.loads((RESULTS/'metrics.json').read_text())['rows']
    survey=[r for r in survey if r['kind']=='full']+json.loads((PARTIAL_RESULTS/'metrics.json').read_text())['rows']
    write(DATA/'survey_metrics.json',dict(rows=survey))
    shutil.copy2(PARTIAL_RESULTS/'confirmation/partial/metrics.json',DATA/'partial_confirmation_metrics.json')
    shutil.copy2(PARTIAL_RESULTS/'verification/check.json',DATA/'native_parity.json')
    shutil.copy2(PARTIAL_RESULTS/'verification/reproduction.json',HERE/'verification/equal_receiver_reproduction.json')
    for name,p in [('survey_protocol',HERE.parent/'equal_receiver_study_20261101/protocol.json'),
         ('full_confirmation_protocol',RESULTS/'confirmation/full/protocol.json'),('partial_confirmation_protocol',PARTIAL_RESULTS/'confirmation/partial/protocol.json'),
         ('selected_protocol',RESULTS/'selected_checks/full/protocol.json'),('mesh_audit',STUDY/'results/saft_workflow/step_03_bandwidth_audit.json')]:
        shutil.copy2(p,DATA/(name+'.json'))
    shutil.copytree(RESULTS/'filter_evidence',DATA/'filter_evidence',dirs_exist_ok=True)
    shutil.copytree(PARTIAL_RESULTS/'filter_evidence',DATA/'filter_evidence',dirs_exist_ok=True)
    old=HERE.parent/'laser_weld_feasibility_20261030'
    shutil.copy2(old/'data/directivity.npz',DATA/'directivity.npz')
    shutil.copytree(old/'data/source_calibration',DATA/'source_calibration',dirs_exist_ok=True)
    # Existing publication artwork has its own audited generators and layout.
    # Collecting numerical evidence must not restore historical figure versions.
    # Preserve original generation source separately; figures are based on the
    # exact model vertices, never hand-redrawn approximations of the bevel.
    prov={str(p.relative_to(DATA)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(DATA.rglob('*')) if p.is_file() and p.name!='sha256.json'}
    write(DATA/'sha256.json',prov)
    print('Collected',len(index),'image configurations',flush=True)

def row(kind,id,stage='confirmation'):
    return next(r for r in INDEX if r['kind']==kind and r['stage']==stage and r['id']==id)
def arrays(r):return np.load(DATA/(r['uid']+'.npz'))
def tips(kind):
    t=np.array(DATASETS[kind]['audit']['defect_benchmark_only']['tips_xy_mm']);t[:,1]=25.4-t[:,1];return t
def flaw_marker(ax,kind,style='outside'):
    t=tips(kind);t=t[np.argsort(t[:,1])]
    if style=='contour':
        line,=ax.plot(*t.T,c='#d2d2d2',lw=.5,ls=(0,(2,2)),zorder=6)
        line.set_path_effects([pe.Stroke(linewidth=1.05,foreground='#202020'),pe.Normal()])
    elif style=='bracket':
        # Leave the complete flaw trace and hot tip pixels unobscured. Short
        # horizontal leaders terminate 0.15 mm to the right of the true tips.
        x=float(t[:,0].max()+.85)
        line,=ax.plot([t[0,0]+.15,x,x,t[1,0]+.15],[t[0,1],t[0,1],t[1,1],t[1,1]],
                     c='#c4c4c4',lw=.55,zorder=6)
        line.set_path_effects([pe.Stroke(linewidth=1.05,foreground='#222222'),pe.Normal()])
    elif style=='outside':
        transform=ax.get_yaxis_transform()
        ax.plot([1.015,1.04,1.04,1.015],[t[0,1],t[0,1],t[1,1],t[1,1]],
                c='#777777',lw=.6,transform=transform,clip_on=False,zorder=6)
        x=float(t[:,0].mean())
        ax.plot([x,x],[1.005,1.035],c='#777777',lw=.7,
                transform=ax.get_xaxis_transform(),clip_on=False,zorder=6)

def outline(ax,kind,marker='outside'):
    d=DATASETS[kind];v=np.array(d['geometry']['weld'])*1000;v=np.vstack([v,v[0]])
    ax.plot(v[:,0],25.4-v[:,1],c='#9b9b9b',lw=.65,zorder=3)
    if kind=='partial':
        p=np.array(d['surface_profile_x_depth_m'])*1000
        ax.fill_between(p[:,0],-10,p[:,1],color='white',lw=0,zorder=4);ax.plot(*p.T,c='#555555',lw=.6,zorder=4)
    flaw_marker(ax,kind,marker)

def save(fig,name):
    verify_layout(fig,name)
    FIG.mkdir(parents=True,exist_ok=True)
    fig.savefig(FIG/(name+'.pdf'),bbox_inches='tight',pad_inches=.05)
    fig.savefig(FIG/(name+'.png'),dpi=600,bbox_inches='tight',pad_inches=.05)
    plt.close(fig)

def heat(ax,r,case='lof',reference=10.,db=45,title=None,wide=False,marker='outside'):
    with arrays(r) as a:
        x,z=a['x_mm'],a['depth_mm'];v=a[case+'_amplitude']*1e6
        value=20*np.log10(np.maximum(v/reference,1e-15))
        dx=np.diff(x).mean();dz=np.diff(z).mean()
        im=ax.imshow(value,cmap=CM,vmin=-db,vmax=0,extent=[x[0]-dx/2,x[-1]+dx/2,z[-1]+dz/2,z[0]-dz/2],interpolation='nearest',aspect='equal')
        ax.set(xlim=(x[0],x[-1]),ylim=(z[-1],z[0]),xlabel='$x$ (mm)')
    outline(ax,r['kind'],marker);ax.set_title(title or label(r['path']),loc='left',pad=10)
    ax.set_xticks([-40,-20,0,20,40] if wide else [-6,0,6] if r['kind']=='full' else [-7.5,0,7.5])
    ax.set_yticks([0,10,20,25.4] if wide else [5,10,15,20,25] if r['kind']=='full' else [17,21,25]);return im

def scale(fig,axes,im,reference,pad=.04):
    fig.colorbar(im,ax=axes,orientation='horizontal',pad=pad,shrink=.9,aspect=45,ticks=np.arange(-40,1,10),
                 label=f'Image amplitude (dB); 0 dB = {reference:g} '+r'$\mu$m/s')

def primary():
    for kind,id,reference in [('full','all',10.),('partial','baseline__1-2_SSS__fine',PARTIAL_REFERENCE)]:
        r=row(kind,id,'selected_checks' if kind=='full' else 'confirmation')
        fig,axs=plt.subplots(1,3,figsize=(6.9,3.9 if kind=='full' else 2.25),layout='constrained')
        for i,(ax,case) in enumerate(zip(axs,['lof','reference','difference'])):
            title=['LOF only','Pristine control','LOF minus pristine'][i]
            im=heat(ax,r,case,reference,title=f'({chr(97+i)}) {title}');
            if i==0:ax.set_ylabel('Depth below OD (mm)')
            else:ax.set_yticklabels([])
        scale(fig,axs,im,reference);save(fig,'primary_'+kind)

def core_full():
    fig,axs=plt.subplots(1,4,figsize=(6.9,3.7),layout='constrained')
    for i,(ax,p) in enumerate(zip(axs,CORE)):
        im=heat(ax,row('full','baseline__'+p+'__fine'),reference=10,title=f'({chr(97+i)}) '+label(p))
        if i==0:ax.set_ylabel('Depth below OD (mm)')
        else:ax.set_yticklabels([])
    scale(fig,axs,im,10);save(fig,'full_common_paths')

def partial_atlas():
    rr=[r for r in INDEX if r['stage']=='numerical' and r['kind']=='partial' and r['variant']=='baseline']
    order=json.loads((DATA/'survey_protocol.json').read_text())['paths']['partial'];rr.sort(key=lambda r:order.index(r['path']))
    fig,axs=plt.subplots(5,4,figsize=(6.9,6.3),layout='constrained')
    for i,ax in enumerate(axs.flat):
        if i>=len(rr):ax.set_visible(False);continue
        im=heat(ax,rr[i],reference=PARTIAL_REFERENCE,db=55,title=label(rr[i]['path']))
        ax.set_xlabel('');ax.set_xticks([-7.5,0,7.5]);ax.set_yticks([17,25])
        if i%4==0:ax.set_ylabel('Depth (mm)')
        else:ax.set_yticklabels([])
    fig.supxlabel('$x$ (mm)',fontsize=8)
    fig.colorbar(im,ax=axs.ravel().tolist(),orientation='horizontal',pad=.02,shrink=.8,aspect=45,label=rf'Image amplitude (dB); common 0 dB = {PARTIAL_REFERENCE:g} $\mu$m/s',ticks=[-50,-40,-30,-20,-10,0])
    save(fig,'partial_path_atlas')

NAMES={'baseline':'0.5-5 MHz / 657 taps','band1_5':'1-5 MHz / 657 taps','band05_3':'0.5-3 MHz / 657 taps',
 'lowpass5_long':'0-5 MHz / 657 taps','bp_short_matched_guard':'0.5-5 MHz / 257 / long guard',
 'bp_short':'0.5-5 MHz / 257 taps','lp_short_matched_guard':'0-5 MHz / 257 / long guard',
 'lp_short':'0-5 MHz / 257 taps','lp_shorter':'0-5 MHz / 129 taps','baseline_initial_rest':'0.5-5 MHz / initial rest',
 'lp_short_initial_rest':'0-5 MHz / 257 / initial rest','butter_lowpass':'5 MHz Butterworth + analytic',
 'broadband_analytic':'Unrestricted analytic','raw_signed':'Unfiltered signed samples'}

def matrices():
    rows=json.loads((DATA/'survey_metrics.json').read_text())['rows'];protocol=json.loads((DATA/'survey_protocol.json').read_text())
    variants=[v['name'] for v in protocol['variants']]
    for kind in ['full','partial']:
        paths=protocol['paths'][kind]
        subsets=[paths] if kind=='full' else [[p for p in paths if '_S' in p], [p for p in paths if '_L' in p]]
        for group,sub in enumerate(subsets):
            values=np.full((14,len(sub)),np.nan);near=values.copy();peak=values.copy()
            for i,v in enumerate(variants):
                for j,p in enumerate(sub):
                    r=next(r for r in rows if r['kind']==kind and r['variant']==v and r['path']==p and r['case']=='lof')
                    if r['cmax_db'] is not None:values[i,j]=r['cmax_db']
                    near[i,j]=r['near_support_pct'];peak[i,j]=r['peak']['distance_mm'] if r['peak'] else np.nan
            fig,ax=plt.subplots(figsize=(6.9,5.1),layout='constrained')
            cm=plt.get_cmap('RdBu').copy();cm.set_bad('#dedede')
            im=ax.imshow(values,vmin=-15,vmax=15,cmap=cm,aspect='auto')
            ax.set_yticks(range(14),[NAMES[v] for v in variants],fontsize=7)
            ax.set_xticks(range(len(sub)),[label(p) for p in sub],rotation=45,ha='right',fontsize=7)
            for i in range(14):
                for j in range(len(sub)):
                    v=values[i,j];t='--' if not np.isfinite(v) else f'{v:.1f}'
                    if np.isfinite(v) and peak[i,j]<=.5:t+='*'
                    ax.text(j,i,t,ha='center',va='center',fontsize=6.6,color='white' if abs(v)>9 else 'black')
            fig.colorbar(im,ax=ax,pad=.02,label=r'$C_{\max}$ (dB)',ticks=[-15,0,15])
            save(fig,f'{kind}_filter_matrix_{group}')

def wide():
    rows=[row('full','wide','selected_checks'),row('partial','baseline__1-2_SSS__wide')]
    fig,axs=plt.subplots(2,1,figsize=(6.9,4.9),layout='constrained')
    for ax,r,ref in zip(axs,rows,[10,PARTIAL_REFERENCE]):
        im=heat(ax,r,reference=ref,wide=True,title=r['kind'].capitalize()+' weld: '+label(r['path']))
        ax.set_ylabel('Depth (mm)')
        fig.colorbar(im,ax=ax,pad=.015,shrink=.7,label=f'dB; 0 dB = {ref} '+r'$\mu$m/s',ticks=[-40,-20,0])
    save(fig,'wide_fields')

def filter_response():
    from scipy import signal
    dt=1e-8;f=np.linspace(0,8e6,4001)
    fig=plt.figure(figsize=(6.9,2.95),layout='constrained')
    gs=fig.add_gridspec(2,2,height_ratios=[1,.12])
    axs=[fig.add_subplot(gs[0,j]) for j in range(2)]
    for lo,hi,half,color,name in [(.5,5,328,'#0072b2','0.5-5 MHz; 657 taps'),(0,5,128,'#d55e00','0-5 MHz; 257 taps')]:
        n=np.arange(-half,half+1);width=(hi-lo)*1e6
        base=2*width*dt*np.sinc(width*n*dt)*np.kaiser(len(n),.1102*(55-8.7))
        h=base*np.exp(2j*np.pi*.5*(lo+hi)*1e6*n*dt)*2/base.sum()
        _,H=signal.freqz(h,worN=f,fs=1/dt)
        axs[0].plot(f/1e6,20*np.log10(np.maximum(abs(H)/2,1e-8)),c=color,label=name,lw=1.2)
        axs[1].plot(n*.01,20*np.log10(np.maximum(abs(h)/max(abs(h)),1e-8)),c=color,label=name,lw=1.1)
    axs[0].set(xlim=(0,8),ylim=(-80,3),xlabel='Frequency (MHz)',ylabel='Gain relative to 2 (dB)',title='(a) Positive-frequency response')
    axs[0].axvline(5,c='#777777',ls=':',lw=.8)
    legend_strip(fig,gs[1,:],*axs[0].get_legend_handles_labels(),ncol=2,fontsize=7)
    axs[1].set(xlim=(-3.4,3.4),ylim=(-80,3),xlabel=r'Time from filter centre ($\mu$s)',ylabel='Normalized kernel magnitude (dB)',title='(b) Temporal support')
    for ax in axs:ax.grid(alpha=.15)
    save(fig,'filter_response')

def marker_comparison():
    fig,axs=plt.subplots(2,3,figsize=(6.9,5.6),layout='constrained')
    for i,(kind,id,ref) in enumerate([('full','all',10),('partial','baseline__1-2_SSS__fine',PARTIAL_REFERENCE)]):
        r=row(kind,id,'selected_checks' if kind=='full' else 'confirmation')
        for j,style in enumerate(['contour','bracket','outside']):
            heat(axs[i,j],r,reference=ref,marker=style,title=['Thin neutral contour','Offset tip bracket','External depth indicator'][j])
            if j==0:axs[i,j].set_ylabel(kind.capitalize()+' weld\nDepth (mm)')
            else:axs[i,j].set_yticklabels([])
    save(fig,'marker_comparison')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--collect',action='store_true');args=parser.parse_args()
    if args.collect:collect()
    global INDEX,DATASETS
    INDEX=json.loads((DATA/'images.json').read_text());DATASETS={k:json.loads((DATA/f'{k}_dataset.json').read_text()) for k in ['full','partial']}
    for function in [primary,core_full,partial_atlas,matrices,wide,filter_response,marker_comparison]:function();print(function.__name__,flush=True)

if __name__=='__main__':main()
