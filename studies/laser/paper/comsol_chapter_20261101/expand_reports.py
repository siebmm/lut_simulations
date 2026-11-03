"""Assemble expanded manuscripts and generate additional figures/tables from evidence.

--assemble is retired: the revised manuscripts are edited directly.
Default only regenerates added numerical artwork/tables from portable data.
"""
from pathlib import Path
import argparse, json
import numpy as np
from scipy.signal import butter, sosfiltfilt
import matplotlib.pyplot as plt
import make_artifacts as f
import make_tables as t
from figure_layout import legend_strip

HERE=Path(__file__).resolve().parent
DATA=HERE/'data'

def write(name,value):
    path=HERE/name;path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(value,encoding='utf-8')

def assemble():
    raise RuntimeError('The historical first-draft assembler is retired. Edit the manuscript .tex files directly; default regeneration updates only figures and tables.')


def calibration():
    a=np.load(DATA/'source_calibration/fluence_100.npz')
    report=json.loads((HERE/'verification/source_calibration_check.json').read_text())['runs']['100']
    times=a['time_s'];sos=butter(4,5e6,fs=1/np.diff(times)[0],output='sos')
    thermal=sosfiltfilt(sos,a['thermal_arc'][:,:,2:4],axis=1)*1e6
    dipole=sosfiltfilt(sos,a['dipole_arc'][:,:,2:4],axis=1)*report['least_squares_dipole_scale']*1e6
    fig=plt.figure(figsize=(6.9,2.5),layout='constrained')
    gs=fig.add_gridspec(2,3,height_ratios=[1,.17])
    axs=[fig.add_subplot(gs[0,j]) for j in range(3)]
    for j in range(2):
        axs[j].plot(times*1e6,thermal[3,:,j],c='#333333',lw=1.5,label='Thermoelastic')
        axs[j].plot(times*1e6,dipole[3,:,j],c='#d55e00',lw=1,ls='--',label='Calibrated dipole')
        axs[j].set(xlim=(.12,.82),xlabel=r'Time ($\mu$s)',ylabel=r'Velocity ($\mu$m/s)',title='('+chr(97+j)+') '+('$v_x$' if j==0 else '$v_y$')+' at 45°')
        axs[j].grid(alpha=.15)
    legend_strip(fig,gs[1,:2],*axs[0].get_legend_handles_labels(),ncol=2,fontsize=7)
    angle=np.array(report['per_angle_deg']);error=np.array(report['per_angle_relative_L2_error'])*100
    axs[2].plot(angle,error,c='#555555',lw=.8)
    axs[2].scatter(angle[1:6],error[1:6],s=22,color='#0072b2',zorder=3,label='Fit directions')
    axs[2].scatter(angle[[0,6]],error[[0,6]],s=26,facecolors='white',edgecolors='#d55e00',zorder=3,label='Withheld')
    axs[2].set(xlim=(-5,95),ylim=(0,1.5),xticks=[0,30,60,90],xlabel='Angle from normal (degrees)',ylabel='Relative waveform error (%)',title='(c) Angular check')
    axs[2].grid(alpha=.15)
    legend_strip(fig,gs[1,2],*axs[2].get_legend_handles_labels(),ncol=1,fontsize=7)
    f.save(fig,'source_calibration')

def atlases():
    fig,axs=plt.subplots(3,4,figsize=(6.9,8.1),layout='constrained')
    cases=['lof','reference','difference'];names=['LOF','Pristine','Difference']
    for j,(case,name) in enumerate(zip(cases,names)):
        for i,p in enumerate(f.CORE):
            ax=axs[j,i];im=f.heat(ax,f.row('full','baseline__'+p+'__fine'),case,10,title=f.label(p))
            if i==0:ax.set_ylabel(name+'\nDepth (mm)')
            else:ax.set_yticklabels([])
    f.scale(fig,axs.ravel().tolist(),im,10);f.save(fig,'full_case_atlas')
    order=json.loads((DATA/'survey_protocol.json').read_text())['paths']['partial']
    for case,name in zip(cases,names):
        fig,axs=plt.subplots(5,4,figsize=(6.9,6.3),layout='constrained')
        for i,ax in enumerate(axs.flat):
            if i>=len(order):ax.set_visible(False);continue
            r=next(r for r in f.INDEX if r['kind']=='partial' and r['stage']=='numerical' and r['path']==order[i] and r['variant']=='baseline')
            im=f.heat(ax,r,case,f.PARTIAL_REFERENCE,55,title=f.label(order[i]));ax.set_xlabel('');ax.set_yticks([17,25])
            if i%4==0:ax.set_ylabel('Depth (mm)')
            else:ax.set_yticklabels([])
        fig.supxlabel('$x$ (mm)',fontsize=8)
        fig.colorbar(im,ax=axs.ravel().tolist(),orientation='horizontal',pad=.02,shrink=.8,aspect=45,
                    label=name+rf': amplitude (dB); common 0 dB = {f.PARTIAL_REFERENCE:g} $\mu$m/s',ticks=[-50,-40,-30,-20,-10,0])
        f.save(fig,'partial_'+case+'_atlas')
    for kind,id,ref in [('full','wide',10),('partial','baseline__1-2_SSS__wide',f.PARTIAL_REFERENCE)]:
        fig,ax=plt.subplots(figsize=(6.9,2.5),layout='constrained')
        r=f.row(kind,id,'selected_checks' if kind=='full' else 'confirmation')
        im=f.heat(ax,r,reference=ref,wide=True,title=f.label(r['path']));ax.set_ylabel('Depth (mm)')
        fig.colorbar(im,ax=ax,pad=.015,label=f'dB; 0 dB = {ref} '+r'$\mu$m/s',ticks=[-40,-20,0])
        f.save(fig,kind+'_wide')

def longtable(name,caption,label,columns,header,lines):
    value=r'{\small\setlength{\tabcolsep}{4pt}'+'\n'
    if name == 'full_confirmation_register':
        value += r'\singlespacing'+'\n'
    value += r'\begin{longtable}{'+columns+'}\n'
    value+=r'\caption{'+caption+r'}\label{'+label+r'}\\'+'\n'+r'\toprule '+header+r'\\\midrule\endfirsthead'+'\n'
    value+=r'\multicolumn{'+str(len(header.split('&')))+r'}{l}{\small Continued from previous page}\\\toprule '+header+r'\\\midrule\endhead'+'\n'
    value+=r'\bottomrule\endfoot'+'\n'+'\n'.join(lines)+'\n'+r'\end{longtable}}'+'\n'
    write(name+'.tex',value)

def tables():
    full=t.rows('full_confirmation_metrics.json');partial=t.rows('partial_confirmation_metrics.json');sel=t.rows('full_selected_metrics.json')
    groups=[('all','All sources and receivers'),('shoulders','Shoulder sources only'),('central','Weld-surface sources only'),('plate','Side-surface receivers only'),('floor','Weld-surface receivers only'),('shoulders_plate','Shoulders and side-surface receivers')]
    rows=[]
    for key,desc in groups:
        a=t.metric(sel,key);b=t.metric(partial,'baseline__1-2_SSS__'+('fine' if key=='all' else key))
        rows.append(desc+' & '+' & '.join([t.num(a['cmax_db']),t.num(a['near_support_pct'],1),t.num(b['cmax_db']),t.num(b['near_support_pct'],1)])+r' \\')
    t.textfile('tables/comsol/main_groups.tex',rows)
    for kind in ['full','partial']:
        rows=[]
        for title,rr,prefix in ([('$S\\mid S$',sel,''),('$L\\mid L$',full,'baseline__1-1_LL__')] if kind=='full' else [('',partial,'baseline__1-2_SSS__')]):
            for key,desc in groups:
                r=t.metric(rr,key if rr is sel else prefix+('fine' if key=='all' else key))
                rows.append((title+' & ' if title else '')+desc+' & '+t.num(r['cmax_db'])+' & '+t.num(r['near_support_pct'],1)+r' \\')
            if kind=='full':rows.append(r'\addlinespace')
        t.textfile(f'tables/comsol/{kind}_groups.tex',rows)
    rows=[]
    for p in f.CORE:
        for case in ['lof','reference','difference']:
            r=t.metric(full,'baseline__'+p+'__fine',case)
            rows.append(t.pathlabel(p)+' & '+{'lof':'LOF','reference':'Pristine','difference':'Difference'}[case]+' & '+
                ' & '.join([t.num(r['peak']['um_s'],3),t.num(r['cmax_db']),t.num(r['c99_db']),t.num(r['near_support_pct'],1),t.num(r['peak']['distance_mm'],3)])+r' \\')
        rows.append(r'\addlinespace')
    t.textfile('tables/comsol/full_controls.tex',rows)
    survey=t.rows('survey_metrics.json');order=json.loads((DATA/'survey_protocol.json').read_text())['paths']['partial'];rows=[]
    for p in order:
        for case in ['lof','reference','difference']:
            r=next(r for r in survey if r['kind']=='partial' and r['variant']=='baseline' and r['path']==p and r['case']==case)
            rows.append(t.pathlabel(p)+' & '+{'lof':'LOF','reference':'Pristine','difference':'Difference'}[case]+' & '+
                ' & '.join([t.num(r['near']['um_s'],2) if r['near'] else '--',t.num(r['cmax_db']),t.num(r['c99_db']),t.num(r['near_support_pct'],1),t.num(r['peak']['distance_mm'],3) if r['peak'] else '--'])+r' \\')
        rows.append(r'\addlinespace')
    longtable('partial_primary_register','All eighteen primary partial-weld paths and three signal cases on the 0.2 mm survey grid. $I_N$ is the near-region maximum in $\\mu$m/s; distances are mm and contrasts dB. N is near support.','comsol:partial:primary_register','@{}llrrrrr@{}',
              r'Path & Case & $I_N$ & $C_{\max}$ & $C_{99}$ & N (\%) & $d_{\rm pk}$',rows)
    variants={'baseline':'BP657','band1_5':'1--5/657','bp_short_matched_guard':'BP257/G328','lp_short_matched_guard':'LP257/G328','lp_short_initial_rest':'LP257/rest'}
    aliases={'fine':'all','finer':'all','all':'all','wide':'wide','central':'weld src','floor':'weld rx','plate':'side rx','shoulders':'shoulder src','shoulders_plate':'shoulder+side','gate15':'gate 15\\%','gate25':'gate 25\\%','lowpass_launch_curve':'LP launch'}
    for kind,rr in [('full',full),('partial',partial)]:
        rows=[]
        configs=[r for r in f.INDEX if r['kind']==kind and r['stage']=='confirmation']
        if kind=='full':configs += [r for r in f.INDEX if r['stage']=='selected_checks']
        for config in configs:
            source=sel if config['stage']=='selected_checks' else rr
            cases=[t.metric(source,config['id'],c) for c in ['lof','reference','difference']]
            suffix=config['id'].split('__')[-1];variant=variants[config['variant']]
            rows.append(variant+' & '+t.pathlabel(config['path'])+' & '+f"{config['selection']['spacing_mm']:g}"+' & '+aliases[suffix]+' & '+
                        ' & '.join(t.num(r['cmax_db']) for r in cases)+' & '+t.num(cases[0]['near_support_pct'],1)+r' \\')
        longtable(kind+'_confirmation_register',
            'Complete '+kind+'-weld confirmations. BP denotes 0.5--5 MHz and LP 0--5 MHz; G328 retains the 3.28 $\\mu$s guard. Rest uses known zero prehistory. L/P/D are LOF, pristine and difference $C_{\\max}$ in dB; N is LOF near support. Src/rx denote source/receiver group reductions; weld and side identify the weld surface and side surfaces.',
            'comsol:'+kind+':register','@{}llrlrrrr@{}',r'Filter & Path & Grid & Change & L & P & D & N (\%)',rows)
    counts={kind:sum(r['kind']==kind for r in f.INDEX) for kind in ['full','partial']}
    assert counts=={'full':104,'partial':280}
    write('verification/expanded_tables.json',json.dumps({'status':'passed','configurations':counts,'partial_primary_rows':54,'full_control_rows':12,'confirmation_rows':{'full':48,'partial':28}},indent=2)+'\n')

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--assemble',action='store_true');args=parser.parse_args()
    if args.assemble:assemble()
    f.INDEX=json.loads((DATA/'images.json').read_text());f.DATASETS={k:json.loads((DATA/f'{k}_dataset.json').read_text()) for k in ['full','partial']}
    tables();calibration();atlases()
    print('Expanded figures and tables generated from existing verified evidence.')

if __name__=='__main__':main()
