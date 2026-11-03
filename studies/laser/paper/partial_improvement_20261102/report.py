"""Create reproducible comparison figures, sensitivity evidence and loadable recipes."""
import copy
import csv
import json
from pathlib import Path
import zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from run import HERE, OUT, ROOT, STUDY, s, protocol, frequency_plan

FIG=OUT/'figures'
SETTINGS=OUT/'settings'
FIG.mkdir(exist_ok=True)
SETTINGS.mkdir(exist_ok=True)
rows=[]
for stage in ('survey','frequency','confirmation'):
    rows.extend(json.loads((OUT/stage/'metrics.json').read_text())['rows'])


def metric(stage,name,case='lof'):
    return next(r for r in rows if r['stage']==stage and r['id']==name and r['case']==case)


def arrays(stage,name):
    with np.load(OUT/stage/name/'arrays.npz') as a:
        return dict(a)


choices=[
    ('survey','all__raw__min2','Current baseline\nEqual mean, 0.5–5 MHz','baseline'),
    ('survey','all__weighted__min2','Earlier improvement\nResponse weighted, 0.5–5 MHz','weighted_05_5'),
    ('frequency','all__weighted__min2__band1_5','Candidate\nResponse weighted, 1–5 MHz','weighted_1_5'),
]
tips=np.array(s.dataset('partial')['audit']['defect_benchmark_only']['tips_xy_mm'])
tips[:,1]=25.4-tips[:,1]
tips=tips[np.argsort(tips[:,1])]
lo,hi=tips[:,1]
cmap=plt.get_cmap('magma').copy()
cmap.set_bad('#cccccc')
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'axes.labelsize':10,'savefig.dpi':180})


def image_panel(ax,a,values,title,wide=False):
    x,z=a['x'],a['z']
    field=values.reshape(len(z),len(x))
    db=20*np.log10(np.maximum(field/40e-6,1e-8))
    im=ax.pcolormesh(x,z,np.ma.masked_invalid(db),cmap=cmap,vmin=-40,vmax=0,shading='nearest',rasterized=True)
    ax.invert_yaxis();ax.set_aspect('equal')
    ax.set_xlabel('x (mm)');ax.set_ylabel('Depth below original surface (mm)')
    ax.set_title(title,pad=9)
    if wide:
        profile=np.array(s.dataset('partial')['surface_profile_x_depth_m'])*1000
        ax.fill_between(profile[:,0],0,profile[:,1],color='white',zorder=2)
        ax.plot(profile[:,0],profile[:,1],color='#8c8c8c',lw=.8)
        ax.set_ylim(25.4,0)
    ax.plot(tips[:,0],tips[:,1],color='#44ddff',ls='--',lw=.8)
    ax.scatter(tips[:,0],tips[:,1],facecolors='none',edgecolors='#44ddff',s=17,lw=.8)
    return im


fig,axes=plt.subplots(2,3,figsize=(15,8.2),layout='constrained')
for col,(stage,name,label,key) in enumerate(choices):
    a=arrays(stage,name)
    for row,case in enumerate(('lof','reference')):
        m=metric(stage,name,case)
        title=(label if row==0 else 'Pristine control, same settings')
        title+=f'\nContrast {m["cmax_db"]:.2f} dB; support {m["near_support_pct"]:.1f}%'
        im=image_panel(axes[row,col],a,a[case+'_2_amplitude'],title)
fig.colorbar(im,ax=axes,shrink=.85,pad=.02,label='Amplitude (dB); 0 dB = 40 µm/s, identical scale in every panel')
fig.suptitle('Partial weld: S | SS, existing COMSOL recordings\nTop: flawed model alone. Bottom: pristine model. Cyan: known defect; grey: unsupported.',fontsize=13)
for ext in ('png','pdf'):fig.savefig(FIG/f'comparison.{ext}')
plt.close(fig)

# Evaluation-only profile: the flaw coordinates are never used in reconstruction.
profile_rows=[]
fig,axes=plt.subplots(1,2,figsize=(12,4.4),layout='constrained')
profile_choices=choices+[('frequency','all__raw__min2__band2_5','Equal mean, 2–5 MHz','equal_2_5')]
colors=['#555555','#d08100','#0077bb','#a43b7d']
for (stage,name,label,key),color in zip(profile_choices,colors):
    a=arrays(stage,name);m=metric(stage,name)
    x,z=a['x'],a['z'];v=a['lof_2_amplitude'].reshape(len(z),len(x))
    centre=tips[0,0]+(z-lo)*(tips[1,0]-tips[0,0])/(hi-lo)
    strip=np.abs(x[None,:]-centre[:,None])<=.35
    masked=np.where(strip&np.isfinite(v),v,-np.inf)
    peaks=masked.max(axis=1);peaks[~np.isfinite(peaks)]=np.nan
    interval=(z>=lo)&(z<=hi)
    relative=20*np.log10(peaks/(m['background']['um_s']*1e-6))
    shown=(z>=17)&(z<=22.5)
    axes[0].plot(z[shown],peaks[shown]*1e6,label=label.replace('\n','; '),color=color)
    axes[1].plot(z[shown],relative[shown],label=label.replace('\n','; '),color=color)
    profile_rows.append(dict(id=name,stage=stage,
        definition='At each depth, maximum amplitude within +/-0.35 mm of the known segment line; evaluation only, not a size estimate.',
        supported_depth_rows_pct=float(100*np.mean(np.isfinite(peaks[interval]))),
        depth_rows_above_strongest_background_pct=float(100*np.mean(relative[interval]>0)),
        depth_rows_within_6dB_of_near_peak_pct=float(100*np.mean(peaks[interval] >= .501187*m['near']['um_s']*1e-6))))
for ax in axes:
    ax.axvspan(lo,hi,color='#b8e5ee',alpha=.3,label='Actual 2.5 mm defect depth')
    ax.axvline(lo,color='#555555',ls=':',lw=.8);ax.axvline(hi,color='#555555',ls=':',lw=.8)
    ax.set(xlim=(17,22.5),xlabel='Depth below original surface (mm)');ax.grid(alpha=.2)
axes[0].set_ylabel('Maximum in a ±0.35 mm strip (µm/s)')
axes[1].set_ylabel('Strip peak / strongest distant background (dB)')
axes[1].axhline(0,color='black',lw=1)
axes[0].legend(fontsize=8)
fig.suptitle('Does the improvement extend along the defect?\nProfiles use the known segment only for assessment; they do not measure defect size.',fontsize=12)
for ext in ('png','pdf'):fig.savefig(FIG/f'defect_profiles.{ext}')
plt.close(fig)
s.write(OUT/'defect_profiles.json',dict(rows=profile_rows))

fig,axes=plt.subplots(2,1,figsize=(14,8.2),layout='constrained')
for ax,name,label in zip(axes,('baseline__wide','weighted_1_5__wide'),('Baseline, 0.5–5 MHz','Candidate: response weighted, 1–5 MHz')):
    a=arrays('confirmation',name);m=metric('confirmation',name)
    im=image_panel(ax,a,a['lof_2_amplitude'],f'{label}; wide-field contrast {m["cmax_db"]:.2f} dB',True)
    ax.set_xlim(-45,45)
fig.colorbar(im,ax=axes,shrink=.9,pad=.02,label='Amplitude (dB); 0 dB = 40 µm/s')
fig.suptitle('Wider-field check: identical settings, expanded image bounds',fontsize=13)
for ext in ('png','pdf'):fig.savefig(FIG/f'wide_field.{ext}')
plt.close(fig)

# Loadable schema-4 files use the production app's existing format.
template=json.loads((HERE.parent/'laser_weld_feasibility_20261030/workbench/saft_settings_partial_S_SS_weighted.json').read_text())
exports=choices+[('frequency','all__raw__min2__band2_5','Equal mean, 2–5 MHz','equal_2_5')]
for stage,name,label,key in exports:
    config=copy.deepcopy(template)
    config['selection']=copy.deepcopy(metric(stage,name)['selection'])
    config['selection']['cases']=['lof']
    config['display'].update(reference_um_s=40,db_range=40)
    s.write(SETTINGS/f'saft_partial_{key}.json',config)

# Independent recombination of every exported case; nested-grid comparison.
checks=[]
plans={
    'survey':protocol()['requests'],
    'frequency':frequency_plan(),
    'confirmation':json.loads((HERE/'confirmation_protocol.json').read_text())['requests'],
}
worst=0.
for stage,jobs in plans.items():
    for job in jobs:
        a=arrays(stage,job['id']);n=len(a['points'])//2;ns=len(job['selection']['sources'])
        for case in ('lof','reference','difference'):
            valid=a[case+'_2_valid'].reshape(ns,n).astype(bool)
            means=a[case+'_2_means'].reshape(ns,n)
            count=valid.sum(axis=0)
            reconstructed=np.sqrt(np.sum(np.where(valid,abs(means)**2,0),axis=0)/np.maximum(count,1))
            reconstructed[count<job['selection']['minimum_sources']]=np.nan
            actual=a[case+'_2_amplitude']
            np.testing.assert_allclose(reconstructed,actual,rtol=5e-14,atol=1e-18,equal_nan=True)
            np.testing.assert_array_equal(count,a[case+'_2_count'])
            worst=max(worst,float(np.nanmax(abs(reconstructed-actual))))
for fine,stage,name in [('baseline__fine','survey','all__raw__min2'),
                        ('weighted_1_5__fine','frequency','all__weighted__min2__band1_5'),
                        ('equal_2_5__fine','frequency','all__raw__min2__band2_5')]:
    a=arrays('confirmation',fine);b=arrays(stage,name)
    maximum=0.
    for case in ('lof','reference','difference'):
        refined=a[case+'_2_amplitude'].reshape(len(a['z']),len(a['x']))[::2,::2]
        original=b[case+'_2_amplitude'].reshape(len(b['z']),len(b['x']))
        np.testing.assert_allclose(refined,original,rtol=2e-8,atol=2e-15,equal_nan=True)
        maximum=max(maximum,float(np.nanmax(abs(refined-original))))
    checks.append(dict(id=fine,passed=True,maximum_absolute_difference=maximum))
s.write(OUT/'confirmation_verification.json',dict(status='passed',configurations=sum(map(len,plans.values())),
    signal_cases=len(rows),independent_source_RMS_maximum_difference=worst,nested_grids=checks))

selected=[]
for stage,name,label,key in exports:
    m=metric(stage,name);ref=metric(stage,name,'reference')
    selected.append(dict(label=label.replace('\n','; '),id=name,stage=stage,
        cmax_db=m['cmax_db'],near_support_pct=m['near_support_pct'],
        peak_um_s=m['peak']['um_s'],peak_distance_mm=m['peak']['distance_mm'],
        one_source_near_pct=m['near_one_source_pct'],
        pristine_peak_um_s=ref['peak']['um_s'],pristine_near_peak_um_s=ref['near']['um_s'],
        lof_near_over_pristine_near_db=float(20*np.log10(m['near']['um_s']/ref['near']['um_s'])),
        settings=f'settings/saft_partial_{key}.json'))
s.write(OUT/'summary.json',dict(configurations=sum(map(len,plans.values())),signal_cases=len(rows),selected=selected,
    confirmation=[{k:r[k] for k in ('id','case','cmax_db','near_support_pct','near_one_source_pct','peak')} for r in rows if r['stage']=='confirmation']))
with (OUT/'summary.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(selected[0]));w.writeheader();w.writerows(selected)
print('Figures, four loadable recipes and numerical checks written.',flush=True)
print(json.dumps(selected,indent=2),flush=True)
