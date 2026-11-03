"""Report the full right-receiver search, including failed wider-field checks."""

# Locate portable scientific inputs without depending on the source checkout.
import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import data_root as _data_root, native_root as _native_root, executable as _executable, native_env as _native_env
import json
import csv
import sys
import zipfile
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from right_paths import DEST, OUT, HERE, STUDY, s, verify_settings

survey=json.loads((DEST/'survey/metrics.json').read_text())['rows']
confirmation=json.loads((DEST/'confirmation/metrics.json').read_text())['rows']
reflection=json.loads((DEST/'reflection/metrics.json').read_text())['rows']
protocol=json.loads((DEST/'confirmation_protocol.json').read_text())


def row(name,case='lof',stage='confirmation'):
    return next(r for r in (confirmation if stage=='confirmation' else reflection) if r['id']==name and r['case']==case)


def arrays(path):
    with np.load(path) as f:return dict(f)


tips=np.array(s.dataset('partial')['audit']['defect_benchmark_only']['tips_xy_mm'])
tips[:,1]=25.4-tips[:,1];tips=tips[np.argsort(tips[:,1])]
figdir=DEST/'figures';figdir.mkdir(exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.titlesize':11,'savefig.dpi':180})
cmap=plt.get_cmap('magma').copy();cmap.set_bad('#cccccc')


def panel(ax,a,value,title,wide=False):
    db=20*np.log10(np.maximum(value.reshape(len(a['z']),len(a['x']))/40e-6,1e-10))
    im=ax.pcolormesh(a['x'],a['z'],np.ma.masked_invalid(db),vmin=-60,vmax=10,
                    shading='nearest',cmap=cmap,rasterized=True)
    ax.set_aspect('equal');ax.invert_yaxis()
    ax.plot(tips[:,0],tips[:,1],'--',color='#44ddff',lw=.9)
    ax.set(xlabel='x (mm)',ylabel='Depth below original surface (mm)',title=title)
    if wide:
        profile=np.array(s.dataset('partial')['surface_profile_x_depth_m'])*1000
        ax.fill_between(profile[:,0],0,profile[:,1],color='white',zorder=2)
        ax.plot(profile[:,0],profile[:,1],color='#888888',lw=.7)
        ax.set(xlim=(-45,45),ylim=(25.4,0))
    return im


comparison=[
    ('Current: S | SS\nAll five sources, all receivers',OUT/'admission_audit/angular_confirmation/fine/arrays.npz'),
    ('Alternative: S | LL\nRight weld source, right receivers',DEST/'confirmation/S_LL_right__fine/arrays.npz'),
    ('Diagnostic: SSL | SS\nLeft weld source, right receivers',DEST/'confirmation/SSL_SS_left__fine/arrays.npz'),
]
fig,axes=plt.subplots(2,3,figsize=(15.5,8),layout='constrained')
for col,(label,path) in enumerate(comparison):
    a=arrays(path);points=a['points'].reshape(-1,2)
    for ri,case in enumerate(('lof','reference')):
        m=s.metrics(a[case+'_2_amplitude'],points,'partial')
        title=(label if ri==0 else 'Pristine control, identical settings')
        title+=f'\nNear peak {m["near"]["um_s"]:.3g} micrometres/s; contrast {m["cmax_db"]:.2f} dB'
        im=panel(axes[ri,col],a,a[case+'_2_amplitude'],title)
fig.colorbar(im,ax=axes,shrink=.85,pad=.015,label='Amplitude (dB relative to 40 micrometres/s); common scale')
fig.suptitle('Partial weld: other paths do contain defect-related information\n0.05 mm grid; weighted 1-5 MHz; 15% angular gate. Cyan: prescribed defect; grey: unsupported.',fontsize=13)
fig.savefig(figdir/'comparison.png');plt.close(fig)

fig,axes=plt.subplots(2,1,figsize=(14,8.3),layout='constrained')
for ax,name,label in zip(axes,('S_LL_right','SSL_SS_left'),('S | LL; right weld source','SSL | SS; left weld source')):
    a=arrays(DEST/f'confirmation/{name}__wide/arrays.npz');m=row(name+'__wide')
    im=panel(ax,a,a['lof_2_amplitude'],f'{label}, right receivers\nWide-field contrast {m["cmax_db"]:.2f} dB; near peak {m["near"]["um_s"]:.3g} micrometres/s',True)
fig.colorbar(im,ax=axes,pad=.015,shrink=.85,label='Amplitude (dB relative to 40 micrometres/s); common scale')
fig.suptitle('Wider-field checks expose false indications outside the weld neighbourhood\n0.25 mm sampling; common amplitude scale, no individual normalisation.',fontsize=13)
fig.savefig(figdir/'wide_field.png');plt.close(fig)

# Illustrate actual, geometry-valid rays. The known segment centre is used only
# for this explanatory diagram, never to generate image delays or select data.
sys.path.insert(0,str(STUDY/'src'))
from verify_saft_datasets import PartialReference
ref=PartialReference(_data_root()/'partial_weld')
focus=tips.mean(axis=0)/1000
receiver=min((i for i,r in enumerate(ref.receivers) if r['group']=='right'),key=lambda i:abs(ref.xy[i,0]-.02))
ray_records=[]
fig,axes=plt.subplots(1,2,figsize=(13,5.2),layout='constrained')
for ax,source,first,last,label in zip(axes,(3,1),('S','SSL'),('LL','SS'),('S | LL: right weld source','SSL | SS: left weld source')):
    inc=ref.ray(ref.sources[source],focus,first);out=ref.ray(focus,ref.xy[receiver],last,True)
    assert inc is not None and out is not None
    profile=ref.profile*1000
    ax.fill_between(profile[:,0],profile[:,1],25.4,color='#e8edf0')
    ax.plot(profile[:,0],profile[:,1],color='#777777',lw=1)
    ax.axhline(25.4,color='#777777',lw=1)
    for route,modes in ((inc,first),(out,last)):
        vertices=route['points']*1000
        for a,b,mode in zip(vertices[:-1],vertices[1:],modes):
            color='#0077bb' if mode=='S' else '#cc6611'
            ax.annotate('',xy=b,xytext=a,arrowprops=dict(arrowstyle='->',color=color,lw=1.8))
            mid=(a+b)/2
            ax.text(mid[0]+.4,mid[1],mode,color=color,fontweight='bold',fontsize=11)
    ax.scatter(* (ref.sources[source]*1000),marker='*',s=100,color='#333333',zorder=5,label='Source')
    ax.scatter(*(ref.xy[receiver]*1000),marker='v',s=45,color='#333333',zorder=5,label='Receiver')
    ax.plot(tips[:,0],tips[:,1],color='#aa2266',lw=2,label='Prescribed defect')
    ax.scatter(*(focus*1000),s=30,color='#aa2266',zorder=5)
    ax.set(xlim=(-5,25),ylim=(26,-1),xlabel='x (mm)',ylabel='Depth (mm)',title=label)
    ax.set_aspect('equal')
    ray_records.append(dict(path=first+'|'+last,source=source,receiver=receiver,focus_m=focus.tolist(),
        incoming_points_m=inc['points'].tolist(),outgoing_points_m=out['points'].tolist(),
        arrival_us=float((inc['time']+out['time']+3e-8)*1e6)))
axes[0].legend(loc='upper left',fontsize=8)
fig.suptitle('Two additional focusing hypotheses\nS = shear; L = longitudinal. The divider marks scattering at the image point.',fontsize=13)
fig.savefig(figdir/'paths.png');plt.close(fig)
s.write(DEST/'illustrative_rays.json',dict(method='Independent geometry implementation; illustrative rays, not unique physical path identification.',rows=ray_records))

# Recipes exactly preserve the native-checked fine-grid selections.
settings=DEST/'settings';settings.mkdir(exist_ok=True)
template=json.loads((OUT/'settings/saft_partial_weighted_1_5.json').read_text())
for name in ('S_LL_right','SSL_SS_left','SSL_SS_two'):
    job=next(j for j in protocol['requests'] if j['id']==name+'__fine')
    session=dict(template,selection=dict(job['selection'],cases=['lof']))
    session['display']=dict(template['display'],reference_um_s=40,db_range=60)
    s.write(settings/('saft_partial_'+name+'.json'),session)
verify_settings()

compact=[]
for stage,rows in (('confirmation',confirmation),('reflection',reflection)):
    for r in rows:
        if r['case']!='lof':continue
        compact.append(dict(stage=stage,id=r['id'],contrast_db=r['cmax_db'],near_support_pct=r['near_support_pct'],
            near_peak_um_s=r['near']['um_s'] if r['near'] else None,
            peak_distance_mm=r['peak']['distance_mm'] if r['peak'] else None,
            lof_over_pristine_near_db=r.get('lof_over_pristine_near_db')))
s.write(DEST/'summary.json',dict(native_configurations=224,preset_replay_checks=3,main_survey_configurations=192,
    note='Survey includes five fixed source-subset views per configuration. Confirmation 22 and reflected-screen follow-up 10. Three signal cases each; no new COMSOL solves.',rows=compact))
with (DEST/'summary.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(compact[0]));w.writeheader();w.writerows(compact)

table='\n'.join(f"  {r['id']:<43} C={r['contrast_db']:7.2f} dB; near support={r['near_support_pct']:6.2f}%; near peak={r['near_peak_um_s']:.4g} um/s"
                 if r['contrast_db'] is not None else f"  {r['id']}: no usable near/background comparison" for r in compact)
body=f'''PARTIAL-WELD RIGHT-RECEIVER PATH STUDY - 2 NOVEMBER 2026

The earlier work did not establish that S|SS is the only useful path.
This study finds a reproducible alternative: S|LL using the rightmost weld
source at x=+1.5 mm and the right-side receiver group. Its fine-grid contrast
is 3.02 dB, near support 96.81%, and near peak 26.97 micrometres/s. The near
peak is 14.37 dB above the matched pristine near peak. The wide-field check
retains positive contrast (4.36 dB) and a peak near the defect. Even/odd
receiver subsets retain positive contrast. Equal receiver averaging loses
the positive contrast; known-response phase/amplitude weighting matters.

S|LL means S from source to image point, then L from image point to backwall
and L from backwall to receiver. Conversion at the image point is a focusing
hypothesis. The ray model does not calibrate the defect scattering coefficient
or prove that a unique physical wave packet follows that route.

SSL|SS using the leftmost weld source also has a defect-related local feature:
fine-grid contrast 3.20 dB with 100% near support, but its near peak is only
0.330 micrometres/s and its wider-field contrast is -28.78 dB. Requiring two
weld sources gives 8.04 dB local contrast but only 43.88% near support and
fails the wider-field check. LSS|S and SLS|LS also fail wider-field checks.
These longer paths are diagnostic views, not replacements for the main image.
Reflected-packet screening and 20% angular controls do not rescue SSL|SS in
the wider field: its contrast remains approximately -27.5 to -28.6 dB.
Enabling all reflected-packet screens removes ALL near support for S|LL.
Its candidate delays therefore overlap rejected reflected-pulse windows;
this screen cannot isolate a surviving near-defect image. The defect/pristine
comparison establishes sensitivity to the defect in the simulation, but does
not uniquely identify an isolated S-to-L defect-scattered packet. Keep this
limitation with the useful alternative rather than treating its path label
as a confirmed physical attribution.

The current all-source/all-receiver S|SS image remains stronger in contrast:
9.10 dB at 0.05 mm with the 15% gate. S|LL supplies another view using a single
weld source and only right receivers. It does not establish improved tip sizing.
In the matched 0.2 mm rightmost-source/right-receiver comparison, S|SS gives
-1.13 dB while S|LL gives +3.03 dB. Thus choosing a path for this acquisition
does matter. Narrowing the right group to fixed near/middle/far windows did
not produce a better robust direct/one-reflection alternative in this scan.

SCOPE
All 132 conventional paths implemented in the native engine were tested with
three weld sources and 157 right receivers (x=4.5534 to 43.5534 mm). The set
includes mode-converted reflections and up to five total legs. The 20 direct
and single-reflection paths were also tested with fixed receiver windows:
  near: x<15 mm, 42 receivers;
  middle: 15<=x<35 mm, 80 receivers;
  far: x>=35 mm, 35 receivers.
There were 192 survey configurations; each retained LOF, pristine and waveform
difference images. Per-source means gave fixed views of each of the three weld
sources and three-source RMS views requiring one or two sources. The source
choice is constant throughout an image; no per-pixel choice of the best source
or path was used. The two-source requirement is not imposed on a single-source
measurement. Another 22 configurations checked selected results; 10 examined
reflected screening and angular thresholds after wider-field failures.

Main settings: weighted 1-5 MHz analytic FIR, 15% angular response threshold,
5% wall threshold, at least 8 real/effective receivers, 25% retained fraction,
CF exponent 1, direct L/S/R pulse screens, no reflected screens initially,
vy, source RMS. Survey grid: 0.2 mm; fine: 0.05 mm; wide: 0.25 mm over
x=+/-45 mm and depth 0.25-25.25 mm. The finer field is x=+/-7.5 mm and
depth 16.65-25.25 mm. Changing spacing and bounds together in the wide run
means its contrast is not a pure field-size sensitivity measurement.

METRICS AND LIMITS
Near support counts valid pixels within 0.5 mm of the known 2.5 mm defect
segment. Contrast is 20 log10(near peak / strongest background more than 1 mm
from the segment). The 0.5-1 mm region is excluded. Coverage is not detection
probability. Unsupported pixels are not tested background. Pristine ratios
use common near support; difference images subtract waveforms before focusing,
not final amplitudes. LOF images do not use pristine subtraction. Measured
amplitudes are simulation outputs, not an experimental energy calibration.
No experimental-noise model, independent defect validation or mesh-convergence
claim is made. Candidate selection is retrospective on one known flaw. This
survey is not exhaustive of all physical scattering mechanisms or settings.

VERIFICATION
All survey three-source RMS images reproduce the exported native image.
Single-source native replay and nested fine-grid checks agree with the survey
source-subset calculations. Protocol files retain input fingerprints and hashes.
All three exported presets reproduce the corresponding LOF image in separate
native LOF-only replays; GUI interaction was not tested.
No production engine, raw data or COMSOL models were changed.

CONFIRMATIONS AND FOLLOW-UP (LOF images)
{table}

FILES
  figures/comparison.png: current image and two alternatives, with pristine
    controls on a common absolute amplitude scale.
  figures/wide_field.png: wider-field S|LL and SSL|SS comparison.
  figures/paths.png: geometry-valid illustrative rays for these two hypotheses.
  settings/saft_partial_S_LL_right.json: useful alternative, single right weld
    source and right receivers. Fine grid exactly matches the native check.
  settings/saft_partial_SSL_SS_left.json: weaker diagnostic with wide-field clutter.
  settings/saft_partial_SSL_SS_two.json: diagnostic with limited coverage.
  summary.csv / summary.json: compact verification comparisons.
  survey/metrics.csv / metrics.json: all surveyed outcomes, including failures
    and plane-1 diagnostics without CF/pulse screening.
  protocol.json, confirmation_protocol.json, reflection_protocol.json: full
    reconstruction requests and the stated reasons for follow-up.
  survey/, confirmation/, reflection/: arrays, timings and native evidence.

REPRODUCE FROM THE REPOSITORY ROOT
Use the existing COMSOL recordings and pinned native SAFTFilterStudy executable.
  python studies/laser/paper/partial_improvement_20261102/right_paths.py run
  python studies/laser/paper/partial_improvement_20261102/right_paths.py confirm
  python studies/laser/paper/partial_improvement_20261102/right_paths.py reflection
  python studies/laser/paper/partial_improvement_20261102/right_paths_report.py
The review ZIP contains figures, presets, protocols and metrics. Full arrays
and original recordings remain in the repository and are required to rerun.
'''
(DEST/'README.txt').write_text(body,encoding='utf-8')
files={}
for path in DEST.iterdir():
    if path.is_file() and path.suffix in ('.json','.txt','.csv'):files[path.name]=path
for sub in ('figures','settings'):
    for path in (DEST/sub).iterdir():files[sub+'/'+path.name]=path
for stage in ('survey','confirmation','reflection'):
    for name in ('metrics.json','verification.json'):
        files[stage+'/'+name]=DEST/stage/name
for name in ('right_paths.py','right_paths_report.py'):
    files['scripts/'+name]=HERE/name
bundle=STUDY/'results/reports/partial_weld_right_paths_20261102.zip'
with zipfile.ZipFile(bundle,'w',zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for name,path in sorted(files.items()):z.write(path,name)
with zipfile.ZipFile(bundle) as z:
    assert z.testzip() is None
    for name,path in files.items():assert z.read(name)==path.read_bytes()
print('Saved figures, presets, summary and verified review bundle:',bundle,flush=True)
