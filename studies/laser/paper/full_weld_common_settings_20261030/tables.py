"""Generate manuscript tables and compact reusable numerical evidence."""

import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import executable as _executable, native_env as _native_env
from analyse import *
import csv

analysis=json.loads((OUT/'analysis.json').read_text());runs=analysis['runs']
index=json.loads((OUT/'index.json').read_text())
def get(k,name='baseline'):return next(r for r in runs if r['path']==PATHS[k] and r['name']==name)
def label(k):return [r'$S\mid S$',r'$L\mid L$',r'$S\mid SS$',r'$SS\mid S$'][k]
def write(name,lines):(HERE/name).write_text('\n'.join(lines)+'\n',encoding='utf-8')
write('table_baseline.tex',[
    label(k)+' & '+f"{m['global_peak_um_s']:.3f} & {m['global_peak_x_mm']:.2f}, {m['global_peak_depth_mm']:.2f} & {m['contrast_peak_db']:.2f} & {m['contrast_p99_db']:.2f} & {m['support_pct']:.1f}"+r' \\'
    for k in range(4) for m in [get(k)['signals']['lof']]])
write('table_reference.tex',[
    label(k)+' & '+' & '.join(f"{get(k)['signals'][c]['near_peak_um_s']:.3f}" for c in ['lof','reference','difference'])+' & '+f"{get(k)['signals']['difference']['contrast_p99_db']:.2f}"+r' \\' for k in range(4)])
write('table_repairs.tex',[
    label(k)+' & '+' & '.join(f"{m[v]:.2f}" for v in ['support_pct','near_peak_um_s','contrast_peak_db','contrast_p99_db'])+r' \\' for k in range(4) for m in [get(k,'clean_fine')['signals']['lof']]])
write('table_grid.tex',[
    label(k)+' & '+f"{b['global_peak_um_s']:.3f} & {f['global_peak_um_s']:.3f} & {20*np.log10(f['global_peak_um_s']/b['global_peak_um_s']):+.3f} & {f['global_peak_x_mm']:.2f}, {f['global_peak_depth_mm']:.2f}"+r' \\' for k in range(4) for b,f in [(get(k)['signals']['lof'],get(k,'fine')['signals']['lof'])]])
source_rows=[]
for k in range(4):
    a=load(index[k]);p=np.array([2.766365013,19.71747506]);j=np.argmin(np.sum((a['points_m']*1000-p)**2,axis=1))
    counts=a['lof_effective'][2,:,j]; valid=a['lof_valid'][2,:,j]
    source_rows.append(label(k)+' & '+' & '.join(str(int(n)) if v else '--' for n,v in zip(counts,valid))+r' \\')
write('table_sources.tex',source_rows)
screen_rows=[]
for screen in ['off','L','S','R','LS','LR','SR','LSR']:
    name='coarse' if screen=='LSR' else 'no_screen' if screen=='off' else 'screen_'+screen
    values=[get(k,name)['signals']['lof'] for k in range(4)]
    screen_rows.append(('None' if screen=='off' else '+'.join(screen))+' & '+' & '.join(f"{m['contrast_peak_db']:.2f}" for m in values)+' & '+str(sum(m['global_peak_distance_mm']<=.5 for m in values))+r' \\')
write('table_screens.tex',screen_rows)
fields=['path','name','signal','supported','support_pct','near_support_pct','global_peak_um_s','global_peak_x_mm','global_peak_depth_mm','near_peak_um_s','background_peak_um_s','background_p99_um_s','contrast_peak_db','contrast_p99_db']
with (HERE/'quantitative_results.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
    for r in runs:
        for case,m in r['signals'].items():w.writerow({**{k:r[k] for k in ['path','name']},'signal':case,**{k:m[k] for k in fields[3:]}})
exe=_executable('SAFTBatch')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
manifest={'scope':'Full weld only; identical L/S/R-screened protocol for all four paths; no partial-weld or reduced-receiver results',
 'dataset_fingerprint':'c3c1a5b2b587a17f6336ccfc20c46ec6c093f084d53c4ee50d6a00013a62a07b',
 'input_files':{str(p.relative_to(ROOT)):sha(p) for p in [HERE/'common_protocol.json',DATA/'dataset.json',exe]},
 'reconstructions':len(runs),'counts_by_grid':{str(s):sum(r['selection']['spacing_mm']==s for r in index) for s in [.05,.1,.2,.25]},
 'common_protocol_verified':analysis['common_protocol_verified'],
 'identical_pulse_windows_verified':analysis['identical_pulse_windows_verified'],
 'numerical_checks':analysis['numerical_checks'],
 'figure_evidence':json.loads((OUT/'figure_evidence.json').read_text())}
(HERE/'reproducibility.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(manifest,indent=2))
