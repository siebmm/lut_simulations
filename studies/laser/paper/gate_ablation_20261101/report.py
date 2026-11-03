"""Summarise the complete, fixed-control gate audit and publish paired checks."""
from pathlib import Path
import csv
import json
import shutil
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import run as audit

OUT, PAPER = audit.OUT, audit.PAPER
sys.path.insert(0, str(PAPER))
import make_artifacts as artwork


def main():
    p=audit.protocol()
    audit.summarise(p)
    rows=json.loads((OUT/'comparison.json').read_text())['rows']
    assert len(rows)==882
    assert len({(r['uid'],r['control'],r['case']) for r in rows})==882
    exact_full=all(r['exactly_equal'] for r in rows if r['kind']=='full' and r['control']=='wall_off')
    assert exact_full
    for kind,uid in audit.PRINCIPAL.items():
        chosen=[r for r in rows if r['uid']==uid]
        assert len(chosen)==9 and all(r['support_changed_pp']==0 for r in chosen)
        assert all(r['exactly_equal'] for r in chosen if r['control']=='wall_off')
    flat=[]
    for r in rows:
        flat.append({k:r[k] for k in ['kind','uid','path','variant','control','case','principal',
                     'exactly_equal','common_relative_l2_pct','common_max_change_pct_peak','support_changed_pp']}
                    | {f'{prefix}_{key}':r[prefix][key] for prefix in ['before','after','delta']
                       for key in ['cmax_db','c99_db','near_support_pct','support_pct']})
    with (OUT/'comparison.csv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    tables=PAPER/'tables/comsol'
    mainrows=[]
    for kind,uid in audit.PRINCIPAL.items():
        selected=[r for r in rows if r['uid']==uid and r['case']=='lof']
        r=next(r for r in selected if r['control']=='wall_off')
        b=r['before']
        lines=[f"5\\% / 1.5 mm & {b['cmax_db']:.2f} & 0.00 & {b['near_support_pct']:.1f}"+r' \\']
        for key,label in [('wall_off',r'0\% / 1.5 mm'),('separation_off',r'5\% / 0 mm'),('both_off',r'0\% / 0 mm')]:
            r=next(r for r in selected if r['control']==key)
            lines.append(label+f" & {r['after']['cmax_db']:.2f} & {r['common_relative_l2_pct']:.2f} & {r['after']['near_support_pct']:.1f}"+r' \\')
        (tables/f'{kind}_gate_controls.tex').write_text('\n'.join(lines)+'\n',encoding='utf-8')
        r=next(r for r in selected if r['control']=='separation_off')
        mainrows.append(f"{kind.capitalize()}: {r['before']['cmax_db']:.2f} to {r['after']['cmax_db']:.2f} dB; relative L2 image change {r['common_relative_l2_pct']:.2f}%; unchanged support.")
    report=[
        'WALL RESPONSE AND SOURCE/RECEIVER EXCLUSION: CONTROLLED AUDIT',
        '',
        'All other controls stay fixed at their published equal-receiver values.',
        '98 baseline configurations: 52 full weld and 46 partial weld. Each is repeated with wall gate off, exclusion off, and both off.',
        '294 ablated configurations / 882 paired signal-case images (LOF, pristine and waveform difference). Principal baselines are also independently rerun in the same native process.',
        '',
        'DECISION: retain the 5% wall gate and 1.5 mm exclusion in the common study protocol. Neither is negligible over the complete tested family.',
        '',
        'WALL GATE: all 52 checked full-weld configurations and the principal partial S|SS image are exactly unchanged when only this gate is disabled. It is not the cause of those principal indications. It changes multiple-reflection partial-weld images substantially, so redundancy for a selected image is not a reason to remove the gate globally.',
        'On the 0.2 mm partial-weld survey, LOF-only relative L2 changes include 49.81% for S|SSS, 7.73% for SS|SS, and 39.84% for SS|SSS. These are image changes, not detection gains.',
        '',
        'SOURCE/RECEIVER EXCLUSION: removing the 1.5 mm exclusion changes both principal images:',
        *mainrows,
        'The exclusion lies outside the prescribed 1.425 mm source footprint. Its exact value is a declared geometric safeguard, not a validated universal near-field boundary or a defect-fitted optimum.',
        '',
        'Relative L2 = 100*norm(test-baseline)/norm(baseline) on common supported pixels, using untruncated amplitudes. Mask changes are reported separately in the complete CSV/JSON. Near-flaw contrast uses the existing fixed evaluation regions; it did not choose which configurations to test.',
        'Disabling both controls tests interaction. On the two principal images it gives exactly the same amplitude arrays as disabling exclusion alone.',
        '',
        'Surface terms: weld surface = exposed top of the weld fill (partial-weld floor at 16.4 mm). Side surfaces = flat top surfaces beside the weld at depth zero. Internal identifiers plate/floor are preserved only for reproducibility.',
    ]
    (OUT/'findings.txt').write_text('\n'.join(report)+'\n',encoding='utf-8')
    verification=dict(status='passed',configurations=98,full=52,partial=46,
        ablated_configurations=294,paired_signal_cases=882,full_wall_gate_images_exactly_equal=156,
        principal_wall_gate_images_exactly_equal=6,primary_support_masks_unchanged=True,
        decisions=dict(wall_gate=.05,min_offset_mm=1.5),
        protocol_sha256=audit.s.sha(audit.HERE/'protocol.json'),
        comparison_sha256=audit.s.sha(OUT/'comparison.json'),
        original_arrays_unchanged=True)
    for cfg in p['configurations']:
        assert audit.s.sha(audit.baseline_folder(cfg)/'arrays.npz')==p['baseline_sha256'][cfg['uid']]
    for kind,uid in audit.PRINCIPAL.items():
        with np.load(OUT/kind/(uid+'__separation_off')/'arrays.npz') as a,np.load(OUT/kind/(uid+'__both_off')/'arrays.npz') as b:
            for case in ('lof','reference','difference'):
                np.testing.assert_array_equal(a[case+'_2_amplitude'],b[case+'_2_amplitude'])
    audit.s.write(OUT/'verification.json',verification)
    destination=PAPER/'data/gate_audit';destination.mkdir(exist_ok=True)
    for name in ['comparison.json','comparison.csv','findings.txt','verification.json']:
        shutil.copy2(OUT/name,destination/name)
    shutil.copy2(audit.HERE/'protocol.json',destination/'protocol.json')
    research=PAPER/'research/studies/laser/paper/gate_ablation_20261101';research.mkdir(parents=True,exist_ok=True)
    for name in ['run.py','report.py','protocol.json']:
        shutil.copy2(audit.HERE/name,research/name)
    artwork.DATASETS={kind:audit.s.dataset(kind) for kind in ('full','partial')}
    selected=[next(r for r in p['configurations'] if r['uid']==audit.PRINCIPAL[k]) for k in ('full','partial')]
    selected.append(next(r for r in p['configurations'] if r['uid']=='numerical_partial_baseline__1-3_SSSS'))
    fig=plt.figure(figsize=(10,9.5),layout='constrained')
    subfigures=fig.subfigures(3,1,height_ratios=[1.8,1,1])
    for sub,cfg in zip(subfigures,selected):
        axs=sub.subplots(1,4,sharex=True,sharey=True)
        kind=cfg['kind'];ref=10 if kind=='full' else 40;db=55 if cfg['stage']=='numerical' else 45
        sub.suptitle(f"{kind.capitalize()} weld: {artwork.label(cfg['path'])}; {cfg['selection']['spacing_mm']:g} mm grid",fontsize=10)
        for ax,key,label in zip(axs,['baseline','wall_off','separation_off','both_off'],['Both controls on','Wall gate off','Exclusion off','Both off']):
            folder=audit.baseline_folder(cfg) if key=='baseline' else OUT/kind/(cfg['uid']+'__'+key)
            with np.load(folder/'arrays.npz') as a:
                # Native grid axes are in mm; flattened ray points are in m.
                x,z=a['x'],a['z']
                xx,zz=np.meshgrid(x,z)
                np.testing.assert_allclose(np.column_stack([xx.ravel(),zz.ravel()])/1000,
                                           a['points'].reshape(-1,2),rtol=0,atol=1e-14)
                amp=a['lof_2_amplitude'].reshape(len(z),len(x))*1e6
                v=20*np.log10(np.maximum(amp/ref,1e-15))
                dx=np.diff(x).mean();dz=np.diff(z).mean()
                im=ax.imshow(v,cmap=artwork.CM,vmin=-db,vmax=0,interpolation='nearest',aspect='equal',
                    extent=[x[0]-dx/2,x[-1]+dx/2,z[-1]+dz/2,z[0]-dz/2])
                ax.set(xlim=(x[0],x[-1]),ylim=(z[-1],z[0]),xlabel='$x$ (mm)',title=label)
                ax.set_xticks([-6,0,6] if kind=='full' else [-7.5,0,7.5])
                ax.set_yticks([5,10,15,20,25] if kind=='full' else [17,21,25])
                artwork.outline(ax,kind)
        axs[0].set_ylabel('Depth below OD (mm)')
        sub.colorbar(im,ax=axs,shrink=.7,pad=.025,ticks=[-db,-20,0],label=f'dB; 0 dB = {ref:g} '+r'$\mu$m/s')
    fig.savefig(OUT/'gate_comparison.png',dpi=300,bbox_inches='tight',pad_inches=.08)
    fig.savefig(OUT/'gate_comparison.svg',bbox_inches='tight',pad_inches=.08)
    plt.close(fig)
    for name in ['gate_comparison.png','gate_comparison.svg']:
        shutil.copy2(OUT/name,destination/name)
    print(json.dumps(verification,indent=2))


if __name__=='__main__':
    main()
