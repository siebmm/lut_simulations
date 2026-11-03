"""Make a portable numerical evidence set without duplicating raw simulations."""

# Locate portable scientific inputs without depending on the source checkout.
import sys as _sys
from pathlib import Path as _Path
for _parent in _Path(__file__).resolve().parents:
    if (_parent / 'studies/laser/src/standalone.py').is_file():
        _sys.path.insert(0, str(_parent / 'studies/laser/src'))
        break
from standalone import data_root as _data_root, native_root as _native_root, executable as _executable, native_env as _native_env
import analyse as a, confirm as c
from analyse import *

def build():
    dest=HERE/'data';dest.mkdir(exist_ok=True);rows=[]
    for kind in ['full','partial']:
        write(dest/f'{kind}_dataset.json',old.dataset(kind))
        entries=[('baseline',r,old.load(kind,r['path'])) for r in old.index(kind) if r['name']=='baseline']
        support=[r for r in old.index(kind) if (kind=='full' and r['path'] in old.CORE and r['name'] in ['overview','coarse','fine','exclude_repairs']) or (kind=='partial' and r['path']=='1-2_SSS' and r['name'] in ['overview','coarse','fine','screen_reflected'])]
        entries += [('previous_check',r,old.load(kind,r['path'],r['name'])) for r in support]
        entries += [('followup',r,a.load(kind,r['path'],r['name'])) for r in a.index(kind)]
        entries += [('confirmation',r,c.load(kind,r['path'],r['name'])) for r in c.index(kind)]
        for stage,r,arr in entries:
            uid=f"{kind}__{r['path']}__{r['name']}";payload={k:arr[k] for k in ['x_mm','depth_mm','points_m']}
            metrics={}
            for case in r['selection']['cases']:
                payload[case]=arr[case+'_maps'][2]
                metrics[case]=old.metrics(payload[case],arr['points_m'],kind)
            # Retain per-source support and complex means for independent source-subset checks.
            if 'lof_maps' in arr:
                payload['lof_means']=arr['lof_means'][2];payload['lof_valid']=arr['lof_valid'][2]
            np.savez_compressed(dest/(uid+'.npz'),**payload)
            rows.append({'uid':uid,'stage':stage,**{k:r[k] for k in ['kind','path','name','selection']},'metrics':metrics,'sha256':sha(dest/(uid+'.npz'))})
    write(dest/'index.json',rows)
    write(dest/'source_subsets.json',json.loads((OUT/'source_subsets.json').read_text()))
    write(dest/'weighted_source_subsets.json',json.loads((c.DEST/'weighted_source_subsets.json').read_text()))
    write(dest/'confirmation_checks.json',json.loads((c.DEST/'checks.json').read_text()))
    write(dest/'locked_protocol.json',json.loads((PREV/'locked_protocol.json').read_text()))
    for n in ['geometry_verification.json','reproducibility.json']:
        p=old.OUT/n
        if p.exists():write(dest/('previous_'+n),json.loads(p.read_text()))
    exe=_executable('SAFTBatch')
    write(dest/'provenance.json',{'native_executable':exe.name,'native_executable_sha256':sha(exe),
        'native_cpp_sha256':sha(_native_root()/'saft.cpp'),
        'datasets':{k:{'input_fingerprint':old.dataset(k)['input_fingerprint'],'manifest_sha256':sha(old.DATA[k]/'dataset.json')} for k in ['full','partial']},
        'scope':'Reconstruction scores from existing deterministic FE records, not new FE solves. Known-case follow-up, not independent optimization validation.',
        'new_configurations':len(a.index('full'))+len(a.index('partial'))+len(c.index('full'))+len(c.index('partial')),
        'new_case_images':sum(len(r['selection']['cases']) for k in ['full','partial'] for r in a.index(k)+c.index(k)),
        'archived_configurations':len(rows),'archived_case_images':sum(len(r['metrics']) for r in rows),
        'units':'x_mm/depth_mm in mm; points_m in m; complex means and amplitude scores in m/s. NaN means unsupported.',
        'raw_data_included':False})
    sys.path.insert(0,str(ROOT/'studies/laser/src'))
    from saft_apertures import directivity,ANGLES
    m=old.dataset('partial')['material'];curves=directivity(('broad',.5,5),m['cl_m_s'],m['cs_m_s'])
    np.savez_compressed(dest/'directivity.npz',angle_deg=ANGLES,**{mode+'_'+field:curves[mode][field] for mode in ['L','S'] for field in ['source','vy']})
    print('Evidence:',len(rows),'configurations;',sum(len(r['metrics']) for r in rows),'case images')
if __name__=='__main__':build()
