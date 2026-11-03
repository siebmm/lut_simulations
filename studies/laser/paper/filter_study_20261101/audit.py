"""Independent filter/extension checks and read-only saved-model metadata audit."""
import json, sys
from zipfile import ZipFile
import numpy as np
from scipy import signal
import study as s


def main():
    evidence={'filters':[], 'extension_checks':[], 'initial_samples':[]}
    frequency=np.linspace(0,15e6,30001)
    for spec in s.variants():
        if spec['method']!='analytic_fir':continue
        h=s.kernel(spec);n=np.arange(-spec['half_samples'],spec['half_samples']+1)
        # scipy freqz uses causal indexing: magnitude is unchanged by centring.
        _,H=signal.freqz(h,worN=frequency,fs=1/s.DT)
        db=20*np.log10(np.maximum(abs(H)/2,1e-16))
        evidence['filters'].append(dict(name=spec['name'],taps=len(h),
            half_support_us=spec['half_samples']*.01,guard_us=spec['support_half_samples']*.01,
            half_power_at_nominal_edges_db=[float(db[np.argmin(abs(frequency-v*1e6))]) for v in [spec['low_MHz'],spec['high_MHz']]],
            magnitude_at_6MHz_relative_to_passband_db=float(db[np.argmin(abs(frequency-6e6))])))
    for kind in ['full','partial']:
        count=len(s.dataset(kind)['audit']['receivers'])
        for case in ['lof','reference']:
            for source in range(5):
                a=np.fromfile(s.LOCAL/'raw'/kind/f'{case}_{source}.f64',dtype='<f8').reshape(count,s.NS)
                evidence['initial_samples'].append(dict(kind=kind,case=case,source=source,max_abs=float(abs(a[:,0]).max())))
                for method in ['butter_lowpass','broadband_analytic']:
                    spec=next(v for v in s.variants() if v['name']==method)
                    b=s.transform(a,spec,32768)[:,328:-329];c=s.transform(a,spec,65536)[:,328:-329]
                    evidence['extension_checks'].append(dict(kind=kind,case=case,source=source,method=method,
                        relative_L2=float(np.linalg.norm(b-c)/np.linalg.norm(c)),
                        max_absolute_change_m_s=float(abs(b-c).max()),
                        max_change_over_global_peak=float(abs(b-c).max()/abs(c).max())))
    sys.path.insert(0,str(s.STUDY/'src'))
    from audit_saft_bandwidth import walk,settings
    evidence['partial_saved_models']=[]
    for case in ['lof','reference']:
        from standalone import ROOT as input_root
        path=input_root/'models/partial_weld'/f'{case}.mph'
        with ZipFile(path) as z:raw=z.read('smodel.json')
        nodes=list(walk(json.loads(raw)));by={n.get('modelEntityPath'):n for n in nodes}
        params={p['name']:p for n in nodes if n.get('apiClass')=='ModelParamGroup' for p in n['settings']}
        row={'case':case,'path':'models/partial_weld/'+case+'.mph', 'metadata_sha256':s.hashlib.sha256(raw).hexdigest(),
             'parameters':{k:params[k] for k in ['mesh_fmax','mesh_nlambda','mesh_factor','steel_cs','steel_cp','dt_out']},
             'physics':{k:v for k,v in settings(by['/physics/elte']).items() if k in ['order_velocity','order_strain','fmax']},
             'sizes':{}}
        for p,n in by.items():
            if p and p.startswith('/mesh/mesh1/feature/'):
                v=settings(n)
                if 'hmax' in v:row['sizes'][p]=v['hmax']
        evidence['partial_saved_models'].append(row)
    s.write(s.OUT/'verification/filter_audit.json',evidence)
    print('Saved filter and partial-mesh audit; worst extension relative L2',max(v['relative_L2'] for v in evidence['extension_checks']))

if __name__=='__main__':main()
