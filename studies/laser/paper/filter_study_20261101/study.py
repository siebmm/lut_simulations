"""Fixed-protocol signal-processing study using the production C++ ray/sum engine.

Commands: prepare, raw, verify, run, summarise. No COMSOL solve is performed.
The external signal banks are temporary, reproducible products of verified CSVs.
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
import argparse, copy, hashlib, json, os, subprocess, sys, time
import numpy as np
from scipy import signal, fft

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[3]; STUDY=ROOT/'studies/laser'
OUT=STUDY/'results/filter_study_20261101'; LOCAL=STUDY/'build/filter_study_20261101'
EXE=_executable('SAFTFilterStudy'); DT=1e-8; NS=5001


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def dataset(kind):return json.loads((_data_root()/f'{kind}_weld/dataset.json').read_text())


PATHS=['1-1_SS','1-1_LL','1-1_LS','1-1_SL','1-2_SSS','2-1_SSS','1-2_LLL','2-1_LLL',
       '1-3_SSSS','2-2_SSSS','3-1_SSSS','1-3_LLLL','2-2_LLLL','3-1_LLLL',
       '2-3_SSSSS','3-2_SSSSS','2-3_LLLLL','3-2_LLLLL']
CORE=['1-1_SS','1-1_LL','1-2_SSS','2-1_SSS']


def variants():
    def fir(name,low=.5,high=5,half=328,guard=None,rest=False):
        return dict(name=name,method='analytic_fir',low_MHz=low,high_MHz=high,
                    half_samples=half,support_half_samples=half if guard is None else guard,
                    known_initial_rest=rest)
    return [fir('baseline'),fir('band1_5',low=1),fir('band05_3',high=3),
            fir('lowpass5_long',low=0),fir('bp_short_matched_guard',half=128,guard=328),
            fir('bp_short',half=128),fir('lp_short_matched_guard',low=0,half=128,guard=328),
            fir('lp_short',low=0,half=128),fir('lp_shorter',low=0,half=64),
            fir('baseline_initial_rest',rest=True),fir('lp_short_initial_rest',low=0,half=128,rest=True),
            dict(name='butter_lowpass',method='butter_analytic',low_MHz=0,high_MHz=5,
                 order=4,support_half_samples=328,known_initial_rest=False),
            dict(name='broadband_analytic',method='broadband_analytic',support_half_samples=328,
                 known_initial_rest=False),
            dict(name='raw_signed',method='raw_signed',support_half_samples=0,known_initial_rest=False)]


def kernel(spec):
    half=spec['half_samples'];n=np.arange(-half,half+1);width=(spec['high_MHz']-spec['low_MHz'])*1e6
    centre=(spec['high_MHz']+spec['low_MHz'])*.5e6
    base=2*width*DT*np.sinc(width*n*DT)*np.kaiser(2*half+1,.1102*(55-8.7))
    return base*np.exp(2j*np.pi*centre*n*DT)*2/base.sum()


def transform(a,spec,fft_size=32768):
    if spec['method']=='analytic_fir':
        h=kernel(spec);half=spec['half_samples']
        # Explicit zero continuation outside the recording. Complete-support
        # runs reject every sample that uses it; initial-rest runs admit only
        # the known pre-excitation zeros, never the unknown post-record tail.
        return signal.fftconvolve(a,h[None,:],mode='full',axes=-1)[:,half:half+NS]
    if spec['method']=='raw_signed':return a.astype(np.complex128)
    padded=np.pad(a,((0,0),(4096,fft_size-NS-4096)))
    if spec['method']=='butter_analytic':
        sos=signal.butter(spec['order'],spec['high_MHz']*1e6,fs=1/DT,output='sos')
        padded=signal.sosfiltfilt(sos,padded,axis=-1,padtype=None)
    return signal.hilbert(padded,axis=-1)[:,4096:4096+NS]


def base(kind,spacing=.2):
    c=json.loads((HERE.parent/'weld_geometry_study_20261030/locked_protocol.json').read_text())['common_controls']
    c=copy.deepcopy(c);c.update(receivers=list(range(len(dataset(kind)['audit']['receivers']))),
        processing='weighted' if kind=='partial' else 'raw',spacing_mm=spacing,
        cases=['lof','reference','difference'])
    if kind=='full':c.update(x_min_mm=-6,x_max_mm=6,depth_min_mm=5,depth_max_mm=25.25)
    return c


def prepare():
    OUT.mkdir(parents=True,exist_ok=True);LOCAL.mkdir(parents=True,exist_ok=True)
    protocol=dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
      intent='Compare bandwidth, filter duration, complete-support policy and no-band-restriction controls before selecting report illustrations.',
      variants=variants(),paths={'full':CORE,'partial':PATHS},
      base={kind:base(kind) for kind in ['full','partial']},
      source_admission='Fixed original 0.5-5 MHz launch/receive response curves and physical thresholds across filter variants, isolating signal processing. Curves do not silently change with the external bank.',
      screen='LSR windows remeasured on the same processed recording; 30 dB tail, 0.1 us margin. No pristine input for LOF windows. Fixed mutes/reflected screens off.',
      edge='FIR support is exact. Initial-rest variants use the model initial condition of zero displacement/velocity for t<0; the post-record tail remains excluded. FFT-Hilbert/IIR controls have noncompact support and retain the original 3.28 us edge guard; extension sensitivity is checked separately.',
      diagnostics='Broadband analytic and raw-signed controls are not calibrated above 5 MHz and cannot establish improved physical resolution. Raw-signed changes the image phase representation and is a diagnostic, not an envelope-equivalent replacement.',
      selection='Known simulated flaws; exploratory development, not blind validation. All filter/path outcomes retained. Main chapter may highlight robust positive examples, but appendices retain competing outcomes, reference controls and expanded fields.',
      confirmation='Follow promising physically defensible candidates with 0.1 mm and wider-field grids, common spatial support comparisons, and all/shoulder/central source groups. No per-pixel tuning.',
      input_fingerprints={k:dataset(k)['input_fingerprint'] for k in ['full','partial']})
    p=HERE/'protocol.json'
    if p.exists():
        previous=json.loads(p.read_text());assert previous['variants']==protocol['variants'];return
    write(p,protocol);print('Frozen',len(protocol['variants']),'variants; 4 full and 18 partial paths',flush=True)


def execute(kind,request,out):
    request_path=LOCAL/f'{kind}_request.json';write(request_path,request)
    env=_native_env()
    out.mkdir(parents=True,exist_ok=True)
    with (out/'compute.log').open('w',encoding='utf-8') as log:
        result=subprocess.run([str(EXE),str(ROOT),str(_data_root()/f'{kind}_weld'),str(request_path),str(out)],
                              env=env,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError(f'Native exit {result.returncode}: '+(out/'compute.log').read_text())


def raw():
    for kind in ['full','partial']:
        folder=LOCAL/'raw'/kind
        if not (folder/'input_verification.json').exists() or len(list(folder.glob('*.f64')))!=10:
            execute(kind,{'dump_raw':True},folder)
        print(kind,'raw recordings verified/exported',flush=True)


def bank(kind,spec):
    folder=LOCAL/'bank'/kind;folder.mkdir(parents=True,exist_ok=True)
    count=len(dataset(kind)['audit']['receivers']);energy=[]
    for case in ['lof','reference']:
        for source in range(5):
            a=np.fromfile(LOCAL/'raw'/kind/f'{case}_{source}.f64',dtype='<f8').reshape(count,NS)
            assert np.max(abs(a[:,0]))<1e-12,'Initial condition is not zero'
            transformed=transform(a,spec)
            transformed.astype('<c16').tofile(folder/f'{case}_{source}.c128')
            f=fft.fftfreq(8192,DT);spectrum=fft.fft(transformed,n=8192,axis=-1);power=np.abs(spectrum)**2
            energy.append(dict(case=case,source=source,above_5MHz_fraction=float(power[:,abs(f)>5e6].sum()/max(power.sum(),1e-300))))
    write(folder/'bank.json',dict(dataset_fingerprint=dataset(kind)['input_fingerprint'],variant=spec,energy=energy))
    write(OUT/'filter_evidence'/f'{kind}_{spec["name"]}.json',json.loads((folder/'bank.json').read_text()))
    return folder


def compact(folder):
    arrays={}
    for suffix,dtype in [('.f64','<f8'),('.c128','<c16'),('.u8','u1'),('.u64','<u8')]:
        for p in folder.glob('*'+suffix):arrays[p.stem]=np.fromfile(p,dtype=dtype)
    np.savez_compressed(folder/'arrays.npz',**arrays)
    for p in folder.iterdir():
        if p.suffix in ('.f64','.c128','.u8','.u64'):p.unlink()


def run():
    protocol=json.loads((HERE/'protocol.json').read_text())
    for kind in ['full','partial']:
        for spec in protocol['variants']:
            output=OUT/'numerical'/kind/spec['name'];paths=protocol['paths'][kind]
            if all((output/p/'arrays.npz').exists() for p in paths):continue
            folder=bank(kind,spec);jobs=[]
            for path in paths:
                c=base(kind);c['path']=path
                jobs.append(dict(id=path,selection=c,signal_processing=dict(spec,bank=str(folder)),
                                 export_windows=path==paths[0]))
            print(kind,spec['name'],'start',flush=True)
            execute(kind,{'requests':jobs},output)
            for path in paths:compact(output/path)
            print(kind,spec['name'],'done',flush=True)


def distances(points,kind):
    t=np.array(dataset(kind)['audit']['defect_benchmark_only']['tips_xy_mm']);t[:,1]=25.4-t[:,1]
    a,b=t;v=b-a;p=points*1000;u=np.clip((p-a)@v/(v@v),0,1)
    return np.linalg.norm(p-a-u[:,None]*v,axis=1)


def metrics(values,points,kind,mask=None):
    d=distances(points,kind);valid=np.isfinite(values)
    if mask is not None:valid &= mask
    near=valid&(d<=.5);back=valid&(d>1)
    def peak(m):
        if not m.any():return None
        j=int(np.argmax(np.where(m,values,-np.inf)))
        return dict(um_s=float(values[j]*1e6),x_mm=float(points[j,0]*1000),depth_mm=float(points[j,1]*1000),distance_mm=float(d[j]))
    n,b=peak(near),peak(back);p99=float(np.percentile(values[back],99)) if back.any() else 0
    ratio=lambda a,b:float(20*np.log10(a/b)) if a>0 and b>0 else None
    return dict(support_pct=float(valid.mean()*100),near_support_pct=float(valid[d<=.5].mean()*100),peak=peak(valid),near=n,background=b,
                cmax_db=ratio(n['um_s'],b['um_s']) if n and b else None,c99_db=ratio(n['um_s']*1e-6,p99) if n else None)


def summarise():
    rows=[];protocol=json.loads((HERE/'protocol.json').read_text())
    for kind in ['full','partial']:
        for spec in protocol['variants']:
            for path in protocol['paths'][kind]:
                folder=OUT/'numerical'/kind/spec['name']/path
                if not (folder/'arrays.npz').exists():continue
                with np.load(folder/'arrays.npz') as a:
                    points=a['points'].reshape(-1,2)
                    baseline=OUT/'numerical'/kind/'baseline'/path/'arrays.npz'
                    with np.load(baseline) as b:mask=np.isfinite(b['lof_2_amplitude'])
                    for case in ['lof','reference','difference']:
                        values=a[case+'_2_amplitude']
                        rows.append(dict(kind=kind,variant=spec['name'],path=path,case=case,
                                         **metrics(values,points,kind),common_baseline_support=metrics(values,points,kind,mask)))
    write(OUT/'metrics.json',dict(rows=rows))
    for kind in ['full','partial']:
        scores=[r for r in rows if r['kind']==kind and r['case']=='lof' and r['cmax_db'] is not None]
        for r in sorted(scores,key=lambda r:r['cmax_db'],reverse=True)[:15]:
            print(kind,r['variant'],r['path'],round(r['cmax_db'],2),round(r['near_support_pct'],1),round(r['peak']['distance_mm'],3),flush=True)


def verify():
    checks=[]
    for kind in ['full','partial']:
        spec=variants()[0];folder=bank(kind,spec);jobs=[]
        for path in CORE:
            c=base(kind);c.update(path=path,spacing_mm=.5)
            points=[[0,.020],[.0028,.019],[.0027,.021],[-.004,.008],[-.003,.024]]
            jobs.append(dict(id=path+'_native',selection=c,points_m=points))
            jobs.append(dict(id=path+'_external',selection=c,points_m=points,
                             signal_processing=dict(spec,bank=str(folder))))
        out=OUT/'verification'/kind;execute(kind,dict(requests=jobs),out)
        worst=0
        for path in CORE:
            for suffix in ['_amplitude.f64','_means.c128','_valid.u8','_count.f64','_screened.f64']:
                dtype='<c16' if 'means' in suffix else 'u1' if 'valid' in suffix else '<f8'
                for case in ['lof','reference','difference']:
                    a=np.fromfile(out/(path+'_native')/(case+'_2'+suffix),dtype=dtype)
                    b=np.fromfile(out/(path+'_external')/(case+'_2'+suffix),dtype=dtype)
                    np.testing.assert_allclose(a,b,rtol=2e-9,atol=1e-15,equal_nan=True)
                    if 'amplitude' in suffix or 'means' in suffix:
                        delta=np.abs(a-b);worst=max(worst,float(np.nanmax(delta)) if np.isfinite(delta).any() else 0)
            for tag in ['native','external']:compact(out/(path+'_'+tag))
        checks.append(dict(kind=kind,paths=4,points=5,signals=3,max_absolute_difference=worst))
    # Independent filter invariants at interior samples (phase and amplitude).
    t=np.arange(NS)*DT;wave=np.cos(2*np.pi*2e6*t)[None,:]
    h=kernel(variants()[0]);gain=np.sum(h*np.exp(-2j*np.pi*2e6*np.arange(-328,329)*DT))
    y=transform(wave,variants()[0]);expected=.5*gain*np.exp(2j*np.pi*2e6*t)
    np.testing.assert_allclose(y[:,500:-500],expected[None,500:-500],atol=5e-5)
    write(OUT/'verification/check.json',dict(status='passed',production_parity=checks,filter_sinusoid='passed'))
    print('Production parity and independent filter check passed',flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['prepare','raw','verify','run','summarise'])
    globals()[parser.parse_args().command]()
