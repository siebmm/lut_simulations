"""Filtering, native execution and scores shared by the fixed paper replays."""
from pathlib import Path
import hashlib,json,subprocess
import numpy as np
from scipy import signal,fft
from standalone import data_root as _data_root, executable as _executable, native_env as _native_env

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[3]; STUDY=ROOT/'studies/laser'
OUT=STUDY/'results/filter_study_20261101'; LOCAL=STUDY/'build/filter_study_20261101'
EXE=_executable('SAFTFilterStudy'); DT=1e-8; NS=5001


def write(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def dataset(kind):return json.loads((_data_root()/f'{kind}_weld/dataset.json').read_text())


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


def execute(kind,request,out):
    request_path=LOCAL/f'{kind}_request.json';write(request_path,request)
    env=_native_env()
    out.mkdir(parents=True,exist_ok=True)
    with (out/'compute.log').open('w',encoding='utf-8') as log:
        result=subprocess.run([str(EXE),str(ROOT),str(_data_root()/f'{kind}_weld'),str(request_path),str(out)],
                              env=env,stdout=log,stderr=subprocess.STDOUT)
    if result.returncode:raise RuntimeError(f'Native exit {result.returncode}: '+(out/'compute.log').read_text())


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
