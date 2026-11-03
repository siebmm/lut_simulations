"""Check all four native export interfaces agree on the same partial-weld pixels."""
import json
from pathlib import Path
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'studies/laser/src'))
from standalone import executable, native_env

def main():
    output = ROOT/'output/interface_check'
    output.mkdir(parents=True, exist_ok=True)
    config = json.loads((ROOT/'saft/examples/partial_weld.json').read_text())['requests'][0]['selection']
    config.update(x_min_mm=2.65,x_max_mm=2.85,depth_min_mm=19.8,depth_max_mm=20.0,spacing_mm=.1,threads=2)
    request = output/'request.json'
    request.write_text(json.dumps({'requests':[{'id':'check','selection':config}]},indent=2)+'\n')
    for name in ['SAFTFilterStudy','SAFTReview','SAFTAdmissionAudit','SAFTBatch']:
        target=output/name
        if name=='SAFTBatch':
            command=[str(executable(name)),'--repository',str(ROOT),'--dataset',str(ROOT/'data/partial_weld'),
                     '--compute',str(request),'--output',str(target)]
        else:
            command=[str(executable(name)),str(ROOT),str(ROOT/'data/partial_weld'),str(request),str(target)]
        subprocess.run(command,env=native_env(),check=True)
    cases=[]
    with np.load(output/'SAFTBatch/0/images.npz') as batch:
        for case in ['lof','reference','difference']:
            expected=np.fromfile(output/f'SAFTFilterStudy/check/{case}_2_amplitude.f64',dtype='<f8')
            if not len(expected) or not np.isfinite(expected).any():
                raise AssertionError(f'No supported pixels for {case}')
            review=np.fromfile(output/f'SAFTReview/check/{case}_2_amplitude.f64',dtype='<f8')
            admission=np.fromfile(output/f'SAFTAdmissionAudit/check/{case}_amplitude.f64',dtype='<f8')
            for other in [review,admission,batch[case+'_maps'][2].ravel()]:
                np.testing.assert_allclose(other,expected,rtol=2e-12,atol=1e-15,equal_nan=True)
            cases.append({'case':case,'pixels':len(expected),'finite':int(np.isfinite(expected).sum())})
    report={'status':'passed','interfaces':4,'checks':cases}
    (output/'check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report))

if __name__=='__main__':
    main()
