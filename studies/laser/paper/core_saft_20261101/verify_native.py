"""Check two principal images against a fresh native LOF-only reconstruction."""
from pathlib import Path
import json
import sys
import numpy as np

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[3]
sys.path.insert(0,str(HERE.parent/'filter_study_20261101'))
import study as s

index=json.loads((HERE/'data/images.json').read_text())
out=ROOT/'studies/laser/results/core_saft_20261101/native_check'
checks=[]
for kind,path in [('full','1-1_LL'),('partial','1-2_SSS')]:
    row=next(r for r in index if r['kind']==kind and r['id']==f'baseline__{path}__fine')
    config=dict(row['selection'],cases=['lof'])
    output=out/kind
    s.execute(kind,{'requests':[{'id':'native','selection':config}]},output)
    s.compact(output/'native')
    with np.load(output/'native/arrays.npz') as native, np.load(HERE/'data'/(row['uid']+'.npz')) as archived:
        expected=archived['amplitude'].ravel();actual=native['lof_2_amplitude']
        np.testing.assert_allclose(actual,expected,rtol=2e-9,atol=1e-15,equal_nan=True)
        np.testing.assert_array_equal(np.isfinite(actual),np.isfinite(expected))
        np.testing.assert_array_equal(native['lof_2_count'],archived['count'].ravel())
        checks.append(dict(kind=kind,path=path,pixels=len(actual),max_absolute_difference_m_s=float(np.nanmax(abs(actual-expected)))))
s.write(HERE/'verification/native_check.json',dict(status='passed',executable_sha256=s.sha(s.EXE),checks=checks,signal_cases=['lof']))
print(checks)
