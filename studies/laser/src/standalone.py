"""Repository-local inputs and native tools, also usable in a replay workspace."""
import os
import json
import hashlib
from pathlib import Path

ROOT = Path(os.environ.get('SAFT_REPOSITORY', Path(__file__).resolve().parents[3])).resolve()

def data_root():
    return ROOT / 'data'

def native_root():
    return ROOT / 'saft/src'

def executable(name):
    suffix = '.exe' if os.name == 'nt' else ''
    candidates = []
    if os.environ.get('SAFT_BIN'):
        candidates.append(Path(os.environ['SAFT_BIN']) / (name + suffix))
    candidates += [ROOT / p / (name + suffix) for p in
                   ['build/bin', 'build/bin/Release', 'build/native/bin', 'build/native/bin/Release']]
    return next((p for p in candidates if p.is_file()), candidates[0])

def native_env():
    env = os.environ.copy()
    if env.get('SAFT_RUNTIME'):
        env['PATH'] = env['SAFT_RUNTIME'] + os.pathsep + env.get('PATH', '')
    return env

def frozen_protocol(path, proposed, hash_key):
    """Keep scientific settings frozen while validating a regression-checked engine."""
    proposed = json.loads(json.dumps(proposed))
    if not path.exists():
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(proposed, indent=2)+'\n')
        return proposed
    archived = json.loads(path.read_text())
    assert {k:v for k,v in archived.items() if k != hash_key} == {k:v for k,v in proposed.items() if k != hash_key}
    mapping = json.loads((ROOT/'verification/engine_portability.json').read_text())
    current = mapping.get('generalized_hashes', mapping['standalone_hashes'])
    assert set(proposed[hash_key]) == set(mapping['standalone_hashes']), 'Unexpected protocol source set'
    assert proposed[hash_key] == {k: current[k] for k in proposed[hash_key]}, 'Native source changed since validation'
    for name, expected in current.items():
        assert hashlib.sha256((native_root()/name).read_bytes()).hexdigest() == expected, 'Native source changed since validation: '+name
    assert archived[hash_key] == mapping['historical_hashes'], 'Unexpected archived engine'
    (path.parent/'replay_source_hashes.json').write_text(json.dumps(current,indent=2)+'\n')
    return archived
