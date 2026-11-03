"""Repository-local inputs and native tools, also usable in a replay workspace."""
import os
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
