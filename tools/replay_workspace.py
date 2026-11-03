"""Materialize the paper replay's inputs from one canonical copy of each file."""
import json
import hashlib
import os
from pathlib import Path
import shutil
import zipfile

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT/'studies/laser/replay_inputs.json'


def contained(root, name):
    path = (root/name).resolve()
    if Path(name).is_absolute() or not path.is_relative_to(root):
        raise ValueError(f'Replay path escapes its root: {name}')
    return path


def prepare_workspace():
    work = Path(os.environ.get('SAFT_WORKSPACE', ROOT/'output/replay')).resolve()
    if ROOT.is_relative_to(work) or (work.is_relative_to(ROOT) and
            not any(work.is_relative_to(ROOT/p) for p in ('output', 'build'))):
        raise ValueError('Use a separate output directory for SAFT_WORKSPACE')
    definition = json.loads(INPUTS.read_text())
    if definition['schema'] != 1:
        raise ValueError('Unsupported replay input schema')
    study = work/'studies/laser'
    sources = [contained(ROOT/'studies/laser', name) for name in definition['modules']]
    code = sources + sorted((ROOT/'saft/src').glob('*.cpp')) + sorted((ROOT/'saft/src').glob('*.hpp'))
    code += sorted((ROOT/'saft/cli').glob('*.cpp'))
    code += [ROOT/'tools'/name for name in ('replay_workspace.py', 'reproduce_main.py', 'reproduce_figures.py')]
    fingerprint = hashlib.sha256(INPUTS.read_bytes() + b''.join(p.read_bytes() for p in code)).hexdigest()
    marker = work/'.replay-source-sha256'
    if work.exists() and any(work.iterdir()) and (not marker.exists() or marker.read_text() != fingerprint):
        raise ValueError('Replay inputs or code changed. Set SAFT_WORKSPACE to a fresh output directory.')
    for source in sources:
        target = study/source.relative_to(ROOT/'studies/laser')
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
    with zipfile.ZipFile(ROOT/'paper/evidence.zip') as evidence:
        for name, origin in definition['inputs'].items():
            target = contained(study, name)
            if 'file' in origin:
                source = contained(ROOT, origin['file'])
                payload = source.read_bytes()
            else:
                payload = evidence.read(origin['member'])
            if hashlib.sha256(payload).hexdigest() != origin['sha256']:
                raise ValueError(f'Canonical replay input checksum mismatch: {name}')
            if target.exists():
                if hashlib.sha256(target.read_bytes()).hexdigest() != origin['sha256']:
                    raise ValueError(f'Replay input changed in the workspace: {name}; use a fresh SAFT_WORKSPACE')
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(payload)
    marker.write_text(fingerprint)
    env = os.environ.copy()
    env['SAFT_REPOSITORY'] = str(ROOT)
    env['PYTHONPATH'] = str(study/'src') + os.pathsep + env.get('PYTHONPATH', '')
    return work, env
