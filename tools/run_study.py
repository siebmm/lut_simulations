"""Run a retained analysis in an isolated workspace, preserving archived evidence."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]

def prepare_workspace():
    work = ROOT / 'output/workspace'
    for p in (ROOT / 'studies/laser').rglob('*'):
        if not p.is_file() or '__pycache__' in p.parts or 'build' in p.relative_to(ROOT).parts:
            continue
        dest = work / p.relative_to(ROOT)
        if not dest.exists() or p.suffix == '.py':
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p, dest)
    env = os.environ.copy()
    env['SAFT_REPOSITORY'] = str(ROOT)
    env['PYTHONPATH'] = str(work / 'studies/laser/src') + os.pathsep + env.get('PYTHONPATH', '')
    return work, env

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('script', help='Path below studies/laser/paper, e.g. filter_study_20261101/study.py')
    parser.add_argument('arguments', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    relative = Path(args.script)
    source = (ROOT / 'studies/laser/paper' / relative).resolve()
    if not source.is_relative_to(ROOT / 'studies/laser/paper') or not source.is_file():
        parser.error('Choose an included analysis script below studies/laser/paper.')
    # Existing generated arrays and protocols survive subsequent commands.
    # Seed files are copied only once so later runs can build on earlier ones.
    work, env = prepare_workspace()
    target = work / source.relative_to(ROOT)
    subprocess.run([sys.executable, str(target), *args.arguments], cwd=work, env=env, check=True)

if __name__ == '__main__':
    main()
