"""Paths within this standalone laser study."""
from pathlib import Path
REPOSITORY = Path(__file__).resolve().parents[3]
STUDY = REPOSITORY / 'studies/laser'

def repository_root(*args):
    return REPOSITORY

def study_root(name='laser', *args):
    if name != 'laser':
        raise ValueError('Only the laser study is included')
    return STUDY
