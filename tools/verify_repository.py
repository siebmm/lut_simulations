"""Check the standalone snapshot, original model/recording hashes and paper assets."""
import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path
import re
import zipfile
from replay_workspace import contained

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT/'verification/standalone_manifest.json'

def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()

def files():
    for base in ['data','models','paper','saft','studies','tools','tests','verification']:
        for p in (ROOT/base).rglob('*'):
            if p.is_file() and not {'build','__pycache__'} & set(p.relative_to(ROOT).parts) and p != MANIFEST:
                yield p
    for name in ['README.md','requirements.txt','CMakeLists.txt','.gitattributes','.gitignore']:
        yield ROOT/name

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write-manifest', action='store_true', help='Record an intentionally updated snapshot.')
    args = parser.parse_args()
    failures = []
    # Fail clearly before attempting to parse an unhydrated archive or array.
    pointers = []
    for p in files():
        if p.is_file() and p.stat().st_size < 1024:
            with p.open('rb') as stream:
                if stream.read(42).startswith(b'version https://git-lfs.github.com/spec/v1'):
                    pointers.append(p.relative_to(ROOT).as_posix())
    if pointers:
        raise SystemExit('Git LFS content is missing. Run git lfs install and git lfs pull, then retry.\n'
                         + '\n'.join(pointers))
    models = json.loads((ROOT/'models/checksums.json').read_text())['models']
    model_hashes = {e['sha256'] for e in models}
    for entry in models:
        p = ROOT/'models'/entry['path']
        if not p.is_file() or sha(p) != entry['sha256']: failures.append(str(p))
    recordings = 0
    for kind in ['full','partial']:
        folder = ROOT/f'data/{kind}_weld'
        d = json.loads((folder/'dataset.json').read_text())
        assert d['schema']==1 and d['kind']=='SAFT acquisition dataset'
        assert d['audit']['status']=='passed' and len(d['audit']['sources'])==5
        assert len(d['audit']['inputs'])==60
        for entry in d['audit']['inputs']:
            p = folder/entry['path'];recordings += 1
            if not p.is_file() or sha(p)!=entry['sha256']: failures.append(str(p))
        for entry in d['model_evidence']:
            if entry['sha256'] not in model_hashes: failures.append('Missing model evidence: '+entry['path'])
    for p in list((ROOT/'studies').rglob('*.py'))+list((ROOT/'models').rglob('*.py'))+list((ROOT/'tools').glob('*.py'))+list((ROOT/'tests').glob('*.py')):
        ast.parse(p.read_text(encoding='utf-8-sig'), filename=str(p))
    study = ROOT/'studies/laser'
    replay = json.loads((study/'replay_inputs.json').read_text())
    assert replay['schema'] == 1
    modules = {contained(study, name) for name in replay['modules']}
    if modules != set(study.rglob('*.py')):
        failures.append('Replay helper modules must match replay_inputs.json exactly')
    canonical = set()
    with zipfile.ZipFile(ROOT/'paper/evidence.zip') as evidence:
        for name, origin in replay['inputs'].items():
            contained(study, name)
            if 'file' in origin:
                source = contained(ROOT, origin['file'])
                canonical.add(source)
                actual = sha(source)
            else:
                actual = hashlib.sha256(evidence.read(origin['member'])).hexdigest()
            if actual != origin['sha256']:
                failures.append('Replay input checksum: '+name)
    allowed_study = modules | canonical | {study/'README.md', study/'replay_inputs.json'}
    retained_study = {p for p in files() if p.is_relative_to(study)}
    if retained_study != {p for p in allowed_study if p.is_relative_to(study)}:
        failures.append('Study tree contains missing or unused replay files')
    content = {}
    for p in files():
        digest = sha(p)
        if digest in content:
            failures.append(f'Duplicate file content: {p.relative_to(ROOT)} and {content[digest]}')
        content[digest] = p.relative_to(ROOT).as_posix()
    graphics = set()
    cites = set()
    for name in ['paper.tex']:
        text = (ROOT/'paper'/name).read_text(encoding='utf-8-sig')
        graphics.update(re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', text))
        cites.update(k.strip() for group in re.findall(r'\\cite\w*\{([^}]+)\}',text) for k in group.split(','))
    for name in graphics:
        if not (ROOT/'paper'/name).is_file(): failures.append('Figure: '+name)
    bib = '\n'.join(p.read_text(encoding='utf-8-sig') for p in (ROOT/'paper').glob('*.bib'))
    missing_cites = sorted(cites-set(re.findall(r'@\w+\s*\{\s*([^,]+)',bib)))
    failures.extend('Missing citation: '+key for key in missing_cites)
    labels = set(re.findall(r'\\label\{([^}]+)\}', text))
    refs = set(re.findall(r'\\(?:ref|eqref|pageref)\{([^}]+)\}', text))
    failures.extend('Missing cross-reference: '+key for key in sorted(refs-labels))
    for extension, expected in [('tex', 'paper.tex'), ('pdf', 'paper.pdf'), ('bib', 'references.bib')]:
        if {p.name for p in (ROOT/'paper').glob('*.'+extension)} != {expected}:
            failures.append('paper/ must contain only the standalone manuscript, PDF and bibliography')
    actual_figures = {p.relative_to(ROOT/'paper').as_posix() for p in (ROOT/'paper/figures').rglob('*') if p.is_file()}
    if actual_figures != graphics:
        failures.append('Paper figure assets must match the manuscript exactly')
    paper_build = json.loads((ROOT/'verification/paper.json').read_text())
    expected_inputs = {'paper/paper.tex', 'paper/references.bib'} | {'paper/'+p for p in graphics}
    if paper_build['status'] != 'passed' or set(paper_build['input_sha256']) != expected_inputs:
        failures.append('Paper build receipt does not cover all manuscript inputs')
    for name, digest in paper_build['input_sha256'].items():
        if sha(ROOT/name) != digest: failures.append('Paper changed since PDF build: '+name)
    if sha(ROOT/'paper/paper.pdf') != paper_build['pdf_sha256']:
        failures.append('Paper PDF does not match its build receipt')
    with zipfile.ZipFile(ROOT/'paper/evidence.zip') as z:
        assert z.testzip() is None
        with z.open('evidence/review/review/atlas_metrics.csv') as stream:
            rows = len(list(csv.DictReader(line.decode('utf-8-sig') for line in stream)))
        assert rows==20412
    if args.write_manifest:
        if failures:
            raise SystemExit('Refusing to record invalid scientific inputs:\n'+'\n'.join(failures))
        entries = [{'path':p.relative_to(ROOT).as_posix(),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(files())]
        MANIFEST.write_text(json.dumps({'schema':1,'files':entries},indent=2)+'\n')
    manifest = json.loads(MANIFEST.read_text())
    actual_names = {p.relative_to(ROOT).as_posix() for p in files()}
    if actual_names != {e['path'] for e in manifest['files']}:
        failures.append('Snapshot file list differs from the repository; review additions/removals before updating it')
    for entry in manifest['files']:
        p = ROOT/entry['path']
        if not p.is_file() or p.stat().st_size!=entry['bytes'] or sha(p)!=entry['sha256']:
            failures.append(entry['path'])
    if failures: raise SystemExit('FAILED:\n'+'\n'.join(failures))
    report = {'status':'passed','models':len(models),'recordings':recordings,
              'figure_assets':len(graphics),'atlas_rows':rows,'snapshot_files':len(manifest['files']),
              'snapshot_bytes':sum(e['bytes'] for e in manifest['files']),
              'paper_pages':paper_build['pages'], 'missing_bibliography_keys':missing_cites}
    (ROOT/'output').mkdir(exist_ok=True)
    (ROOT/'output/integrity_check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__ == '__main__':
    main()
