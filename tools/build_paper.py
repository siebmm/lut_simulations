"""Compile the sole manuscript and record the source/PDF correspondence."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tectonic', default='tectonic', help='Tectonic executable or path')
    args = parser.parse_args()
    compiler = shutil.which(args.tectonic)
    if not compiler:
        parser.error('Install Tectonic or supply --tectonic PATH; see paper/README.md.')
    paper = ROOT/'paper'
    work = ROOT/'build/paper'
    work.mkdir(parents=True, exist_ok=True)
    result = subprocess.run([compiler, '--keep-logs', '--outdir', str(work), 'paper.tex'],
                            cwd=paper, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    (work/'build.log').write_text(result.stdout, encoding='utf-8')
    if result.returncode:
        raise SystemExit(result.stdout)
    log = (work/'paper.log').read_text(encoding='utf-8', errors='replace')
    if re.search(r'undefined|Overfull \\[hv]box|Missing character:|multiply defined', log, re.I):
        raise SystemExit('Paper has unresolved references or layout warnings; inspect build/paper/paper.log.')
    pages = re.search(r'Output written on .*?\((\d+) pages?', log)
    if not pages:
        raise SystemExit('Could not confirm the compiled page count.')
    source = paper/'paper.tex'
    text = source.read_text(encoding='utf-8')
    figures = set(re.findall(r'\\includegraphics(?:\[[^]]*\])?\{([^}]+)\}', text))
    inputs = [source, paper/'references.bib', *(paper/name for name in sorted(figures))]
    shutil.copy2(work/'paper.pdf', paper/'paper.pdf')
    report = dict(status='passed', pages=int(pages[1]), figure_assets=len(figures),
                  compiler=subprocess.check_output([compiler, '--version'], text=True).strip(),
                  input_sha256={p.relative_to(ROOT).as_posix(): sha(p) for p in inputs},
                  pdf_sha256=sha(paper/'paper.pdf'), unresolved_references=False,
                  overfull_boxes=False)
    (ROOT/'verification/paper.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(f'Built paper/paper.pdf ({pages[1]} pages); source/PDF hashes recorded in verification/paper.json.')


if __name__ == '__main__':
    main()
