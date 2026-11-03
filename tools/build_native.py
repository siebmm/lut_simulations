"""Build the four headless tools with an existing Qt qmake/compiler installation."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--qmake', default='qmake')
    parser.add_argument('--make', default='mingw32-make' if os.name == 'nt' else 'make')
    parser.add_argument('--jobs', type=int, default=4)
    args = parser.parse_args()
    qmake = Path(shutil.which(args.qmake) or args.qmake).resolve()
    make = Path(shutil.which(args.make) or args.make).resolve()
    if not qmake.is_file() or not make.is_file():
        parser.error('Supply --qmake and --make from compatible Qt/compiler installations.')
    env = os.environ.copy()
    env['PATH'] = os.pathsep.join([str(qmake.parent), str(make.parent), env.get('PATH', '')])
    build = ROOT / 'build/qmake'
    build.mkdir(parents=True, exist_ok=True)
    common = ['saft.cpp', 'screening.cpp', 'dataset_geometry.cpp', 'groove_paths.cpp']
    for target, cli in [('SAFTFilterStudy', 'filter'), ('SAFTReview', 'review'),
                        ('SAFTAdmissionAudit', 'admission'), ('SAFTBatch', 'batch')]:
        sources = [ROOT / 'saft/src' / name for name in common] + [ROOT / f'saft/cli/{cli}_cli.cpp']
        quote = lambda p: '"' + p.as_posix() + '"'
        project = '\n'.join(['QT += core', 'QT -= gui', 'CONFIG += console c++17 release',
            'CONFIG -= app_bundle debug debug_and_release', 'DEFINES += QT_NO_KEYWORDS',
            'TARGET = ' + target, 'DESTDIR = ../bin', 'OBJECTS_DIR = objects',
            'INCLUDEPATH += ' + quote(ROOT / 'saft/src'),
            'SOURCES += ' + ' '.join(map(quote, sources))]) + '\n'
        (build / 'tool.pro').write_text(project)
        subprocess.run([str(qmake), 'tool.pro'], cwd=build, env=env, check=True)
        subprocess.run([str(make), f'-j{args.jobs}'], cwd=build, env=env, check=True)
    if os.name == 'nt':
        # Local build products only. No compiler/runtime binaries are versioned.
        qtbin = Path(subprocess.check_output([str(qmake), '-query', 'QT_INSTALL_BINS'], env=env, text=True).strip())
        for name in ['Qt6Core.dll', 'libgcc_s_seh-1.dll', 'libstdc++-6.dll', 'libwinpthread-1.dll']:
            source = next((p / name for p in [qtbin, make.parent] if (p / name).is_file()), None)
            if source is None:
                raise FileNotFoundError(f'Runtime library missing: {name}')
            shutil.copy2(source, ROOT / 'build/bin' / name)
    print('Built all four tools in build/bin')

if __name__ == '__main__':
    main()
