"""Compute the local thermoelastic and dipole sources, then compare their waveforms."""
import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'models/laser_generation'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--comsol', help='COMSOL bin directory (containing comsolbatch and comsolcompile)')
    parser.add_argument('--fluence', type=float, default=10., help='Peak absorbed fluence in J/m^2 (default: 10)')
    parser.add_argument('--output', type=Path, default=ROOT/'output/laser_generation', help='Empty directory for this run')
    parser.add_argument('--rebuild', action='store_true', help='Construct the geometry and physics from Java instead of the saved original model')
    parser.add_argument('--prepare-only', action='store_true', help='Prepare a portable COMSOL model for manual inspection, without solving')
    parser.add_argument('--cores', type=int, default=4)
    args = parser.parse_args()
    if not math.isfinite(args.fluence) or args.fluence <= 0 or args.cores < 1:
        parser.error('Use positive fluence and core count.')
    work = args.output.resolve()
    if ROOT.is_relative_to(work) or (work.is_relative_to(ROOT) and not work.is_relative_to(ROOT/'output')):
        parser.error('Choose an output directory under output/ or outside this repository.')
    if work.exists() and any(work.iterdir()):
        parser.error('Choose an empty --output directory for each simulation.')
    if os.name == 'nt':
        folder = (Path(args.comsol) if args.comsol else Path(r'C:\Program Files\COMSOL\COMSOL64\Multiphysics\bin\win64')).resolve()
        compile_cmd = [str(folder/'comsolcompile.exe')]
        batch_cmd = [str(folder/'comsolbatch.exe')]
    else:
        executable = str(Path(args.comsol)/'comsol') if args.comsol else shutil.which('comsol')
        if not executable:
            parser.error('COMSOL was not found; supply --comsol with its bin directory.')
        executable = str(Path(executable).resolve())
        compile_cmd, batch_cmd = [executable, 'compile'], [executable, 'batch']
    for command in (compile_cmd, batch_cmd):
        if not Path(command[0]).is_file():
            parser.error('COMSOL was not found; supply --comsol with its bin directory.')
    (work/'models').mkdir(parents=True)
    (work/'exports').mkdir()
    (work/'run.properties').write_text(f'fluence_J_m2={args.fluence:.17g}\nprepare_only={str(args.prepare_only).lower()}\n')
    shutil.copy2(SOURCE/'laser_source_calibration.mph', work/'models/laser_source_calibration.mph')

    def execute(name, marker):
        shutil.copy2(SOURCE/(name+'.java'), work/(name+'.java'))
        print('COMSOL:', name, flush=True)
        with (work/(name+'_compile.log')).open('w') as log:
            subprocess.run(compile_cmd+[name+'.java'], cwd=work, stdout=log, stderr=subprocess.STDOUT, check=True)
        if not (work/(name+'.class')).is_file():
            raise RuntimeError(f'COMSOL did not compile {name}; see {work}/{name}_compile.log')
        with (work/(name+'_console.log')).open('w') as log:
            result = subprocess.run(batch_cmd+['-inputfile', name+'.class', '-nosave', '-batchlog', name+'.log', '-np', str(args.cores)],
                                    cwd=work, stdout=log, stderr=subprocess.STDOUT)
        logs = '\n'.join(p.read_text(errors='replace') for p in (work/(name+'.log'),work/(name+'_console.log')) if p.exists())
        if 'License error:' in logs:
            raise SystemExit(f'COMSOL could not obtain a licence. Check the licence-server connection, then use a fresh output folder. Details: {work}/{name}.log')
        if result.returncode != 0 or marker not in logs:
            raise RuntimeError(f'COMSOL did not finish {name}; inspect the logs in {work}')

    if args.rebuild:
        execute('BuildLaserCalibration', 'CALIBRATION_BUILD_OK')
    execute('RunLaserCalibration', 'LASER_MODEL_PREPARED' if args.prepare_only else 'LASER_CALIBRATION_COMPLETE')
    if args.prepare_only:
        print('Prepared model:', work/'models/laser_source_calibration.mph')
    else:
        subprocess.run([sys.executable, str(SOURCE/'analyze.py'), '--exports', str(work/'exports'), '--output', str(work)], check=True)
        print('Solved model, waveforms, temperature comparison and calibration:', work)


if __name__ == '__main__':
    main()
