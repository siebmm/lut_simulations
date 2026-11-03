"""Regenerate final-paper artwork in output/workspace; published PDFs are preserved."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
from run_study import ROOT, prepare_workspace

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--group', choices=['baseline','review','all'], default='all')
    args = parser.parse_args()
    work, env = prepare_workspace()
    paper = work/'studies/laser/paper'
    def run(relative):
        print('Figures:', relative, flush=True)
        subprocess.run([sys.executable, str(paper/relative)], cwd=work, env=env, check=True)
    if args.group in ('baseline','all'):
        for name in ['verify_data.py','source_calibration.py','make_artifacts.py','expand_reports.py',
                     'geometry_artwork.py','recording_artwork.py','joint_acquisition_artwork.py']:
            run('comsol_chapter_20261101/'+name)
        run('explanatory_saft_20261101/make_figures.py')
        run('explanatory_saft_20261101/source_requirement.py')
        run('explanatory_saft_20261101/source_requirement_report.py')
        # These explanatory panels reuse the archived, numerically checked data.
        source = paper/'comsol_chapter_20261101/figures/comsol'
        dest = paper/'explanatory_saft_20261101/figures/explanatory'
        for name in ['primary_full','partial_joint','full_common_paths','source_calibration']:
            p = source/(name+'.pdf')
            if p.exists(): shutil.copy2(p, dest/p.name)
        run('combined_explanatory_saft_20261101/prepare.py')
    if args.group in ('review','all'):
        run('partial_improvement_20261102/replay_figure_inputs.py')
    print('Generated artwork is under output/workspace/studies/laser; paper/ is unchanged.')

if __name__ == '__main__':
    main()
