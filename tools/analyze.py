"""Validate, reconstruct, plot and summarize any supported SAFT dataset."""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from dataset_io import in_loop, load_dataset, read_json

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'studies/laser/src'))
from standalone import executable, native_env


def metrics(amplitude, points, regions):
    finite = np.isfinite(amplitude)
    result = dict(pixels=len(points), supported_pixels=int(finite.sum()),
                  support_pct=100*float(finite.mean()), peak_m_s=None, peak_position_m=None)
    if finite.any():
        index = int(np.nanargmax(amplitude))
        result.update(peak_m_s=float(amplitude[index]), peak_position_m=points[index].tolist())
    for name, polygon in regions.items():
        polygon = np.asarray(polygon, dtype=float)
        if polygon.ndim != 2 or polygon.shape[1] != 2 or len(polygon) < 3 or not np.isfinite(polygon).all():
            raise ValueError(f'Region {name} must be a finite polygon in x,z metres')
        mask = np.array([bool(in_loop(p, polygon)) for p in points])
        supported = mask & finite
        result[name] = dict(pixels=int(mask.sum()), supported_pixels=int(supported.sum()),
                            peak_m_s=float(amplitude[supported].max()) if supported.any() else None)
    if 'target' in result and 'background' in result:
        target, background = (result[k]['peak_m_s'] for k in ('target', 'background'))
        result['target_background_db'] = float(20*np.log10(target/background)) if target and background else None
    return result


def analyze(dataset, request_path, output):
    dataset, request_path, output = map(lambda p: Path(p).resolve(), (dataset, request_path, output))
    manifest, geometry, acquisition, validation = load_dataset(dataset)
    request = read_json(request_path)
    if not isinstance(request.get('requests'), list) or not request['requests']:
        raise ValueError('Analysis request needs a nonempty requests array')
    for job in request['requests']:
        if 'signal_processing' in job:
            raise ValueError('External processed banks use SAFTFilterStudy; analyze.py uses the native analytic FIR')
    if output.exists() and any(output.iterdir()):
        raise ValueError('Choose an empty analysis output directory')
    program = executable('SAFTBatch')
    if not program.is_file():
        raise FileNotFoundError('Build the native tools first; see README.md')
    output.mkdir(parents=True, exist_ok=True)
    request['verify_inputs'] = True
    saved_request = output / 'request.json'
    saved_request.write_text(json.dumps(request, indent=2)+'\n')
    (output / 'validation.json').write_text(json.dumps(validation, indent=2)+'\n')
    command = [str(program), '--repository', str(ROOT), '--dataset', str(dataset),
               '--compute', str(saved_request), '--output', str(output / 'native')]
    subprocess.run(command, env=native_env(), check=True)
    if geometry is not None:
        polygons = [np.array(loop['vertices_m'])*1000 for loop in geometry['loops']]
        sources = np.array([s['position_m'] for s in acquisition['sources']])*1000
        receivers = np.array([r['position_m'] for r in acquisition['receivers']])*1000
    else:
        h = manifest['material']['thickness_m']
        polygons = [np.column_stack([np.array(p)[:, 0], h-np.array(p)[:, 1]])*1000
                    for p in manifest['geometry'].values()]
        sources = np.array([[s['xy_m'][0], h-s['xy_m'][1]] for s in manifest['audit']['sources']])*1000
        receivers = np.array([[r['x_m'], h-r['y_m']] for r in manifest['audit']['receivers']])*1000
    rows = []
    for i, job in enumerate(request['requests']):
        folder = output / 'native' / str(i)
        selection = read_json(folder / 'selection.json')['selection']
        with np.load(folder / 'images.npz') as arrays:
            points = arrays['points_m']
            for case in selection['cases']:
                amplitude = arrays[case+'_maps'][2]
                row = dict(request=job.get('id', str(i)), case=case, path=selection['path'],
                           **metrics(amplitude, points, request.get('regions', {})))
                rows.append(row)
                fig, ax = plt.subplots(figsize=(8, 5.4), constrained_layout=True)
                x, z = arrays['x_mm'], arrays['depth_mm']
                if len(points) == len(x)*len(z) and len(x) > 1 and len(z) > 1:
                    artist = ax.pcolormesh(x, z, amplitude.reshape(len(z), len(x))*1e6, shading='nearest', cmap='magma')
                else:
                    artist = ax.scatter(points[:, 0]*1000, points[:, 1]*1000, c=amplitude*1e6, cmap='magma')
                for polygon in polygons:
                    closed = np.vstack([polygon, polygon[0]])
                    ax.plot(closed[:, 0], closed[:, 1], color='0.4', lw=.9)
                ax.scatter(sources[:, 0], sources[:, 1], marker='*', s=35, c='#126e82', label='Sources')
                ax.scatter(receivers[:, 0], receivers[:, 1], marker='v', s=12, c='#22a699', label='Receivers')
                ax.set(xlabel='x (mm)', ylabel='z / depth (mm)', aspect='equal',
                       title=f'{manifest["name"]}\n{row["request"]}: {case} / {selection["path"]}')
                ax.invert_yaxis()
                ax.legend(loc='best', fontsize=8)
                fig.colorbar(artist, ax=ax, label='Reconstruction amplitude (µm/s)')
                fig.savefig(output / f'{i}_{case}.png', dpi=160)
                plt.close(fig)
    report = dict(status='passed', dataset=str(dataset), dataset_name=manifest['name'],
                  created_utc=datetime.now(timezone.utc).isoformat(), results=rows,
                  interpretation='SAFT amplitudes and optional user-defined region contrasts; no automatic defect or experimental SNR claim.')
    (output / 'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    fields = ['request', 'case', 'path', 'supported_pixels', 'pixels', 'support_pct', 'peak_m_s', 'target_background_db']
    with (output / 'summary.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', type=Path, required=True)
    parser.add_argument('--request', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    try:
        report = analyze(args.dataset, args.request, args.output)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f'Analysis failed: {error}\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
