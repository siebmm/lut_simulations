"""Generate a small, deterministic non-weld dataset with point-scatterer echoes."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
from dataset_io import read_json, validate_geometry, visible

ROOT = Path(__file__).resolve().parents[1]


def generate(output, source_count=7, receiver_count=13, samples=2401, dt=20e-9, start=-2e-6):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError('Choose an empty output directory to preserve existing data')
    example = ROOT / 'saft/examples/general_2d'
    for name in ('dataset.json', 'geometry.json', 'analysis.json'):
        shutil.copy2(example / name, output / name)
    geometry = read_json(output / 'geometry.json')
    polygons, _ = validate_geometry(geometry)
    sources = [dict(id=f'shot{i}', position_m=[float(x), 0.], trigger_s=i*30e-9, inward_normal=[0., 1.])
               for i, x in enumerate(np.linspace(-.014, .016, source_count))]
    receivers = [dict(id=f'sensor{i}', position_m=[float(x), 0.], group='array', inward_normal=[0., 1.])
                 for i, x in enumerate(np.linspace(-.017, .020, receiver_count))]
    acquisition = dict(schema=2, velocity_units='m/s',
                       sampling=dict(samples=samples, dt_s=dt, start_s=start), filter_half_samples=128,
                       response=dict(model='isotropic'), sources=sources, receivers=receivers,
                       cases=['measurement', 'baseline'], components=['vx', 'vz'],
                       difference=dict(minuend='measurement', subtrahend='baseline'), recordings=[])
    times = start + np.arange(samples)*dt
    target = np.array([.003, .012])
    speed = geometry['material']['cl_m_s']
    def packet(arrival):
        tau = times-arrival
        return np.exp(-.5*(tau/.18e-6)**2)*np.cos(2*np.pi*2e6*tau)
    for source in sources:
        a = np.array(source['position_m'])
        echo = np.zeros((samples, receiver_count))
        baseline = np.zeros_like(echo)
        for j, receiver in enumerate(receivers):
            b = np.array(receiver['position_m'])
            baseline[:, j] = .05e-3*packet(source['trigger_s']+np.linalg.norm(a-b)/speed)
            if visible(a, target, polygons) and visible(target, b, polygons):
                arrival = source['trigger_s']+(np.linalg.norm(a-target)+np.linalg.norm(target-b))/speed
                echo[:, j] = 1e-3*packet(arrival)
        for case, values in (('measurement', echo+baseline), ('baseline', baseline)):
            for component, scale in (('vx', .5), ('vz', 1.)):
                relative = f'histories/{source["id"]}_{case}_{component}.csv'
                path = output / relative
                path.parent.mkdir(exist_ok=True)
                np.savetxt(path, np.column_stack([times, scale*values]), delimiter=',', fmt='%.17g',
                           header='time_s,'+','.join(r['id'] for r in receivers), comments='')
                acquisition['recordings'].append(dict(source=source['id'], case=case, component=component, path=relative))
    (output / 'acquisition.json').write_text(json.dumps(acquisition, indent=2)+'\n', encoding='utf-8')
    return output


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=ROOT / 'output/example_dataset')
    args = parser.parse_args()
    print(generate(args.output).resolve())


if __name__ == '__main__':
    main()
