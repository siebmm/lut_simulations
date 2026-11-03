"""Check the imported model pair and the short COMSOL source/receiver runs."""
from study_paths import STUDY
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import re
import xml.etree.ElementTree as ET
import numpy as np

ORIGINAL_HASHES = {
    'reference': '79615dc1c806db3c19629a3280df65d1c4540401f5380e29d9e65cac8ffd639e',
    'lof': '27af244a294cc423d1ebf1c8bef055de4a7e8b5a4cfbd30dd9126a4a6737ceb2',
}


def export(path, nvar=1):
    lines = path.read_text(encoding='utf-8-sig').splitlines()
    header = next(x for x in lines if x.startswith('%') and '@ t=' in x)
    t = np.array([float(x) for x in re.findall(r'@\s*t=([+\-\d.eE]+)', header)])
    data = np.atleast_2d(np.loadtxt(path, comments='%'))
    assert data.shape[1] == len(t)+2 and np.isfinite(data).all(), str(path)
    units = next(x.split(':', 1)[1].strip() for x in lines if x.startswith('% Length unit:'))
    factor = {'mm': 1e-3, 'm': 1.0}[units]
    coords = data[:, :2]*factor
    idx = np.argsort(coords[:, 0])
    times = t.reshape(-1, nvar)
    assert np.all(times == times[:, :1])
    return times[:, 0], coords[idx], data[idx, 2:].reshape(len(idx), -1, nvar)


def main():
    stats = {'scope': '0-0.2 us source/receiver smoke test; no defect echo or SAFT validation',
             'models': {}, 'receivers': {}}
    all_fields = {}
    for kind in ('reference', 'lof'):
        original = Path(STUDY / f'models/original_weld_geometry/baseline/weld_{kind}_v3.mph')
        assert hashlib.sha256(original.read_bytes()).hexdigest() == ORIGINAL_HASHES[kind]
        with ZipFile(original) as old, ZipFile(STUDY / f'models/weld_{kind}_laser.mph') as new:
            assert old.read('mesh1.mphbin') == new.read('mesh1.mphbin')
            for name in old.namelist():
                if name.startswith('geometry') and name.endswith('.mphbin'):
                    assert old.read(name) == new.read(name), name
            root = ET.fromstring(new.read('dmodel.xml'))
            params = {x.get('name'): x.get('expr') for x in root.find('ModelParam').findall('expressions')}
            assert params['t_end'] == '12[us]' and params['dt_out'] == '10[ns]'
            assert params['laser_F_abs'] == '10[J/m^2]'
            load = next(x for x in root.iter('PhysicsFeature') if x.get('tag') == 'load1')
            force = next(x.get('value') for x in load.findall('param') if x.get('param') == 'forceReferenceArea')
            assert "'laser_tx'" in force and "'pulse(t)'" not in force
            for tag in ('left', 'right', 'weld'):
                ds = next(x for x in root.iter('DatasetFeature') if x.get('tag') == tag)
                assert any(x.get('name') == 'p:snapping' and x.get('value') == 'boundary'
                           for x in ds.findall('propertyValue'))
            stats['models'][kind] = {'original_unchanged': True, 'geometry_and_mesh_unchanged': True,
                                     'full_record_s': 12e-6, 'source_force': force}
        all_fields[kind] = {}
        coords_all = []
        for group, count in [('left', 45), ('right', 45), ('weld', 23)]:
            fields = {}
            path = Path(STUDY / 'results/short_run_validation/receivers')/f'{kind}_{group}.txt'
            t, coords, all_values = export(path, 4)
            for i, component in enumerate(('vx', 'vy', 'ux', 'uy')):
                values = all_values[:, :, i]
                assert values.shape == (count, 101), (path, values.shape)
                assert np.allclose(t, np.arange(101)*2e-9, atol=1e-15, rtol=1e-10)
                assert np.unique(coords, axis=0).shape[0] == count
                expected_x = {'left': np.arange(-15, -3.99, .25),
                              'right': np.arange(4, 15.01, .25),
                              'weld': np.arange(-2.75, 2.751, .25)}[group]*1e-3
                expected_y = .01033930387414418 if group == 'weld' else .0198
                assert np.allclose(coords[:, 0], expected_x, atol=1e-12, rtol=0)
                assert np.allclose(coords[:, 1], expected_y, atol=1e-12, rtol=0)
                fields[component] = values
            coords_all.extend(coords.tolist())
            stats['receivers'][f'{kind}_{group}'] = {
                'count': count, 'time_samples': len(t),
                'peak_vertical_displacement_m': float(np.max(np.abs(fields['uy']))),
                'peak_vertical_velocity_m_per_s': float(np.max(np.abs(fields['vy']))),
            }
            all_fields[kind][group] = fields
        assert np.unique(np.array(coords_all), axis=0).shape[0] == 113
        assert stats['receivers'][f'{kind}_weld']['peak_vertical_velocity_m_per_s'] > 1e-9
    # The defect cannot yet have been reached; nearby early waveforms should agree.
    ref = all_fields['reference']['weld']['vy']
    lof = all_fields['lof']['weld']['vy']
    stats['early_reference_defect_velocity_relative_L2_difference'] = float(np.linalg.norm(lof-ref)/np.linalg.norm(ref))
    Path(STUDY / 'results/weld_wave_propagation/weld_laser_validation.json').write_text(json.dumps(stats, indent=2)+'\n')
    print(json.dumps(stats, indent=2))


if __name__ == '__main__':
    main()
