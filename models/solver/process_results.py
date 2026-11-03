"""Validate native-step probes, resample to 10 ns, and render stored field grids."""
import csv
import gzip
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.path import Path as PolygonPath

TIMES_US = (0.5, 1, 2, 4, 8, 12, 20, 30, 40)
GROUPS = ('left', 'right', 'weld')


def canonical_raw(folder, group, names):
    target = folder / (group + '_raw.csv')
    if target.exists():
        return target
    logs = '\n'.join(p.read_text(errors='replace') for p in (folder / 'comsol.log', folder / 'console.log') if p.exists())
    line = next((line for line in logs.splitlines() if line.startswith('PROBE_HEADERS ' + group + '|')), None)
    if line is None:
        raise ValueError('Missing COMSOL probe headers for ' + group)
    headers = line.split('|')[1:]
    order = [0]
    for name in names:
        matches = [i for i, h in enumerate(headers) if name in h]
        if len(matches) != 1 or matches[0] == 0:
            raise ValueError('Missing or ambiguous probe header: ' + name)
        order.append(matches[0])
    raw = np.loadtxt(folder / (group + '_native.txt'), comments='%', ndmin=2)
    if raw.shape[1] != len(headers) or len(headers) != len(names) + 1:
        raise ValueError('Incomplete native probe table: ' + group)
    temp = target.with_suffix('.tmp')
    np.savetxt(temp, raw[:, order], delimiter=',', header='time_s,' + ','.join(names), comments='', fmt='%.17g')
    temp.replace(target)
    return target


def load_raw(path, expected):
    with path.open() as f:
        header = next(csv.reader(f))
    if header != ['time_s'] + expected:
        raise ValueError('Unexpected columns: ' + str(path))
    a = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2)
    if a.shape[1] != len(header) or len(a) < 2 or not np.isfinite(a).all():
        raise ValueError('Incomplete or nonfinite probe export: ' + str(path))
    dt = np.diff(a[:, 0])
    # Accept repeated output at exactly the same time only if the values agree.
    same = dt == 0
    if (dt < 0).any() or not np.allclose(a[1:][same], a[:-1][same], rtol=1e-10, atol=1e-15):
        raise ValueError('Nonmonotonic or inconsistent probe times: ' + str(path))
    a = a[np.r_[True, ~same]]
    if abs(a[0, 0]) > 1e-14 or a[-1, 0] < 40e-6 - 1e-13 or np.max(np.diff(a[:, 0])) > 10e-9 * 1.001:
        raise ValueError('Probe record is truncated or too coarsely sampled: ' + str(path))
    return a


def comsol_export(path, times, xy, values, expressions, unit):
    """Same receiver-row layout consumed by verify_weld_laser.export(path, nvar)."""
    with path.open('w') as f:
        f.write('% Portable laser: linearly resampled native-step probes\n% Length unit: m\n')
        f.write('% x y ' + ' '.join(f'{e} ({unit}) @ t={t:.17g}' for t in times for e in expressions) + '\n')
        np.savetxt(f, np.column_stack([xy, values.reshape(len(xy), -1)]), fmt='%.17g')


def process_case(folder, root):
    folder, root = Path(folder), Path(root)
    if not (folder / 'EXPORT_COMPLETE').is_file():
        raise ValueError('COMSOL has not finished all exports: ' + str(folder))
    with (root / 'metadata/receivers.csv').open() as f:
        receivers = list(csv.DictReader(f))
    times = np.arange(4001, dtype=float) * 1e-8
    processed = folder / 'receivers'
    processed.mkdir(exist_ok=True)
    report = {'sampling': 'native solver steps; linear interpolation onto 10 ns grid',
              'time_samples': len(times), 'units': 'm, s, m/s', 'groups': {}}
    for group in GROUPS:
        rows = [r for r in receivers if r['group'] == group]
        ids = [r['receiver_id'] for r in rows]
        xy = np.array([[float(r['x_m']), float(r['y_m'])] for r in rows])
        names = [r + '_' + e for r in ids for e in ('vx', 'vy')]
        raw = load_raw(canonical_raw(folder, group, names), names)
        v = np.column_stack([np.interp(times, raw[:, 0], raw[:, k]) for k in range(1, raw.shape[1])])
        v = v.reshape(len(times), len(rows), 2).transpose(1, 0, 2)
        if not np.any(v):
            raise ValueError('All-zero receiver group: ' + group)
        values = np.concatenate([v, v[:, :, 1:2]], axis=2)  # vertical LDV assumption
        np.savez_compressed(processed / (group + '.npz'), time_s=times, xy_m=xy,
                            receiver_id=ids, vx_m_s=v[:, :, 0], vy_m_s=v[:, :, 1], v_ldv_m_s=v[:, :, 1])
        comsol_export(processed / (group + '.txt'), times, xy, values, ('vx', 'vy', 'v_ldv'), 'm/s')
        for k, expr in enumerate(('vx', 'vy', 'v_ldv')):
            np.savetxt(processed / f'{group}_{expr}.csv', np.column_stack([times, values[:, :, k].T]),
                       delimiter=',', header='time_s,' + ','.join(ids), comments='', fmt='%.17g')
        report['groups'][group] = {'receivers': len(rows), 'native_samples': len(raw),
                                   'max_native_step_s': float(np.diff(raw[:, 0]).max()),
                                   'peak_abs_vy_m_s': float(abs(v[:, :, 1]).max())}
    names = ['Sminus_laser_tx', 'Splus_laser_tx']
    source = load_raw(canonical_raw(folder, 'source', names), names)
    if not np.allclose(source[:, 1], -source[:, 2], rtol=1e-9, atol=1e-9) or source[-1, 2] <= 0:
        raise ValueError('Source traction dipole check failed')
    if np.max(abs(source[:, 1:3]), axis=0).min() == 0:
        raise ValueError('Source is inactive')
    (folder / 'validation.json').write_text(json.dumps(report, indent=2) + '\n')
    for time in TIMES_US:
        native = folder / f'field_{time:04.1f}us.txt'
        packed = folder / f'field_{time:04.1f}us.csv.gz'
        if native.exists():
            a = np.loadtxt(native, comments='%')
            temp = packed.with_suffix('.tmp')
            with gzip.open(temp, 'wt') as f:
                np.savetxt(f, a, delimiter=',', header='x_m,y_m,vx_m_s,vy_m_s', comments='', fmt='%.17g')
            temp.replace(packed)
            native.unlink()
    # Native-step CSVs retain every sample at full precision. Remove only redundant
    # copies created by this runner, after successful validation; failed runs keep them.
    for group in (*GROUPS, 'source'):
        for suffix in ('_stream.txt', '_native.txt'):
            (folder / (group + suffix)).unlink(missing_ok=True)
    return report


def read_field(path, polygons):
    a = np.loadtxt(path, delimiter=',', skiprows=1)
    if a.shape != (1501 * 255, 4):
        raise ValueError('Incomplete snapshot: ' + str(path))
    xx, yy = np.meshgrid(np.linspace(-.075, .075, 1501), np.linspace(0, .0254, 255))
    if not np.allclose(a[:, :2], np.column_stack([xx.ravel(), yy.ravel()]), rtol=0, atol=1e-13):
        raise ValueError('Unexpected snapshot grid: ' + str(path))
    inside = np.zeros(len(a), dtype=bool)
    interior = np.zeros(len(a), dtype=bool)
    for polygon in polygons.values():
        p = np.asarray(polygon)
        # contains_points radius assumes a counterclockwise polygon.
        area = np.sum(p[:, 0] * np.roll(p[:, 1], -1) - p[:, 1] * np.roll(p[:, 0], -1))
        if area < 0:
            p = p[::-1]
        shape = PolygonPath(np.vstack([p, p[0]]))
        inside |= shape.contains_points(a[:, :2], radius=1e-10)
        interior |= shape.contains_points(a[:, :2], radius=-1e-8)
    if not np.isfinite(a[interior, 2:]).all():
        raise ValueError('Missing interior field values: ' + str(path))
    a[~inside, 2:] = np.nan
    vy = a[:, 3].reshape(255, 1501)
    if not np.nanmax(abs(vy)) > 0:
        raise ValueError('All-zero field: ' + str(path))
    return a[:, 2:].reshape(255, 1501, 2)


def render_snapshots(run, root):
    run, root = Path(run), Path(root)
    manifest = json.loads((root / 'metadata/setup_manifest.json').read_text())
    polygons = manifest['polygons_xy_m']
    cases = [c for c in ('reference', 'lof') if (run / c / 'EXPORT_COMPLETE').exists()]
    if not cases:
        raise ValueError('No complete exports in ' + str(run))
    scales = {}
    for time in TIMES_US:
        fields = {case: read_field(run / case / f'field_{time:04.1f}us.csv.gz', polygons) for case in cases}
        vmax = max(float(np.nanmax(abs(v[:, :, 1]))) for v in fields.values()) * 1e6
        speed_max = max(float(np.nanmax(np.hypot(v[:, :, 0], v[:, :, 1]))) for v in fields.values()) * 1e6
        scales[str(time)] = {'vy_symmetric': vmax, 'speed_maximum': speed_max}
        for case, v in fields.items():
            target = run / case / 'snapshots'
            target.mkdir(exist_ok=True)
            fig = plt.figure(figsize=(12, 7), constrained_layout=True)
            grid = fig.add_gridspec(2, 2, height_ratios=(1, 2))
            axes = [fig.add_subplot(grid[0, :]), fig.add_subplot(grid[1, 0]), fig.add_subplot(grid[1, 1])]
            for index, ax in enumerate(axes):
                speed = index == 2
                data = np.hypot(v[:, :, 0], v[:, :, 1]) if speed else v[:, :, 1]
                im = ax.imshow(data * 1e6, origin='lower', extent=(-75.05, 75.05, -.05, 25.45),
                               cmap='inferno' if speed else 'RdBu_r', vmin=0 if speed else -vmax,
                               vmax=speed_max if speed else vmax, interpolation='none', aspect='equal')
                for polygon in polygons.values():
                    p = np.asarray(polygon) * 1e3
                    p = np.vstack([p, p[0]])
                    ax.plot(p[:, 0], p[:, 1], color='0.25', linewidth=.55)
                ax.plot(0, 9, marker='*', color='#c98e00', markersize=7)
                ax.set(xlabel='x (mm)', ylabel='y above backwall (mm)', ylim=(-.5, 26.5))
                ax.set_title('Velocity magnitude |v|' if speed else ('Full field: vertical velocity' if index == 0 else 'Weld region: vertical velocity'), fontsize=10)
                fig.colorbar(im, ax=ax, label='μm/s', shrink=.85, pad=.02, panchor=False)
            axes[0].set_xlim(-75, 75)
            for x in (-60, 60):
                axes[0].axvline(x, color='0.45', linestyle='--', linewidth=.7)
            axes[1].set_xlim(-15, 15)
            axes[2].set_xlim(-15, 15)
            fig.suptitle(f'{case.upper()} | {time:g} μs | unfiltered ultrasonic field')
            fig.savefig(target / f'wavefield_{time:04.1f}us.png', dpi=200)
            plt.close(fig)
    (run / 'snapshot_scales.json').write_text(json.dumps({'scale': 'shared between available cases at each time',
                                                         'symmetric_limits_um_s': scales}, indent=2) + '\n')


def process_run(run, root):
    run = Path(run)
    for case in ('reference', 'lof'):
        if (run / case / 'EXPORT_COMPLETE').exists():
            process_case(run / case, root)
    render_snapshots(run, root)
    state_file = run / 'run.json'
    if state_file.exists():
        state = json.loads(state_file.read_text())
        for case in ('reference', 'lof'):
            if (run / case / 'validation.json').exists():
                state['cases'][case] = 'validated'
        if all(state['cases'].get(c) == 'validated' for c in state.get('requested_cases', ['reference', 'lof'])):
            state['status'] = 'complete'
            state.pop('error', None)
        state_file.write_text(json.dumps(state, indent=2) + '\n')
    print('Processed outputs:', run)
