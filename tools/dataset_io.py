"""Read and validate the public geometry/acquisition contract without a solver licence."""
import csv
import hashlib
import json
from pathlib import Path
import re

import numpy as np

IDENTIFIER = re.compile(r'^[A-Za-z][A-Za-z0-9_-]*$')
TOL = 1e-10


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'),
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError(f'Nonfinite JSON value: {x}')))


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def local_file(root, relative):
    require(isinstance(relative, str) and relative and '\\' not in relative,
            'Input paths must be relative and use forward slashes')
    path = (root / relative).resolve()
    require(not Path(relative).is_absolute() and path.is_relative_to(root) and path.is_file(),
            f'Input path is missing or outside the dataset: {relative}')
    return path


def identifier(value):
    require(isinstance(value, str) and IDENTIFIER.fullmatch(value), f'Invalid identifier: {value!r}')
    return value


def finite(value, name):
    require(isinstance(value, (int, float)) and not isinstance(value, bool) and np.isfinite(value),
            f'{name} must be a finite number')
    return float(value)


def point(value):
    require(isinstance(value, list) and len(value) == 2, 'Coordinates must be [x,z] in metres')
    return np.array([finite(v, 'coordinate') for v in value])


def cross(a, b):
    return a[0]*b[1] - a[1]*b[0]


def on_segment(p, a, b):
    return abs(cross(b-a, p-a)) <= TOL*np.linalg.norm(b-a) and np.dot(p-a, p-b) <= TOL**2


def in_loop(p, polygon):
    """Return 0 outside, 1 inside, 2 on the boundary."""
    result = False
    for a, b in zip(polygon, np.roll(polygon, -1, axis=0)):
        if on_segment(p, a, b):
            return 2
        if (a[1] > p[1]) != (b[1] > p[1]) and p[0] < (b[0]-a[0])*(p[1]-a[1])/(b[1]-a[1])+a[0]:
            result = not result
    return int(result)


def contains(p, polygons):
    return bool(in_loop(p, polygons[0])) and not any(in_loop(p, hole) == 1 for hole in polygons[1:])


def visible(a, b, polygons):
    if not contains(a, polygons) or not contains(b, polygons):
        return False
    direction = b-a
    if np.linalg.norm(direction) < TOL:
        return True
    cuts = [0., 1.]
    for polygon in polygons:
        for c, d in zip(polygon, np.roll(polygon, -1, axis=0)):
            det = cross(direction, d-c)
            if abs(det) > 1e-20:
                t, u = cross(c-a, d-c)/det, cross(c-a, direction)/det
                if 0 < t < 1 and -TOL <= u <= 1+TOL:
                    cuts.append(t)
            elif abs(cross(direction, c-a)) < TOL*np.linalg.norm(direction):
                cuts.extend(float(np.dot(p-a, direction)/np.dot(direction, direction)) for p in (c, d)
                            if 0 < np.dot(p-a, direction)/np.dot(direction, direction) < 1)
    cuts.sort()
    return all(contains(a+direction*(lo+hi)/2, polygons) for lo, hi in zip(cuts[:-1], cuts[1:]))


def validate_geometry(geometry):
    require(geometry.get('schema') == 2 and geometry.get('coordinate_system') == 'x_z_m',
            'Geometry needs schema 2 and coordinate_system x_z_m')
    material = geometry['material']
    cl, cs, cr, rho = [finite(material[k], k) for k in ('cl_m_s', 'cs_m_s', 'cr_m_s', 'rho_kg_m3')]
    require(cl > cs > cr > 0 and rho > 0, 'Require cl > cs > cr > 0 and positive density')
    polygons, edges, names = [], [], set()
    for index, loop in enumerate(geometry['loops']):
        require(loop['role'] == ('outer' if index == 0 else 'hole'), 'First loop must be outer; subsequent loops must be holes')
        vertices = np.array([point(p) for p in loop['vertices_m']])
        ids = loop['boundary_ids']
        require(len(vertices) >= 3 and len(ids) == len(vertices), 'One boundary id per edge; omit repeated closing vertex')
        area = sum(cross(a, b) for a, b in zip(vertices, np.roll(vertices, -1, axis=0)))
        require(area > 1e-20 if index == 0 else area < -1e-20,
                'Outer loop needs positive signed area; holes negative signed area in x,z')
        for edge_index, (name, a, b) in enumerate(zip(ids, vertices, np.roll(vertices, -1, axis=0))):
            identifier(name)
            require(name not in names, f'Duplicate boundary id: {name}')
            require(np.linalg.norm(b-a) > TOL, 'Zero-length polygon edge')
            names.add(name)
            edges.append((index, edge_index, name, a, b))
        polygons.append(vertices)
    require(polygons, 'At least one polygon is required')
    for i, (loop, ei, _, a, b) in enumerate(edges):
        for other, ej, _, c, d in edges[i+1:]:
            adjacent = loop == other and ((ei-ej) % len(polygons[loop]) in (1, len(polygons[loop])-1))
            u, v, w = b-a, d-c, c-a
            det = cross(u, v)
            if abs(det) > 1e-20:
                t, s = cross(w, v)/det, cross(w, u)/det
                require(adjacent or not (-TOL <= t <= 1+TOL and -TOL <= s <= 1+TOL), 'Polygon boundaries cross or touch')
            elif abs(cross(w, u)) < 1e-20:
                t0, t1 = np.dot(w, u)/np.dot(u, u), np.dot(d-a, u)/np.dot(u, u)
                overlap = min(1., max(t0, t1))-max(0., min(t0, t1))
                require(overlap < TOL and (adjacent or overlap < -TOL), 'Polygon boundaries overlap or touch')
    for i, hole in enumerate(polygons[1:], 1):
        require(all(in_loop(p, polygons[0]) == 1 for p in hole), 'Hole is outside or touching the outer boundary')
        require(not any(in_loop(hole[0], other) for j, other in enumerate(polygons[1:], 1) if i != j), 'Nested holes are unsupported')
    path_ids = {'1-1_LL', '1-1_LS', '1-1_SL', '1-1_SS'}
    for path in geometry.get('paths', []):
        name = identifier(path['id'])
        require(name not in path_ids, f'Duplicate path id: {name}')
        path_ids.add(name)
        for branch in ('source', 'receiver'):
            modes, faces = path[branch]['modes'], path[branch]['boundaries']
            require(re.fullmatch('[LS]{1,5}', modes) and len(modes) == len(faces)+1,
                    'Each branch needs one L/S mode per leg, and one boundary per reflection (up to four)')
            require(all(f in names for f in faces), 'Unknown path boundary')
            require(all(a != b for a, b in zip(faces[:-1], faces[1:])), 'Adjacent reflections cannot use the same boundary')
    require(all(name in names for name in geometry.get('screening_boundaries', {}).values()), 'Unknown screening boundary')
    return polygons, edges


def load_dataset(folder, check_histories=True):
    """Return the validated public definitions and a compact integrity report."""
    root = Path(folder).resolve()
    manifest = read_json(root / 'dataset.json')
    require(manifest.get('kind') == 'SAFT acquisition dataset', 'Not a SAFT acquisition dataset')
    if manifest.get('schema') == 1:
        if check_histories:
            for entry in manifest['audit']['inputs']:
                path = local_file(root, entry['path'])
                require(sha(path) == entry['sha256'], f'Legacy input checksum mismatch: {entry["path"]}')
        report = dict(status='passed', schema=1, name=manifest['name'],
                      sources=len(manifest['audit']['sources']), receivers=len(manifest['audit']['receivers']),
                      recordings=len(manifest['audit']['inputs']), samples=5001, dt_s=1e-8)
        return manifest, None, None, report
    require(manifest.get('schema') == 2, 'Unsupported dataset schema')
    require(isinstance(manifest.get('name'), str) and manifest['name'], 'Dataset name is required')
    geometry_path, acquisition_path = (local_file(root, manifest[k]) for k in ('geometry', 'acquisition'))
    geometry, acquisition = read_json(geometry_path), read_json(acquisition_path)
    polygons, edges = validate_geometry(geometry)
    require(acquisition.get('schema') == 2 and acquisition.get('velocity_units') == 'm/s',
            'Acquisition needs schema 2 and velocity_units m/s')
    sampling = acquisition['sampling']
    samples, dt, start = [finite(sampling[k], k) for k in ('samples', 'dt_s', 'start_s')]
    require(samples == int(samples) and 3 <= samples <= 100000000 and dt > 0, 'Invalid uniform time grid')
    samples = int(samples)
    half = finite(acquisition['filter_half_samples'], 'filter_half_samples')
    require(half == int(half) and half >= 1 and 2*half+2 <= samples and half <= 1000000, 'Filter support must fit in the recording')
    response = acquisition['response']
    require(response['model'] in ('isotropic', 'laser_free_surface'), 'Unknown response model')
    laser = response['model'] == 'laser_free_surface'
    if laser:
        require(finite(response['line_half_width_m'], 'line_half_width_m') > 0 and
                finite(response['pulse_fwhm_s'], 'pulse_fwhm_s') > 0, 'Laser width and pulse duration must be positive')
    source_ids, receiver_ids = [], []
    for key, ids in (('sources', source_ids), ('receivers', receiver_ids)):
        for entry in acquisition[key]:
            name = identifier(entry['id'])
            require(name not in ids, f'Duplicate {key} id: {name}')
            ids.append(name)
            p = point(entry['position_m'])
            require(contains(p, polygons), f'{name} is outside the specimen')
            require(not laser or 'inward_normal' in entry, 'Laser response requires endpoint inward normals')
            n = point(entry.get('inward_normal', [0., 1.]))
            require(abs(np.linalg.norm(n)-1) < 1e-8, 'Endpoint inward_normal must be a unit vector')
            if laser:
                matches = []
                for _, _, _, a, b in edges:
                    face_normal = np.array([-(b-a)[1], (b-a)[0]])/np.linalg.norm(b-a)
                    matches.append(on_segment(p, a, b) and np.dot(p-a, b-a) > 1e-16 and
                                   np.dot(p-b, a-b) > 1e-16 and np.linalg.norm(n-face_normal) < 1e-8)
                require(any(matches), f'{name} requires a face-interior position and its inward normal')
            if key == 'sources':
                finite(entry['trigger_s'], 'trigger_s')
            else:
                identifier(entry['group'])
        require(ids, f'At least one {key} entry is required')
    cases, components = acquisition['cases'], acquisition['components']
    require(cases and len(cases) == len(set(cases)) and 'difference' not in cases, 'Cases must be nonempty, unique and exclude reserved difference')
    for case in cases:
        identifier(case)
    require(components and len(components) == len(set(components)) and set(components) <= {'vx', 'vz'}, 'Components must be unique vx and/or vz')
    if 'difference' in acquisition:
        pair = acquisition['difference']
        require(pair['minuend'] in cases and pair['subtrahend'] in cases and pair['minuend'] != pair['subtrahend'], 'Difference needs two distinct available cases')
    expected = {(s, c, v) for s in source_ids for c in cases for v in components}
    seen, paths, hashes = set(), set(), {manifest['geometry']: sha(geometry_path), manifest['acquisition']: sha(acquisition_path)}
    times = start + np.arange(samples)*dt
    require(np.isfinite(times).all() and (np.diff(times) > 0).all(), 'Time grid must be finite and strictly increasing at float64 precision')
    for recording in acquisition['recordings']:
        key = tuple(recording[k] for k in ('source', 'case', 'component'))
        require(key in expected and key not in seen, f'Unknown or duplicate recording: {key}')
        seen.add(key)
        path = local_file(root, recording['path'])
        require(path not in paths and path not in (geometry_path, acquisition_path, root/'dataset.json'), 'Recording paths must be unique')
        paths.add(path)
        if check_histories:
            with path.open(encoding='utf-8-sig', newline='') as stream:
                rows = csv.reader(stream)
                require(next(rows) == ['time_s', *receiver_ids], f'CSV receiver header/order mismatch: {path}')
                row_count = 0
                for row in rows:
                    require(len(row) == len(receiver_ids)+1, f'Blank row or incorrect CSV column count: {path}')
                    row_count += 1
                require(row_count == samples, f'CSV sample count mismatch: {path}')
            values = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2, encoding='utf-8-sig')
            require(values.shape == (samples, len(receiver_ids)+1), f'CSV shape mismatch: {path}')
            require(np.isfinite(values).all(), f'Nonfinite velocity history: {path}')
            require(np.all(np.abs(values[:, 0]-times) <= max(1e-15, dt*1e-7)), f'CSV time grid mismatch: {path}')
            require((np.diff(values[:, 0]) > 0).all(), f'Duplicate or decreasing CSV times: {path}')
            hashes[recording['path']] = sha(path)
    require(seen == expected, f'Missing recordings: {sorted(expected-seen)}')
    report = dict(status='passed', schema=2, name=manifest['name'], sources=len(source_ids),
                  receivers=len(receiver_ids), recordings=len(seen), samples=samples, dt_s=dt,
                  start_s=start, histories_checked=check_histories, input_sha256=hashes)
    return manifest, geometry, acquisition, report
