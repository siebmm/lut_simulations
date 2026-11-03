# Core COMSOL models

These are the original model inputs for the completed full weld and the
partial J-bevel weld with 9 mm fill. All twelve `.mph` files retain their
original bytes; `checksums.json` records their SHA-256 hashes with paths relative
to this directory. They require COMSOL Multiphysics 6.4 and the relevant physics
licences.

## Full weld

`full_weld/<source>/reference.mph` and `full_weld/<source>/lof.mph` contain the
configured pristine and lack-of-fusion cases for each of five laser positions:

| Source directory | x position (mm) |
|---|---:|
| `shoulder_left` | -12 |
| `bead_left` | -1.5 |
| `bead_centre` | 0 |
| `bead_right` | 1.5 |
| `shoulder_right` | 12 |

Each file already contains the source position, geometry, mesh, and solver
settings for that shot. Open the desired file and run its time-dependent
study (`sol1`) with output paths set for your installation.

## Partial weld

`partial_weld/reference.mph` and `partial_weld/lof.mph` are the original
reference/defect model pair. The five source definitions are in
`partial_weld/metadata/multisource.json`. To run another shot, configure
`laser_x`, `laser_y`, and the `sourceFace` boundary selection for that source;
the source and receiver definitions needed for this are retained locally.

## Metadata

Each weld's `metadata/` directory retains the model setup, source definitions,
and receiver coordinates. The shared flaw definition is stored once in
`defect_geometry.json`. Coordinates use the original model's
x/y convention; the paper and SAFT use depth = plate thickness - y. Historical
paths inside the original setup metadata describe the source study and are not
dependencies of this repository.

The canonical receiver histories and repair information are included in
`../data/`. The saved `.mph` inputs are the authoritative models. Legacy solver
automation and export templates are outside this replay-focused repository.
Export histories for other models using the contract in `../saft/FORMAT.md`.
No new COMSOL solve was performed during packaging.

Solved field snapshots, original native-time probe tables and machine-specific
batch-control infrastructure are omitted. They are not required to reconstruct
the published images from the included acquisitions.
