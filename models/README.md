# The physical models

The simulation separates **wave generation** from **wave propagation**. A small
model resolves the heated layer under the laser. Its calibrated equivalent
force then excites the larger weld models without resolving heat conduction
throughout the weld calculation.

## Local laser generation

Open [laser_generation/laser_source_calibration.mph](laser_generation/laser_source_calibration.mph)
in COMSOL. This is the original model copied from the NDT study, with its
geometry, mesh and two source descriptions intact. Follow the
[laser-generation guide](laser_generation/README.md) for the governing physics,
parameters, solve procedure and comparison with the recorded calibration.

## Full weld

`full_weld/` contains a sound model (`reference.mph`) and a lack-of-fusion model
(`lof.mph`) for each of five laser positions:

| Folder | Laser position x |
|---|---:|
| `shoulder_left` | -12 mm |
| `bead_left` | -1.5 mm |
| `bead_centre` | 0 mm |
| `bead_right` | 1.5 mm |
| `shoulder_right` | 12 mm |

These ten files already contain the source position, geometry, mesh and solver
settings for their respective shots. Open the desired file, save a working
copy under `output/`, and compute its time-dependent study (`sol1`). Set export
destinations to your working folder before exporting results.

## Partial weld

`partial_weld/reference.mph` and `partial_weld/lof.mph` describe the sound and
defective J-bevel weld with **9 mm fill**. The same pair is reused for five
source positions; this is why there are two partial-weld files and ten
full-weld files.

The positions and source boundaries are listed in
`partial_weld/metadata/multisource.json`. For another shot, set `laser_x`,
`laser_y` and the `sourceFace` boundary selection accordingly in a working copy.
The other physical parameters remain those of the chosen specimen.

## Geometry, coordinates and recordings

Each weld's `metadata/` directory contains the model parameters, source
positions and receiver coordinates. `defect_geometry.json` defines the shared
lack-of-fusion flaw. The COMSOL models use x and upward y; the paper's image
depth is plate thickness minus y. Apply the same sign change to vertical
velocity when converting coordinates.

The [supplied receiver histories](../data/README.md) are the results used for
the paper's images. To image a different specimen, export its geometry and
recordings using the [SAFT data format](../saft/FORMAT.md).

All thirteen supplied `.mph` files are original model inputs, with checksums
in `checksums.json`. New calculations need COMSOL 6.4 and licences for the
physics interfaces in the model. Solved fields and new exports belong in a
separate working directory.
