# Surface velocity histories

These files contain the receiver signals from the full-weld and partial-weld
simulations. Each laser position has a **sound reference** and a **lack-of-fusion
case**. Their difference helps identify the part of the wavefield caused by
the defect.

| Quantity | Full weld | Partial weld |
|---|---:|---:|
| Laser positions | 5 | 5 |
| Receivers | 337 | 331 |
| Recorded components | vx and vy | vx and vy |
| Time interval | 0-50 microseconds | 0-50 microseconds |
| Sample spacing | 10 ns | 10 ns |
| Samples per trace | 5,001 | 5,001 |

## Read a recording

Within either specimen folder, the path identifies the physical measurement:

```text
sources/<laser position>/<reference or lof>/receivers/<left, right or weld>_<vx or vy>.csv
```

`vx` is horizontal velocity and `vy` is upward vertical velocity, both in m/s.
The first CSV column, `time_s`, is time in seconds. The remaining columns are
receiver traces. `dataset.json` gives their coordinates, ordering, source
positions and material properties. The receiver labels restart within each
group; use the associated coordinates when combining groups.

The files contain the resampled histories used in the paper. Their metadata
also records any repairs to the original exports. Original recording bytes
are preserved and can be checked with `python tools/verify_repository.py`.

## Use the signals for imaging

Follow the [SAFT guide](../saft/README.md) to reconstruct the specimen from
these traces, or run `python tools/reproduce_main.py` to repeat the main images.
The recorded interval includes propagation and scattering through the weld;
the short local laser-generation comparison is described
[separately](../models/laser_generation/README.md).

For another specimen, supply its own geometry, source/receiver coordinates and
velocity histories using the [data-format guide](../saft/FORMAT.md). Its counts
and sampling do not need to match these examples. If depth points downward,
convert the vertical velocity consistently: `vz = -vy`.

The CSV files total about 1.57 GB and are stored with Git LFS. Follow the
[setup instructions](../tools/README.md#software-setup) to download their contents.
