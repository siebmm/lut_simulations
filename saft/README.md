# C++ SAFT tool

`src/` contains the reconstruction engine, weld geometry and ray paths,
wall response, and packet screening. `cli/` contains four
command-line entry points:

| Executable | Purpose |
|---|---|
| `SAFTFilterStudy` | Reconstruction grids or selected points, optional external processed-signal banks, and raw-history export |
| `SAFTReview` | Reconstruction grids with per-path statistics and the groove-path catalogue |
| `SAFTAdmissionAudit` | Per-source admission statistics for the appendix's threshold checks |
| `SAFTBatch` | Headless version of the original desktop batch exporter, including geometry checks and `images.npz` |

Build from the repository root using the commands in the top-level README.
Qt 6 Core must be available at runtime as well as when compiling.

For a new model, use [the schema-2 geometry and acquisition format](FORMAT.md).
All four executables share this reader and support variable source/receiver
counts, timing, arbitrary case names, and finite polygon boundaries with voids.
`tools/analyze.py` provides validation, reconstruction, plots and metrics.

## Inputs

The filter, review and admission programs take four positional arguments:

```text
SAFTReview ROOT DATASET REQUEST OUTPUT
SAFTFilterStudy ROOT DATASET REQUEST OUTPUT
```

- `ROOT`: a writable working directory; caches go to `ROOT/build/saft_cache/`.
- `DATASET`: `data/full_weld` or `data/partial_weld`, each containing
  `dataset.json` and its receiver histories. Paths are resolved from the current
  working directory.
- `REQUEST`: a JSON object with a `requests` array. Each job has an `id` and a
  `selection` object parsed by `Config::fromJson` in `src/saft.cpp`.
- `OUTPUT`: the directory for generated arrays and JSON provenance.

The archived weld acquisition manifest has schema 1 and kind `SAFT acquisition dataset`.
It supplies material, geometry, source/receiver coordinates, audit hashes, and
repair information. Receiver histories are referenced under
`sources/<source>/<reference|lof>/receivers/<left|right|weld>_<vx|vy>.csv`.
The legacy schema-1 adapter expects five sources, at most 337 receivers, and
5,001 samples at 10 ns intervals. Schema 2 reads these dimensions and its
recording paths from `acquisition.json`. A COMSOL `.mph` file is not a SAFT dataset.

## Partial-weld example

`examples/partial_weld.json` retains the final 0.025 mm muted reconstruction
and its matched unmuted control, including their original selections. Run from
the repository root:

```powershell
.\build\bin\SAFTReview.exe . data/partial_weld saft/examples/partial_weld.json output/partial
```

For multi-configuration builds, use `build/bin/Release/SAFTReview.exe`.
On Linux/macOS, use `./build/bin/SAFTReview`. These fine grids require several
GB of working memory/cache space and can take time to compute.

Outputs include grid coordinates, binary float64/complex128 arrays, request
settings, timing, and provenance. `python tools/reproduce_main.py` regenerates
the full-weld processed bank, runs the published grids, and verifies metrics
against the included archive. `python tools/reproduce_figures.py` uses the
retained figure generators and reconstructs the additional review grids needed.

The older geometry studies use `SAFTBatch --dataset DATASET --compute REQUEST
--output OUTPUT`, with optional `--repository ROOT`. It exports the same
numbered folders and arrays as the original desktop application's compute mode.
All four tools share the reconstruction and physics source in `src/`.
