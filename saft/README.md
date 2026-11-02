# C++ SAFT tool

`src/` contains the reconstruction engine, weld geometry and ray paths,
wall response, and packet screening. `cli/` contains the two existing
command-line entry points:

| Executable | Purpose |
|---|---|
| `SAFTFilterStudy` | Reconstruction grids or selected points, optional external processed-signal banks, and raw-history export |
| `SAFTReview` | Reconstruction grids with per-path statistics and the groove-path catalogue |

Build from the repository root using the commands in the top-level README.
Qt 6 Core must be available at runtime as well as when compiling.

## Inputs

Both programs take four positional arguments:

```text
SAFTReview ROOT DATASET REQUEST OUTPUT
SAFTFilterStudy ROOT DATASET REQUEST OUTPUT
```

- `ROOT`: a writable working directory; caches go to `ROOT/build/saft_cache/`.
- `DATASET`: an external SAFT acquisition folder with `dataset.json` and its
  receiver CSV files. Supply the full-weld or partial-weld dataset from the full
  AllseasNDT study. Dataset paths are resolved from the current working directory.
- `REQUEST`: a JSON object with a `requests` array. Each job has an `id` and a
  `selection` object parsed by `Config::fromJson` in `src/saft.cpp`.
- `OUTPUT`: the directory for generated arrays and JSON provenance.

The acquisition manifest has schema 1 and kind `SAFT acquisition dataset`.
It supplies material, geometry, source/receiver coordinates, audit hashes, and
repair information. Receiver histories are referenced under
`sources/<source>/<reference|lof>/receivers/<left|right|weld>_<vx|vy>.csv`.
The existing engine expects five sources, at most 337 receivers, and 5,001
samples at 10 ns intervals. A COMSOL `.mph` file is not a SAFT dataset.

## Partial-weld example

`examples/partial_weld.json` retains the final 0.025 mm muted reconstruction
and its matched unmuted control, including their original selections. With
the external partial-weld acquisition available, run from the repository root:

```powershell
.\build\bin\SAFTReview.exe . "C:\path\to\partial_weld" saft/examples/partial_weld.json output/partial
```

For multi-configuration builds, use `build/bin/Release/SAFTReview.exe`.
On Linux/macOS, use `./build/bin/SAFTReview`. These fine grids require several
GB of working memory/cache space and can take time to compute.

Outputs include grid coordinates, binary float64/complex128 arrays, request
settings, timing, and provenance. The analysis scripts that combine these
arrays into the paper's metrics and figures remain in the full study. Its
full-weld result also uses a custom processed-signal bank; the default engine
filter alone does not reproduce that result.
