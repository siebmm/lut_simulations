# Running the scientific calculations

Choose the calculation you need:

| Aim | Command or guide |
|---|---|
| Understand and repeat laser generation | [Local-source model and workflow](../models/laser_generation/README.md) |
| Compare the recorded thermal and dipole waves | `python models/laser_generation/analyze.py --stored-fluence 10` |
| Reconstruct the paper's main weld images | `python tools/reproduce_main.py` |
| Regenerate the figures | `python tools/reproduce_figures.py` |
| Image another specimen | [Geometry and velocity-history format](../saft/FORMAT.md) |
| Check travel times and focusing independently | `python studies/laser/src/verify_saft_datasets.py` |

All commands below run from the repository root. Results go under `output/`.

## Software setup

**For waveform comparisons:** install Python 3.12 and the packages listed in
`requirements.txt`. **For SAFT imaging:** also build the reconstruction tools
as described below. **For new finite-element solves:** use COMSOL 6.4 with the
relevant physics licences. The recorded-data analyses need no COMSOL licence.

Download the large recordings and reference arrays with Git LFS:

```sh
git lfs install
git lfs pull
python -m venv .venv
```

Activate the environment using `source .venv/bin/activate` on Linux/macOS, or
`.venv\Scripts\activate.bat` in Windows Command Prompt. In PowerShell, use
`.venv\Scripts\python.exe` in place of `python` if activation is restricted.

```sh
python -m pip install -r requirements.txt
python tools/verify_repository.py
```

The last command checks that the models, recordings, paper and numerical inputs
are complete and match the supplied reference files.

### Build the imaging tools

Install CMake 3.21 or newer, a C++17 compiler and matching Qt 6 Core development
files, then run:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="<Qt installation>"
cmake --build build --config Release
```

The executables appear in `build/bin/` or `build/bin/Release/`. Make the matching
Qt and compiler runtime libraries available through `PATH`.

Alternatively, use an installed Qt qmake and its matching compiler:

```sh
python tools/build_native.py --qmake <path-to-qmake> --make <path-to-make>
```

On Windows, use Qt's matching MinGW `mingw32-make.exe`; this helper also copies
the required runtime libraries into `build/bin/`. `SAFT_BIN` selects another
executable directory and `SAFT_RUNTIME` adds runtime-library directories.

## Reproduce the paper

```sh
python tools/reproduce_main.py
python tools/reproduce_figures.py
```

The first command regenerates the full-weld image, partial-weld image and
matched partial-weld control without early-arrival suppression. It checks
nine signal-case metric sets against the paper. Use `--case full` or
`--case partial` to compute just one specimen.

The second command regenerates the baseline comparisons, explanatory diagrams
and five final partial-weld review figures. Use `--group baseline` or
`--group review` for just that part. Fine-grid calculations need several GB of
working space and memory.

The supplied recordings and reference results stay unchanged. Figure work is
assembled under `output/replay/`. Set `SAFT_WORKSPACE` to another empty directory
for a fresh calculation. If the analysis code or input definitions change, use
a fresh workspace so old computed results cannot enter the comparison. Run the
two reproduction commands sequentially because they share intermediate files.

## Analyze another specimen

A small synthetic example demonstrates the full sequence:

```sh
python tools/make_example.py
python tools/validate_dataset.py output/example_dataset
python tools/analyze.py --dataset output/example_dataset --request output/example_dataset/analysis.json --output output/example_analysis
```

Then replace the example with your own geometry, source/receiver positions and
recorded velocities using [the input guide](../saft/FORMAT.md). The outputs
include image arrays, plots, peak locations and optional defect-region contrasts.

## Check the calculations

| Check | Command |
|---|---|
| Original model/recording contents and complete paper | `python tools/verify_repository.py` |
| Agreement between the four SAFT programs | `python tools/verify_interfaces.py` |
| Independent ray times, focusing, receiver mapping and boundary response | `python studies/laser/src/verify_saft_datasets.py` |
| Numerical and data-format tests | `python -m unittest discover -s tests -v` |

Run Python normally, without `-O`, so scientific assertions remain active.
The checks produce reports under `output/`. Reference records are kept in
`verification/`; `replay_inputs.json` identifies the inputs used for the paper.

## Rebuild the paper PDF

Install Tectonic, then run `python tools/build_paper.py`. Use `--tectonic PATH`
if needed. Its first run downloads standard TeX resources. It writes
`paper/paper.pdf` and checks citations, cross-references and text overflow.

After intentional changes to the retained scientific files, review the changes
and use `python tools/verify_repository.py --write-manifest` to update the
recorded file checksums. Temporary builds, simulations and plots belong under
`build/` or `output/`, which are excluded from Git.
