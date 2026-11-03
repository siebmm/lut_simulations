# 2D elastic-wave analysis and laser-weld reproductions

Standalone scientific inputs and analysis for the full weld and the 9 mm
partial J-bevel weld: original COMSOL models, verified receiver histories,
headless C++ SAFT, a complete standalone paper, and its numerical evidence.

**Read the [full paper (PDF)](paper/paper.pdf): Laser-ultrasonic imaging of full
and partial welds.** It includes the detailed methods, results and complete
comparison registers.

The shared SAFT engine also accepts other homogeneous, isotropic 2D specimens.
Supply separate geometry and acquisition JSON files plus velocity-history CSVs.
Polygonal boundaries and voids, material speeds, source/receiver layouts, case
names, sample counts and time intervals are defined by the dataset.

## Analyze your own model

After installing the dependencies and building the native tools below:

```sh
python tools/validate_dataset.py path/to/my_dataset
python tools/analyze.py --dataset path/to/my_dataset --request path/to/analysis.json --output output/my_analysis
```

See [the input format and examples](saft/FORMAT.md) for the geometry, CSV contract,
supported physics, reflected paths and analysis settings. Start with a runnable
non-weld example using `python tools/make_example.py`; then use its
`output/example_dataset/analysis.json` with `tools/analyze.py`. Results include
numerical arrays, images, support/peak summaries and optional region contrasts.
The frozen weld studies below retain their original protocols and metric checks.

| Directory | Contents |
|---|---|
| [models/](models/README.md) | Full-weld and partial-weld COMSOL model inputs, acquisition metadata, and model checksums |
| [data/](data/README.md) | Both acquisition manifests and all 120 velocity-history CSVs |
| [saft/](saft/README.md) | C++17 reconstruction engine and four headless scientific exporters |
| [paper/](paper/README.md) | One complete standalone paper, its rendered PDF, source, cited figures, bibliography and electronic supplement |
| [studies/laser/](studies/laser/README.md) | Required replay helpers, canonical reference arrays and an explicit input index |
| [tools/](tools/README.md) | General analysis, integrity checks, main-result and figure replay, and paper build |
| [tests/](tests/test_general_datasets.py) | General-dataset contract and independent numerical regressions |
| `verification/` | Snapshot hashes and current numerical/PDF validation receipts |

## Verify and reproduce

Requirements: Python 3.12, Git LFS, CMake 3.21 or newer, a C++17 compiler,
and matching Qt 6 Core development files. COMSOL is only needed for new
finite-element solves; all replay and verification below use the supplied data.

From a fresh checkout, retrieve the large scientific inputs first:

```sh
git lfs install
git lfs pull
```

Create a Python environment from the repository root:

```sh
python -m venv .venv
```

Activate it with `source .venv/bin/activate` on Linux/macOS or
`.venv\Scripts\activate.bat` in Windows Command Prompt. In PowerShell, you can
use `.venv\Scripts\python.exe` instead of `python` in the commands below,
without changing your script execution policy.

```sh
python -m pip install -r requirements.txt
python tools/verify_repository.py
```

Build the native tools as below, then:

```sh
python tools/verify_interfaces.py
python studies/laser/src/verify_saft_datasets.py
python tools/reproduce_main.py
python tools/reproduce_figures.py
```

The interface check reconstructs nine pixels through all four exporters. The
physics check independently verifies ray travel times and focusing with NumPy
and SciPy, plus receiver mapping, air masks, muting and isolated dataset use.
Reports go to `output/interface_check/check.json` and
`output/physics_check/check.json`. These checks need no desktop application.

The main replay checks all three signal cases against the archived metrics.
Expected LOF results are **12.79250 dB** for the full weld and **14.44694 dB**
for the partial weld, with **9.08623 dB** for the matched unmuted partial control.
These are retrospective contrasts for the prescribed simulated flaw.
Replays write to `output/` and caches to `build/`; original evidence is preserved.
Allow several GB of temporary disk space and RAM for fine-grid reconstruction.

No input is read from the original NDT checkout. The pinned Python dependencies
are the versions used for the archived replay; no Python package installation
of this repository is needed.

## Build the SAFT tool

Install CMake 3.21 or newer, a C++17 compiler, and Qt 6 Core development files
for that compiler. From this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="<Qt installation>"
cmake --build build --config Release
```

The executables are `SAFTFilterStudy`, `SAFTReview`, `SAFTAdmissionAudit` and
`SAFTBatch`. They are written to
`build/bin/` (or a configuration subdirectory such as `build/bin/Release/`).
See [saft/README.md](saft/README.md) for inputs and invocation.

Alternatively, with an installed compatible Qt qmake and compiler:

```sh
python tools/build_native.py --qmake <path-to-qmake> --make <path-to-make>
```

On Windows use Qt's matching MinGW `mingw32-make.exe`. This helper places the
executables and required Windows runtime libraries in the ignored `build/bin/`.
For a CMake build, make Qt/compiler runtime libraries available through `PATH`.
`SAFT_BIN` can select another executable directory; `SAFT_RUNTIME` can add
runtime-library directories to the replay environment.

## Scope

The retained subset supports the standalone paper. It includes
all 20,412 partial-review comparison rows and the 384 baseline image
configurations used for 1,152 signal-case comparisons. The full-weld processed
bank is regenerated from the recordings using the exact saved filter recipe.
Replay inputs are stored once, with duplicate logical paths resolved by
`studies/laser/replay_inputs.json`. Matching comparison records are read directly
from `paper/evidence.zip`. The replay tools materialize those inputs under
`output/replay/`; `SAFT_WORKSPACE` selects another empty output directory.

Only helpers and reference data used by the documented replays remain in the
study tree. Historical packaging, broad survey runners, generated logs and
unused snapshots are omitted. The complete paper comparison record stays in
the electronic supplement; new parameter studies use the generic analysis
workflow. The integrity check rejects unused study files, duplicate contents,
missing inputs and additions absent from the checked snapshot.

The `.mph` files are model inputs rather than saved time-dependent solutions.
Open them with COMSOL Multiphysics 6.4 and the relevant physics licences.
The model files and detailed scientific evidence are preserved from the
original package. Repeating finite-element solves requires COMSOL; replaying
SAFT and regenerating the figures does not. COMSOL solves were not repeated
during packaging. See [paper/README.md](paper/README.md) to rebuild the paper
with `python tools/build_paper.py`. All cited bibliography entries are included.
