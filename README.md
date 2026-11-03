# Laser-weld models and SAFT

Standalone scientific inputs and analysis for the full weld and the 9 mm
partial J-bevel weld: original COMSOL models, verified receiver histories,
headless C++ SAFT, final chapter/appendix, and their numerical evidence.

| Directory | Contents |
|---|---|
| [models/](models/README.md) | Full-weld and partial-weld COMSOL model inputs, acquisition metadata, and model checksums |
| [data/](data/README.md) | Both acquisition manifests and all 120 velocity-history CSVs |
| [saft/](saft/README.md) | C++17 reconstruction engine and four headless scientific exporters |
| [paper/](paper/README.md) | Final chapter, detailed appendix, figures, bibliography, appendix PDF and electronic supplement |
| [studies/laser/](studies/laser/README.md) | Analysis/figure scripts, frozen protocols, comparison records and compact figure arrays |
| [tools/](tools/README.md) | Integrity checks, main-result replay and isolated study/figure runners |
| `verification/` | Import provenance, snapshot hashes and validation receipts |

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

The retained subset covers the final chapter and detailed appendix. It includes
all 20,412 partial-review comparison rows and the 384 baseline image
configurations used for 1,152 signal-case comparisons. The full-weld processed
bank is regenerated from the recordings using the exact saved filter recipe.
Large intermediate survey arrays, caches, duplicate native-time exports, old
paper versions, GUI builds, toolchains and unrelated experiments are omitted.
Their required scientific protocols and reconstruction code are retained.

The `.mph` files are model inputs rather than saved time-dependent solutions.
Open them with COMSOL Multiphysics 6.4 and the relevant physics licences.
The model files and scientific document contents are preserved from the
original package. Repeating finite-element solves requires COMSOL; replaying
SAFT and regenerating the figures does not. COMSOL solves were not repeated
during packaging. See [paper/README.md](paper/README.md) for report integration
and the standalone LaTeX driver. All cited bibliography entries are included.
