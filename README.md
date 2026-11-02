# Laser-weld models and SAFT

The compact companion to the full AllseasNDT study: the core COMSOL models,
the C++ SAFT tool, and the final COMSOL chapter with its detailed SAFT appendix.

| Directory | Contents |
|---|---|
| [models/](models/README.md) | Full-weld and partial-weld COMSOL model inputs, acquisition metadata, and model checksums |
| [saft/](saft/README.md) | C++17 reconstruction engine, two command-line interfaces, and a partial-weld request example |
| [paper/](paper/README.md) | Final chapter, detailed appendix, required figures, bibliography entries, and the existing appendix PDF |

## Build the SAFT tool

Install CMake 3.21 or newer, a C++17 compiler, and Qt 6 Core development files
for that compiler. From this directory:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="<Qt installation>"
cmake --build build --config Release
```

The executables are `SAFTFilterStudy` and `SAFTReview`. They are written to
`build/bin/` (or a configuration subdirectory such as `build/bin/Release/`).
See [saft/README.md](saft/README.md) for inputs and invocation.

## Scope

Receiver histories, reconstruction outputs, parameter sweeps, analysis notebooks
and scripts, machine-specific batch launchers, packaged Windows runtimes, and
release-verification bundles belong in the full AllseasNDT study. They are not
part of this source-and-paper repository. Running SAFT requires an external
acquisition dataset; this repository alone does not replay the published results.

The `.mph` files are model inputs rather than saved time-dependent solutions.
Open them with COMSOL Multiphysics 6.4 and the relevant physics licences.
The model files and scientific document contents are preserved from the
original package.
