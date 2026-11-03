# Local laser generation and source calibration

**Open [laser_source_calibration.mph](laser_source_calibration.mph) in COMSOL 6.4.**
This is the original, unsolved model from
`NDT/studies/laser/models/laser_source_calibration.mph`, copied without changing
its contents. It includes the geometry, mesh, material, two studies and receiver
definitions. It is a local source model, not a weld specimen.

## What the model answers

A laser heats a thin surface layer. Thermal expansion launches longitudinal,
shear and surface waves into the steel. Here we compare its outgoing waves
with those produced by an equivalent tangential surface force, which can then
be used in a larger weld calculation without resolving the heated layer.

**Study 1: heat input and thermoelastic wave generation.** Fourier heat
conduction supplies the temperature field for linear thermal expansion and
dynamic solid mechanics. The mechanical motion does not feed heat back into
the temperature calculation.

**Study 2: equivalent dipole, elastic only.** A tangential force distribution
replaces the heated layer. Its magnitude is fitted to the outgoing velocity
waveforms from Study 1. Both calculations use the same geometry and mesh.

The model is plane strain: the laser is a line extending out of the section.
It uses constant material properties and contains no ablation, plasticity,
optical-field calculation or weld defect. It is a numerical calibration,
not a measurement of a particular laser.

## Geometry and parameters

The simulated half-patch extends from x = 0 to 4 mm and y = -4 to 0 mm.
The top surface is irradiated; symmetry at x = 0 represents the other half.
The remote right and bottom boundaries are free. The short observation window
precedes strong reflections from those artificial boundaries.

| Quantity | Value |
|---|---:|
| Peak absorbed fluence in the original model | 10 J/m² |
| Fluence used for the weld source | 100 J/m² |
| Gaussian pulse duration, FWHM | 6.5 ns |
| Pulse centre | 30 ns |
| Gaussian line half-width, at 1/e² intensity | 0.5 mm |
| Young's modulus / Poisson ratio | 212.443 GPa / 0.288978 |
| Density | 7859.98 kg/m³ |
| Conductivity / heat capacity | 52 W/(m K) / 480 J/(kg K) |
| Thermal expansion coefficient | 11 × 10⁻⁶ K⁻¹ |
| Initial temperature | 293.15 K |
| Recorded interval | 0–0.9 µs, every 2 ns |

Nine graded surface layers resolve the thermal skin: thicknesses range from
50 nm to 12.8 µm, totalling 25.55 µm. The bulk mesh target is 0.08 mm.
Seven receivers lie on a 2 mm radius arc, at 0–90° from the downward vertical;
five more lie on the surface at x = 0, 0.5, 1, 1.5 and 2 mm.

The absorbed surface heat flux is

```text
q(x,t) = F_abs exp[-2(x/w0)^2] g(t)
g(t)   = exp[-(t-t0)^2/(2 sigma^2)] / (sqrt(2 pi) sigma)
sigma  = pulse_fwhm / (2 sqrt(2 ln 2))
```

Its full-width energy per unit line length is `F_abs w0 sqrt(pi/2)`.
Absorbed fluence is the heat deposited in the specimen; incident laser energy
also depends on absorption and illuminated line length.

The equivalent surface traction is

```text
tx(x,t) = s A0 (4x/w0) exp[-2(x/w0)^2] H(t),     ty = 0
A0      = E alpha F_abs / [(1-nu) rho Cp w0]
H(t)    = integral of g(t) from -infinity to t
```

The two halves exert opposite tangential forces. The multiplier `s` is the
calibration result. Both velocity components, angles 15–75° and times
0.12–0.82 µs enter the least-squares fit after identical 5 MHz low-pass filtering.
The 0° and 90° receivers provide comparisons outside the fitted angles.

## Examine the existing comparison

After installing the [Python dependencies](../../tools/README.md#software-setup),
run from the repository root:

```sh
python models/laser_generation/analyze.py --stored-fluence 10
```

This writes `calibration_waveform_comparison.png`, `calibration_temperature.png`
and numerical results to `output/laser_comparison/`. Use `--stored-fluence 100`
for the retained higher-fluence comparison. No COMSOL solve is needed.

At 10 J/m², the recorded peak temperature rise is **7.834 K**, the fitted
dipole multiplier is **1.084434**, and the fitted bulk-wave difference is
**0.931%** in relative L2 norm. The independent one-dimensional heat-conduction
comparison differs by **0.0693%**. Changes to the spot, material, bandwidth or
nearby boundaries need reassessment.

## Repeat the COMSOL calculation

The sequence is: load a working copy, solve Study 1, solve Study 2,
export both sets of receiver histories, integrate thermal energy, then fit the
dipole and plot the comparisons. From the repository root:

```sh
python tools/laser_generation.py --fluence 10 --output output/laser_10
```

Use `--fluence 100 --output output/laser_100` for the weld-source fluence.
Each calculation needs an empty output folder. The command uses the original
saved model by default. Add `--rebuild` to construct the benchmark from
`BuildLaserCalibration.java` before solving. The requested fluence is
applied explicitly in either case.

For a different COMSOL installation, add `--comsol "path/to/COMSOL/bin"`
(`bin/win64` on Windows). COMSOL must have access to its licence server and the
model's heat-transfer and solid-mechanics licences.

For **manual inspection in COMSOL**, add `--prepare-only`. Open the resulting
`output/laser_10/models/laser_source_calibration.mph`; it has the requested
fluence and local export paths. You can inspect the parameters and compute the
two named studies in the GUI. The command-line sequence also performs the
energy integral and writes the physical parameters needed by the comparison.

The original file contains export settings from its first run. The runner
updates export locations and uses COMSOL's velocity variables for cut-point
exports. It also applies the displacement scaling and early time-step limits
used in the historical workflow. The supplied original `.mph` stays unchanged.

## Use the calibration in a weld

`calibration_comparison.json` reports `least_squares_dipole_scale`. For the
included weld models, the corresponding parameter is `laser_fit`; the retained
value is 1.0844342935972646. Keep fluence, width, pulse duration and material
consistent with the calibration. Apply the signed tangential traction in the
local surface direction. A normal pressure pulse describes a different source.

The local model's 0.9 µs exports diagnose source generation. The
[weld models](../README.md) produce the longer propagation histories used for
[SAFT imaging](../../saft/README.md).

## Files in this folder

| File | Role |
|---|---|
| `laser_source_calibration.mph` | Original local COMSOL model |
| `BuildLaserCalibration.java` | Geometry, material, heat flux, dipole, mesh and studies |
| `RunLaserCalibration.java` | Both solves, receiver exports and thermal-energy integral |
| `analyze.py` | Waveform fit, temperature/energy checks and comparison plots |

Reference waveforms are stored once in
`studies/laser/paper/comsol_chapter_20261101/data/source_calibration/`.
The [paper](../../paper/paper.pdf) contains the detailed derivation and references.
