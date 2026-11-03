PARTIAL-WELD RIGHT-RECEIVER PATH STUDY - 2 NOVEMBER 2026

The earlier work did not establish that S|SS is the only useful path.
This study finds a reproducible alternative: S|LL using the rightmost weld
source at x=+1.5 mm and the right-side receiver group. Its fine-grid contrast
is 3.02 dB, near support 96.81%, and near peak 26.97 micrometres/s. The near
peak is 14.37 dB above the matched pristine near peak. The wide-field check
retains positive contrast (4.36 dB) and a peak near the defect. Even/odd
receiver subsets retain positive contrast. Equal receiver averaging loses
the positive contrast; known-response phase/amplitude weighting matters.

S|LL means S from source to image point, then L from image point to backwall
and L from backwall to receiver. Conversion at the image point is a focusing
hypothesis. The ray model does not calibrate the defect scattering coefficient
or prove that a unique physical wave packet follows that route.

SSL|SS using the leftmost weld source also has a defect-related local feature:
fine-grid contrast 3.20 dB with 100% near support, but its near peak is only
0.330 micrometres/s and its wider-field contrast is -28.78 dB. Requiring two
weld sources gives 8.04 dB local contrast but only 43.88% near support and
fails the wider-field check. LSS|S and SLS|LS also fail wider-field checks.
These longer paths are diagnostic views, not replacements for the main image.
Reflected-packet screening and 20% angular controls do not rescue SSL|SS in
the wider field: its contrast remains approximately -27.5 to -28.6 dB.
Enabling all reflected-packet screens removes ALL near support for S|LL.
Its candidate delays therefore overlap rejected reflected-pulse windows;
this screen cannot isolate a surviving near-defect image. The defect/pristine
comparison establishes sensitivity to the defect in the simulation, but does
not uniquely identify an isolated S-to-L defect-scattered packet. Keep this
limitation with the useful alternative rather than treating its path label
as a confirmed physical attribution.

The current all-source/all-receiver S|SS image remains stronger in contrast:
9.10 dB at 0.05 mm with the 15% gate. S|LL supplies another view using a single
weld source and only right receivers. It does not establish improved tip sizing.
In the matched 0.2 mm rightmost-source/right-receiver comparison, S|SS gives
-1.13 dB while S|LL gives +3.03 dB. Thus choosing a path for this acquisition
does matter. Narrowing the right group to fixed near/middle/far windows did
not produce a better robust direct/one-reflection alternative in this scan.

SCOPE
All 132 conventional paths implemented in the native engine were tested with
three weld sources and 157 right receivers (x=4.5534 to 43.5534 mm). The set
includes mode-converted reflections and up to five total legs. The 20 direct
and single-reflection paths were also tested with fixed receiver windows:
  near: x<15 mm, 42 receivers;
  middle: 15<=x<35 mm, 80 receivers;
  far: x>=35 mm, 35 receivers.
There were 192 survey configurations; each retained LOF, pristine and waveform
difference images. Per-source means gave fixed views of each of the three weld
sources and three-source RMS views requiring one or two sources. The source
choice is constant throughout an image; no per-pixel choice of the best source
or path was used. The two-source requirement is not imposed on a single-source
measurement. Another 22 configurations checked selected results; 10 examined
reflected screening and angular thresholds after wider-field failures.

Main settings: weighted 1-5 MHz analytic FIR, 15% angular response threshold,
5% wall threshold, at least 8 real/effective receivers, 25% retained fraction,
CF exponent 1, direct L/S/R pulse screens, no reflected screens initially,
vy, source RMS. Survey grid: 0.2 mm; fine: 0.05 mm; wide: 0.25 mm over
x=+/-45 mm and depth 0.25-25.25 mm. The finer field is x=+/-7.5 mm and
depth 16.65-25.25 mm. Changing spacing and bounds together in the wide run
means its contrast is not a pure field-size sensitivity measurement.

METRICS AND LIMITS
Near support counts valid pixels within 0.5 mm of the known 2.5 mm defect
segment. Contrast is 20 log10(near peak / strongest background more than 1 mm
from the segment). The 0.5-1 mm region is excluded. Coverage is not detection
probability. Unsupported pixels are not tested background. Pristine ratios
use common near support; difference images subtract waveforms before focusing,
not final amplitudes. LOF images do not use pristine subtraction. Measured
amplitudes are simulation outputs, not an experimental energy calibration.
No experimental-noise model, independent defect validation or mesh-convergence
claim is made. Candidate selection is retrospective on one known flaw. This
survey is not exhaustive of all physical scattering mechanisms or settings.

VERIFICATION
All survey three-source RMS images reproduce the exported native image.
Single-source native replay and nested fine-grid checks agree with the survey
source-subset calculations. Protocol files retain input fingerprints and hashes.
All three exported presets reproduce the corresponding LOF image in separate
native LOF-only replays; GUI interaction was not tested.
No production engine, raw data or COMSOL models were changed.

CONFIRMATIONS AND FOLLOW-UP (LOF images)
  S_LL_right__replay                          C=   3.03 dB; near support= 95.18%; near peak=26.84 um/s
  S_LL_right__fine                            C=   3.02 dB; near support= 96.81%; near peak=26.97 um/s
  S_LL_right__wide                            C=   4.36 dB; near support= 96.23%; near peak=26.88 um/s
  S_LL_right__gate20                          C=   3.07 dB; near support= 95.76%; near peak=26.94 um/s
  S_LL_right__band05_gate20                   C=   2.49 dB; near support= 95.76%; near peak=26.12 um/s
  S_LL_right__equal                           C=  -0.42 dB; near support= 96.67%; near peak=20.31 um/s
  S_LL_right__even                            C=   2.93 dB; near support= 96.67%; near peak=29.18 um/s
  S_LL_right__odd                             C=   3.16 dB; near support= 96.67%; near peak=25.21 um/s
  SSL_SS_left__replay                         C=   3.15 dB; near support=100.00%; near peak=0.3256 um/s
  SSL_SS_left__fine                           C=   3.20 dB; near support=100.00%; near peak=0.3297 um/s
  SSL_SS_left__wide                           C= -28.78 dB; near support=100.00%; near peak=0.3179 um/s
  SSL_SS_left__gate20                         C=   7.90 dB; near support= 76.67%; near peak=0.3355 um/s
  SSL_SS_left__band05_gate20                  C=   4.51 dB; near support= 70.91%; near peak=0.3715 um/s
  SSL_SS_left__equal                          C=  -9.22 dB; near support=100.00%; near peak=0.07759 um/s
  SSL_SS_left__even                           C=   5.70 dB; near support=100.00%; near peak=0.4308 um/s
  SSL_SS_left__odd                            C=   4.72 dB; near support=100.00%; near peak=0.4041 um/s
  SSL_SS_two__fine                            C=   8.04 dB; near support= 43.88%; near peak=0.2304 um/s
  SSL_SS_two__wide                            C= -35.52 dB; near support= 45.28%; near peak=0.21 um/s
  LSS_S_left__fine                            C=   4.15 dB; near support= 36.96%; near peak=1.949 um/s
  LSS_S_left__wide                            C= -28.35 dB; near support= 41.51%; near peak=1.942 um/s
  SLS_LS_left__fine                           C=   0.71 dB; near support= 99.54%; near peak=0.01122 um/s
  SLS_LS_left__wide                           C= -63.83 dB; near support=100.00%; near peak=0.01111 um/s
  SSL_SS_left__all_reflected__fine            C=   3.59 dB; near support=100.00%; near peak=0.3819 um/s
  SSL_SS_left__repeated_shear__fine           C=   3.59 dB; near support=100.00%; near peak=0.3819 um/s
  SSL_SS_left__gate20__fine                   C=   7.96 dB; near support= 76.73%; near peak=0.34 um/s
  SSL_SS_left__gate20_all_reflected__fine     C=   7.75 dB; near support= 76.73%; near peak=0.3839 um/s
  SSL_SS_left__all_reflected__wide            C= -27.47 dB; near support=100.00%; near peak=0.3695 um/s
  SSL_SS_left__repeated_shear__wide           C= -27.47 dB; near support=100.00%; near peak=0.3695 um/s
  SSL_SS_left__gate20__wide                   C= -28.60 dB; near support= 79.25%; near peak=0.3246 um/s
  SSL_SS_left__gate20_all_reflected__wide     C= -27.47 dB; near support= 79.25%; near peak=0.3695 um/s
  S_LL_right__all_reflected__fine: no usable near/background comparison
  S_LL_right__all_reflected__wide: no usable near/background comparison

FILES
  figures/comparison.png: current image and two alternatives, with pristine
    controls on a common absolute amplitude scale.
  figures/wide_field.png: wider-field S|LL and SSL|SS comparison.
  figures/paths.png: geometry-valid illustrative rays for these two hypotheses.
  settings/saft_partial_S_LL_right.json: useful alternative, single right weld
    source and right receivers. Fine grid exactly matches the native check.
  settings/saft_partial_SSL_SS_left.json: weaker diagnostic with wide-field clutter.
  settings/saft_partial_SSL_SS_two.json: diagnostic with limited coverage.
  summary.csv / summary.json: compact verification comparisons.
  survey/metrics.csv / metrics.json: all surveyed outcomes, including failures
    and plane-1 diagnostics without CF/pulse screening.
  protocol.json, confirmation_protocol.json, reflection_protocol.json: full
    reconstruction requests and the stated reasons for follow-up.
  survey/, confirmation/, reflection/: arrays, timings and native evidence.

REPRODUCE FROM THE REPOSITORY ROOT
Use the existing COMSOL recordings and pinned native SAFTFilterStudy executable.
  python studies/laser/paper/partial_improvement_20261102/right_paths.py run
  python studies/laser/paper/partial_improvement_20261102/right_paths.py confirm
  python studies/laser/paper/partial_improvement_20261102/right_paths.py reflection
  python studies/laser/paper/partial_improvement_20261102/right_paths_report.py
The review ZIP contains figures, presets, protocols and metrics. Full arrays
and original recordings remain in the repository and are required to rerun.
