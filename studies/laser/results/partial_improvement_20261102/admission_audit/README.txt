PARTIAL-WELD SAFT ADMISSION AUDIT - 2 NOVEMBER 2026

The 25% receiver-retention rule and the receiver-count floors do not restrict
the current S | SS reconstruction. Lowering them does not recover any pixels.
The angular response threshold does affect the missing near-defect coverage.

WHAT THE 25% RULE MEANS
A source contributes at a pixel when:
  accepted >= max(minimum_receivers, ceil(fraction * eligible))
  effective >= minimum_effective_receivers (native tolerance 1e-10)
The image then requires at least two contributing sources in the main comparison.
Eligible receivers have already passed geometry, angular response, wall,
source-receiver offset and aperture gates. Accepted receivers additionally have
complete filter/interpolation support and pass pulse screening. Thus 25% refers
to this eligible set for one source/pixel, not 25% of all 331 receivers.
Effective count is (sum weights)^2 / sum(weights squared); it measures weight
concentration, not a count of statistically independent observations.

TESTED RANGE
Fractions: 0, 10, 25, 50, 100%; minimum receivers: 1, 4, 8, 16, 32;
minimum effective receivers: 0, 4, 8, 16; minimum sources: 1, 2, 3.
All combinations were replayed for the equal-mean 0.5-5 MHz baseline and the
weighted 1-5 MHz candidate: 600 combinations, each in LOF, pristine and
waveform-difference cases. These are statistic replays, not new COMSOL solves.
Eight separate native reconstructions checked the replay (24 signal cases).
No qualifying source changed across receiver-guard choices. Where eligible
receivers exist, all survive; the minimum accepted count is 224. The weighted
effective count is at least 189.09. These conclusions apply to these fields
and recipes; they do not establish that the guards are always unnecessary.

THE ACTIVE LIMIT
On the 0.1 mm grid, 14 of 330 pixels within 0.5 mm of the defect lack two-source
support. The left shoulder's incoming ray is blocked by geometry, while all
three weld sources fall below the 20% launch-response threshold. Only the
right shoulder is admitted. The known defect line itself already had two-source
support at all 101 sampled locations in the preceding support audit.

Fresh native angular-threshold comparisons, weighted 1-5 MHz, two sources:
  Angular gate  Near support    Peak-to-background contrast
    20%          95.76%             8.59 dB
    15%          96.67%             8.70 dB
    10%         100.00%             0.90 dB
     5%         100.00%             2.19 dB
The Workbench angular threshold applies jointly to launch and receive response;
these tests change both. The 10% run with the fraction guard disabled is also
unchanged. Relaxing to 15% adds three near pixels on the 0.1 mm grid. On the
0.05 mm grid, coverage changes 95.82 -> 96.81%, contrast 9.08 -> 9.10 dB.
On the wider 0.25 mm grid, coverage changes 94.34 -> 96.23%, contrast
9.20 -> 9.12 dB. This is a small coverage gain with essentially unchanged
contrast, not an established contrast improvement. The 10% and 5% gates obtain
100% neighbourhood support but introduce substantially stronger background.

METRICS AND LIMITS
Near support is the fraction of sampled pixels within 0.5 mm of the prescribed
defect segment with a valid image value. It is not detection probability.
Contrast is 20 log10(near peak / distant background peak), with background
more than 1 mm from the segment. Unsupported pixels are excluded. Coverage
changes can expose new background, so contrast is not evaluated on a fixed mask.
The 15% choice is retrospective on one known simulated flaw; it is an optional
sensitivity setting, not a validated optimum or improved defect-tip sizing.
Geometry checks remain active. No COMSOL model, raw data or production engine
was changed. These are reconstruction changes using existing recordings.

FILES (relative to this README)
  saft_partial_weighted_1_5_angular15.json: optional Workbench preset.
  summary.json: compact verified comparison and receiver-count ranges.
  all_results.csv / all_results.json: all threshold replay outcomes.
  protocol.json: declared matrix, anchor settings and engine source hashes.
  statistics/: native per-source arrays and provenance for both anchors.
  qualification_diagnostics.json: source counts at missing near pixels.
  missing_pixel_launch_audit.json: independent ray/launch-gain diagnosis.
  baseline_verification.json / confirmation_verification.json: native parity.
  angular/: requests, arrays and metrics for 15%, 10%, 5% and fraction-off control.
  angular_confirmation/: finer/wider requests, arrays, metrics and grid check.

REPRODUCE FROM THE REPOSITORY ROOT
Prerequisites: existing partial-weld recordings, previous improvement study
settings/arrays, pinned Qt/MinGW toolchain and SAFTFilterStudy executable.
Build the exporter, which links the unchanged production engine:
  powershell -NoProfile -ExecutionPolicy Bypass -File studies/laser/paper/partial_improvement_20261102/build_admission.ps1
Run each command in this order with:
  python studies/laser/paper/partial_improvement_20261102/admission_audit.py COMMAND
Commands: prepare, dump, analyse, confirm, angular, angular_confirm, report.
Completed array exports are reused; frozen protocols must agree before reuse.
