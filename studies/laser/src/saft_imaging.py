"""Step 4: finite-support quadrature bandpass and direct bistatic DAS.

Read-only float64 input. Geometry and masks never use the known defect.
Complex sums retain phase within each shot; multi-shot RMS is the default.
"""
import base64
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np

DT = 1e-8
FS = 1 / DT
TAPS = 657
HALF = TAPS // 2
BETA = .1102 * (55 - 8.7)
TRANSITION_HZ = .5e6
BANDS = [("broad", .5, 5.), ("low", .5, 1.5), ("middle", 1.5, 3.), ("high", 3., 5.)]
MODES = ["LL", "LS", "SL", "SS"]  # Source-to-pixel, then pixel-to-receiver.
GROUPS = ["left", "weld", "right"]
COHORTS = ["all", "clean"]
X_MM = np.arange(-15, 15.001, .25)
DEPTH_MM = np.arange(2, 24.001, .25)
MIN_OFFSET_M = .002


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def kernel(low_MHz, high_MHz):
    """Finite complex FIR: real part even, imaginary part odd; positive-band gain 2.

    Nominal edges are approximately half-amplitude. The real part is a
    symmetric bandpass; its quadrature partner shares its exact support.
    """
    if not 0 < low_MHz < high_MHz < FS / 2e6:
        raise ValueError("Invalid band edges")
    n = np.arange(-HALF, HALF + 1)
    width, centre = (high_MHz - low_MHz) * 1e6, (high_MHz + low_MHz) * .5e6
    h = 2 * width / FS * np.sinc(width * n / FS) * np.exp(2j * np.pi * centre * n / FS)
    h *= np.kaiser(TAPS, BETA)
    h *= 2 / np.sum(h * np.exp(-2j * np.pi * centre * n / FS)).real
    return h


def filter_record(values, h):
    """Linear convolution, centred without delay; padded edges are later excluded."""
    values = np.asarray(values, dtype=np.float64)
    h = np.asarray(h, dtype=np.complex128)
    if values.ndim != 2 or h.ndim != 1 or len(h) % 2 != 1 or not np.isfinite(values).all():
        raise ValueError("Expected finite receiver-by-time data and an odd FIR")
    length = 1 << (values.shape[-1] + len(h) - 2).bit_length()
    result = np.fft.ifft(np.fft.fft(values, n=length, axis=-1) * np.fft.fft(h, n=length), axis=-1)
    half = len(h) // 2
    return result[:, half:half + values.shape[-1]]


def filter_evidence(band, h):
    ident, low, high = band
    nfft = 131072
    f = np.fft.fftfreq(nfft, DT)
    H = np.fft.fft(h, nfft) * np.exp(2j * np.pi * f * HALF * DT) / 2
    passed = (f >= (low + .25) * 1e6) & (f <= (high - .25) * 1e6)
    stopped = (f <= (low - .25) * 1e6) | (f >= (high + .25) * 1e6)
    deviation = float(np.max(np.abs(H[passed] - 1)))
    stop_db = float(20 * np.log10(np.max(np.abs(H[stopped]))))
    if deviation > .005 or stop_db > -50:
        raise ValueError("Filter fails the measured passband/stopband specification")
    shown = (f >= 0) & (f <= 8e6) & (np.arange(nfft) % 8 == 0)
    return {"id": ident, "nominal_edges_MHz": [low, high],
            "passband_MHz": [low + .25, high - .25],
            "stopband_edges_MHz": [low - .25, high + .25],
            "maximum_complex_passband_error": deviation, "maximum_stopband_dB": stop_db,
            "positive_frequency_gain": 2, "coefficient_real": h.real.tolist(),
            "coefficient_imag": h.imag.tolist(), "frequency_MHz": (f[shown] / 1e6).tolist(),
            "normalized_response_dB": (20 * np.log10(np.maximum(np.abs(H[shown]), 1e-12))).tolist()}


def grid_points(thickness_m):
    xx, depth = np.meshgrid(X_MM / 1000, DEPTH_MM / 1000)
    return np.column_stack((xx.ravel(), thickness_m - depth.ravel()))


def direct_times(source_xy, receiver_xy, points, speeds, trigger_s):
    """No reflecting surfaces: mode conversion, when hypothesized, is at the pixel."""
    source_xy, receiver_xy, points = map(lambda x: np.asarray(x, float), (source_xy, receiver_xy, points))
    cs, cr = speeds
    if cs <= 0 or cr <= 0:
        raise ValueError("Wave speeds must be positive")
    source_leg = np.linalg.norm(points - source_xy, axis=1)
    receive_leg = np.linalg.norm(points[:, None, :] - receiver_xy[None, :, :], axis=-1)
    return trigger_s + source_leg[:, None] / cs + receive_leg / cr


def interpolation_plan(times, repairs, samples=5001, half=HALF, dt=DT):
    """Reject every output whose two interpolation stencils touch repairs or edges.

    Repair bounds receive an additional native-export neighbour guard of dt.
    Both bracketing samples are required, even at an exactly integer delay.
    """
    times = np.asarray(times, dtype=float)
    if not np.isfinite(times).all():
        raise ValueError("Nonfinite delays")
    pos = times / dt
    lower = np.floor(pos).astype(np.int32)
    alpha = pos - lower
    first, last = lower - half, lower + 1 + half
    valid = (first >= 0) & (last < samples)
    clean = valid.copy()
    for gap in repairs:
        start, end = gap["start_us"] * 1e-6 - dt, gap["end_us"] * 1e-6 + dt
        clean &= ~((first * dt <= end) & (last * dt >= start))
    return {"lower": np.clip(lower, 0, samples - 2), "alpha": alpha, "all": valid, "clean": clean}


def sample_sum(filtered, plan, mask):
    """Signed complex interpolation followed by receiver summation, never magnitude first."""
    r = np.arange(filtered.shape[0])[None, :]
    lower, a = plan["lower"], plan["alpha"]
    values = (1 - a) * filtered[r, lower] + a * filtered[r, lower + 1]
    return np.sum(np.where(mask, values, 0), axis=1)


def paired_repairs(audit, source, group):
    return [gap for row in audit["records"] if row["source"] == source
            for g in row["groups"] if g["group"] == group for gap in g["repair_intervals"]]


def prepare_geometry(audit, arrivals, shot, points):
    speed = {"L": arrivals["material"]["cl_m_s"], "S": arrivals["material"]["cs_m_s"]}
    plans, counts, eligible = {}, {}, {}
    for group in GROUPS:
        rr = [r for r in audit["receivers"] if r["group"] == group]
        xy = np.array([[r["x_m"], r["y_m"]] for r in rr])
        admitted = np.abs(xy[:, 0] - shot["xy_m"][0]) >= MIN_OFFSET_M - 1e-14
        eligible[group] = int(admitted.sum())
        for mode in MODES:
            times = direct_times(shot["xy_m"], xy, points, [speed[m] for m in mode], shot["source_trigger_s"])
            plan = interpolation_plan(times, paired_repairs(audit, shot["id"], group))
            for cohort in COHORTS:
                plan[cohort] &= admitted[None, :]
                counts[f"{mode}/{group}/{cohort}"] = plan[cohort].sum(axis=1).astype(np.uint16)
            plans[(mode, group)] = plan
    return plans, counts, eligible


def encode(a, dtype):
    a = np.asarray(a)
    rounded = np.asarray(a, dtype=dtype)
    error = float(np.max(np.abs(a - rounded.astype(a.dtype))))
    peak = float(np.max(np.abs(a)))
    return base64.b64encode(rounded.tobytes()).decode("ascii"), {
        "max_absolute_error": error, "peak_relative_error": error / peak if peak else 0}


def write_packet(path, callback, payload):
    path.write_text(callback + "(" + json.dumps(payload, separators=(",", ":"), allow_nan=False) + ");\n", encoding="utf-8")
    return {"path": "imaging_data/" + path.name, "sha256": sha(path), "bytes": path.stat().st_size}


def source_mean(sums, counts, eligible, mode, cohort, case="difference", groups=GROUPS, floor=.25):
    n = sum(counts[f"{mode}/{g}/{cohort}"] for g in groups).astype(float)
    threshold = max(8, int(np.ceil(floor * sum(eligible[g] for g in groups))))
    total = sum(sums[f"{mode}/{g}/{cohort}/{case}"] for g in groups)
    return np.divide(total, n, out=np.zeros_like(total), where=n > 0), n >= threshold, n


def combine_sources(means, valid, coherent=False, minimum_sources=2):
    means, valid = np.asarray(means), np.asarray(valid, bool)
    count = valid.sum(axis=0)
    if coherent:
        value = np.abs(np.sum(np.where(valid, means, 0), axis=0) / np.maximum(count, 1))
    else:
        value = np.sqrt(np.sum(np.where(valid, np.abs(means)**2, 0), axis=0) / np.maximum(count, 1))
    value[count < minimum_sources] = np.nan
    return value, count


def generate(audit, arrays, arrivals, out, evidence, paper):
    points = grid_points(arrivals["material"]["thickness_m"])
    folder = out / "imaging_data"
    folder.mkdir(exist_ok=True)
    filters = {b[0]: kernel(*b[1:]) for b in BANDS}
    report = {"schema_version": 1, "step": 4, "status": "passed",
              "created_utc": datetime.now(timezone.utc).isoformat(),
              "baseline_input_fingerprint": audit["input_fingerprint"],
              "physics_input_fingerprint": arrivals["physics_input_fingerprint"],
              "implementation_sha256": sha(__file__), "modes": MODES, "groups": GROUPS,
              "grid": {"x_mm": X_MM.tolist(), "depth_mm": DEPTH_MM.tolist(),
                       "order": "depth-major then x; depth measured downward from OD", "pixels": len(points)},
              "filters": [filter_evidence(b, filters[b[0]]) for b in BANDS],
              "method": {"taps": TAPS, "half_support_us": HALF * DT * 1e6, "dt_s": DT,
                         "kaiser_beta": BETA, "transition_width_MHz": TRANSITION_HZ / 1e6,
                         "repair_neighbour_guard_us": DT * 1e6, "minimum_source_receiver_offset_mm": 2,
                         "default_coverage_fraction": .25, "minimum_receivers_per_source": 8,
                         "minimum_sources_for_multisource_image": 2,
                         "default_combination": "RMS of absolute coherent receiver means per source",
                         "delays": "source trigger + source-to-pixel / mode1 speed + pixel-to-receiver / mode2 speed",
                         "quadrature": "Finite complex FIR with approximately analytic positive-band response; no global Hilbert transform",
                         "repair_policy": "Same paired-union support mask for reference, LOF and difference; no zeroing of gaps",
                         "precision": "Filtering and DAS in float64/complex128; separate reference and pre-subtracted difference sums published as float32"},
              "sources": [], "summary": [], "limitations": [
                  "Direct-path hypotheses only. LL, LS, SL, SS describe the source and receiver legs, with conversion at the trial pixel.",
                  "No reflected-path reconstruction, fitted time shifts, directivity correction, gain compensation or source-polarity correction.",
                  "Reference subtraction is an ideal simulation diagnostic, not demonstrated field-inspection performance.",
                  "The 0.5-5 MHz nominal band has finite transitions; the upper transition extends to 5.25 MHz and does not validate that content.",
                  "The quadrature FIR is only approximately analytic, with its negative-frequency residual measured in the stopband checks.",
                  "Per-source receiver means preserve velocity units but are not calibrated defect amplitudes. Coverage and aperture vary spatially.",
                  "RMS source combination avoids cross-shot phase cancellation; it has no coherent transmit-aperture resolution claim.",
                  "Repair-included versus repair-excluded changes also change receiver coverage and cannot isolate repair error.",
                  "A sharp maximum is not automatically a defect tip. The known outline is a hidden benchmark overlay and never enters reconstruction.",
                  "Model convergence above or below 5 MHz, experimental noise and physical sizing accuracy remain unverified."]}
    aggregates = {}
    for shot in audit["sources"]:
        source = shot["id"]
        print("Direct SAFT geometry: " + source, flush=True)
        plans, counts, eligible = prepare_geometry(audit, arrivals, shot, points)
        payload = {"key": source, "fingerprint": audit["input_fingerprint"],
                   "implementation": report["implementation_sha256"], "counts": {}, "eligible": eligible}
        for key, a in counts.items():
            payload["counts"][key] = encode(a, "<u2")[0]
        entry = {"id": source, "eligible": eligible, "geometry": write_packet(folder / (source + "_geometry.js"), "window.SAFT_IMAGE_GEOMETRY", payload), "packets": []}
        for band, _, _ in BANDS:
            for component in ["vy", "vx"]:
                print(f"Direct SAFT: {source}/{band}/{component}", flush=True)
                sums = {}
                for group in GROUPS:
                    suffix = ("vx",) if component == "vx" else ()
                    ref, lof = [arrays[(source, case, group) + suffix] for case in ["reference", "lof"]]
                    for case, values in [("reference", ref), ("difference", lof - ref)]:
                        filtered = filter_record(values, filters[band])
                        for mode in MODES:
                            plan = plans[(mode, group)]
                            for cohort in COHORTS:
                                sums[f"{mode}/{group}/{cohort}/{case}"] = sample_sum(filtered, plan, plan[cohort])
                packet = {"key": f"{source}/{band}/{component}", "fingerprint": audit["input_fingerprint"],
                          "implementation": report["implementation_sha256"], "sums": {}}
                errors = {}
                for key, a in sums.items():
                    # Each pixel: real then imaginary. Difference is never recovered by subtracting rounded large fields.
                    packet["sums"][key], errors[key] = encode(np.column_stack((a.real, a.imag)), "<f4")
                entry["packets"].append({"band": band, "component": component, "precision": errors,
                                         **write_packet(folder / f"{source}_{band}_{component}.js", "window.SAFT_IMAGE_DATA", packet)})
                if component == "vy":
                    for mode in MODES:
                        for cohort in COHORTS:
                            for case in ["reference", "difference"]:
                                mean, valid, n = source_mean(sums, counts, eligible, mode, cohort, case)
                                aggregates.setdefault((band, mode, cohort, case), []).append((mean, valid, n))
        report["sources"].append(entry)
    maps = {}
    for (band, mode, cohort, case), rows in aggregates.items():
        value, count = combine_sources([r[0] for r in rows], [r[1] for r in rows])
        finite = np.isfinite(value)
        peak = int(np.nanargmax(value)) if finite.any() else None
        result = {"band": band, "mode": mode, "cohort": cohort, "case": case, "component": "vy",
                  "source_selection": "all", "combination": "rms", "coverage_fraction": .25,
                  "valid_pixels": int(finite.sum()), "total_pixels": len(points),
                  "maximum_m_s": float(value[peak]) if peak is not None else None,
                  "maximum_x_mm": float(points[peak, 0] * 1000) if peak is not None else None,
                  "maximum_depth_mm": float(arrivals["material"]["thickness_m"] * 1000 - points[peak, 1] * 1000) if peak is not None else None,
                  "source_count_at_maximum": int(count[peak]) if peak is not None else 0}
        report["summary"].append(result)
        maps[(band, mode, cohort, case)] = value.reshape(len(DEPTH_MM), len(X_MM))
    report["figures"] = make_figures(report, maps, evidence / "figures")
    (evidence / "step_04_imaging.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (out / "step_04_imaging.json").write_bytes((evidence / "step_04_imaging.json").read_bytes())
    (out / "imaging.js").write_text("window.SAFT_IMAGING=" + json.dumps(report, separators=(",", ":"), allow_nan=False) + ";\n", encoding="utf-8")
    write_latex(report, paper / "full_weld_saft_step_04.tex")
    return report


def make_figures(report, maps, folder):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    folder.mkdir(exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(8, 7), constrained_layout=True)
    for ax, mode in zip(axes.flat, MODES):
        a = maps[("middle", mode, "clean", "difference")]
        peak = np.nanmax(a) if np.isfinite(a).any() else 1
        db = 20 * np.log10(np.maximum(a / max(peak, 1e-30), 1e-3))
        im = ax.imshow(db, origin="upper", extent=[X_MM[0] - .125, X_MM[-1] + .125, DEPTH_MM[-1] + .125, DEPTH_MM[0] - .125],
                       vmin=-40, vmax=0, cmap="magma", aspect="equal", interpolation="nearest")
        ax.set(title=mode + " | peak " + f"{peak:.2e} m/s", xlabel="x (mm)", ylabel="OD depth (mm)")
    fig.colorbar(im, ax=axes.ravel().tolist(), label="dB relative to each panel maximum")
    fig.suptitle("Direct paths | nominal 1.5-3 MHz | vertical LOF-reference\nPaired repair exclusion; RMS across sources; no defect overlay")
    fig.savefig(folder / "direct_saft_baseline.png", dpi=180)
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 3.5), constrained_layout=True)
    for f in report["filters"]:
        ax.plot(f["frequency_MHz"], f["normalized_response_dB"], label=f["id"])
    ax.set(xlim=(0, 6), ylim=(-90, 3), xlabel="Frequency (MHz)", ylabel="Positive-band response / 2 (dB)", title="Finite quadrature filters: 657 taps; half-support 3.28 microseconds")
    ax.axvline(5, color="grey", ls="--", lw=.8)
    ax.grid(alpha=.2); ax.legend(ncol=4)
    fig.savefig(folder / "saft_filter_response.png", dpi=180)
    plt.close(fig)
    return [{"path": n, "sha256": sha(folder / n)} for n in ["direct_saft_baseline.png", "saft_filter_response.png"]]


def write_latex(report, target):
    lines = ["% Generated by saft_imaging.py; do not edit numerical results.",
             r"\subsection{Step 4 numerical results}",
             "The filters use 657 taps at 10 ns spacing: both quadratures have half-support "
             r"$3.28\us$. Values below come from a 131072-point frequency-response check.",
             r"\begin{center}\begin{tabular}{lrrr}\toprule Band (MHz) & Passband (MHz) & Max. error & Stopband (dB)\\\midrule"]
    for f in report["filters"]:
        lo, hi = f["nominal_edges_MHz"]; a, b = f["passband_MHz"]
        lines.append(f"{lo:g}--{hi:g} & {a:g}--{b:g} & {100*f['maximum_complex_passband_error']:.3f}\\% & {f['maximum_stopband_dB']:.2f}" + r"\\")
    lines += [r"\bottomrule\end{tabular}\end{center}",
              "The following whole-grid maxima use vertical paired differences, all source positions, "
              "RMS source combination and paired repair exclusion. A source needs at least "
              "25 percent of its eligible receivers (minimum eight); a pixel needs at least two sources. "
              "These are image maxima, not identified defect tips or measured defect heights.",
              r"\begin{longtable}{llrrrr}\toprule Band & Mode & Valid pixels & Peak (m/s) & $x$ (mm) & Depth (mm)\\\midrule\endhead"]
    for r in report["summary"]:
        if r["case"] != "difference" or r["cohort"] != "clean":
            continue
        if r["maximum_m_s"] is None:
            lines.append(f"{r['band']} & {r['mode']} & 0 & --- & --- & ---" + r"\\")
        else:
            lines.append(f"{r['band']} & {r['mode']} & {r['valid_pixels']} & {r['maximum_m_s']:.3g} & {r['maximum_x_mm']:.2f} & {r['maximum_depth_mm']:.2f}" + r"\\")
    lines += [r"\bottomrule\end{longtable}",
              r"Full evidence: \path{studies/laser/results/saft_workflow/step_04_imaging.json}. "
              "Each lazy-loaded map packet has a checksum and a recorded float32 display-rounding bound."]
    target.write_text("\n\n".join(lines) + "\n", encoding="utf-8")
