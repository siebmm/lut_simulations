"""Read saved MPH metadata and mesh hashes; do not solve or modify models.

The derived wavelength/size ratios use configured size limits, not measured
element edges. They are resolution diagnostics, not convergence certificates.
"""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
from zipfile import ZipFile

STUDY = Path(__file__).resolve().parents[1]
SOURCES = ["shoulder_left", "bead_left", "bead_centre", "bead_right", "shoulder_right"]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def walk(node):
    yield node
    for child in node.get("nodes", []):
        yield from walk(child)


def settings(node):
    return {s["name"]: s for s in node["settings"]}


def inspect(path):
    raw = path.read_bytes()
    with ZipFile(path) as z:
        smodel = z.read("smodel.json")
        dmodel = z.read("dmodel.xml")
        mesh_hash = sha(z.read("mesh1.mphbin"))
    nodes = list(walk(json.loads(smodel)))
    bypath = {n.get("modelEntityPath"): n for n in nodes}
    params = {s["name"]: s for n in nodes if n.get("apiClass") == "ModelParamGroup"
              for s in n["settings"]}
    parameter_names = ["mesh_fmax", "mesh_nlambda", "mesh_factor", "steel_cs", "steel_cp",
                       "dt_out", "laser_dt_cap"]
    values = {k: float(params[k]["scalarReal"]) for k in parameter_names}
    physics = settings(bypath["/physics/elte"])
    mesh = settings(bypath["/mesh/mesh1"])
    solver = settings(bypath["/sol/sol1/feature/t1"])
    if mesh["meshtype"]["apiValue"] != "manual":
        raise ValueError(f"Expected the user-controlled mesh: {path}")
    limits = {}
    xml = ET.fromstring(dmodel)
    for region, tag, expression, factor in [
        ("bulk", "size", "mesh_factor*steel_cs/mesh_fmax/mesh_nlambda", values["mesh_factor"]),
        ("weld", "weldSize", "steel_cs/mesh_fmax/mesh_nlambda", 1),
    ]:
        size = settings(bypath[f"/mesh/mesh1/feature/{tag}"])
        if size["hmax"]["value"] != expression:
            raise ValueError(f"Unrecognized size expression in {path}: {tag}")
        feature = xml.find(f".//MeshFeature[@tag='{tag}']")
        if feature is None or feature.findtext("status") != "BUILT":
            raise ValueError(f"No built size feature in {path}: {tag}")
        key = feature.findtext("key", "")
        match = re.search(r"hmax([0-9.Ee+\-]+)h", key)
        if match is None:
            raise ValueError(f"Missing stored maximum size in {path}: {tag}")
        stored = float(match.group(1))
        calculated = factor * values["steel_cs"] / values["mesh_fmax"] / values["mesh_nlambda"]
        if abs(stored - calculated) > 1e-14:
            raise ValueError(f"Stored size and parameter expression disagree: {path}")
        limits[region] = {"expression": expression, "configured_hmax_m": calculated,
                          "stored_built_hmax_m": stored,
                          "named_selection": feature.findtext("selection/named")}
    if limits["weld"]["named_selection"] != "/selection/geom1_weld_dom":
        raise ValueError(f"Unexpected weld-domain selection: {path}")
    validation_path = path.with_name("validation.json")
    validation_raw = validation_path.read_bytes()
    validation = json.loads(validation_raw)
    return {
        "path": path.relative_to(STUDY).as_posix(), "mph_sha256": sha(raw),
        "smodel_sha256": sha(smodel), "dmodel_sha256": sha(dmodel),
        "mesh_binary_sha256": mesh_hash, "parameters_SI": values,
        "mesh_type": mesh["meshtype"]["apiValue"], "mesh_size_limits": limits,
        "velocity_order": int(physics["order_velocity"]["apiValue"]),
        "strain_order": int(physics["order_strain"]["apiValue"]),
        "physics_node_fmax_display": physics["fmax"]["value"],
        "solver_method": solver["timemethodexp"]["apiValue"],
        "solver_order": int(solver["erkorder"]["apiValue"]),
        "validation_path": validation_path.relative_to(STUDY).as_posix(),
        "validation_sha256": sha(validation_raw),
        "historically_recorded_max_native_step_s": {
            k: v["max_native_step_s"] for k, v in validation["groups"].items()},
    }


def main():
    root = STUDY / "results/full_weld_complete/sources"
    models = [inspect(root / source / case / "configured.mph")
              for source in SOURCES for case in ["reference", "lof"]]
    common_keys = ["parameters_SI", "mesh_type", "mesh_size_limits", "velocity_order",
                   "strain_order", "physics_node_fmax_display", "solver_method", "solver_order"]
    common = {k: models[0][k] for k in common_keys}
    if any({k: m[k] for k in common_keys} != common for m in models):
        raise ValueError("Bandwidth-related settings differ between acquisitions")
    mesh_hashes = sorted({m["mesh_binary_sha256"] for m in models})
    if len(mesh_hashes) != 1:
        raise ValueError("Saved meshes differ between acquisitions")
    p = common["parameters_SI"]
    rows = []
    for region, limits in common["mesh_size_limits"].items():
        h = limits["configured_hmax_m"]
        for mode, speed in [("S", p["steel_cs"]), ("L", p["steel_cp"])]:
            rows.append({"region": region, "mode": mode, "configured_hmax_mm": h * 1e3,
                         "wavelength_over_hmax": {str(f): speed / (f * 1e6 * h) for f in [5, 8, 10]},
                         "frequency_at_1_5_elements_per_wavelength_MHz": speed / (1.5 * h) / 1e6})
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "Saved-model bandwidth audit accompanying Step 3; no new COMSOL solve",
        "script_sha256": sha(Path(__file__).read_bytes()),
        "model_count": len(models), "common_settings": common,
        "all_ten_mesh_binaries_identical": len(mesh_hashes) == 1,
        "mesh_binary_sha256": mesh_hashes[0], "models": models,
        "derived_resolution_diagnostics": rows,
        "export_sampling_MHz": 1 / p["dt_out"] / 1e6,
        "export_temporal_nyquist_MHz": .5 / p["dt_out"] / 1e6,
        "plate_receiver_pitch_m": .00025,
        "single_leg_grazing_spatial_nyquist_MHz": {
            "S": p["steel_cs"] / (.0005 * 1e6), "L": p["steel_cp"] / (.0005 * 1e6)},
        "limitations": [
            "Configured hmax is not a measurement of the largest actual mesh edge.",
            "The 1.5-element quartic guideline is not a validated frequency cutoff.",
            "Physics-node fmax is not the driver of this manual mesh; the size expressions use mesh_fmax.",
            "Native time-step values are read from historical validation files, not remeasured.",
            "Identical meshes do not prove cancellation of numerical errors in independently evolved fields.",
            "Output Nyquist and local receiver-pitch bounds do not establish propagation or SAFT accuracy.",
        ],
        "guidance": [
            "https://doc.comsol.com/6.3/doc/com.comsol.help.aco/aco_ug_elastic_waves.06.58.html",
            "https://www.comsol.com/blogs/introduction-to-the-elastic-waves-time-explicit-interface/",
        ],
    }
    output = STUDY / "results/saft_workflow/step_03_bandwidth_audit.json"
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"output": str(output), "models": len(models), "common_mesh": mesh_hashes[0],
                      "resolution_diagnostics": rows}, indent=2))


if __name__ == "__main__":
    main()
