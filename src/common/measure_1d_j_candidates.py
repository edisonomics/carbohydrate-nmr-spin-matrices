#!/usr/bin/env python3
"""Estimate candidate scalar couplings from processed 1-D proton spectra.

The processed Bruker ``1r`` spectrum is the frequency-domain representation of
the raw ``fid``.  This tool records both paths, detects resolved positive lines
near provisional proton shifts, and reports line spacings in Hz.  It is a
screening/measurement aid: crowded or second-order multiplets remain REVIEW
and no matrix or J value is modified automatically.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import numpy as np
from scipy.signal import find_peaks

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bruker_metadata import parse_jcamp, read_processed_1r, required_number  # noqa: E402
from carbohydrate_config import load_config  # noqa: E402


def load_spectrum(experiment: Path, procno: str) -> dict[str, Any]:
    processed = experiment / "pdata" / procno
    acqus_path, procs_path, spectrum_path = experiment / "acqus", processed / "procs", processed / "1r"
    for path in (acqus_path, procs_path, spectrum_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    acqus, procs = parse_jcamp(acqus_path), parse_jcamp(procs_path)
    si = required_number(procs, "SI", int)
    sf = required_number(procs, "SF", float)
    sw_p = required_number(procs, "SW_p", float)
    offset = required_number(procs, "OFFSET", float)
    intensity = np.asarray(read_processed_1r(
        spectrum_path, point_count=si,
        bytordp=required_number(procs, "BYTORDP", int),
        dtypp=required_number(procs, "DTYPP", int),
        nc_proc=required_number(procs, "NC_proc", int),
    ), dtype=float)
    step_ppm = (sw_p / sf) / max(si - 1, 1)
    ppm = offset - np.arange(si, dtype=float) * step_ppm
    raw_fid = experiment / "fid"
    # A processed Bruker dataset may have been copied to a new experiment
    # number (for example 1100_MHz/11) while the original raw FID remains in
    # the sibling experiment directory (1100_MHz/1).
    if not raw_fid.is_file():
        sibling_fids = sorted(experiment.parent.glob("*/fid"))
        if len(sibling_fids) == 1:
            raw_fid = sibling_fids[0]
    return {
        "experiment": str(experiment),
        "raw_fid": str(raw_fid) if raw_fid.is_file() else None,
        "processed": str(spectrum_path),
        "field_mhz": required_number(acqus, "SFO1", float),
        "sf_mhz": sf,
        "resolution_hz": sw_p / si,
        "ppm": ppm,
        "intensity": intensity,
        "pulse_program": acqus.get("PULPROG", "").strip("<>") or "unknown",
    }


def detect_lines(spectrum: dict[str, Any], target: float, window: float,
                 prominence_fraction: float, minimum_distance_hz: float) -> dict[str, Any]:
    ppm = spectrum["ppm"]
    intensity = spectrum["intensity"]
    mask = (ppm >= target - window) & (ppm <= target + window)
    if int(mask.sum()) < 8:
        return {"target_ppm": target, "lines": [], "status": "INSUFFICIENT_POINTS"}
    x = ppm[mask]
    y = intensity[mask]
    baseline = float(np.median(y))
    z = y - baseline
    # Absorption-mode lines are positive.  Robust scale limits weak noise
    # peaks without requiring a molecule-specific absolute intensity.
    noise = float(np.median(np.abs(z - np.median(z)))) or 1.0
    prominence = max(noise * 6.0, float(np.max(z)) * prominence_fraction)
    distance = max(1, int(round(minimum_distance_hz / spectrum["resolution_hz"])))
    peaks, props = find_peaks(z, prominence=prominence, distance=distance)
    lines = [
        {"ppm": float(x[index]), "intensity": float(z[index]),
         "prominence": float(props["prominences"][n])}
        for n, index in enumerate(peaks)
    ]
    lines.sort(key=lambda item: item["intensity"], reverse=True)
    return {
        "target_ppm": target,
        "window_ppm": window,
        "noise_proxy": noise,
        "prominence_threshold": prominence,
        "lines": lines,
        "status": "OK" if lines else "NO_RESOLVED_LINES",
    }


def spacings(lines: list[dict[str, Any]], sf_mhz: float) -> list[float]:
    values = sorted(float(line["ppm"]) for line in lines)
    return sorted((values[j] - values[i]) * sf_mhz for i in range(len(values)) for j in range(i + 1, len(values)))


def screen_coupling(coupling: dict[str, Any], field: dict[str, Any]) -> dict[str, Any]:
    expected = abs(float(coupling.get("j_hz", 0.0)))
    atoms = {str(item["atom"]): item for item in field["atoms"]}
    checks = []
    for atom_id in (str(coupling["i"]), str(coupling["j"])):
        values = [value for value in atoms.get(atom_id, {}).get("pairwise_spacings_hz", []) if 0.5 <= value <= 20.0]
        if values:
            nearest = min(values, key=lambda value: abs(value - expected))
            checks.append({"atom": atom_id, "nearest_spacing_hz": nearest,
                           "absolute_difference_hz": abs(nearest - expected)})
        else:
            checks.append({"atom": atom_id, "nearest_spacing_hz": None,
                           "absolute_difference_hz": None})
    return {"i": coupling["i"], "j": coupling["j"], "provisional_j_hz": coupling.get("j_hz"),
            "field_mhz": field["field_mhz"], "screening": checks,
            "interpretation": "screening only; nearest line spacing is not an automatic J assignment"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    parser.add_argument("--fields", nargs="+", default=None)
    parser.add_argument("--window-ppm", type=float, default=0.035)
    parser.add_argument("--prominence-fraction", type=float, default=0.08)
    parser.add_argument("--minimum-distance-hz", type=float, default=0.8)
    args = parser.parse_args()

    root = args.repo_root.resolve()
    config = load_config(root, args.molecule)
    seed_path = root / "data" / args.molecule / f"{args.molecule}_provisional_observations.json"
    if not seed_path.is_file():
        raise FileNotFoundError(seed_path)
    seed = json.loads(seed_path.read_text(encoding="utf-8"))
    atoms = list(seed.get("atoms", []))
    wanted = {str(field) for field in args.fields} if args.fields else None
    datasets = [item for item in config.get("datasets", [])
                if (wanted is None or str(item.get("key")) in wanted)
                and str(item.get("nucleus", "1H")).upper() == "1H"]
    if not datasets:
        raise SystemExit("No selected 1H datasets found")

    output_rows: list[dict[str, Any]] = []
    field_reports: list[dict[str, Any]] = []
    for dataset in datasets:
        experiment = root / "data" / args.molecule / str(dataset["relative_dir"])
        spectrum = load_spectrum(experiment, str(dataset.get("procno", "1")))
        atom_reports = []
        for atom in atoms:
            result = detect_lines(spectrum, float(atom["shift_ppm"]), args.window_ppm,
                                  args.prominence_fraction, args.minimum_distance_hz)
            result.update({"atom": atom["id"], "shift_ppm": atom["shift_ppm"]})
            result["pairwise_spacings_hz"] = spacings(result["lines"], spectrum["sf_mhz"])
            atom_reports.append(result)
        field_reports.append({
            "field_mhz": dataset.get("field_mhz"),
            "acquisition": dataset.get("acquisition"),
            "experiment": spectrum["experiment"],
            "raw_fid": spectrum["raw_fid"],
            "processed_spectrum": spectrum["processed"],
            "pulse_program": spectrum["pulse_program"],
            "sf_mhz": spectrum["sf_mhz"],
            "resolution_hz": spectrum["resolution_hz"],
            "atoms": atom_reports,
        })
        for atom_report in atom_reports:
            for spacing in atom_report["pairwise_spacings_hz"]:
                output_rows.append({
                    "field_mhz": dataset.get("field_mhz"),
                    "atom": atom_report["atom"],
                    "target_shift_ppm": atom_report["shift_ppm"],
                    "spacing_hz": spacing,
                    "raw_fid": spectrum["raw_fid"],
                    "processed_spectrum": spectrum["processed"],
                })

    output_dir = root / "outputs" / args.molecule / "j_measurements"
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "j_candidate_measurements.json"
    csv_path = output_dir / "j_candidate_spacings.csv"
    report = {
        "molecule": args.molecule,
        "status": "REVIEW",
        "seed_file": str(seed_path),
        "method": "resolved positive line spacing in processed 1r spectra; raw FID provenance recorded",
        "parameters": {
            "window_ppm": args.window_ppm,
            "prominence_fraction": args.prominence_fraction,
            "minimum_distance_hz": args.minimum_distance_hz,
        },
        "fields": field_reports,
        "coupling_screening": [screen_coupling(coupling, field)
                               for field in field_reports for coupling in seed.get("couplings", [])],
        "interpretation": "Candidate spacings require manual multiplet assignment and cross-field agreement before entering the spin matrix.",
        "matrix_updated": False,
    }
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["field_mhz", "atom", "target_shift_ppm", "spacing_hz", "raw_fid", "processed_spectrum"])
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"Wrote {json_path}")
    print(f"Wrote {csv_path}")
    print("STATUS: REVIEW — candidate spacings only; matrix unchanged.")
    for field in field_reports:
        print(f"{field['field_mhz']} MHz: raw FID={'yes' if field['raw_fid'] else 'no'}, resolution={field['resolution_hz']:.3f} Hz/point")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
