#!/usr/bin/env python3
"""Prepare imported Bruker spectra for a carbohydrate molecule.

The molecule name is used only to load its configuration. All acquisition
values and the DSS correction are read from the imported Bruker files.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bruker_metadata import (  # noqa: E402
    find_dss_peak,
    indices_in,
    parse_jcamp,
    read_processed_1r,
    required_number,
)
from carbohydrate_config import load_config  # noqa: E402


def reason_for(ppm: float, processing: dict[str, object]) -> str:
    fit = tuple(processing.get("fit_region_ppm", [0.0, 10.0]))
    water = tuple(processing.get("water_region_ppm", [4.65, 4.90]))
    artifact = tuple(processing.get("artifact_region_ppm", [5.15, 5.30]))
    if not fit[0] <= ppm <= fit[1]:
        return "outside_fit_region"
    if water[0] <= ppm <= water[1]:
        return "water"
    # Some carbohydrates have a real anomeric resonance in the generic
    # sucrose artifact window.  An empty range explicitly disables masking.
    if len(artifact) >= 2 and artifact[0] <= ppm <= artifact[1]:
        return "artifact"
    return "included"


def write_csv(path: Path, header: tuple[str, ...], rows: list[tuple[object, ...]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def dataset_nucleus(repo_root: Path, molecule: str, item: dict[str, object]) -> str:
    """Return the configured nucleus, falling back to Bruker acqus metadata."""
    configured = str(item.get("nucleus", "")).strip()
    if configured:
        return configured
    experiment = repo_root / "data" / molecule / str(item.get("relative_dir", ""))
    acqus_path = experiment / "acqus"
    if acqus_path.is_file():
        acqus = parse_jcamp(acqus_path)
        return acqus.get("NUC1", "").strip().strip("<>")
    return ""


def prepare_dataset(repo_root: Path, molecule: str, item: dict[str, object], output_dir: Path,
                    processing: dict[str, object]) -> dict[str, object]:
    experiment = repo_root / "data" / molecule / str(item["relative_dir"])
    procno = str(item.get("procno", "1"))
    processed = experiment / "pdata" / procno
    acqus_path, procs_path, spectrum_path = experiment / "acqus", processed / "procs", processed / "1r"
    for path in (acqus_path, procs_path, spectrum_path):
        if not path.is_file():
            raise FileNotFoundError(f"Required Bruker input is missing: {path}")

    acqus, procs = parse_jcamp(acqus_path), parse_jcamp(procs_path)
    nucleus = acqus.get("NUC1", "").strip().strip("<>") or "unknown"
    if nucleus.upper() != "1H":
        raise ValueError(f"{experiment} is {nucleus}, not 1H")
    si = required_number(procs, "SI", int)
    sf = required_number(procs, "SF", float)
    sw_p = required_number(procs, "SW_p", float)
    offset = required_number(procs, "OFFSET", float)
    intensity = read_processed_1r(
        spectrum_path,
        point_count=si,
        bytordp=required_number(procs, "BYTORDP", int),
        dtypp=required_number(procs, "DTYPP", int),
        nc_proc=required_number(procs, "NC_proc", int),
    )
    ppm_step = (sw_p / sf) / max(si - 1, 1)
    ppm_raw = [offset - index * ppm_step for index in range(si)]
    dss_region = tuple(processing.get("dss_search_region_ppm", [-0.2, 0.2]))
    dss_raw, dss_snr = find_dss_peak(ppm_raw, intensity, dss_region)
    ppm_ref = [value - dss_raw for value in ppm_raw]
    reasons = [reason_for(value, processing) for value in ppm_ref]
    included = [index for index, reason in enumerate(reasons) if reason == "included"]
    if not included:
        raise ValueError(f"No fit points remain for {experiment}")
    baseline = statistics.median(intensity[index] for index in included)
    corrected = [value - baseline for value in intensity]
    key = str(item["key"])
    stem = f"{key}_MHz_{item.get('acquisition', 'experiment')}"
    full_name, fit_name = f"{stem}_full.csv", f"{stem}_fit.csv"
    full_rows = [
        (f"{raw:.12g}", f"{ref:.12g}", f"{raw_y:.12g}", f"{corr:.12g}",
         int(reason == "included"), reason)
        for raw, ref, raw_y, corr, reason in zip(ppm_raw, ppm_ref, intensity, corrected, reasons)
    ]
    write_csv(
        output_dir / full_name,
        ("ppm_raw", "ppm_dss", "intensity_scaled", "intensity_baseline_corrected", "included_in_fit", "exclusion_reason"),
        full_rows,
    )
    fit_rows = [
        (f"{ppm_ref[index]:.12g}", f"{intensity[index]:.12g}", f"{corrected[index]:.12g}")
        for index in included
    ]
    write_csv(output_dir / fit_name, ("ppm_dss", "intensity_scaled", "intensity_baseline_corrected"), fit_rows)
    return {
        "key": key,
        "field_mhz": item.get("field_mhz"),
        "acquisition": item.get("acquisition"),
        "relative_dir": item.get("relative_dir"),
        "procno": procno,
        "nucleus": nucleus,
        "pulse_program": acqus.get("PULPROG", "").strip().strip("<>"),
        "concentration_mM": item.get("concentration_mM"),
        "sfo1_mhz": required_number(acqus, "SFO1", float),
        "o1_hz": required_number(acqus, "O1", float),
        "sw_acq_hz": required_number(acqus, "SW_h", float),
        "td_acquired": required_number(acqus, "TD", int),
        "ncomplex": required_number(acqus, "TD", int) // 2,
        "sf_mhz": sf,
        "sw_hz": sw_p,
        "offset_dss_ppm": offset - dss_raw,
        "bytordp": required_number(procs, "BYTORDP", int),
        "dtypp": required_number(procs, "DTYPP", int),
        "nc_proc": required_number(procs, "NC_proc", int),
        "dss_raw_ppm": dss_raw,
        "dss_snr_proxy": dss_snr,
        "points": si,
        "fit_points": len(included),
        "full_spectrum": full_name,
        "fit_spectrum": fit_name,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--dataset", action="append", help="advanced: dataset key; repeat to select fields")
    parser.add_argument("--all", action="store_true", help="prepare every imported 1H dataset separately (default)")
    parser.add_argument("--reference-only", action="store_true", help="advanced: prepare only the automatic single-field reference")
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    args = parser.parse_args()
    config = load_config(args.repo_root, args.molecule)
    processing = config.get("processing", {})
    output_dir = args.output_dir or args.repo_root / "outputs" / args.molecule / "prepared"
    if not output_dir.is_absolute():
        output_dir = args.repo_root / output_dir
    if args.all or (not args.dataset and not args.reference_only):
        selected_keys = {
            str(item["key"])
            for item in config.get("datasets", [])
            if dataset_nucleus(args.repo_root, args.molecule, item).upper() == "1H"
        }
        selection_note = "all imported 1H datasets (separate-condition diagnostics)"
        print(f"AUTO-SELECTION: {selection_note}")
    elif args.dataset:
        selected_keys = set(args.dataset)
        selection_note = "user-selected dataset"
    else:
        plan_path = args.repo_root / "outputs" / args.molecule / "multifield_plan.json"
        plan = {}
        if plan_path.is_file():
            try:
                import json
                plan = json.loads(plan_path.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                plan = {}
        selected_keys = set(str(key) for key in plan.get("training_fields", []))
        selection_note = "planned training field(s)"
        if not selected_keys:
            proton = [
                item for item in config.get("datasets", [])
                if dataset_nucleus(args.repo_root, args.molecule, item).upper() == "1H"
            ]
            if not proton:
                raise SystemExit("No 1H datasets are available for automatic preparation")
            reference = min(proton, key=lambda item: float(item.get("field_mhz", 1e9)))
            selected_keys = {str(reference["key"])}
            selection_note = "single-field reference (no matched multifield set was available)"
            print(f"AUTO-SELECTION: {reference['field_mhz']} MHz — {selection_note}")
    selected = [item for item in config.get("datasets", []) if str(item["key"]) in selected_keys]
    if not selected:
        raise SystemExit("No matching datasets were found in the molecule configuration")
    rows = []
    for item in selected:
        row = prepare_dataset(args.repo_root, args.molecule, item, output_dir, processing)
        rows.append(row)
        print(f"Prepared {row['field_mhz']} MHz/{row['nucleus']}: "
              f"DSS {row['dss_raw_ppm']:+.6f} ppm; {row['fit_points']} fit points")
    summary = output_dir / "preparation_summary.csv"
    write_csv(summary, tuple(rows[0].keys()), [tuple(row.values()) for row in rows])
    print(f"Wrote {summary}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
