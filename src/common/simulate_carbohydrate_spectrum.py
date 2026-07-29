#!/usr/bin/env python3
"""Forward-simulate one prepared carbohydrate spectrum from its matrix."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src" / "common"))
sys.path.insert(0, str(REPO_ROOT / "src" / "sucrose" / "bayes_astro"))
from bruker_metadata import parse_jcamp, required_number  # noqa: E402
from carbohydrate_config import load_config  # noqa: E402
from carbohydrate_model import component_specs, mixture_spectrum  # noqa: E402


def simulate_row(args, config, prep_dir: Path, row: dict[str, str]):
    dataset = next(item for item in config["datasets"] if str(item["key"]) == row["key"])
    experiment = args.repo_root / "data" / args.molecule / str(dataset["relative_dir"])
    acqus = parse_jcamp(experiment / "acqus")
    sfo1 = required_number(acqus, "SFO1", float)
    carrier = required_number(acqus, "O1", float) / sfo1
    fit_path = prep_dir / row["fit_spectrum"]
    with fit_path.open(newline="", encoding="utf-8") as handle:
        experimental = list(csv.DictReader(handle))
    ppm = np.array([float(item["ppm_dss"]) for item in experimental])
    exp = np.array([float(item["intensity_baseline_corrected"]) for item in experimental])
    lb_hz = float(config.get("independent_model", {}).get("lb_hz", 1.0))
    sim = mixture_spectrum(config, args.repo_root, ppm, sfo1, carrier, lb_hz)
    exp = exp - np.median(exp)
    if np.max(exp) > 0:
        exp = exp / np.max(exp)
    if np.max(sim) > 0:
        sim = sim / np.max(sim)
    correlation = float(np.corrcoef(exp, sim)[0, 1]) if np.std(exp) and np.std(sim) else float("nan")
    rmse = float(np.sqrt(np.mean((exp - sim) ** 2)))
    output_dir = args.repo_root / "outputs" / args.molecule / "simulation"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_csv = output_dir / f"{row['key']}_MHz_matrix_simulation.csv"
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("ppm_dss", "experimental_normalized", "matrix_simulation_normalized"))
        writer.writerows((f"{x:.12g}", f"{y:.12g}", f"{z:.12g}") for x, y, z in zip(ppm, exp, sim))
    result = {
        "molecule": args.molecule, "dataset": row["key"], "field_mhz": row["field_mhz"],
        "matrix_file": str(config.get("matrix_file") or "component_matrices"),
        "model_type": config.get("model_type", "single"),
        "component_linewidths_hz": [
            {"name": spec["name"], "linewidth_hz": spec.get("linewidth_hz"),
             "provenance": spec.get("linewidth_provenance")}
            for spec in component_specs(config, args.repo_root)
            if spec.get("linewidth_hz") is not None
        ],
        "linewidth_hz": lb_hz,
        "correlation_r": correlation, "rmse": rmse, "simulation_csv": str(output_csv),
    }
    result_path = output_dir / f"{row['key']}_MHz_matrix_simulation.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"Simulated {args.molecule} {row['field_mhz']} MHz: r={correlation:.4f}, RMSE={rmse:.4f}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--dataset", help="advanced: dataset key; default simulates all prepared fields")
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    args = parser.parse_args()
    config = load_config(args.repo_root, args.molecule)
    prep_dir = args.repo_root / "outputs" / args.molecule / "prepared"
    summary_path = prep_dir / "preparation_summary.csv"
    if not summary_path.is_file():
        raise SystemExit("Prepare the spectra first with prepare_carbohydrate_spectra.py")
    with summary_path.open(newline="", encoding="utf-8") as handle:
        summaries = list(csv.DictReader(handle))
    selected = [row for row in summaries if not args.dataset or row["key"] == str(args.dataset)]
    if not selected:
        raise SystemExit("Requested dataset is not present in the preparation summary")
    for row in selected:
        simulate_row(args, config, prep_dir, row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
