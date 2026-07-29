#!/usr/bin/env python3
"""Build and test a non-destructive HSQC-constrained matrix candidate.

This is intentionally a conservative first 2-D refinement: it moves only
diagonal proton shifts to the nearest deposited HSQC proton coordinate and
leaves all J couplings unchanged.  The original configuration is never
overwritten.  The candidate is then nuisance-fitted against every prepared
1-D field so the report answers whether the 2-D constraint helped.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "common"))
from bruker_metadata import parse_jcamp, required_number  # noqa: E402
from carbohydrate_config import load_config  # noqa: E402
from carbohydrate_model import component_specs, mixture_spectrum  # noqa: E402


def _read_hsqc_proton_targets(observations_path: Path) -> list[float]:
    observations = json.loads(observations_path.read_text(encoding="utf-8"))
    targets: list[float] = []
    for peak_list in observations.get("peak_data", []):
        if str(peak_list.get("dimensions")) != "2":
            continue
        if "HSQC" not in str(peak_list.get("experiment_name", "")).upper():
            continue
        for peak in peak_list.get("peaks", []):
            for coordinate in peak.get("coordinates", []):
                if str(coordinate.get("atom_type", "")).upper() != "H":
                    continue
                try:
                    targets.append(float(coordinate["value"]))
                except (KeyError, TypeError, ValueError):
                    continue
    # Preserve every deposited coordinate for provenance, but use unique ppm
    # values for nearest-neighbour matching of unresolved multiplets.
    return sorted({round(value, 6) for value in targets})


def _write_matrix(path: Path, matrix: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(" ".join(f"{value:.9f}" for value in row) for row in matrix) + "\n",
        encoding="utf-8",
    )


def _candidate_config(config: dict, repo_root: Path, output_dir: Path, targets: list[float]):
    candidate = copy.deepcopy(config)
    assignments = []
    for component in candidate.get("components", []):
        source_path = repo_root / str(component["matrix_file"])
        matrix = np.loadtxt(source_path, dtype=float)
        original = np.diag(matrix).copy()
        shifted = original.copy()
        for index, value in enumerate(original):
            target = min(targets, key=lambda target_value: abs(target_value - value))
            shifted[index] = target
            assignments.append({
                "component": component.get("name", "component"),
                "spin_index": index,
                "original_shift_ppm": float(value),
                "hsqc_target_ppm": float(target),
                "delta_ppm": float(target - value),
                "match_distance_ppm": float(abs(target - value)),
            })
        matrix[np.diag_indices_from(matrix)] = shifted
        output_path = output_dir / f"{component.get('name', 'component')}_hsqc_matrix.txt"
        _write_matrix(output_path, matrix)
        component["matrix_file"] = str(output_path.resolve())
        component["provenance"] = (
            "HSQC-constrained candidate; diagonal shifts moved to deposited BMRB HSQC coordinates; "
            "J couplings unchanged"
        )
    return candidate, assignments


def _fit_one(config: dict, repo_root: Path, prep_dir: Path, row: dict[str, str]) -> dict:
    with (prep_dir / row["fit_spectrum"]).open(newline="", encoding="utf-8") as handle:
        experimental = list(csv.DictReader(handle))
    ppm = np.array([float(item["ppm_dss"]) for item in experimental])
    exp = np.array([float(item["intensity_baseline_corrected"]) for item in experimental])
    exp -= np.median(exp)
    if np.max(exp) > 0:
        exp /= np.max(exp)

    dataset = next(item for item in config["datasets"] if str(item["key"]) == row["key"])
    acqus = parse_jcamp(repo_root / "data" / config["name"] / str(dataset["relative_dir"]) / "acqus")
    sfo1 = required_number(acqus, "SFO1", float)
    carrier = required_number(acqus, "O1", float) / sfo1

    def model(params: np.ndarray) -> np.ndarray:
        lb_hz, offset_ppm, scale, baseline = params
        raw = mixture_spectrum(config, repo_root, ppm, sfo1, carrier, lb_hz, offset_ppm)
        if np.max(raw) > 0:
            raw = raw / np.max(raw)
        return scale * raw + baseline

    result = minimize(
        lambda params: float(np.mean((model(params) - exp) ** 2)),
        x0=np.array([1.0, 0.0, 1.0, 0.0]),
        method="L-BFGS-B",
        bounds=[(0.1, 10.0), (-0.05, 0.05), (0.01, 10.0), (-1.0, 1.0)],
        options={"maxiter": 2000},
    )
    fitted = model(result.x)
    correlation = float(np.corrcoef(exp, fitted)[0, 1]) if np.std(exp) and np.std(fitted) else float("nan")
    return {
        "field_mhz": row["field_mhz"],
        "correlation_r": correlation,
        "rmse": float(np.sqrt(np.mean((exp - fitted) ** 2))),
        "fit_parameters": {
            "linewidth_hz": float(result.x[0]),
            "offset_ppm": float(result.x[1]),
            "scale": float(result.x[2]),
            "baseline": float(result.x[3]),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    config = load_config(repo_root, args.molecule)
    observations_path = repo_root / "data" / args.molecule / "bmrb"
    observation_files = sorted(observations_path.glob("*/spectral_observations.json"))
    if not observation_files:
        raise SystemExit("No BMRB spectral_observations.json found; run query_bmrb_entry.py first")
    targets = _read_hsqc_proton_targets(observation_files[-1])
    if not targets:
        raise SystemExit("No numeric 2-D HSQC proton coordinates were deposited for this molecule")
    if not config.get("components"):
        raise SystemExit("HSQC candidate refinement currently requires an explicit component model")

    output_dir = repo_root / "outputs" / args.molecule / "hsqc_refinement"
    output_dir.mkdir(parents=True, exist_ok=True)
    candidate, assignments = _candidate_config(config, repo_root, output_dir, targets)
    (output_dir / "hsqc_targets_ppm.json").write_text(json.dumps(targets, indent=2) + "\n", encoding="utf-8")
    (output_dir / "assignments.json").write_text(json.dumps(assignments, indent=2) + "\n", encoding="utf-8")

    prep_dir = repo_root / "outputs" / args.molecule / "prepared"
    with (prep_dir / "preparation_summary.csv").open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    baseline_rows = [_fit_one(config, repo_root, prep_dir, row) for row in rows]
    candidate_rows = [_fit_one(candidate, repo_root, prep_dir, row) for row in rows]
    comparison = []
    for baseline, refined in zip(baseline_rows, candidate_rows):
        comparison.append({
            "field_mhz": baseline["field_mhz"],
            "baseline": baseline,
            "hsqc_constrained": refined,
            "delta_r": refined["correlation_r"] - baseline["correlation_r"],
            "delta_rmse": baseline["rmse"] - refined["rmse"],
        })
    report = {
        "molecule": args.molecule,
        "source_observations": str(observation_files[-1]),
        "hsqc_target_count": len(targets),
        "constraint_policy": "nearest deposited HSQC proton coordinate; J couplings unchanged",
        "publication_ready": False,
        "comparison": comparison,
        "note": "This is a candidate diagnostic, not a final matrix. Resolve ambiguous HSQC assignments with the student's 2-D data before publication.",
    }
    report_path = output_dir / "hsqc_refinement_comparison.json"
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"HSQC proton targets: {len(targets)}")
    for item in comparison:
        print(
            f"{item['field_mhz']} MHz: baseline r={item['baseline']['correlation_r']:.4f}, "
            f"HSQC r={item['hsqc_constrained']['correlation_r']:.4f}; "
            f"delta r={item['delta_r']:+.4f}"
        )
    print(f"Wrote {report_path}")
    print("Status: REVIEW (candidate only; not publication-ready)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
