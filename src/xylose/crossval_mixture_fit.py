#!/usr/bin/env python3
"""Three-fold multifield cross-validation for the xylose mixture model.

Each fold refines the shared alpha/beta matrices and alpha fraction on two
fields. The held-out field receives only nuisance fitting before the provisional
seed and transferred candidate are compared.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

import joint_mixture_fit as jm


FOLDS = [
    (("600", "900"), "1100"),
    (("600", "1100"), "900"),
    (("900", "1100"), "600"),
]
OUTPUT_DIR = jm.OUTPUT_DIR / "cross_validation"


def _parameter_names() -> list[str]:
    names: list[str] = []
    for spec, index in zip(jm.CONFIG["components"], jm.COUPLING_INDEX):
        component = str(spec["name"])
        atoms = [str(atom) for atom in spec["atom_ids"]]
        names.extend(f"{component}.shift.{atom}" for atom in atoms)
        names.extend(
            f"{component}.J.{atoms[i]}__{atoms[j]}" for i, j in index
        )
    names.append("alpha_fraction")
    return names


def _fit_fold(training_keys: tuple[str, str], fields: dict[str, dict]):
    physical0 = jm._initial_physical()
    physical_lower, physical_upper = jm._physical_bounds()
    nuisance0 = jm._initial_nuisance()
    nuisance_lower, nuisance_upper = jm._nuisance_bounds()
    theta0 = np.concatenate([physical0] + [nuisance0.copy() for _ in training_keys])
    lower = np.concatenate([physical_lower] + [nuisance_lower.copy() for _ in training_keys])
    upper = np.concatenate([physical_upper] + [nuisance_upper.copy() for _ in training_keys])

    def residual(theta: np.ndarray) -> np.ndarray:
        matrices, alpha_fraction, cursor = jm._unpack_physical(theta)
        pieces = []
        for index, key in enumerate(training_keys):
            nuisance = theta[cursor + 5 * index:cursor + 5 * (index + 1)]
            model = jm._field_model(matrices, alpha_fraction, fields[key], nuisance)
            pieces.append(
                (model - fields[key]["y"]) / float(fields[key]["sigma"])
            )
        return np.concatenate(pieces)

    solution = least_squares(
        residual, theta0, bounds=(lower, upper), method="trf", x_scale="jac",
        tr_solver="lsmr", max_nfev=int(jm.MODEL.get("max_nfev", 1000)),
        verbose=0,
    )
    matrices, alpha_fraction, _ = jm._unpack_physical(solution.x)
    return solution, matrices, alpha_fraction, solution.x[:jm.N_PHYSICAL]


def _score(matrices: list[np.ndarray], alpha_fraction: float, field: dict):
    nuisance, curve, correlation, rmse = jm._fit_nuisance(
        matrices, alpha_fraction, field
    )
    return {
        "r": correlation,
        "rmse": rmse,
        "nuisance": {
            "alpha_lb_hz": float(nuisance[0]),
            "beta_lb_hz": float(nuisance[1]),
            "offset_ppm": float(nuisance[2]),
            "scale": float(nuisance[3]),
            "baseline": float(nuisance[4]),
        },
        "curve": curve,
    }


def main() -> int:
    summaries = jm._summary_rows()
    fields = {key: jm._load_field(summaries[key]) for key in ("600", "900", "1100")}
    seed_vector = jm._initial_physical()
    seed_matrices, seed_fraction, _ = jm._unpack_physical(seed_vector)
    names = _parameter_names()
    if len(names) != jm.N_PHYSICAL:
        raise ValueError("Parameter labels do not match the physical vector")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fold_results = []
    fold_vectors = []
    table_rows = []

    for training, heldout in FOLDS:
        print(f"Training on {'+'.join(training)} MHz; holding out {heldout} MHz", flush=True)
        solution, matrices, alpha_fraction, physical = _fit_fold(training, fields)
        fold_vectors.append(physical)
        baseline = _score(seed_matrices, seed_fraction, fields[heldout])
        candidate = _score(matrices, alpha_fraction, fields[heldout])
        delta_r = candidate["r"] - baseline["r"]
        delta_rmse = baseline["rmse"] - candidate["rmse"]
        generalizes = delta_r >= 0.0 and delta_rmse >= 0.0
        verdict = "IMPROVED" if generalizes else "NO_TRANSFER_GAIN"

        fold_dir = OUTPUT_DIR / f"holdout_{heldout}"
        fold_dir.mkdir(parents=True, exist_ok=True)
        alpha_path = fold_dir / "xylose_alpha_matrix.txt"
        beta_path = fold_dir / "xylose_beta_matrix.txt"
        np.savetxt(alpha_path, matrices[0], fmt="%.9f")
        np.savetxt(beta_path, matrices[1], fmt="%.9f")
        overlay_path = fold_dir / f"xylose_{heldout}_blind_prediction.csv"
        with overlay_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(("ppm_dss", "experimental_normalized", "seed", "candidate"))
            writer.writerows(zip(
                fields[heldout]["ppm"], fields[heldout]["y"],
                baseline["curve"], candidate["curve"],
            ))

        row = {
            "training_fields": "+".join(training),
            "heldout_field": heldout,
            "baseline_r": baseline["r"],
            "candidate_r": candidate["r"],
            "delta_r": delta_r,
            "baseline_rmse": baseline["rmse"],
            "candidate_rmse": candidate["rmse"],
            "delta_rmse": delta_rmse,
            "alpha_fraction": alpha_fraction,
            "optimizer_success": bool(solution.success),
            "nfev": int(solution.nfev),
            "verdict": verdict,
        }
        table_rows.append(row)
        fold_results.append({
            **row,
            "optimizer_message": str(solution.message),
            "optimizer_cost": float(solution.cost),
            "baseline_nuisance": baseline["nuisance"],
            "candidate_nuisance": candidate["nuisance"],
            "matrix_files": {"alpha": str(alpha_path), "beta": str(beta_path)},
            "overlay_csv": str(overlay_path),
            "physical_parameters": {
                name: float(value) for name, value in zip(names, physical)
            },
        })
        print(
            f"  blind {heldout}: r {baseline['r']:.4f} -> {candidate['r']:.4f} "
            f"({delta_r:+.4f}); RMSE {baseline['rmse']:.4f} -> "
            f"{candidate['rmse']:.4f} ({delta_rmse:+.4f} improvement); {verdict}",
            flush=True,
        )

    metrics_path = OUTPUT_DIR / "three_fold_metrics.csv"
    with metrics_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table_rows[0]))
        writer.writeheader()
        writer.writerows(table_rows)

    vector_array = np.asarray(fold_vectors)
    stability_rows = []
    for index, name in enumerate(names):
        values = vector_array[:, index]
        stability_rows.append({
            "parameter": name,
            "seed": float(seed_vector[index]),
            "fold_holdout_1100": float(values[0]),
            "fold_holdout_900": float(values[1]),
            "fold_holdout_600": float(values[2]),
            "mean": float(np.mean(values)),
            "std": float(np.std(values, ddof=1)),
            "range": float(np.ptp(values)),
        })
    stability_path = OUTPUT_DIR / "parameter_stability.csv"
    with stability_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(stability_rows[0]))
        writer.writeheader()
        writer.writerows(stability_rows)

    passed = sum(row["verdict"] == "IMPROVED" for row in table_rows)
    status = "PASS" if passed == len(FOLDS) else "REVIEW"
    summary = {
        "status": status,
        "criterion": "candidate must increase r and decrease RMSE on every blind field",
        "passed_folds": passed,
        "total_folds": len(FOLDS),
        "folds": fold_results,
        "parameter_stability_csv": str(stability_path),
        "metrics_csv": str(metrics_path),
    }
    summary_path = OUTPUT_DIR / "three_fold_cross_validation.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"THREE-FOLD STATUS: {status} ({passed}/{len(FOLDS)} blind fields improved)")
    print(f"Wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
