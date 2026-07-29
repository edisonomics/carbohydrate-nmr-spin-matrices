#!/usr/bin/env python3
"""Jointly refine the provisional alpha/beta xylose spin matrices.

The configured training fields share all physical parameters: 12 chemical
shifts, 12 nonzero J couplings, and one alpha fraction. Each field has its own
alpha/beta linewidths, calibration offset, scale, and baseline. Validation
fields never change the physical matrices; only their nuisance parameters are
fit before scoring baseline and transferred candidates.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


REPO_ROOT = Path(__file__).resolve().parents[2]
COMMON = REPO_ROOT / "src" / "common"
SIM = REPO_ROOT / "src" / "sucrose" / "bayes_astro"
sys.path.insert(0, str(COMMON))
sys.path.insert(0, str(SIM))

from carbohydrate_config import load_config  # noqa: E402
from carbohydrate_model import component_specs  # noqa: E402
from multifield_quality_gate import evaluate_rows  # noqa: E402
from sucrose_sim import lorentzian_spectrum, sucrose_sticks  # noqa: E402


CONFIG = load_config(REPO_ROOT, "xylose")
MODEL = CONFIG["independent_model"]
OUTPUT_DIR = REPO_ROOT / "outputs" / "xylose"
PREP_DIR = OUTPUT_DIR / "prepared"


def _summary_rows() -> dict[str, dict[str, str]]:
    path = PREP_DIR / "preparation_summary.csv"
    if not path.is_file():
        raise FileNotFoundError(
            f"Missing {path}; run prepare_carbohydrate_spectra.py first."
        )
    with path.open(newline="", encoding="utf-8") as handle:
        return {str(row["key"]): row for row in csv.DictReader(handle)}


def _load_field(row: dict[str, str]) -> dict[str, np.ndarray | float | str]:
    with (PREP_DIR / row["fit_spectrum"]).open(newline="", encoding="utf-8") as handle:
        points = list(csv.DictReader(handle))
    stride = int(MODEL.get("fit_stride", 1))
    ppm = np.asarray([float(item["ppm_dss"]) for item in points], dtype=float)[::stride]
    intensity = np.asarray(
        [float(item["intensity_baseline_corrected"]) for item in points], dtype=float
    )[::stride]
    intensity -= np.median(intensity)
    maximum = np.max(intensity)
    if not np.isfinite(maximum) or maximum <= 0:
        raise ValueError(f"Cannot normalize prepared spectrum {row['fit_spectrum']}")
    intensity /= maximum
    sfo1 = float(row["sfo1_mhz"])
    return {
        "key": str(row["key"]),
        "ppm": ppm,
        "y": intensity,
        "sfo1": sfo1,
        "carrier": float(row["o1_hz"]) / sfo1,
        "sigma": float(MODEL["noise_sigma"]),
    }


SPECS = component_specs(CONFIG, REPO_ROOT)
if [spec["name"] for spec in SPECS] != ["alpha", "beta"]:
    raise ValueError("The xylose joint fit requires configured alpha and beta components")

MATRIX0 = [np.asarray(spec["matrix"], dtype=float) for spec in SPECS]
COUPLING_INDEX = [
    [(i, j) for i in range(matrix.shape[0]) for j in range(i + 1, matrix.shape[1])
     if matrix[i, j] != 0.0]
    for matrix in MATRIX0
]
N_SPINS = [matrix.shape[0] for matrix in MATRIX0]
N_PHYSICAL = sum(N_SPINS) + sum(len(index) for index in COUPLING_INDEX) + 1


def _build_matrix(shifts: np.ndarray, couplings: np.ndarray,
                  coupling_index: list[tuple[int, int]]) -> np.ndarray:
    matrix = np.diag(shifts).astype(float)
    for (i, j), value in zip(coupling_index, couplings):
        matrix[i, j] = value
        matrix[j, i] = value
    return matrix


def _unpack_physical(theta: np.ndarray) -> tuple[list[np.ndarray], float, int]:
    cursor = 0
    matrices = []
    for nspins, index in zip(N_SPINS, COUPLING_INDEX):
        shifts = theta[cursor:cursor + nspins]
        cursor += nspins
        couplings = theta[cursor:cursor + len(index)]
        cursor += len(index)
        matrices.append(_build_matrix(shifts, couplings, index))
    alpha_fraction = float(theta[cursor])
    return matrices, alpha_fraction, cursor + 1


def _component_spectrum(matrix: np.ndarray, field: dict, linewidth: float,
                        offset_ppm: float) -> np.ndarray:
    frequencies, intensities = sucrose_sticks(
        matrix, float(field["sfo1"]), float(field["carrier"]),
        blocks=[list(range(matrix.shape[0]))],
    )
    return lorentzian_spectrum(
        field["ppm"], frequencies, intensities, linewidth,
        float(field["sfo1"]), offset_ppm,
    )


def _field_model(matrices: list[np.ndarray], alpha_fraction: float, field: dict,
                 nuisance: np.ndarray) -> np.ndarray:
    alpha_lb, beta_lb, offset_ppm, scale, baseline = nuisance
    alpha = _component_spectrum(matrices[0], field, alpha_lb, offset_ppm)
    beta = _component_spectrum(matrices[1], field, beta_lb, offset_ppm)
    spectrum = alpha_fraction * alpha + (1.0 - alpha_fraction) * beta
    maximum = np.max(spectrum)
    if not np.isfinite(maximum) or maximum <= 0:
        return np.full_like(field["ppm"], 1e3)
    return scale * (spectrum / maximum) + baseline


def _initial_physical() -> np.ndarray:
    values: list[float] = []
    for matrix, index in zip(MATRIX0, COUPLING_INDEX):
        values.extend(np.diag(matrix))
        values.extend(matrix[i, j] for i, j in index)
    values.append(float(SPECS[0]["fraction"]))
    return np.asarray(values, dtype=float)


def _physical_bounds() -> tuple[np.ndarray, np.ndarray]:
    shift_bound = float(MODEL["shift_bound_ppm"])
    coupling_bound = float(MODEL["coupling_bound_hz"])
    lower: list[float] = []
    upper: list[float] = []
    for matrix, index in zip(MATRIX0, COUPLING_INDEX):
        shifts = np.diag(matrix)
        couplings = np.asarray([matrix[i, j] for i, j in index])
        lower.extend(shifts - shift_bound)
        upper.extend(shifts + shift_bound)
        lower.extend(couplings - coupling_bound)
        upper.extend(couplings + coupling_bound)
    fraction_bounds = MODEL.get("fraction_bounds", [0.05, 0.95])
    lower.append(float(fraction_bounds[0]))
    upper.append(float(fraction_bounds[1]))
    return np.asarray(lower), np.asarray(upper)


def _initial_nuisance() -> np.ndarray:
    return np.asarray([
        float(SPECS[0].get("linewidth_hz") or 1.5),
        float(SPECS[1].get("linewidth_hz") or 1.5),
        0.0, 1.0, 0.0,
    ])


def _nuisance_bounds() -> tuple[np.ndarray, np.ndarray]:
    linewidth = MODEL["linewidth_bounds_hz"]
    offset = float(MODEL["offset_bound_ppm"])
    return (
        np.asarray([linewidth[0], linewidth[0], -offset, 0.25, -0.10]),
        np.asarray([linewidth[1], linewidth[1], offset, 2.00, 0.10]),
    )


def _fit_nuisance(matrices: list[np.ndarray], alpha_fraction: float,
                  field: dict) -> tuple[np.ndarray, np.ndarray, float, float]:
    lower, upper = _nuisance_bounds()

    def residual(nuisance: np.ndarray) -> np.ndarray:
        return (_field_model(matrices, alpha_fraction, field, nuisance) - field["y"]) / float(field["sigma"])

    solution = least_squares(
        residual, _initial_nuisance(), bounds=(lower, upper), method="trf",
        max_nfev=int(MODEL.get("nuisance_max_nfev", 800)),
    )
    fitted = _field_model(matrices, alpha_fraction, field, solution.x)
    good = np.isfinite(fitted) & np.isfinite(field["y"])
    correlation = float(np.corrcoef(fitted[good], field["y"][good])[0, 1])
    rmse = float(np.sqrt(np.mean((fitted[good] - field["y"][good]) ** 2)))
    return solution.x, fitted, correlation, rmse


def _write_overlay(field: dict, baseline: np.ndarray, fitted: np.ndarray) -> Path:
    path = OUTPUT_DIR / f"xylose_{field['key']}_baseline_vs_candidate.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("ppm_dss", "experimental_normalized", "baseline", "candidate"))
        writer.writerows(zip(field["ppm"], field["y"], baseline, fitted))
    return path


def main() -> int:
    summaries = _summary_rows()
    training_keys = [str(key) for key in MODEL["fields_for_joint_fit"]]
    validation_keys = [str(key) for key in MODEL["validation_fields"]]
    fields = {key: _load_field(summaries[key]) for key in training_keys + validation_keys}

    physical0 = _initial_physical()
    physical_lower, physical_upper = _physical_bounds()
    nuisance0 = _initial_nuisance()
    nuisance_lower, nuisance_upper = _nuisance_bounds()
    theta0 = np.concatenate([physical0] + [nuisance0.copy() for _ in training_keys])
    lower = np.concatenate([physical_lower] + [nuisance_lower.copy() for _ in training_keys])
    upper = np.concatenate([physical_upper] + [nuisance_upper.copy() for _ in training_keys])

    def residual(theta: np.ndarray) -> np.ndarray:
        matrices, alpha_fraction, cursor = _unpack_physical(theta)
        pieces = []
        for index, key in enumerate(training_keys):
            nuisance = theta[cursor + 5 * index:cursor + 5 * (index + 1)]
            model = _field_model(matrices, alpha_fraction, fields[key], nuisance)
            pieces.append((model - fields[key]["y"]) / float(fields[key]["sigma"]))
        return np.concatenate(pieces)

    print(
        f"Fitting {N_PHYSICAL} shared physical parameters on "
        f"{', '.join(training_keys)} MHz ({len(theta0)} total parameters)."
    )
    solution = least_squares(
        residual, theta0, bounds=(lower, upper), method="trf", x_scale="jac",
        tr_solver="lsmr", max_nfev=int(MODEL.get("max_nfev", 1000)), verbose=2,
    )
    matrices_fit, fraction_fit, _ = _unpack_physical(solution.x)
    matrices_seed, fraction_seed, _ = _unpack_physical(physical0)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    alpha_out = OUTPUT_DIR / "xylose_alpha_matrix_fit_600_900.txt"
    beta_out = OUTPUT_DIR / "xylose_beta_matrix_fit_600_900.txt"
    np.savetxt(alpha_out, matrices_fit[0], fmt="%.9f")
    np.savetxt(beta_out, matrices_fit[1], fmt="%.9f")

    report_rows: list[dict[str, str | float]] = []
    field_details = []
    for key in training_keys + validation_keys:
        role = "training" if key in training_keys else "validation"
        seed_nuisance, seed_curve, seed_r, seed_rmse = _fit_nuisance(
            matrices_seed, fraction_seed, fields[key]
        )
        fit_nuisance, fit_curve, fit_r, fit_rmse = _fit_nuisance(
            matrices_fit, fraction_fit, fields[key]
        )
        overlay = _write_overlay(fields[key], seed_curve, fit_curve)
        for model_name, nuisance, r_value, rmse_value in (
            ("baseline", seed_nuisance, seed_r, seed_rmse),
            ("fit", fit_nuisance, fit_r, fit_rmse),
        ):
            report_rows.append({
                "field_mhz": key,
                "role": role,
                "model": model_name,
                "r": r_value,
                "rmse": rmse_value,
                "alpha_lb_hz": float(nuisance[0]),
                "beta_lb_hz": float(nuisance[1]),
                "offset_ppm": float(nuisance[2]),
                "scale": float(nuisance[3]),
                "baseline": float(nuisance[4]),
            })
        field_details.append({
            "field_mhz": key,
            "role": role,
            "baseline": {"r": seed_r, "rmse": seed_rmse},
            "candidate": {"r": fit_r, "rmse": fit_rmse},
            "overlay_csv": str(overlay),
        })
        print(
            f"{key} MHz {role}: r {seed_r:.4f} -> {fit_r:.4f}; "
            f"RMSE {seed_rmse:.4f} -> {fit_rmse:.4f}"
        )

    report_path = OUTPUT_DIR / "xylose_transfer_600_900_to_1100.csv"
    with report_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(report_rows[0]))
        writer.writeheader()
        writer.writerows(report_rows)

    quality = evaluate_rows(report_rows, CONFIG)
    quality_path = OUTPUT_DIR / "multifield_quality_gate.json"
    quality_path.write_text(json.dumps(quality, indent=2) + "\n", encoding="utf-8")

    changes = []
    for name, seed, fitted, index in zip(
        ("alpha", "beta"), matrices_seed, matrices_fit, COUPLING_INDEX
    ):
        changes.append({
            "component": name,
            "shift_changes": [
                {"spin_index": i + 1, "seed_ppm": float(seed[i, i]),
                 "fitted_ppm": float(fitted[i, i]),
                 "delta_ppm": float(fitted[i, i] - seed[i, i])}
                for i in range(seed.shape[0])
            ],
            "coupling_changes": [
                {"spin_i": i + 1, "spin_j": j + 1, "seed_hz": float(seed[i, j]),
                 "fitted_hz": float(fitted[i, j]),
                 "delta_hz": float(fitted[i, j] - seed[i, j])}
                for i, j in index
            ],
        })

    summary = {
        "optimizer": {
            "success": bool(solution.success),
            "message": str(solution.message),
            "nfev": int(solution.nfev),
            "cost": float(solution.cost),
        },
        "training_fields": training_keys,
        "validation_fields": validation_keys,
        "shared_physical_parameter_count": N_PHYSICAL,
        "total_training_parameter_count": len(theta0),
        "alpha_fraction": {"seed": fraction_seed, "fitted": fraction_fit},
        "matrix_files": {"alpha": str(alpha_out), "beta": str(beta_out)},
        "parameter_changes": changes,
        "fields": field_details,
        "quality_gate": quality,
    }
    summary_path = OUTPUT_DIR / "xylose_multifield_fit_600_900.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"QUALITY GATE: {quality['status']}")
    print(f"Fitted alpha fraction: {fraction_seed:.4f} -> {fraction_fit:.4f}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {report_path}")
    return 0 if quality["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
