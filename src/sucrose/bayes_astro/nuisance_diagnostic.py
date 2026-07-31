#!/usr/bin/env python3
"""Test additive nuisance subtraction without changing the sucrose matrix.

The physical matrix is held fixed.  For each field this diagnostic compares a
constant-baseline fit with a weak quadratic-baseline fit, then reports the
agreement before and after subtracting the fitted nuisance baseline.  It is a
diagnostic only: no matrix, prepared spectrum, or quality-gate result is
overwritten.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[2]
sys.path.insert(0, str(HERE))
import joint_fit as jf  # noqa: E402


def correlation(a: np.ndarray, b: np.ndarray) -> float:
    if np.std(a) == 0 or np.std(b) == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def fit_with_baseline(matrix: np.ndarray, data: dict[str, np.ndarray], degree: int) -> dict[str, object]:
    """Fit linewidth, offset, scale, and a bounded low-order baseline."""
    x = np.asarray(data["ppm"], dtype=float)
    y = np.asarray(data["y"], dtype=float)
    xmid = float(np.mean(x))
    xhalf = max(float(np.ptp(x)) / 2.0, 1.0)
    xn = (x - xmid) / xhalf
    nbaseline = degree + 1 if degree >= 0 else 0

    def components(params: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        lb, offset, scale = params[:3]
        physical = scale * jf.field_model(matrix, data, lb, offset)
        nuisance = np.zeros_like(x)
        for power, coefficient in enumerate(params[3:]):
            nuisance += coefficient * xn**power
        return physical, nuisance, physical + nuisance

    # The penalty keeps a flexible polynomial from absorbing real resonances.
    # It is deliberately small compared with the spectrum residual, but makes
    # the quadratic test conservative and reproducible.
    ridge = 1.0e-3

    def residual(params: np.ndarray) -> np.ndarray:
        _, _, prediction = components(params)
        penalty = np.sqrt(ridge) * np.asarray(params[3:], dtype=float)
        return np.concatenate((prediction - y, penalty))

    x0 = np.zeros(3 + nbaseline, dtype=float)
    x0[:3] = (1.5, 0.0, 1.0)
    lower = np.concatenate(([0.3, -0.03, 0.01], np.full(nbaseline, -0.5)))
    upper = np.concatenate(([5.0, 0.03, 10.0], np.full(nbaseline, 0.5)))
    solution = least_squares(residual, x0, bounds=(lower, upper), max_nfev=1200)
    physical, nuisance, prediction = components(solution.x)
    corrected = y - nuisance
    return {
        "degree": degree,
        "parameters": solution.x,
        "physical": physical,
        "nuisance": nuisance,
        "prediction": prediction,
        "corrected": corrected,
        "r_raw": correlation(y, prediction),
        "rmse_raw": float(np.sqrt(np.mean((y - prediction) ** 2))),
        "r_corrected": correlation(corrected, physical),
        "rmse_corrected": float(np.sqrt(np.mean((corrected - physical) ** 2))),
        "nuisance_range": float(np.ptp(nuisance)),
        "optimizer_success": bool(solution.success),
        "optimizer_message": str(solution.message),
    }


def load_roles(repo_root: Path) -> dict[str, str]:
    plan_path = repo_root / "outputs" / "sucrose" / "multifield_plan.json"
    roles: dict[str, str] = {}
    if plan_path.is_file():
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        for key in plan.get("training_fields", []):
            roles[str(key)] = "training"
        for key in plan.get("validation_fields", []):
            roles[str(key)] = "validation"
    if not roles:
        for item in jf.CONFIG.get("datasets", []):
            key = str(item.get("key", ""))
            if key:
                roles[key] = str(item.get("role", "unassigned"))
    return roles


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument(
        "--candidate-matrix",
        type=Path,
        help="optional fixed refined matrix; defaults to outputs/sucrose/sucrose_matrix_fit_600_900.txt",
    )
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    output_dir = repo_root / "outputs" / "sucrose" / "nuisance_test"
    output_dir.mkdir(parents=True, exist_ok=True)

    gis_matrix = jf.MATRIX0
    candidate_path = args.candidate_matrix or repo_root / "outputs" / "sucrose" / "sucrose_matrix_fit_600_900.txt"
    if not candidate_path.is_absolute():
        candidate_path = repo_root / candidate_path
    candidate_matrix = np.loadtxt(candidate_path) if candidate_path.is_file() else gis_matrix
    candidate_label = "600/900 candidate" if candidate_path.is_file() else "GISSMO seed (candidate unavailable)"
    roles = load_roles(repo_root)
    summary = jf.load_preparation_summary(jf.SUMMARY_FILE)

    all_rows: list[dict[str, object]] = []
    plot_data: list[tuple[str, str, dict[str, object], dict[str, object]]] = []
    print("Nuisance diagnostic: fixed matrix, no matrix refinement")
    print(f"Candidate matrix: {candidate_label}")
    for field_key in sorted(summary, key=lambda value: float(value)):
        field = jf.field_from_metadata(field_key)
        data = jf.load_field(field)
        matrices = [("GISSMO", gis_matrix)]
        if candidate_path.is_file():
            matrices.append((candidate_label, candidate_matrix))
        field_results: dict[str, dict[str, object]] = {}
        for matrix_label, matrix in matrices:
            for degree in (-1, 0, 2):
                result = fit_with_baseline(matrix, data, degree=degree)
                pars = np.asarray(result["parameters"])
                degree_label = "none" if degree < 0 else str(degree)
                all_rows.append({
                    "field_mhz": field_key,
                    "role": roles.get(field_key, "unassigned"),
                    "matrix": matrix_label,
                    "baseline_degree": degree_label,
                    "r_raw": result["r_raw"],
                    "rmse_raw": result["rmse_raw"],
                    "r_corrected": result["r_corrected"],
                    "rmse_corrected": result["rmse_corrected"],
                    "lb_hz": float(pars[0]),
                    "offset_ppm": float(pars[1]),
                    "scale": float(pars[2]),
                    "nuisance_range": result["nuisance_range"],
                    "optimizer_success": result["optimizer_success"],
                })
                if degree in (-1, 2):
                    print(
                        f"{field_key} MHz ({roles.get(field_key, 'unassigned'):10s}) {matrix_label:24s} "
                        f"baseline={degree_label:4s} r={result['r_raw']:.4f} -> corrected r={result['r_corrected']:.4f}; "
                        f"RMSE={result['rmse_raw']:.4f} -> {result['rmse_corrected']:.4f}; "
                        f"nuisance range={result['nuisance_range']:.4g}"
                    )
                if degree == 2:
                    field_results[matrix_label] = result
        trace_label = candidate_label if candidate_path.is_file() else "GISSMO"
        plot_data.append((field_key, roles.get(field_key, "unassigned"), data, field_results[trace_label]))

        # Save the experimental trace and both model diagnostics for auditability.
        candidate_result = field_results[trace_label]
        path = output_dir / f"{field_key}_MHz_nuisance_trace.csv"
        rows = []
        for index, ppm in enumerate(data["ppm"]):
            rows.append({
                "ppm": float(ppm),
                "experimental": float(data["y"][index]),
                "candidate_physical": float(candidate_result["physical"][index]),
                "candidate_nuisance": float(candidate_result["nuisance"][index]),
                "candidate_corrected_experimental": float(candidate_result["corrected"][index]),
                "candidate_full_prediction": float(candidate_result["prediction"][index]),
            })
        write_csv(path, rows, list(rows[0]))

    fields = [
        "field_mhz", "role", "matrix", "baseline_degree", "r_raw", "rmse_raw",
        "r_corrected", "rmse_corrected", "lb_hz", "offset_ppm", "scale",
        "nuisance_range", "optimizer_success",
    ]
    write_csv(output_dir / "nuisance_summary.csv", all_rows, fields)
    report = {
        "status": "DIAGNOSTIC_ONLY",
        "candidate_matrix": str(candidate_path) if candidate_path.is_file() else "GISSMO seed",
        "baseline_model": "bounded quadratic in normalized ppm",
        "matrix_refinement": False,
        "official_outputs_modified": False,
        "fields": all_rows,
    }
    (output_dir / "nuisance_summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(4, 1, figsize=(13, 10), sharex=True)
        for ax, (field_key, role, data, result) in zip(axes, plot_data):
            ax.plot(data["ppm"], data["y"], color="#b4bd00", lw=0.8, label="experiment")
            ax.plot(data["ppm"], result["physical"], color="#ba0c2f", lw=0.8, ls="--", label="candidate physical")
            ax.plot(data["ppm"], result["prediction"], color="#004e60", lw=0.7, ls=":", label="candidate + nuisance")
            ax.set_ylabel(f"{field_key} MHz")
            ax.text(0.01, 0.88, role, transform=ax.transAxes, fontsize=9,
                    bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "0.5"})
            ax.grid(alpha=0.25)
        axes[0].legend(loc="upper right", fontsize=8)
        axes[-1].set_xlabel("1H chemical shift (ppm)")
        axes[-1].invert_xaxis()
        fig.suptitle("Sucrose additive-nuisance diagnostic (candidate matrix fixed)")
        fig.tight_layout()
        fig.savefig(output_dir / "nuisance_diagnostic.png", dpi=220, facecolor="white")
        plt.close(fig)
    except ImportError:
        print("matplotlib unavailable; skipped diagnostic plot")

    print(f"Wrote {output_dir / 'nuisance_summary.csv'}")
    print(f"Wrote {output_dir / 'nuisance_summary.json'}")
    print(f"Wrote {output_dir / 'nuisance_diagnostic.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
