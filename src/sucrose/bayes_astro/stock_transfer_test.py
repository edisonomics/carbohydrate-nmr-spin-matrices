#!/usr/bin/env python3
"""Fit the same-stock fields, then test transfer to the 1100 MHz stock.

The 600/800/900 MHz spectra are treated as the training set. The 1100 MHz
spectrum is held out because it came from a different stock solution and tube.
Only linewidth and calibration offset are optimized for the held-out spectrum;
the physical matrix is never changed while scoring it.
"""

import csv
import os
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

import joint_fit as jf

sys.path.insert(0, str(Path(jf.REPO_ROOT) / "src" / "common"))
from multifield_quality_gate import evaluate_rows


REPO_ROOT = Path(jf.REPO_ROOT)
SUMMARY_FILE = REPO_ROOT / "outputs" / "sucrose" / "prepared" / "preparation_summary.csv"
PLAN_FILE = REPO_ROOT / "outputs" / "sucrose" / "multifield_plan.json"


def summary_rows():
    with SUMMARY_FILE.open(newline="", encoding="utf-8") as handle:
        return {row["field_mhz"]: row for row in csv.DictReader(handle)}


def field_from_summary(row):
    return dict(
        SFO1=float(row["sfo1_mhz"]),
        O1=float(row["o1_hz"]),
        exp=os.path.join(jf.DATA_ROOT, row["relative_dir"], "pdata", row["procno"], "1r"),
        SI=int(row["points"]),
        SF=float(row["sf_mhz"]),
        SW_p=float(row["sw_hz"]),
        OFFSET=float(row["offset_dss_ppm"]),
        NC=int(row["nc_proc"]),
    )


def best_nuisance(matrix, data):
    lo, hi = jf.MODEL_CONFIG["linewidth_bounds_hz"]
    off_bound = float(jf.MODEL_CONFIG["offset_bound_ppm"])

    def residual(pars):
        return (jf.field_model(matrix, data, pars[0], pars[1]) - data["y"]) / data["sig"]

    sol = least_squares(
        residual,
        [float(np.mean(jf.MODEL_CONFIG["initial_linewidth_hz"])), 0.0],
        bounds=([lo, -off_bound], [hi, off_bound]),
        max_nfev=int(jf.MODEL_CONFIG.get("nuisance_max_nfev", 800)),
    )
    simulated = jf.field_model(matrix, data, sol.x[0], sol.x[1])
    good = np.isfinite(simulated) & np.isfinite(data["y"])
    return {
        "r": float(np.corrcoef(simulated[good], data["y"][good])[0, 1]),
        "rmse": float(np.sqrt(np.mean((simulated[good] - data["y"][good]) ** 2))),
        "lb_hz": float(sol.x[0]),
        "offset_ppm": float(sol.x[1]),
    }


def main():
    fields = summary_rows()
    theta0, bounds = jf.theta0_and_bounds()
    if PLAN_FILE.is_file():
        plan = __import__("json").loads(PLAN_FILE.read_text(encoding="utf-8"))
        validation_keys = [str(key) for key in plan.get("validation_fields", [])]
    else:
        validation_keys = [str(key) for key in jf.CONFIG["independent_model"]["validation_fields"]]
    print("Fitting shared matrix on clean training fields:", ", ".join(jf.FIELD_KEYS.values()))
    solution = least_squares(
        jf.residuals,
        theta0,
        bounds=bounds,
        method="trf",
        x_scale="jac",
        tr_solver="lsmr",
        verbose=2,
        max_nfev=int(jf.MODEL_CONFIG.get("transfer_max_nfev", 500)),
    )
    fitted_matrix = jf.unpack(solution.x)[0]

    results = []
    for name in jf.NAMES:
        training = best_nuisance(jf.unpack(solution.x)[0], jf.DATA[name])
        training_key = jf.FIELD_KEYS[name]
        print(f"\n===== training fit: {training_key} MHz =====")
        print("600/900 matrix      " + "  ".join(f"{k}={v:.6g}" for k, v in training.items()))
        results.append({"field_mhz": training_key, "role": "training", "model": "fit", **training})

    for validation_key in validation_keys:
        heldout = jf.load_field(field_from_summary(fields[validation_key]))
        baseline = best_nuisance(jf.MATRIX0, heldout)
        transferred = best_nuisance(fitted_matrix, heldout)
        print(f"\n===== held-out transfer: {validation_key} MHz =====")
        print("GISSMO baseline     " + "  ".join(f"{k}={v:.6g}" for k, v in baseline.items()))
        print("600/900 matrix      " + "  ".join(f"{k}={v:.6g}" for k, v in transferred.items()))
        results.extend([
            {"field_mhz": validation_key, "role": "validation", "model": "baseline", **baseline},
            {"field_mhz": validation_key, "role": "validation", "model": "fit", **transferred},
        ])
    print("\nInterpretation: improvement on 800 tests tube transfer; improvement on")
    print("1100 tests stock/tube transfer. A failure does not invalidate the matrix;")
    print("it identifies a sample-specific mismatch outside the shared physics.")

    out_dir = REPO_ROOT / "outputs" / "sucrose"
    out_dir.mkdir(parents=True, exist_ok=True)
    matrix_out = out_dir / "sucrose_matrix_fit_600_900.txt"
    np.savetxt(matrix_out, fitted_matrix, fmt="%.9f", delimiter="\t")
    report_out = out_dir / "stock_transfer_600_900_to_800_1100.csv"
    with report_out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["field_mhz", "role", "model", "r", "rmse", "lb_hz", "offset_ppm"])
        writer.writeheader()
        writer.writerows(results)
    print(f"\nWrote {matrix_out}")
    print(f"Wrote {report_out}")

    quality = evaluate_rows(results, jf.CONFIG)
    quality_out = out_dir / "multifield_quality_gate.json"
    quality_out.write_text(__import__("json").dumps(quality, indent=2) + "\n", encoding="utf-8")
    print(f"\nQUALITY GATE: {quality['status']}")
    for reason in quality["reasons"]:
        print(f"  - {reason}")
    print(f"Wrote {quality_out}")
    return 0 if quality["status"] == "PASS" else (2 if quality["status"] == "REDO" else 1)


if __name__ == "__main__":
    raise SystemExit(main())
