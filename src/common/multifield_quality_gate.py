#!/usr/bin/env python3
"""Quality gate for multi-field carbohydrate transfer reports.

The gate is intentionally independent of Spinach or SciPy. It consumes the
standard CSV written by a molecule-specific transfer test and returns a
machine-readable PASS/REVIEW/REDO decision for teaching and automation.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


def evaluate_rows(rows: list[dict[str, Any]], config: dict[str, Any]) -> dict[str, Any]:
    gate = config["quality_gate"]
    validation = [row for row in rows if row.get("role") == "validation" and row.get("model") == "fit"]
    training = [row for row in rows if row.get("role") == "training" and row.get("model") == "fit"]
    baselines = {
        row["field_mhz"]: row for row in rows
        if row.get("role") == "validation" and row.get("model") == "baseline"
    }

    reasons: list[str] = []
    field_results: list[dict[str, Any]] = []

    if len(validation) < int(gate["min_validation_fields"]):
        reasons.append(
            f"only {len(validation)} validation field(s) present; "
            f"need {gate['min_validation_fields']}"
        )

    for row in training + validation:
        r = float(row["r"])
        rmse = float(row["rmse"])
        offset = abs(float(row["offset_ppm"]))
        is_validation = row.get("role") == "validation"
        min_r = float(gate["validation_min_r"] if is_validation else gate["training_min_r"])
        max_rmse = float(gate["validation_max_rmse"] if is_validation else gate["training_max_rmse"])
        delta_r = None
        delta_rmse = None
        if is_validation and row["field_mhz"] in baselines:
            baseline = baselines[row["field_mhz"]]
            delta_r = r - float(baseline["r"])
            delta_rmse = float(baseline["rmse"]) - rmse

        checks = {
            "r": r >= min_r,
            "rmse": rmse <= max_rmse,
            "offset": offset <= float(gate["max_abs_offset_ppm"]),
            "transfer": (not is_validation) or (
                delta_r is not None
                and delta_r >= float(gate["minimum_delta_r"])
                and delta_rmse >= float(gate["minimum_delta_rmse"])
            ),
        }
        field_pass = all(checks.values())
        field_results.append({
            "field_mhz": row["field_mhz"],
            "role": row["role"],
            "r": r,
            "rmse": rmse,
            "offset_ppm": float(row["offset_ppm"]),
            "delta_r": delta_r,
            "delta_rmse": delta_rmse,
            "checks": checks,
            "pass": field_pass,
        })
        if not field_pass:
            reasons.append(f"{row['role']} field {row['field_mhz']} failed one or more quality checks")

    hard_fail = any(not item["pass"] for item in field_results) or len(validation) < int(gate["min_validation_fields"])
    if hard_fail:
        # A borderline result is useful feedback, but a clearly poor result
        # should tell the student to reacquire/reprocess rather than tweak the
        # matrix indefinitely.
        clearly_poor = any(
            item["r"] < float(gate["borderline_r"])
            or item["rmse"] > float(gate["borderline_rmse"])
            for item in field_results
        )
        status = "REDO" if clearly_poor else "REVIEW"
    else:
        status = "PASS"

    return {
        "status": status,
        "reasons": reasons,
        "thresholds": gate,
        "fields": field_results,
    }


def evaluate_csv(report_path: Path, config_path: Path) -> dict[str, Any]:
    with report_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    with config_path.open(encoding="utf-8") as handle:
        config = json.load(handle)
    return evaluate_rows(rows, config)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--json", type=Path, help="write the machine-readable decision")
    args = parser.parse_args()
    result = evaluate_csv(args.report, args.config)
    print(f"QUALITY GATE: {result['status']}")
    for reason in result["reasons"]:
        print(f"  - {reason}")
    for field in result["fields"]:
        print(f"  {field['role']:10s} {field['field_mhz']:>4s} MHz: "
              f"r={field['r']:.4f}, RMSE={field['rmse']:.4f}, "
              f"{'PASS' if field['pass'] else 'FAIL'}")
    if args.json:
        args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(f"Wrote {args.json}")
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
