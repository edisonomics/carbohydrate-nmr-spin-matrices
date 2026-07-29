#!/usr/bin/env python3
"""Compare baseline GISSMO and candidate Spinach multifield summaries."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_summary(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        return {str(row["field_mhz"]): row for row in rows}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path,
                        default=Path(__file__).resolve().parents[2])
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    out_dir = args.repo_root / "outputs" / args.molecule
    baseline_path = args.baseline or (out_dir / "spinach_multifield_summary.csv")
    candidate_path = args.candidate or (out_dir / "spinach_multifield_summary_candidate.csv")
    output_path = args.output or (out_dir / "spinach_candidate_vs_gissmo.csv")

    if not baseline_path.is_file():
        raise FileNotFoundError(f"Missing baseline summary: {baseline_path}")
    if not candidate_path.is_file():
        raise FileNotFoundError(
            f"Missing candidate summary: {candidate_path}. "
            "Run the candidate Spinach workflow first."
        )

    baseline = read_summary(baseline_path)
    candidate = read_summary(candidate_path)
    common = sorted(set(baseline) & set(candidate), key=lambda value: float(value))
    if not common:
        raise ValueError("The baseline and candidate summaries share no fields.")

    config_path = args.repo_root / "data" / args.molecule / f"{args.molecule}_config.json"
    config = json.loads(config_path.read_text(encoding="utf-8")) if config_path.is_file() else {}
    roles = {str(item.get("key")): str(item.get("role", "unassigned"))
             for item in config.get("datasets", [])}

    rows = []
    print("field   role        baseline_r  candidate_r  delta_r   baseline_rmse  candidate_rmse  delta_rmse  result")
    for field in common:
        b = baseline[field]
        c = candidate[field]
        br = float(b["r_spinach_vs_expt"])
        cr = float(c["r_spinach_vs_expt"])
        be = float(b["rmse_spinach_vs_expt"])
        ce = float(c["rmse_spinach_vs_expt"])
        delta_r = cr - br
        delta_rmse = ce - be
        result = "IMPROVED" if delta_r >= 0 and delta_rmse <= 0 else (
            "WORSE" if delta_r < 0 or delta_rmse > 0 else "UNCHANGED")
        role = roles.get(field, "unassigned")
        print(f"{field:>4}   {role:<10}  {br:10.4f}  {cr:11.4f}  {delta_r:+.4f}"
              f"   {be:13.4f}  {ce:14.4f}  {delta_rmse:+.4f}  {result}")
        rows.append({
            "field_mhz": field,
            "role": role,
            "baseline_r": br,
            "candidate_r": cr,
            "delta_r": delta_r,
            "baseline_rmse": be,
            "candidate_rmse": ce,
            "delta_rmse": delta_rmse,
            "result": result,
        })

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
