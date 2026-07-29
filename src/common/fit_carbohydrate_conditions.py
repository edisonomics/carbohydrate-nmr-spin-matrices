#!/usr/bin/env python3
"""Fit every prepared carbohydrate field independently and summarize results."""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    prep_dir = args.repo_root / "outputs" / args.molecule / "prepared"
    summary_path = prep_dir / "preparation_summary.csv"
    if not summary_path.is_file():
        raise SystemExit("Prepare the spectra first with prepare_carbohydrate_spectra.py")
    with summary_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    fit_script = Path(__file__).with_name("fit_carbohydrate_spectrum.py")
    results = []
    for row in rows:
        command = [
            sys.executable, str(fit_script),
            "--molecule", args.molecule,
            "--dataset", str(row["key"]),
            "--repo-root", str(args.repo_root),
        ]
        subprocess.run(command, check=True)
        result_path = args.repo_root / "outputs" / args.molecule / "fit" / f"{row['key']}_MHz_nuisance_fit.json"
        result = json.loads(result_path.read_text(encoding="utf-8"))
        results.append(result)
        print(f"{result['field_mhz']} MHz: r={result['correlation_r']:.4f}, RMSE={result['rmse']:.4f}")
    output = args.repo_root / "outputs" / args.molecule / "fit" / "condition_summary.json"
    output.write_text(json.dumps({"molecule": args.molecule, "fields": results}, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output}")
    print("These are independent condition fits; no shared multifield fit was performed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
