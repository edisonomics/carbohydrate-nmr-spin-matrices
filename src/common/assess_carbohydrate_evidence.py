#!/usr/bin/env python3
"""Record Bubb-guided chemistry checks alongside BMRB/GISSMO provenance."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "common"))
from bubb_rules import assess_config  # noqa: E402
from carbohydrate_config import load_config  # noqa: E402


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=ROOT)
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    config = load_config(repo_root, args.molecule)
    bmrb: dict = {}
    provenance_dir = repo_root / "data" / args.molecule / "bmrb"
    if provenance_dir.is_dir():
        candidates = sorted(provenance_dir.glob("*/provenance.json"))
        if candidates:
            bmrb = json.loads(candidates[-1].read_text(encoding="utf-8"))
            bmrb["gissmo_matrix_file"] = bmrb.get("gissmo", {}).get("matrix_file")
    report = assess_config(args.molecule, config, bmrb=bmrb)
    report["bmrb_artifacts"] = bmrb
    assignment_report = repo_root / "outputs" / args.molecule / "bmrb_2d_assignment_evidence.json"
    if assignment_report.is_file():
        report["bmrb_2d_assignment_evidence"] = json.loads(
            assignment_report.read_text(encoding="utf-8")
        )
    local_2d_report = repo_root / "outputs" / args.molecule / "local_2d_evidence.json"
    if local_2d_report.is_file():
        report["local_2d_evidence"] = json.loads(
            local_2d_report.read_text(encoding="utf-8")
        )
    assignment_review = repo_root / "outputs" / args.molecule / "2d_assignment_review.json"
    if assignment_review.is_file():
        report["local_2d_assignment_review"] = json.loads(
            assignment_review.read_text(encoding="utf-8")
        )
    # J candidates are deliberately a separate evidence product.  Their
    # presence should make the review state visible without silently promoting
    # a line spacing to a signed scalar coupling or changing the matrix.
    j_report = repo_root / "outputs" / args.molecule / "j_measurements" / "j_candidate_measurements.json"
    if j_report.is_file():
        report["j_measurement_evidence"] = json.loads(
            j_report.read_text(encoding="utf-8")
        )
    else:
        report["j_measurement_evidence"] = {
            "status": "NOT_RUN",
            "required_for_new_matrix": True,
            "interpretation": "Run the independent J screen when resolved 1-D splittings or a J-resolved experiment are available; otherwise document the Bubb/BMRB prior."
        }
    output = repo_root / "outputs" / args.molecule / "chemistry_evidence.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Bubb profile: {report['profile']['bubb_profile']}")
    print(f"BMRB shifts: {'available' if report['checks']['bmrb_shifts_available'] else 'not recorded'}")
    print(f"Model check: {report['status']}")
    print(f"J measurement stage: {report['j_measurement_evidence']['status']}")
    for warning in report["warnings"]:
        print(f"REVIEW: {warning}")
    print(f"Wrote {output}")
    return 0 if report["status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
