#!/usr/bin/env python3
"""Summarize deposited 2-D assignment evidence without inventing couplings.

The report distinguishes an experiment being listed in a BMRB entry from a
numeric, assigned peak table being deposited.  Proton coordinates from HSQC
can support shift/one-bond H-C assignments; they do not prove H-H couplings or
COSY/TOCSY connectivity.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def _nearest(values: list[float], target: float) -> tuple[float | None, float | None]:
    if not values:
        return None, None
    value = min(values, key=lambda item: abs(item - target))
    return value, abs(value - target)


def _load_candidate_shifts(path: Path) -> list[dict[str, Any]]:
    if not path.is_file():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    return [dict(item) for item in data.get("atoms", [])]


def build_report(observations: dict[str, Any], candidates: list[dict[str, Any]],
                 source: str, tolerance_ppm: float) -> dict[str, Any]:
    peak_lists = observations.get("peak_data", [])
    two_d = [item for item in peak_lists if str(item.get("dimensions")) == "2"]
    numeric_two_d = [item for item in two_d if int(item.get("numeric_peak_count", 0)) > 0]
    inventory_names = [
        str(item.get("name", ""))
        for item in observations.get("experiments", [])
        if "2D" in str(item.get("name", "")).upper()
    ]
    names = inventory_names or [str(item.get("experiment_name", "")) for item in two_d]
    numeric_names = [str(item.get("experiment_name", "")) for item in numeric_two_d]
    missing_numeric = [name for name in names if name not in numeric_names]

    hsqc_rows: list[dict[str, Any]] = []
    hsqc_proton_values: list[float] = []
    for peak_list in numeric_two_d:
        dimensions = peak_list.get("dimension_metadata", [])
        is_hsqc = "HSQC" in str(peak_list.get("experiment_name", "")).upper()
        if not is_hsqc:
            continue
        h_dim_ids = {
            str(dim.get("id")) for dim in dimensions
            if str(dim.get("atom_type", "")).upper() == "H"
        }
        for peak in peak_list.get("peaks", []):
            proton = None
            carbon = None
            assignments: list[str] = []
            for coordinate in peak.get("coordinates", []):
                if str(coordinate.get("dimension_id")) in h_dim_ids:
                    proton = float(coordinate["value"])
                elif str(coordinate.get("atom_type", "")).upper() == "C":
                    carbon = float(coordinate["value"])
            for group in peak.get("assigned_atoms", []):
                if str(group.get("dimension_id")) in h_dim_ids:
                    assignments.extend(str(atom) for atom in group.get("atom_ids", []))
            if proton is None:
                continue
            hsqc_proton_values.append(proton)
            hsqc_rows.append({
                "peak_id": str(peak.get("id")),
                "proton_ppm": proton,
                "carbon_ppm": carbon,
                "assigned_proton_atoms": sorted(set(assignments)),
            })

    candidate_matches: list[dict[str, Any]] = []
    for candidate in candidates:
        if "shift_ppm" not in candidate:
            continue
        target = float(candidate["shift_ppm"])
        nearest, distance = _nearest(hsqc_proton_values, target)
        candidate_matches.append({
            "candidate_atom": str(candidate.get("id", "")),
            "candidate_shift_ppm": target,
            "nearest_hsqc_proton_ppm": nearest,
            "difference_ppm": distance,
            "within_diagnostic_tolerance": bool(distance is not None and distance <= tolerance_ppm),
            "interpretation": "shift-supported only; alpha/beta identity is not confirmed by this match",
        })

    return {
        "source": source,
        "status": "REVIEW",
        "interpretation": {
            "numeric_hsqc_support": bool(hsqc_rows),
            "cosy_tocsy_numeric_tables_available": any(
                "COSY" in name.upper() or "TOCSY" in name.upper()
                for name in numeric_names
            ),
            "h_h_couplings_confirmed": False,
            "alpha_beta_assignments_confirmed": False,
            "reason": "HSQC supports proton/carbon shift correlations; COSY/TOCSY numeric coordinates or assigned cross-peaks are required to confirm H-H connectivity and couplings.",
        },
        "two_d_experiments_listed": names,
        "numeric_two_d_experiments": numeric_names,
        "two_d_experiments_without_numeric_peak_table": missing_numeric,
        "hsqc_peak_count": len(hsqc_rows),
        "hsqc_peaks": hsqc_rows,
        "candidate_shift_matches": candidate_matches,
        "diagnostic_tolerance_ppm": tolerance_ppm,
        "candidate_couplings_are_not_confirmed": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--entry", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--tolerance-ppm", type=float, default=0.02)
    args = parser.parse_args()

    repo = args.repo_root.resolve()
    entry_dir = repo / "data" / args.molecule / "bmrb" / args.entry
    observation_path = entry_dir / "spectral_observations.json"
    if not observation_path.is_file():
        raise FileNotFoundError(f"Missing {observation_path}; query the BMRB entry first.")
    observations = json.loads(observation_path.read_text(encoding="utf-8"))
    candidate_path = repo / "data" / args.molecule / f"{args.molecule}_provisional_observations.json"
    report = build_report(
        observations,
        _load_candidate_shifts(candidate_path),
        source=f"BMRB {args.entry}",
        tolerance_ppm=float(args.tolerance_ppm),
    )
    report["observation_file"] = str(observation_path)
    report["candidate_observation_file"] = str(candidate_path) if candidate_path.is_file() else None
    output = repo / "outputs" / args.molecule / "bmrb_2d_assignment_evidence.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"BMRB 2-D experiments listed: {len(report['two_d_experiments_listed'])}")
    print(f"Numeric 2-D experiments: {', '.join(report['numeric_two_d_experiments']) or 'none'}")
    print(f"Numeric HSQC peaks: {report['hsqc_peak_count']}")
    print("H-H couplings confirmed: no")
    print("Alpha/beta assignments confirmed: no")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
