#!/usr/bin/env python3
"""Compare exported TopSpin 2-D peak lists with a provisional seed.

The output is deliberately a review artifact.  It identifies candidate H-H
connectivities and H-C shift matches, but never changes a spin matrix or
turns a 2-D cross-peak into a scalar J coupling.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def read_peaklist(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not path.is_file():
        return rows
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        tokens = line.split()
        if len(tokens) < 6 or not tokens[0].lstrip("-").isdigit():
            continue
        try:
            rows.append({
                "index": int(tokens[0]),
                "f2_ppm": float(tokens[3]),
                "f1_ppm": float(tokens[4]),
                "intensity": float(tokens[5]),
                "annotation": " ".join(tokens[6:]) if len(tokens) > 6 else "",
            })
        except ValueError:
            continue
    return rows


def candidate_matches(value: float, atoms: list[dict[str, Any]], tolerance: float) -> list[dict[str, Any]]:
    return [
        {"id": str(atom["id"]), "shift_ppm": float(atom["shift_ppm"]),
         "difference_ppm": abs(float(atom["shift_ppm"]) - value)}
        for atom in atoms
        if abs(float(atom["shift_ppm"]) - value) <= tolerance
    ]


def h_h_review(rows: list[dict[str, Any]], atoms: list[dict[str, Any]], tolerance: float) -> dict[str, Any]:
    matched: list[dict[str, Any]] = []
    for row in rows:
        if abs(row["f1_ppm"] - row["f2_ppm"]) <= tolerance:
            continue
        f2 = candidate_matches(row["f2_ppm"], atoms, tolerance)
        f1 = candidate_matches(row["f1_ppm"], atoms, tolerance)
        if not f1 or not f2:
            continue
        matched.append({**row, "f2_candidates": f2, "f1_candidates": f1})
    return {
        "input_peak_count": len(rows),
        "off_diagonal_candidate_count": len(matched),
        "candidate_crosspeaks": matched,
    }


def h_c_review(rows: list[dict[str, Any]], atoms: list[dict[str, Any]], tolerance: float) -> dict[str, Any]:
    matched = []
    for row in rows:
        candidates = candidate_matches(row["f2_ppm"], atoms, tolerance)
        if candidates:
            matched.append({**row, "proton_candidates": candidates})
    return {
        "input_peak_count": len(rows),
        "proton_shift_supported_count": len(matched),
        "candidate_hc_peaks": matched,
    }


def coupling_support(couplings: list[dict[str, Any]], hh_reviews: dict[str, dict[str, Any]]) -> list[dict[str, Any]]:
    result = []
    for coupling in couplings:
        i, j = str(coupling["i"]), str(coupling["j"])
        observations = []
        for experiment, review in hh_reviews.items():
            for peak in review["candidate_crosspeaks"]:
                f2_ids = {item["id"] for item in peak["f2_candidates"]}
                f1_ids = {item["id"] for item in peak["f1_candidates"]}
                if (i in f2_ids and j in f1_ids) or (j in f2_ids and i in f1_ids):
                    observations.append({"experiment": experiment, "index": peak["index"],
                                         "f2_ppm": peak["f2_ppm"], "f1_ppm": peak["f1_ppm"]})
        result.append({
            "i": i, "j": j, "provisional_j_hz": coupling.get("j_hz"),
            "hh_crosspeak_support": observations,
            "connectivity_supported": bool(observations),
            "interpretation": "supports H-H connectivity only; provisional J remains unchanged",
        })
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--tolerance-ppm", type=float, default=0.05)
    args = parser.parse_args()
    root = args.repo_root.resolve()
    seed_path = root / "data" / args.molecule / f"{args.molecule}_provisional_observations.json"
    if not seed_path.is_file():
        raise FileNotFoundError(seed_path)
    seed = json.loads(seed_path.read_text(encoding="utf-8"))
    atoms = list(seed.get("atoms", []))
    base = root / "data" / args.molecule / "2d"
    files = {
        "cosy": base / "bmse000026_10" / "cosy_pp2dml_peaklist.txt",
        "tocsy": base / "bmse000026_4" / "tocsy_pp2d_peaklist.txt",
        "hsqc": base / "bmse000026_8" / "hsqc_pp2dml_peaklist.txt",
        "hmbc": base / "bmse000026_9" / "hmbc_pp2dml_peaklist.txt",
    }
    rows = {name: read_peaklist(path) for name, path in files.items()}
    hh = {name: h_h_review(rows[name], atoms, args.tolerance_ppm) for name in ("cosy", "tocsy")}
    hc = {name: h_c_review(rows[name], atoms, args.tolerance_ppm) for name in ("hsqc", "hmbc")}
    report = {
        "molecule": args.molecule,
        "status": "REVIEW",
        "seed_file": str(seed_path),
        "peak_list_files": {name: str(path) for name, path in files.items()},
        "tolerance_ppm": args.tolerance_ppm,
        "h_h_connectivity_review": hh,
        "h_c_shift_review": hc,
        "provisional_coupling_support": coupling_support(seed.get("couplings", []), hh),
        "interpretation": {
            "cosy_tocsy": "A matched off-diagonal peak supports a candidate H-H connectivity, subject to manual review.",
            "hsqc_hmbc": "A matched proton coordinate supports a shift assignment; carbon values remain evidence for review.",
            "matrix_policy": "No shifts, couplings, or matrix entries were modified by this report.",
        },
    }
    output = root / "outputs" / args.molecule / "2d_assignment_review.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"COSY candidate crosspeaks: {hh['cosy']['off_diagonal_candidate_count']}")
    print(f"TOCSY candidate crosspeaks: {hh['tocsy']['off_diagonal_candidate_count']}")
    print(f"HSQC proton-shift matches: {hc['hsqc']['proton_shift_supported_count']}")
    print(f"HMBC proton-shift matches: {hc['hmbc']['proton_shift_supported_count']}")
    print(f"Wrote {output}")
    print("Status: REVIEW — evidence summarized; matrix unchanged.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
