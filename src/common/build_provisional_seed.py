#!/usr/bin/env python3
"""Build a provisional carbohydrate spin matrix from Bubb-guided observations.

This does not guess a structure or invent shifts.  The student supplies the
spin labels and the values extracted from the spectra (COSY/TOCSY/HSQC,
anomeric assignments, and resolved multiplet splittings).  Bubb-style
interpretation supplies the connectivity/stereochemical rationale; this tool
records that provenance and emits the numeric matrix used by the fit.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def build_matrix(observations: dict[str, Any]) -> tuple[list[str], list[list[float]]]:
    atoms = observations.get("atoms", [])
    couplings = observations.get("couplings", [])
    if not atoms:
        raise ValueError("observations.atoms is empty; provide the assigned spins")
    ids = [str(atom["id"]) for atom in atoms]
    if len(set(ids)) != len(ids):
        raise ValueError("atom ids must be unique")
    shifts = {}
    for atom in atoms:
        if "shift_ppm" not in atom:
            raise ValueError(f"missing shift_ppm for atom {atom.get('id')}")
        shifts[str(atom["id"])] = float(atom["shift_ppm"])
    index = {atom_id: i for i, atom_id in enumerate(ids)}
    matrix = [[0.0 for _ in ids] for _ in ids]
    for atom_id, shift in shifts.items():
        matrix[index[atom_id]][index[atom_id]] = shift
    for coupling in couplings:
        left, right = str(coupling["i"]), str(coupling["j"])
        if left not in index or right not in index:
            raise ValueError(f"coupling references an unknown atom: {left}, {right}")
        value = float(coupling["j_hz"])
        matrix[index[left]][index[right]] = value
        matrix[index[right]][index[left]] = value
    return ids, matrix


def write_matrix(path: Path, matrix: list[list[float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(" ".join(f"{value:.9f}" for value in row) for row in matrix) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-matrix", type=Path)
    parser.add_argument("--output-manifest", type=Path)
    args = parser.parse_args()

    observations = json.loads(args.observations.read_text(encoding="utf-8"))
    atom_ids, matrix = build_matrix(observations)
    repo_root = args.repo_root.resolve()
    matrix_path = args.output_matrix or repo_root / "data" / args.molecule / "matrix" / f"{args.molecule}_spin_matrix_provisional.txt"
    manifest_path = args.output_manifest or repo_root / "outputs" / args.molecule / "seed_selection.json"
    write_matrix(matrix_path, matrix)
    try:
        relative_matrix = str(matrix_path.resolve().relative_to(repo_root))
    except ValueError:
        relative_matrix = str(matrix_path.resolve())
    manifest = {
        "molecule": args.molecule,
        "requested_source": "bubb_rules",
        "selected_source": "bubb_rules_plus_spectra",
        "status": "PROVISIONAL_SEED",
        "confidence": "provisional",
        "matrix_file": relative_matrix,
        "atom_ids": atom_ids,
        "provenance_required": True,
        "provenance": {
            "guidance": "Bubb 2003 structural reporter and coupling interpretation",
            "spectral_observations": str(args.observations.resolve()),
            "observation_count": len(observations.get("atoms", [])),
            "coupling_count": len(observations.get("couplings", [])),
        },
        "next_step": "Measure or document independent J evidence, manually review candidate splittings, then refine against training fields and pass held-out transfer/identifiability checks before publishing as a GISSMO matrix.",
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote provisional matrix: {matrix_path}")
    print(f"Wrote seed manifest: {manifest_path}")
    print("Status: PROVISIONAL_SEED (fit allowed, publication not yet allowed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
