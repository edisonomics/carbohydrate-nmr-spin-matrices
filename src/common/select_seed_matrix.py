#!/usr/bin/env python3
"""Select and document a carbohydrate spin-matrix seed.

The endpoint is a publishable GISSMO matrix.  Bubb-rule selection is therefore
the start of a provisional-seed workflow, not a dead end; use
``build_provisional_seed.py`` to combine Bubb-guided connectivity with values
measured from COSY/TOCSY/HSQC and resolved multiplets.
"""

from __future__ import annotations

import argparse
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SOURCES = ("auto", "gissmo", "bmrb", "curated_library", "bubb_rules", "ensemble")


def load_config(repo_root: Path, molecule: str) -> dict[str, Any]:
    path = repo_root / "data" / molecule / f"{molecule}_config.json"
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def default_candidates(repo_root: Path, molecule: str, config: dict[str, Any]) -> dict[str, Path | None]:
    matrix = config.get("matrix_file")
    gissmo = repo_root / matrix if matrix else None
    library = repo_root / "data" / molecule / "seeds" / "curated_library.json"
    return {
        "gissmo": gissmo if gissmo and gissmo.is_file() else None,
        "bmrb": None,
        "curated_library": library if library.is_file() else None,
        "bubb_rules": None,
    }


def choose_source(requested: str, candidates: dict[str, Path | None], priority: list[str]) -> tuple[str, Path | None]:
    if requested != "auto":
        return requested, candidates.get(requested)
    for source in priority:
        if source == "bubb_rules":
            return source, None
        if candidates.get(source) is not None:
            return source, candidates[source]
    return "bubb_rules", None


def make_manifest(molecule: str, source: str, path: Path | None, requested: str) -> dict[str, Any]:
    verified = source == "gissmo" and path is not None
    provisional = source in {"bmrb", "curated_library", "bubb_rules", "ensemble"}
    if verified:
        status = "READY"
        next_step = "Use this matrix as the physical seed and run multifield transfer validation."
    elif source in {"bmrb", "curated_library"} and path is not None:
        status = "PROVISIONAL"
        next_step = "Construct/verify the numeric matrix from this source before fitting."
    else:
        status = "NEEDS_SPECTRAL_ASSIGNMENT"
        next_step = "Use build_provisional_seed.py with Bubb-guided spectral observations, run the independent J-measurement stage, then refine and publish a GISSMO matrix."
    return {
        "molecule": molecule,
        "requested_source": requested,
        "selected_source": source,
        "status": status,
        "confidence": "verified" if verified else ("provisional" if provisional else "unknown"),
        "matrix_file": str(path) if path is not None else None,
        "provenance_required": True,
        "next_step": next_step,
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--seed-source", choices=SOURCES, default="auto")
    parser.add_argument("--bmrb", type=Path, help="BMRB-derived numeric seed or source manifest")
    parser.add_argument("--curated-library", type=Path)
    parser.add_argument("--interactive", action="store_true", help="ask for a source before selecting")
    parser.add_argument("--yes", action="store_true", help="accept automatic selection without prompting")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--copy-matrix", type=Path, help="copy the selected numeric matrix to this path")
    args = parser.parse_args()

    repo_root = args.repo_root.resolve()
    config = load_config(repo_root, args.molecule)
    candidates = default_candidates(repo_root, args.molecule, config)
    if args.bmrb:
        candidates["bmrb"] = args.bmrb.resolve()
    if args.curated_library:
        candidates["curated_library"] = args.curated_library.resolve()

    requested = args.seed_source
    if args.interactive and not args.yes:
        print("Seed source options: " + ", ".join(SOURCES))
        answer = input(f"Choose seed source [{requested}]: ").strip().lower()
        if answer:
            if answer not in SOURCES:
                raise SystemExit(f"Unknown seed source {answer!r}")
            requested = answer

    priority = config.get("seed_selection", {}).get(
        "priority", ["gissmo", "bmrb", "curated_library", "bubb_rules"]
    )
    source, path = choose_source(requested, candidates, priority)
    manifest = make_manifest(args.molecule, source, path, requested)
    if path is not None:
        try:
            manifest["matrix_file"] = str(path.resolve().relative_to(repo_root))
        except ValueError:
            manifest["matrix_file"] = str(path.resolve())

    output = args.output or repo_root / "outputs" / args.molecule / "seed_selection.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    if args.copy_matrix:
        if path is None or not path.is_file():
            raise SystemExit("The selected source has no numeric matrix to copy.")
        args.copy_matrix.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, args.copy_matrix)
        try:
            manifest["matrix_file"] = str(args.copy_matrix.resolve().relative_to(repo_root))
        except ValueError:
            manifest["matrix_file"] = str(args.copy_matrix.resolve())
        output.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print(f"Seed source selected: {source}")
    print(f"Status: {manifest['status']}")
    if path is not None:
        print(f"Source: {path}")
    print(f"Next step: {manifest['next_step']}")
    print(f"Wrote {output}")
    return 0 if manifest["status"] == "READY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
