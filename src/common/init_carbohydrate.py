#!/usr/bin/env python3
"""Create a new carbohydrate configuration from an optional seed manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def create_config(root: Path, molecule: str, *, overwrite: bool = False) -> Path:
    root = root.resolve()
    data_dir = root / "data" / molecule
    config_path = data_dir / f"{molecule}_config.json"
    if config_path.exists() and not overwrite:
        return config_path
    manifest_path = root / "outputs" / molecule / "seed_selection.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    atom_ids = manifest.get("atom_ids", [])
    matrix_file = manifest.get("matrix_file")
    blocks = [list(range(1, len(atom_ids) + 1))] if matrix_file and atom_ids else []
    config = {
        "name": molecule,
        "chemistry": {"bubb_profile": molecule},
        "seed_selection": {"default_source": "auto", "priority": ["gissmo", "bmrb", "curated_library", "bubb_rules"], "require_provenance": True, "allow_provisional": True},
        # Numeric J values are an evidence stage, not an automatic matrix
        # edit.  A new carbohydrate may use Bubb/BMRB values as a starting
        # prior, but publication requires review of independently resolved
        # splittings (or a J-resolved experiment) where available.
        "j_measurement": {
            "required_for_new_matrix": True,
            "allowed_methods": ["resolved_1h_multiplet", "j_resolved_2d"],
            "fallback": "bubb_or_bmrb_prior",
            "candidate_report": "outputs/{molecule}/j_measurements/j_candidate_measurements.json",
            "matrix_update_requires_manual_review": True
        },
        "matrix_file": matrix_file,
        "atom_ids": atom_ids,
        "blocks": blocks,
        "datasets": [],
        # Do not assume the sucrose artifact band for a new carbohydrate: it
        # may contain a real anomeric resonance.  Add an artifact range only
        # after inspecting that molecule's spectra.
        "processing": {"fit_region_ppm": [2.50, 6.20], "water_region_ppm": [4.65, 4.90], "artifact_region_ppm": [], "crowded_region_ppm": [3.40, 4.00], "anomeric_region_ppm": [4.30, 5.80], "dss_search_region_ppm": [-0.20, 0.20], "anomeric_reference_ppm": 5.40, "baseline": "median_fit_region", "normalization": "anomeric_peak"},
        "independent_model": {"grid_points": 20000, "lb_hz": 1.0, "fit_stride": 3, "noise_sigma": 0.05, "shift_bound_ppm": 0.05, "coupling_bound_hz": 3.0, "linewidth_bounds_hz": [0.3, 5.0], "offset_bound_ppm": 0.03, "initial_linewidth_hz": [1.5, 1.5], "fields_for_joint_fit": [], "validation_fields": []},
        "quality_gate": {"training_min_r": 0.90, "training_max_rmse": 0.10, "validation_min_r": 0.90, "validation_max_rmse": 0.10, "minimum_delta_r": 0.0, "minimum_delta_rmse": 0.0, "borderline_r": 0.85, "borderline_rmse": 0.15, "max_abs_offset_ppm": 0.03, "min_validation_fields": 1},
        "publication": {"target": "GISSMO", "required_seed_status": "READY", "required_quality_status": "PASS"}
    }
    data_dir.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return config_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    root = args.repo_root.resolve()
    config_path = root / "data" / args.molecule / f"{args.molecule}_config.json"
    if config_path.exists() and not args.force:
        raise SystemExit(f"Configuration already exists: {config_path} (use --force to replace it)")
    create_config(root, args.molecule, overwrite=args.force)
    manifest = json.loads((root / "outputs" / args.molecule / "seed_selection.json").read_text()) if (root / "outputs" / args.molecule / "seed_selection.json").is_file() else {}
    matrix_file = manifest.get("matrix_file")
    print(f"Created {config_path}")
    if matrix_file:
        print(f"Loaded seed matrix from {matrix_file}")
    else:
        print("No seed matrix yet; add GISSMO/BMRB/Bubb spectral evidence next.")
    print("Now import datasets with import_bruker_dataset.py.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
