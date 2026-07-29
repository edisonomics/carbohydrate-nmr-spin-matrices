"""Resolve the selected matrix seed for downstream runners."""

from __future__ import annotations

import json
from pathlib import Path


def resolve_seed_matrix(repo_root: Path, molecule: str, config: dict) -> Path:
    manifest_path = repo_root / "outputs" / molecule / "seed_selection.json"
    if manifest_path.is_file():
        with manifest_path.open(encoding="utf-8") as handle:
            manifest = json.load(handle)
        status = manifest.get("status")
        allowed = {"READY"}
        if config.get("seed_selection", {}).get("allow_provisional", False):
            allowed.add("PROVISIONAL_SEED")
        if status not in allowed:
            raise RuntimeError(
                f"Seed selection is {status!r}; provide a numeric matrix and "
                f"spectral provenance before fitting. See {manifest_path}."
            )
        raw_path = manifest.get("matrix_file")
        if not raw_path:
            raise RuntimeError(f"READY seed manifest has no matrix_file: {manifest_path}")
        path = Path(raw_path)
        if not path.is_absolute():
            path = repo_root / path
    else:
        raw_path = config.get("matrix_file")
        if not raw_path:
            raise RuntimeError(f"No matrix_file is configured for {molecule}")
        path = repo_root / raw_path
    if not path.is_file():
        raise FileNotFoundError(f"Selected matrix seed does not exist: {path}")
    return path.resolve()
