"""Load the repository's molecule/carbohydrate workflow configuration."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_config(repo_root: Path, name: str = "sucrose") -> dict[str, Any]:
    path = repo_root / "data" / name / f"{name}_config.json"
    if not path.is_file():
        raise FileNotFoundError(f"Missing carbohydrate configuration: {path}")
    with path.open("r", encoding="utf-8") as handle:
        config = json.load(handle)
    if config.get("name") != name:
        raise ValueError(f"{path} declares name={config.get('name')!r}, expected {name!r}")
    return config


def dataset_by_key(config: dict[str, Any], key: str) -> dict[str, Any]:
    for dataset in config.get("datasets", []):
        if str(dataset.get("key")) == str(key):
            return dataset
    raise KeyError(f"Dataset {key!r} is not listed in the carbohydrate configuration")
