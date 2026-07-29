#!/usr/bin/env python3
"""Import a processed Bruker 1D experiment into a carbohydrate project."""

from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path
from typing import Any


def acqus_value(path: Path, key: str, default: str = "unknown") -> str:
    pattern = re.compile(rf"^##\${re.escape(key)}=\s*<?([^>\r\n]+)", re.MULTILINE)
    match = pattern.search(path.read_text(encoding="utf-8", errors="ignore"))
    return match.group(1).strip() if match else default


def locate_source(source: Path, procno: str) -> tuple[Path, Path]:
    """Return (experiment_dir, processed_dir) from an experiment or pdata path."""
    source = source.resolve()
    if (source / "acqus").is_file():
        experiment = source
    elif (source.parent.parent / "acqus").is_file() and source.parent.name == "pdata":
        experiment = source.parent.parent
    elif (source.parent / "acqus").is_file() and source.parent.name == "pdata":
        experiment = source.parent
    else:
        raise FileNotFoundError(f"Could not find acqus above source path: {source}")
    processed = experiment / "pdata" / procno
    if not processed.is_dir() and source.is_dir() and (source / "1r").is_file():
        processed = source
    return experiment, processed


def import_dataset(args: argparse.Namespace) -> dict[str, Any]:
    repo_root = args.repo_root.resolve()
    source_experiment, source_processed = locate_source(Path(args.source), args.procno)
    required = [source_experiment / "acqus", source_processed / "1r", source_processed / "procs"]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required Bruker files:\n  " + "\n  ".join(missing))

    config_path = repo_root / "data" / args.molecule / f"{args.molecule}_config.json"
    if not config_path.is_file():
        raise FileNotFoundError(f"Missing molecule configuration: {config_path}")
    config = json.loads(config_path.read_text(encoding="utf-8"))
    field_value = float(args.field_mhz)
    key = str(int(field_value)) if field_value.is_integer() else str(field_value)
    existing = [item for item in config.get("datasets", []) if str(item.get("key")) == key]
    if existing and not args.replace:
        raise ValueError(f"Dataset key {key!r} already exists; use --replace to overwrite it")

    relative_dir = f"{key}_MHz/{args.experiment}"
    destination_experiment = repo_root / "data" / args.molecule / relative_dir
    destination_processed = destination_experiment / "pdata" / args.procno
    destination_processed.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_experiment / "acqus", destination_experiment / "acqus")
    shutil.copy2(source_processed / "1r", destination_processed / "1r")
    shutil.copy2(source_processed / "procs", destination_processed / "procs")

    dataset = {
        "key": key,
        "field_mhz": int(field_value) if field_value.is_integer() else field_value,
        "acquisition": str(args.experiment),
        "relative_dir": relative_dir,
        "procno": str(args.procno),
        "sample_id": args.sample_id,
        "tube_id": args.tube_id,
        "concentration_mM": args.concentration_mM,
        "role": args.role,
        "nucleus": acqus_value(source_experiment / "acqus", "NUC1"),
        "pulse_program": acqus_value(source_experiment / "acqus", "PULPROG"),
    }
    config["datasets"] = [item for item in config.get("datasets", []) if str(item.get("key")) != key]
    config["datasets"].append(dataset)
    config["datasets"].sort(key=lambda item: float(item["field_mhz"]))
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return {"dataset": dataset, "destination": str(destination_experiment), "config": str(config_path)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--source", type=Path, required=True, help="Bruker experiment directory or pdata/<procno> directory")
    parser.add_argument("--field-mhz", required=True, type=float)
    parser.add_argument("--experiment", required=True, help="Experiment/acquisition label, e.g. 12")
    parser.add_argument("--procno", default="1")
    parser.add_argument("--sample-id", default="unassigned_sample")
    parser.add_argument("--tube-id", default="unassigned_tube")
    parser.add_argument("--concentration-mM", type=float)
    parser.add_argument("--role", choices=("training", "validation", "unassigned"), default="unassigned")
    parser.add_argument("--replace", action="store_true", help="replace an existing dataset with the same field key")
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    result = import_dataset(args)
    print(f"Imported {result['dataset']['field_mhz']} MHz dataset")
    print(f"Copied Bruker files to {result['destination']}")
    print(f"Updated {result['config']}")
    print("Next: run plan_multifield_split.py, prepare_sucrose_spectra.py, then the unit tests.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
