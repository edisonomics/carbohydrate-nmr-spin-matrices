#!/usr/bin/env python3
"""Plan a metadata-aware multi-field matrix fit for any carbohydrate."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any


def plan(config: dict[str, Any]) -> dict[str, Any]:
    datasets = config.get("datasets", [])
    n = len(datasets)
    if n < 2:
        return {"status": "NO_MATRIX_VALIDATION", "reason": "Only one field is available.", "training_fields": [str(datasets[0]["key"])] if datasets else [], "validation_fields": [], "groups": {}}

    groups: dict[str, list[str]] = {}
    condition_groups: dict[str, list[str]] = {}
    for item in datasets:
        sample = str(item.get("sample_id", "unspecified_sample"))
        tube = str(item.get("tube_id", "unspecified_tube"))
        group = f"{sample}::{tube}"
        groups.setdefault(group, []).append(str(item["key"]))
        concentration = str(item.get("concentration_mM", "unknown"))
        nucleus = str(item.get("nucleus", "unknown")).lower()
        pulse = str(item.get("pulse_program", "unknown")).lower()
        condition = f"{sample}::{tube}::{concentration}::{nucleus}::{pulse}"
        condition_groups.setdefault(condition, []).append(str(item["key"]))

    explicit_training = [str(item["key"]) for item in datasets if item.get("role") == "training"]
    explicit_validation = [str(item["key"]) for item in datasets if item.get("role") == "validation"]
    if explicit_training and explicit_validation:
        training, validation, method = explicit_training, explicit_validation, "config_roles"
    else:
        compatible = {key: keys for key, keys in condition_groups.items() if len(keys) >= 2}
        if not compatible:
            return {
                "status": "NEEDS_ASSIGNMENT",
                "reason": "No two fields share the same sample, tube, concentration, and observed nucleus.",
                "training_fields": [],
                "validation_fields": [],
                "unassigned_fields": [str(item["key"]) for item in datasets],
                "groups": groups,
                "condition_groups": condition_groups,
                "field_count": n,
            }
        largest = max(compatible, key=lambda key: len(compatible[key]))
        compatible_fields = compatible[largest]
        training = compatible_fields[:-1]
        validation = compatible_fields[-1:]
        method = "same_condition_holdout"
        unassigned = [
            str(item["key"])
            for item in datasets
            if str(item["key"]) not in compatible_fields
        ]
    status = "READY_FOR_VALIDATION" if training and validation else "NO_MATRIX_VALIDATION"
    return {"status": status, "method": method, "training_fields": training, "validation_fields": validation, "unassigned_fields": unassigned if "unassigned" in locals() else [], "groups": groups, "condition_groups": condition_groups, "field_count": n}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--json", type=Path)
    args = parser.parse_args()
    config_path = args.repo_root / "data" / args.molecule / f"{args.molecule}_config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    result = plan(config)
    print(f"STATUS: {result['status']}")
    print(f"METHOD: {result.get('method', 'none')}")
    print("TRAINING: " + ", ".join(result["training_fields"]) or "none")
    print("VALIDATION: " + ", ".join(result["validation_fields"]) or "none")
    if result.get("unassigned_fields"):
        print("UNASSIGNED: " + ", ".join(result["unassigned_fields"]))
    if result.get("reason"):
        print("REASON: " + result["reason"])
    for group, fields in result.get("groups", {}).items():
        print(f"GROUP {group}: {', '.join(fields)}")
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(f"Wrote {args.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
