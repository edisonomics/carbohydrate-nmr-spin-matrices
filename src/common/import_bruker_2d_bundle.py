#!/usr/bin/env python3
"""Import a folder of Bruker 2-D experiments as assignment evidence.

Two-dimensional experiments are deliberately stored outside the 1-D dataset
configuration.  They can support atom assignments and connectivity review,
but they are not interchangeable with 1-D spectra used by the Spinach fit.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
from pathlib import Path
from typing import Any


def value(path: Path, key: str, default: str = "unknown") -> str:
    if not path.is_file():
        return default
    match = re.search(rf"^##\${re.escape(key)}=\s*<?([^>\r\n]+)",
                      path.read_text(encoding="utf-8", errors="ignore"), re.M)
    return match.group(1).strip() if match else default


def number(path: Path, key: str) -> float | None:
    raw = value(path, key, "")
    try:
        return float(raw)
    except (TypeError, ValueError):
        return None


def label(experiment: Path) -> str:
    name = (experiment / "experiment_name.txt").read_text(
        encoding="utf-8", errors="ignore").strip() if (experiment / "experiment_name.txt").is_file() else experiment.name
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_").lower() or experiment.name


def is_2d(experiment: Path) -> bool:
    pdata = experiment / "pdata"
    processed = any((p / "2rr").is_file() for p in pdata.iterdir() if p.is_dir()) if pdata.is_dir() else False
    # Older Bruker exports may provide an FT2/TAB assignment file without a
    # pdata/1/2rr matrix.  It is still valid 2-D evidence and must be retained.
    legacy = any(experiment.glob(pattern) for pattern in ("*.ft2", "*.ft", "*.tab"))
    return processed or legacy


def metadata(experiment: Path) -> dict[str, Any]:
    acqus = experiment / "acqus"
    acqu2 = experiment / "acqu2"
    proc = experiment / "pdata" / "1" / "proc"
    proc2 = experiment / "pdata" / "1" / "proc2"
    return {
        "source_name": experiment.name,
        "experiment_name": (experiment / "experiment_name.txt").read_text(encoding="utf-8", errors="ignore").strip()
            if (experiment / "experiment_name.txt").is_file() else experiment.name,
        "label": label(experiment),
        "field_mhz": number(acqus, "BF1") or number(acqus, "SFO1"),
        "nucleus_1": value(acqus, "NUC1"),
        "nucleus_2": value(acqus, "NUC2"),
        "indirect_nucleus": value(acqu2, "NUC1"),
        "pulse_program": value(acqus, "PULPROG"),
        "acquisition_points_direct": number(acqus, "TD"),
        "acquisition_points_indirect": number(acqu2, "TD"),
        "processed_points_direct": number(proc, "SI"),
        "processed_points_indirect": number(proc2, "SI"),
        "direct_offset_ppm": number(proc, "OFFSET"),
        "indirect_offset_ppm": number(proc2, "OFFSET"),
        "processed_files": sorted(p.name for p in (experiment / "pdata" / "1").iterdir())
            if (experiment / "pdata" / "1").is_dir() else [],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--replace", action="store_true")
    args = parser.parse_args()

    root = args.repo_root.resolve()
    source_root = args.source_root.resolve()
    if not source_root.is_dir():
        raise FileNotFoundError(source_root)
    candidates = sorted({p.parent for p in source_root.rglob("acqus") if is_2d(p.parent)})
    if not candidates:
        raise SystemExit(f"No processed Bruker 2-D experiments found below {source_root}")

    destination_root = root / "data" / args.molecule / "2d"
    index_path = destination_root / "index.json"
    index: dict[str, Any] = {"molecule": args.molecule, "experiments": []}
    if index_path.is_file() and not args.replace:
        index = json.loads(index_path.read_text(encoding="utf-8"))
    prior = {str(item.get("source_name")): item for item in index.get("experiments", [])}

    for source in candidates:
        item = metadata(source)
        key = str(item["source_name"])
        destination = destination_root / key
        if destination.exists() and not args.replace:
            raise FileExistsError(f"{destination} exists; use --replace")
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(source, destination)
        item["relative_dir"] = str(destination.relative_to(root))
        prior[key] = item
        print(f"Imported {item['experiment_name']} -> {destination}")

    index["experiments"] = sorted(prior.values(), key=lambda x: str(x.get("source_name")))
    index["source_root"] = str(source_root)
    destination_root.mkdir(parents=True, exist_ok=True)
    index_path.write_text(json.dumps(index, indent=2) + "\n", encoding="utf-8")
    print(f"Imported {len(candidates)} 2-D experiment(s)")
    print(f"Wrote {index_path}")
    print("These are assignment evidence, not 1-D fit fields.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
