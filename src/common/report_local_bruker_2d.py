#!/usr/bin/env python3
"""Report locally deposited Bruker 2-D matrices and peak tables.

This is intentionally an evidence report, not an automatic assignment tool.
COSY/TOCSY matrices are retained for review; a peak table or a human-checked
cross-peak list is required before changing the spin topology or J values.
"""

from __future__ import annotations

import argparse
import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any


def read_peak_table(path: Path) -> dict[str, Any] | None:
    text = path.read_text(encoding="utf-8", errors="ignore")
    variables: list[str] = []
    for line in text.splitlines():
        if line.startswith("VARS"):
            variables = line.split()[1:]
            break
    if not variables:
        return None
    rows: list[dict[str, Any]] = []
    start = False
    for line in text.splitlines():
        if line.startswith("FORMAT"):
            start = True
            continue
        if not start or not line.strip() or line.startswith(("REMARK", "DATA", "VARS")):
            continue
        tokens = line.split()
        if len(tokens) < len(variables):
            continue
        try:
            row: dict[str, Any] = {}
            for name, token in zip(variables, tokens):
                try:
                    row[name.lower()] = float(token)
                except ValueError:
                    row[name.lower()] = token
            rows.append(row)
        except (TypeError, ValueError):
            continue
    return {"file": str(path), "variables": variables, "peak_count": len(rows), "peaks": rows}


def matrix_shape(proc: Path, proc2: Path) -> list[int] | None:
    def si(path: Path) -> int | None:
        if not path.is_file():
            return None
        match = re.search(r"^##\$SI=\s*<?([^>\r\n]+)", path.read_text(encoding="utf-8", errors="ignore"), re.M)
        try:
            return int(float(match.group(1))) if match else None
        except (TypeError, ValueError):
            return None
    direct, indirect = si(proc), si(proc2)
    return [indirect, direct] if direct and indirect else None


def read_topspin_peaklist(path: Path) -> dict[str, Any] | None:
    if not path.is_file():
        return None
    try:
        root = ET.parse(path).getroot()
    except ET.ParseError:
        return None
    peaks = []
    for item in root.findall(".//Peak2D"):
        try:
            f1, f2 = float(item.attrib["F1"]), float(item.attrib["F2"])
        except (KeyError, ValueError):
            continue
        peaks.append({
            "f1_ppm": f1,
            "f2_ppm": f2,
            "intensity": float(item.attrib.get("intensity", "nan")),
            "annotation": item.attrib.get("annotation", ""),
            "type": item.attrib.get("type", ""),
        })
    return {
        "file": str(path),
        "peak_count": len(peaks),
        "peaks": peaks,
        "off_diagonal_count": sum(abs(p["f1_ppm"] - p["f2_ppm"]) > 0.02 for p in peaks),
        "note": "Machine-picked coordinates; off-diagonal peaks support connectivity review but do not directly provide scalar J values.",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[2])
    args = parser.parse_args()
    root = args.repo_root.resolve()
    index_path = root / "data" / args.molecule / "2d" / "index.json"
    if not index_path.is_file():
        raise FileNotFoundError(f"Missing {index_path}; import the 2-D bundle first.")
    index = json.loads(index_path.read_text(encoding="utf-8"))
    experiments: list[dict[str, Any]] = []
    for item in index.get("experiments", []):
        directory = root / str(item["relative_dir"])
        proc = directory / "pdata" / "1" / "proc"
        proc2 = directory / "pdata" / "1" / "proc2"
        tables = [read_peak_table(path) for path in sorted(directory.glob("*.tab"))]
        tables = [table for table in tables if table is not None]
        # Bruker exports often include both a named peak table and a generic
        # ``test2.tab`` copy.  Keep both provenance files, but report unique
        # coordinates so the same peaks are not counted twice.
        coordinates = {
            (round(float(peak["x_ppm"]), 6), round(float(peak["y_ppm"]), 6))
            for table in tables for peak in table["peaks"]
            if "x_ppm" in peak and "y_ppm" in peak
        }
        topspin_path = directory / "pdata" / "1" / "peaklist.xml"
        if not topspin_path.is_file():
            source_name = str(item.get("source_name", ""))
            match = re.match(r"(?P<dataset>.+)_(?P<expno>\d+)$", source_name)
            if match:
                candidates = (root / "data" / args.molecule / "2d" / "topspin").glob(
                    f"*/{match.group('expno')}/pdata/1/peaklist.xml"
                )
                topspin_path = next(iter(candidates), topspin_path)
        topspin_peaklist = read_topspin_peaklist(topspin_path)
        if topspin_peaklist is not None and str(item.get("nucleus_1", "")).lower() != str(item.get("indirect_nucleus", "")).lower():
            # F1/F2 are different nuclei for HSQC/HMBC; an F1/F2 difference
            # is not a COSY/TOCSY off-diagonal connectivity measure.
            topspin_peaklist["off_diagonal_count"] = None
        experiments.append({
            "source_name": item.get("source_name"),
            "experiment_name": item.get("experiment_name"),
            "pulse_program": item.get("pulse_program"),
            "field_mhz": item.get("field_mhz"),
            "relative_dir": item.get("relative_dir"),
            "processed_matrix_shape": matrix_shape(proc, proc2),
            "has_2rr_matrix": (directory / "pdata" / "1" / "2rr").is_file(),
            "peak_tables": tables,
            "unique_peak_count": len(coordinates),
            "topspin_peaklist": topspin_peaklist,
            "interpretation": (
                "Numeric peak table available for review; assignments are not changed automatically."
                if tables else
                "Processed 2-D matrix archived; peak picking/cross-peak review is still required."
            ),
        })
    report = {
        "molecule": args.molecule,
        "source": str(index.get("source_root", "")),
        "status": "REVIEW",
        "experiments": experiments,
        "assignment_policy": {
            "cosy_tocsy": "Use checked cross-peaks to support H-H connectivity; do not convert peak coordinates directly into scalar J values.",
            "hsqc_hmbc": "Use H-C correlations to support atom/shift assignments; do not infer H-H topology from HSQC alone.",
            "matrix_update": "No spin matrix is changed by this report. A student or mentor must approve assignments and measured couplings.",
        },
    }
    output = root / "outputs" / args.molecule / "local_2d_evidence.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"Local 2-D experiments: {len(experiments)}")
    for item in experiments:
        count = int(item["unique_peak_count"])
        print(f"  {item['experiment_name']}: matrix={item['processed_matrix_shape']}, numeric peaks={count}")
    print(f"Wrote {output}")
    print("Status: REVIEW — evidence imported; no assignments or J couplings were invented.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
