#!/usr/bin/env python3
"""Extract reusable acquisition metadata from a processed Bruker 1D dataset.

The acquisition numbers come from ``acqus`` and ``pdata/<procno>/procs``.
The DSS correction is measured from the processed ``1r`` spectrum, rather
than typed into a MATLAB driver.  The module uses only the Python standard
library so it can run before NumPy, SciPy, MATLAB, or Spinach.
"""

from __future__ import annotations

import array
import csv
import math
import re
import statistics
import sys
from pathlib import Path
from typing import Iterable


def parse_jcamp(path: Path) -> dict[str, str]:
    """Read scalar Bruker JCAMP parameters from a text file."""
    values: dict[str, str] = {}
    pattern = re.compile(r"^##\$(?P<key>[^=]+)=\s*(?P<value>.*)$")
    with path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            match = pattern.match(line.rstrip("\r\n"))
            if match:
                values[match.group("key").strip()] = match.group("value").strip()
    return values


def required_number(params: dict[str, str], name: str, cast):
    if name not in params:
        raise ValueError(f"Required Bruker parameter {name!r} is missing")
    raw = params[name].strip().strip("<>")
    try:
        return cast(float(raw)) if cast is int else cast(raw)
    except ValueError as exc:
        raise ValueError(f"Invalid Bruker parameter {name}={params[name]!r}") from exc


def read_processed_1r(path: Path, *, point_count: int, bytordp: int,
                      dtypp: int, nc_proc: int) -> list[float]:
    """Read a processed Bruker int32 ``1r`` spectrum."""
    if dtypp != 0:
        raise ValueError(f"Unsupported DTYPP={dtypp} in {path}; expected int32")
    if bytordp not in (0, 1):
        raise ValueError(f"Unsupported BYTORDP={bytordp} in {path}")
    raw = path.read_bytes()
    expected_bytes = point_count * 4
    if len(raw) != expected_bytes:
        raise ValueError(
            f"{path} contains {len(raw)} bytes; SI={point_count} expects {expected_bytes}"
        )
    values = array.array("i")
    values.frombytes(raw)
    file_is_little_endian = bytordp == 0
    if file_is_little_endian != (sys.byteorder == "little"):
        values.byteswap()
    scale = math.ldexp(1.0, nc_proc)
    return [float(value) * scale for value in values]


def indices_in(values: Iterable[float], bounds: tuple[float, float]) -> list[int]:
    low, high = bounds
    return [index for index, value in enumerate(values) if low <= value <= high]


def find_dss_peak(ppm: list[float], intensity: list[float],
                  search_region: tuple[float, float] = (-0.20, 0.20)) -> tuple[float, float]:
    """Find the DSS peak and return its ppm position plus an SNR proxy."""
    candidates = indices_in(ppm, search_region)
    if len(candidates) < 3:
        raise ValueError("The processed spectrum does not cover the DSS search region")
    local = [intensity[index] for index in candidates]
    baseline = statistics.median(local)
    scores = [abs(intensity[index] - baseline) for index in candidates]
    best_local = max(range(len(scores)), key=scores.__getitem__)
    best = candidates[best_local]
    if best_local in (0, len(candidates) - 1):
        raise ValueError("The strongest DSS candidate lies on the search-window boundary")
    left, center, right = scores[best_local - 1 : best_local + 2]
    denominator = left - 2.0 * center + right
    fractional_index = 0.0
    if denominator != 0.0:
        candidate = 0.5 * (left - right) / denominator
        if abs(candidate) <= 1.0:
            fractional_index = candidate
    point_step = ppm[best + 1] - ppm[best]
    dss_ppm = ppm[best] + fractional_index * point_step
    noise = statistics.median(abs(value - baseline) for value in local) or 1.0
    snr_proxy = center / noise
    if snr_proxy < 20.0:
        raise ValueError(f"DSS candidate is not sufficiently prominent (SNR={snr_proxy:.1f})")
    return dss_ppm, snr_proxy


def extract_dataset(experiment_dir: Path, procno: str = "1") -> dict[str, object]:
    """Extract acquisition, processing, and DSS-reference metadata."""
    experiment_dir = experiment_dir.resolve()
    processed_dir = experiment_dir / "pdata" / str(procno)
    acqus_path = experiment_dir / "acqus"
    procs_path = processed_dir / "procs"
    spectrum_path = processed_dir / "1r"
    for path in (acqus_path, procs_path, spectrum_path):
        if not path.is_file():
            raise FileNotFoundError(f"Required Bruker input is missing: {path}")

    acqus = parse_jcamp(acqus_path)
    procs = parse_jcamp(procs_path)
    si = required_number(procs, "SI", int)
    sf = required_number(procs, "SF", float)
    sw_p = required_number(procs, "SW_p", float)
    offset_raw = required_number(procs, "OFFSET", float)
    bytordp = required_number(procs, "BYTORDP", int)
    dtypp = required_number(procs, "DTYPP", int)
    nc_proc = required_number(procs, "NC_proc", int)
    intensity = read_processed_1r(
        spectrum_path, point_count=si, bytordp=bytordp,
        dtypp=dtypp, nc_proc=nc_proc,
    )
    ppm_step = (sw_p / sf) / max(si - 1, 1)
    ppm_raw = [offset_raw - index * ppm_step for index in range(si)]
    dss_raw, dss_snr = find_dss_peak(ppm_raw, intensity)

    sfo1 = required_number(acqus, "SFO1", float)
    nucleus = acqus.get("NUC1", "").strip().strip("<>") or "unknown"
    o1 = required_number(acqus, "O1", float)
    sw_h = required_number(acqus, "SW_h", float)
    td = required_number(acqus, "TD", int)
    ncomplex = td // 2 if td % 2 == 0 else td
    return {
        "experiment_dir": str(experiment_dir),
        "processed_dir": str(processed_dir),
        "processed_file": str(spectrum_path),
        "procno": str(procno),
        "field_mhz": f"{sfo1:.12g}",
        "nucleus": nucleus,
        "sfo1_mhz": f"{sfo1:.12g}",
        "o1_hz": f"{o1:.12g}",
        "sw_acq_hz": f"{sw_h:.12g}",
        "td_acquired": td,
        "ncomplex": ncomplex,
        "pulse_program": acqus.get("PULPROG", "").strip("<>"),
        "scans": required_number(acqus, "NS", int),
        "si": si,
        "sf_mhz": f"{sf:.12g}",
        "sw_p_hz": f"{sw_p:.12g}",
        "offset_raw_ppm": f"{offset_raw:.12g}",
        "dss_raw_ppm": f"{dss_raw:.9f}",
        "dss_shift_applied_ppm": f"{-dss_raw:.9f}",
        "offset_dss_ppm": f"{offset_raw - dss_raw:.12g}",
        "dss_snr_proxy": f"{dss_snr:.1f}",
        "bytordp": bytordp,
        "dtypp": dtypp,
        "nc_proc": nc_proc,
    }


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--dataset", action="append", metavar="EXPERIMENT_DIR",
        help="Bruker experiment directory; repeat for multiple fields",
    )
    source.add_argument(
        "--molecule", metavar="NAME",
        help="Use every imported dataset listed in data/NAME/NAME_config.json",
    )
    parser.add_argument(
        "--procno", action="append", metavar="N",
        help="pdata number for each dataset (one value applies to all; default 1)",
    )
    parser.add_argument("--output", type=Path, help="Output CSV path")
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    if args.molecule:
        from carbohydrate_config import load_config

        config = load_config(repo_root, args.molecule)
        datasets = [
            repo_root / "data" / args.molecule / item["relative_dir"]
            for item in config.get("datasets", [])
        ]
        procno = [str(item.get("procno", "1")) for item in config.get("datasets", [])]
        if not datasets:
            parser.error(f"no datasets listed for molecule {args.molecule!r}")
        output = args.output or repo_root / "outputs" / args.molecule / "bruker_metadata.csv"
    else:
        datasets = [Path(path) for path in args.dataset]
        procno = args.procno or ["1"]
        if len(procno) not in (1, len(datasets)):
            parser.error("provide either one --procno or one --procno per --dataset")
        if len(procno) == 1:
            procno *= len(datasets)
        output = args.output
        if output is None:
            parser.error("--output is required when using --dataset")

    rows = [extract_dataset(path, p) for path, p in zip(datasets, procno)]
    if not output.is_absolute():
        output = repo_root / output
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {output}")
    for row in rows:
        print(
            f"{float(row['field_mhz']):.3f} MHz ({row['nucleus']}): "
            f"raw OFFSET {row['offset_raw_ppm']} ppm; "
            f"DSS {row['dss_raw_ppm']} ppm; DSS-referenced OFFSET {row['offset_dss_ppm']} ppm"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
