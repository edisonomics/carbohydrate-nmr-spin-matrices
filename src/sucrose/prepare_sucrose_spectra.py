#!/usr/bin/env python3
"""Prepare the official multifield sucrose spectra for fitting.

For each selected Bruker acquisition this script:

1. reads the processed ``pdata/<procno>/1r`` spectrum and its ``procs`` file;
2. reconstructs the processed ppm axis;
3. finds the DSS peak in -0.20..0.20 ppm and moves it to 0 ppm;
4. marks the 3.00..5.80 ppm sucrose fitting region while excluding water
   (4.65..4.90 ppm) and the unexplained artifact band (5.15..5.30 ppm);
5. writes a full annotated spectrum, a fit-only spectrum, and a summary CSV.

The implementation intentionally uses only the Python standard library so the
preparation step does not depend on NumPy, SciPy, MATLAB, or Spinach.
"""

from __future__ import annotations

import argparse
import array
import csv
import math
import re
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from carbohydrate_config import load_config  # noqa: E402


REPO_ROOT = Path(__file__).resolve().parents[2]
CONFIG = load_config(REPO_ROOT, "sucrose")
PROCESSING = CONFIG["processing"]
SUCROSE_REGION = tuple(PROCESSING["fit_region_ppm"])
WATER_REGION = tuple(PROCESSING["water_region_ppm"])
ARTIFACT_REGION = tuple(PROCESSING["artifact_region_ppm"])
DSS_SEARCH_REGION = tuple(PROCESSING["dss_search_region_ppm"])
ANOMERIC_REGION = tuple(PROCESSING["anomeric_region_ppm"])


@dataclass(frozen=True)
class Dataset:
    key: str
    field_mhz: int
    acquisition: str
    relative_dir: str
    procno: str = "1"

    @property
    def experiment_dir(self) -> Path:
        return REPO_ROOT / "data" / "sucrose" / self.relative_dir

    @property
    def processed_dir(self) -> Path:
        return self.experiment_dir / "pdata" / self.procno

    @property
    def output_stem(self) -> str:
        return f"{self.field_mhz}_MHz_{self.acquisition}"


DATASETS = tuple(
    Dataset(
        str(item["key"]), int(item["field_mhz"]), str(item["acquisition"]),
        str(item["relative_dir"]), procno=str(item.get("procno", "1")),
    )
    for item in CONFIG["datasets"]
)


def parse_jcamp(path: Path) -> dict[str, str]:
    """Read single-line Bruker JCAMP parameters from a text file."""
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
    """Read an official dataset's processed Bruker 1r spectrum."""
    if dtypp != 0:
        raise ValueError(
            f"Unsupported DTYPP={dtypp} in {path}; official inputs are int32 (DTYPP=0)"
        )
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
    host_is_little_endian = sys.byteorder == "little"
    if file_is_little_endian != host_is_little_endian:
        values.byteswap()

    scale = math.ldexp(1.0, nc_proc)  # exactly 2**NC_proc, including negative values
    return [float(value) * scale for value in values]


def indices_in(values: Iterable[float], bounds: tuple[float, float]) -> list[int]:
    low, high = bounds
    return [index for index, value in enumerate(values) if low <= value <= high]


def find_dss_peak(ppm: list[float], intensity: list[float]) -> tuple[float, float]:
    """Return a sub-point DSS position and a simple signal-to-noise proxy."""
    candidates = indices_in(ppm, DSS_SEARCH_REGION)
    if len(candidates) < 3:
        raise ValueError("The processed spectrum does not cover the DSS search region")

    local = [intensity[index] for index in candidates]
    baseline = statistics.median(local)
    scores = [abs(intensity[index] - baseline) for index in candidates]
    best_local = max(range(len(scores)), key=scores.__getitem__)
    best = candidates[best_local]

    if best_local in (0, len(candidates) - 1):
        raise ValueError("The strongest DSS candidate lies on the search-window boundary")

    # Quadratic interpolation of the absolute peak amplitude supplies a stable
    # sub-point reference while retaining the same peak-picking rule as the
    # historical dss_reference_600_900.py workflow.
    left, center, right = scores[best_local - 1 : best_local + 2]
    denominator = left - 2.0 * center + right
    fractional_index = 0.0
    if denominator != 0.0:
        candidate = 0.5 * (left - right) / denominator
        if abs(candidate) <= 1.0:
            fractional_index = candidate
    point_step = ppm[best + 1] - ppm[best]
    dss_ppm = ppm[best] + fractional_index * point_step

    absolute_deviations = [abs(value - baseline) for value in local]
    noise = statistics.median(absolute_deviations) or 1.0
    snr_proxy = center / noise
    if snr_proxy < 20.0:
        raise ValueError(
            f"DSS candidate is not sufficiently prominent (SNR proxy={snr_proxy:.1f})"
        )
    return dss_ppm, snr_proxy


def exclusion_reason(ppm: float) -> str:
    if not SUCROSE_REGION[0] <= ppm <= SUCROSE_REGION[1]:
        return "outside_sucrose_region"
    if WATER_REGION[0] <= ppm <= WATER_REGION[1]:
        return "water"
    if ARTIFACT_REGION[0] <= ppm <= ARTIFACT_REGION[1]:
        return "artifact_5p15_5p30"
    return "included"


def write_full_spectrum(path: Path, ppm_raw: list[float], ppm_dss: list[float],
                        intensity: list[float], corrected: list[float],
                        reasons: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow((
            "ppm_raw", "ppm_dss", "intensity_scaled",
            "intensity_baseline_corrected", "included_in_fit", "exclusion_reason",
        ))
        for row in zip(ppm_raw, ppm_dss, intensity, corrected, reasons):
            raw_ppm, ref_ppm, raw_y, corrected_y, reason = row
            writer.writerow((
                f"{raw_ppm:.12g}", f"{ref_ppm:.12g}", f"{raw_y:.12g}",
                f"{corrected_y:.12g}", int(reason == "included"), reason,
            ))


def write_fit_spectrum(path: Path, ppm_dss: list[float], intensity: list[float],
                       corrected: list[float], reasons: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(("ppm_dss", "intensity_scaled", "intensity_baseline_corrected"))
        for ppm, raw_y, corrected_y, reason in zip(ppm_dss, intensity, corrected, reasons):
            if reason == "included":
                writer.writerow((f"{ppm:.12g}", f"{raw_y:.12g}", f"{corrected_y:.12g}"))


def prepare_dataset(dataset: Dataset, output_dir: Path) -> dict[str, object]:
    procs_path = dataset.processed_dir / "procs"
    spectrum_path = dataset.processed_dir / "1r"
    acqus_path = dataset.experiment_dir / "acqus"
    for path in (procs_path, spectrum_path, acqus_path):
        if not path.is_file():
            raise FileNotFoundError(f"Required input is missing: {path}")

    procs = parse_jcamp(procs_path)
    acqus = parse_jcamp(acqus_path)
    sfo1_mhz = required_number(acqus, "SFO1", float)
    o1_hz = required_number(acqus, "O1", float)
    sw_acq_hz = required_number(acqus, "SW_h", float)
    td_acquired = required_number(acqus, "TD", int)
    ncomplex = td_acquired // 2 if td_acquired % 2 == 0 else td_acquired
    si = required_number(procs, "SI", int)
    sf = required_number(procs, "SF", float)
    sw_p = required_number(procs, "SW_p", float)
    offset = required_number(procs, "OFFSET", float)
    bytordp = required_number(procs, "BYTORDP", int)
    dtypp = required_number(procs, "DTYPP", int)
    nc_proc = required_number(procs, "NC_proc", int)

    intensity = read_processed_1r(
        spectrum_path,
        point_count=si,
        bytordp=bytordp,
        dtypp=dtypp,
        nc_proc=nc_proc,
    )
    spectral_width_ppm = sw_p / sf
    denominator = max(si - 1, 1)
    ppm_raw = [offset - index * spectral_width_ppm / denominator for index in range(si)]

    dss_raw_ppm, dss_snr_proxy = find_dss_peak(ppm_raw, intensity)
    ppm_dss = [ppm - dss_raw_ppm for ppm in ppm_raw]
    reasons = [exclusion_reason(ppm) for ppm in ppm_dss]
    fit_indices = [index for index, reason in enumerate(reasons) if reason == "included"]
    if not fit_indices:
        raise ValueError(f"No points remain in the fit mask for {dataset.relative_dir}")

    baseline = statistics.median(intensity[index] for index in fit_indices)
    corrected = [value - baseline for value in intensity]
    anomeric = indices_in(ppm_dss, ANOMERIC_REGION)
    if not anomeric:
        raise ValueError(f"No points cover the anomeric region for {dataset.relative_dir}")
    anomeric_index = max(anomeric, key=lambda index: corrected[index])

    full_name = f"{dataset.output_stem}_full.csv"
    fit_name = f"{dataset.output_stem}_fit.csv"
    write_full_spectrum(
        output_dir / full_name, ppm_raw, ppm_dss, intensity, corrected, reasons
    )
    write_fit_spectrum(output_dir / fit_name, ppm_dss, intensity, corrected, reasons)

    pulse_program = acqus.get("PULPROG", "").strip("<>")
    scans = required_number(acqus, "NS", int)
    counts = {reason: reasons.count(reason) for reason in set(reasons)}
    return {
        "field_mhz": dataset.field_mhz,
        "acquisition": dataset.acquisition,
        "relative_dir": dataset.relative_dir,
        "procno": dataset.procno,
        "pulse_program": pulse_program,
        "scans": scans,
        "sfo1_mhz": f"{sfo1_mhz:.12g}",
        "o1_hz": f"{o1_hz:.12g}",
        "sw_acq_hz": f"{sw_acq_hz:.12g}",
        "td_acquired": td_acquired,
        "ncomplex": ncomplex,
        "points": si,
        "sf_mhz": f"{sf:.12g}",
        "sw_hz": f"{sw_p:.12g}",
        "offset_ppm": f"{offset:.12g}",
        "offset_dss_ppm": f"{offset - dss_raw_ppm:.12g}",
        "dss_raw_ppm": f"{dss_raw_ppm:.9f}",
        "dss_shift_applied_ppm": f"{-dss_raw_ppm:.9f}",
        "bytordp": bytordp,
        "dtypp": dtypp,
        "nc_proc": nc_proc,
        "dss_snr_proxy": f"{dss_snr_proxy:.1f}",
        "anomeric_peak_ppm_dss": f"{ppm_dss[anomeric_index]:.6f}",
        "fit_points": counts.get("included", 0),
        "water_points_excluded": counts.get("water", 0),
        "artifact_points_excluded": counts.get("artifact_5p15_5p30", 0),
        "full_spectrum": full_name,
        "fit_spectrum": fit_name,
    }


def write_summary(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset",
        action="append",
        choices=[dataset.key for dataset in DATASETS],
        help="prepare only this field strength; repeat for multiple fields",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "outputs" / "sucrose" / "prepared",
        help="output directory (default: outputs/sucrose/prepared)",
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    if not output_dir.is_absolute():
        output_dir = REPO_ROOT / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    selected = [
        dataset for dataset in DATASETS
        if not args.dataset or dataset.key in set(args.dataset)
    ]
    rows = []
    for dataset in selected:
        row = prepare_dataset(dataset, output_dir)
        rows.append(row)
        print(
            f"{dataset.field_mhz:4d} MHz/{dataset.acquisition}: "
            f"DSS {row['dss_raw_ppm']} ppm -> 0; "
            f"anomeric {row['anomeric_peak_ppm_dss']} ppm; "
            f"fit points {row['fit_points']}"
        )

    summary_path = output_dir / "preparation_summary.csv"
    write_summary(summary_path, rows)
    print(f"Wrote {summary_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
