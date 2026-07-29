#!/usr/bin/env python3
"""Create a poster-style stacked four-field sucrose comparison.

Inputs are the latest candidate Spinach curve exports.  The curves are moved
to the common DSS axis using the fitted per-field ppm offset, then normalized
within each row for display.  This is presentation-only; it does not refit
the matrix or alter any data.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path("/Users/cece/Desktop/edison_lab/final_repo")
OUT = ROOT / "outputs/sucrose/common_reference_plots"
FIELDS = (("600", "training"), ("800", "validation"),
          ("900", "training"), ("1100", "validation"))
HIGH, LOW = 5.80, 3.00

HEDGES = "#B4BD00"   # experiment
GLORY = "#E4002B"    # candidate Spinach
OLYMPIC = "#004E60"  # deposited GISSMO
BLACK = "#000000"
GRID = "#D9E2E5"


def load_field(field: str):
    folder = ROOT / f"outputs/sucrose/{field}MHz_spinach_candidate"
    curve_file = next(folder.glob(f"*curves_{field}MHz.csv"))
    summary_file = next(folder.glob("*summary.json"))
    summary = json.loads(summary_file.read_text())
    offset = float(summary["ppm_offset_fitted"])
    rows = []
    with curve_file.open(newline="") as fh:
        for row in csv.DictReader(fh):
            def number(key: str) -> float:
                value = row[key].strip()
                return float(value) if value and value.lower() != "nan" else np.nan
            # model coordinate -> common DSS-referenced coordinate
            rows.append([float(row["ppm"]) - offset,
                         number("experiment"),
                         number("candidate_spinach"),
                         number("gissmo")])
    data = np.asarray(sorted(rows, key=lambda x: x[0]), dtype=float)
    mask = (data[:, 0] >= LOW) & (data[:, 0] <= HIGH)
    data = data[mask]
    # Equalize each row for a readable poster comparison.
    scale = np.nanmax(data[:, 1:])
    data[:, 1:] /= scale
    return data, offset


def corr(a: np.ndarray, b: np.ndarray) -> float:
    good = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[good], b[good])[0, 1]) if good.sum() > 2 else float("nan")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    # Share the ppm axis, but keep independent y-limits so each vertically
    # offset trace remains visible.
    fig, axes = plt.subplots(4, 1, figsize=(16, 11.5), sharex=True, sharey=False)
    fig.patch.set_facecolor("white")
    baseline_step = 1.35
    traces = []

    for index, ((field, role), ax) in enumerate(zip(FIELDS, axes)):
        data, offset = load_field(field)
        base = (len(FIELDS) - 1 - index) * baseline_step
        traces.append((data, base))
        ax.set_facecolor("white")
        ax.plot(data[:, 0], data[:, 1] + base, color=HEDGES, lw=2.0, zorder=3)
        ax.plot(data[:, 0], data[:, 2] + base, color=GLORY, lw=1.9, ls="--", zorder=4)
        ax.plot(data[:, 0], data[:, 3] + base, color=OLYMPIC, lw=1.8, ls=":", zorder=2)
        r_se = corr(data[:, 2], data[:, 1])
        r_ge = corr(data[:, 3], data[:, 1])
        ax.text(0.012, 0.80,
                f"{field} MHz — {role}  |  Spinach r = {r_se:.4f}  |  GISSMO r = {r_ge:.4f}",
                transform=ax.transAxes, fontsize=13, weight="bold", color=BLACK,
                bbox=dict(facecolor="white", edgecolor="#777777", alpha=0.93, pad=3.5))
        ax.text(0.965, 0.80, f"DSS axis correction {(-offset):+.5f} ppm",
                transform=ax.transAxes, ha="right", fontsize=9.5, color="#333333")
        ax.set_ylim(-0.04 + base, 1.10 + base)
        ax.set_yticks([base, base + 0.5, base + 1.0])
        ax.set_yticklabels(["0", "0.5", "1.0"], fontsize=10)
        ax.grid(True, color=GRID, linewidth=0.8, alpha=0.85)
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", labelsize=11, width=1.2)
        for spine in ax.spines.values():
            spine.set_color(BLACK)
            spine.set_linewidth(1.3)

    axes[-1].set_xlim(HIGH, LOW)
    axes[-1].set_xlabel(r"$^1$H chemical shift (ppm), common DSS-referenced axis", fontsize=15, weight="bold")
    fig.supylabel("Row-normalized intensity (vertical offsets for display)", fontsize=15, weight="bold")
    fig.suptitle("Sucrose multi-field validation — shared 14-spin matrix",
                 fontsize=21, weight="bold", y=0.988)
    fig.text(0.5, 0.952, "600/900 MHz training; 800/1100 MHz blind validation",
             ha="center", fontsize=13, color="#333333")
    handles = [
        Line2D([0], [0], color=HEDGES, lw=2.5, label="Experimental spectrum"),
        Line2D([0], [0], color=GLORY, lw=2.4, ls="--", label="Candidate Spinach fit"),
        Line2D([0], [0], color=OLYMPIC, lw=2.3, ls=":", label="GISSMO published simulation"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.917),
               ncol=3, fontsize=12, frameon=True, facecolor="white", edgecolor="#555555")
    fig.text(0.5, 0.012, "Display normalization and vertical offsets do not change the fitted spectra.",
             ha="center", fontsize=10.5, color="#333333")
    fig.tight_layout(rect=(0.065, 0.04, 0.985, 0.875), h_pad=0.45)

    png = OUT / "sucrose_four_field_stacked_common_DSS.png"
    pdf = OUT / "sucrose_four_field_stacked_common_DSS.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
