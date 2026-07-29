#!/usr/bin/env python3
"""Create the sucrose four-field stacked plot with frequency axes in Hz.

The common DSS-referenced ppm coordinate is converted separately for each
field using its measured proton SFO1.  Thus each row has its own Hz scale,
which is the physically correct representation across magnetic fields.
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
HIGH_PPM, LOW_PPM = 5.80, 3.00

HEDGES = "#B4BD00"
GLORY = "#E4002B"
OLYMPIC = "#004E60"
BLACK = "#000000"
GRID = "#D9E2E5"


def sfo1_for(field: str) -> float:
    with (ROOT / "outputs/sucrose/bruker_metadata.csv").open(newline="") as fh:
        for row in csv.DictReader(fh):
            if round(float(row["field_mhz"])) == int(field):
                return float(row["sfo1_mhz"])
    raise FileNotFoundError(f"No SFO1 metadata for {field} MHz")


def load_field(field: str):
    folder = ROOT / f"outputs/sucrose/{field}MHz_spinach_candidate"
    curve_file = next(folder.glob(f"*curves_{field}MHz.csv"))
    summary_file = next(folder.glob("*summary.json"))
    offset = float(json.loads(summary_file.read_text())["ppm_offset_fitted"])
    rows = []
    with curve_file.open(newline="") as fh:
        for row in csv.DictReader(fh):
            def number(key: str) -> float:
                value = row[key].strip()
                return float(value) if value and value.lower() != "nan" else np.nan
            ppm_dss = float(row["ppm"]) - offset
            rows.append([ppm_dss, number("experiment"), number("candidate_spinach"), number("gissmo")])
    data = np.asarray(sorted(rows, key=lambda x: x[0]), dtype=float)
    data = data[(data[:, 0] >= LOW_PPM) & (data[:, 0] <= HIGH_PPM)]
    data[:, 1:] /= np.nanmax(data[:, 1:])
    sfo1 = sfo1_for(field)
    data[:, 0] *= sfo1  # Hz relative to DSS (DSS = 0 ppm)
    return data, sfo1, offset


def corr(a: np.ndarray, b: np.ndarray) -> float:
    good = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[good], b[good])[0, 1]) if good.sum() > 2 else float("nan")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(4, 1, figsize=(16, 11.5), sharex=True, sharey=False)
    fig.patch.set_facecolor("white")
    sfo_values = [sfo1_for(field) for field, _ in FIELDS]
    for index, ((field, role), ax) in enumerate(zip(FIELDS, axes)):
        data, sfo1, offset = load_field(field)
        ax.set_facecolor("white")
        ax.plot(data[:, 0], data[:, 1], color=HEDGES, lw=2.0, label="Experimental spectrum")
        ax.plot(data[:, 0], data[:, 2], color=GLORY, lw=1.9, ls="--", label="Candidate Spinach fit")
        ax.plot(data[:, 0], data[:, 3], color=OLYMPIC, lw=1.8, ls=":", label="GISSMO published simulation")
        r_se = corr(data[:, 2], data[:, 1])
        r_ge = corr(data[:, 3], data[:, 1])
        ax.set_ylim(-0.04, 1.10)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_yticklabels(["0", "0.5", "1.0"], fontsize=10)
        ax.grid(True, color=GRID, linewidth=0.8, alpha=0.85)
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", labelsize=11, width=1.2)
        for spine in ax.spines.values():
            spine.set_color(BLACK)
            spine.set_linewidth(1.3)
        ax.text(0.012, 0.80,
                f"{field} MHz — {role}  |  Spinach r = {r_se:.4f}  |  GISSMO r = {r_ge:.4f}",
                transform=ax.transAxes, fontsize=13, weight="bold", color=BLACK,
                bbox=dict(facecolor="white", edgecolor="#777777", alpha=0.93, pad=3.5))
        ax.text(0.965, 0.80, f"SFO1 = {sfo1:.3f} MHz  |  DSS correction {(-offset):+.5f} ppm",
                transform=ax.transAxes, ha="right", fontsize=9.5, color="#333333")
        if index < len(axes) - 1:
            ax.set_xlabel("")
        else:
            ax.set_xlabel(r"$^1$H frequency relative to DSS (Hz)", fontsize=15, weight="bold")

    # One shared Hz scale makes the field-dependent spreading visible: the
    # same ppm interval occupies more Hz at 1100 MHz than at 600 MHz.
    axes[-1].set_xlim(HIGH_PPM * max(sfo_values), LOW_PPM * min(sfo_values))
    fig.supylabel("Row-normalized intensity", fontsize=15, weight="bold")
    fig.suptitle("Sucrose multi-field validation — frequency-domain view",
                 fontsize=21, weight="bold", y=0.988)
    fig.text(0.5, 0.952, "All rows share one Hz scale; DSS is 0 Hz and 1100 MHz is most expanded",
             ha="center", fontsize=13, color="#333333")
    handles = [
        Line2D([0], [0], color=HEDGES, lw=2.5, label="Experimental spectrum"),
        Line2D([0], [0], color=GLORY, lw=2.4, ls="--", label="Candidate Spinach fit"),
        Line2D([0], [0], color=OLYMPIC, lw=2.3, ls=":", label="GISSMO published simulation"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.917),
               ncol=3, fontsize=12, frameon=True, facecolor="white", edgecolor="#555555")
    fig.text(0.5, 0.012, "Conversion is display-only; the fitted matrix and ppm values are unchanged.",
             ha="center", fontsize=10.5, color="#333333")
    fig.tight_layout(rect=(0.065, 0.04, 0.985, 0.875), h_pad=0.45)
    png = OUT / "sucrose_four_field_stacked_common_DSS_Hz_shared_scale.png"
    pdf = OUT / "sucrose_four_field_stacked_common_DSS_Hz_shared_scale.pdf"
    fig.savefig(png, dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(png)
    print(pdf)


if __name__ == "__main__":
    main()
