#!/usr/bin/env python3
"""Create four-field panels for two diagnostic sucrose regions."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


ROOT = Path("/Users/cece/Desktop/edison_lab/final_repo")
OUT = ROOT / "outputs/sucrose/common_reference_plots"
FIELDS = (("600", "training"), ("800", "validation"),
          ("900", "training"), ("1100", "validation"))
HEDGES, GLORY, OLYMPIC = "#B4BD00", "#E4002B", "#004E60"
GRID, BLACK = "#D9E2E5", "#000000"
WINDOWS = {
    "crowded_3p67ppm": (3.69, 3.65),
    "anomeric_5p40ppm": (5.43, 5.38),
}


def read_field(field: str, limits: tuple[float, float]):
    high, low = limits
    folder = ROOT / f"outputs/sucrose/{field}MHz_spinach_candidate"
    curve = next(folder.glob(f"*curves_{field}MHz.csv"))
    summary = json.loads(next(folder.glob("*summary.json")).read_text())
    offset = float(summary["ppm_offset_fitted"])
    rows = []
    with curve.open(newline="") as fh:
        for row in csv.DictReader(fh):
            ppm = float(row["ppm"]) - offset
            if low <= ppm <= high:
                def val(key: str) -> float:
                    x = row[key].strip()
                    return float(x) if x and x.lower() != "nan" else np.nan
                rows.append((ppm, val("experiment"), val("candidate_spinach"), val("gissmo")))
    return np.asarray(sorted(rows, key=lambda x: x[0]), dtype=float)


def corr(a: np.ndarray, b: np.ndarray) -> float:
    good = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(a[good], b[good])[0, 1]) if good.sum() > 2 else float("nan")


def panel(name: str, limits: tuple[float, float]) -> tuple[Path, Path]:
    high, low = limits
    fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharex=True, sharey=False)
    fig.patch.set_facecolor("white")
    for ax, (field, role) in zip(axes.flat, FIELDS):
        data = read_field(field, limits)
        ax.set_facecolor("white")
        ax.plot(data[:, 0], data[:, 1], color=HEDGES, lw=2.5, label="Experiment")
        ax.plot(data[:, 0], data[:, 2], color=GLORY, lw=2.3, ls="--", label="Candidate Spinach")
        ax.plot(data[:, 0], data[:, 3], color=OLYMPIC, lw=2.3, ls=":", label="GISSMO")
        ymax = max(0.08, float(np.nanmax(data[:, 1:])) * 1.14)
        ax.set_xlim(high, low)
        ax.set_ylim(-0.04 * ymax, ymax)
        ax.grid(True, color=GRID, linewidth=0.9, alpha=0.9)
        ax.set_axisbelow(True)
        ax.tick_params(labelsize=12, width=1.3)
        for spine in ax.spines.values():
            spine.set_color(BLACK)
            spine.set_linewidth(1.6)
        ax.set_title(f"{field} MHz — {role}", fontsize=16, weight="bold", pad=9)
        ax.text(0.02, 0.88,
                f"Spinach r = {corr(data[:,2], data[:,1]):.4f}\nGISSMO r = {corr(data[:,3], data[:,1]):.4f}",
                transform=ax.transAxes, fontsize=11.5,
                bbox=dict(facecolor="white", edgecolor="#777777", alpha=0.94))
    for ax in axes[0, :]: ax.set_xlabel("")
    for ax in axes[:, 1]: ax.set_ylabel("")
    axes[1, 0].set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=14, weight="bold")
    axes[1, 1].set_xlabel(r"$^1$H chemical shift (ppm)", fontsize=14, weight="bold")
    axes[0, 0].set_ylabel("Normalized intensity", fontsize=14, weight="bold")
    axes[1, 0].set_ylabel("Normalized intensity", fontsize=14, weight="bold")
    handles = [
        plt.Line2D([0], [0], color=HEDGES, lw=2.5, label="Experiment"),
        plt.Line2D([0], [0], color=GLORY, lw=2.4, ls="--", label="Candidate Spinach"),
        plt.Line2D([0], [0], color=OLYMPIC, lw=2.4, ls=":", label="GISSMO"),
    ]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=3,
               fontsize=12, frameon=True, facecolor="white", edgecolor="#555555")
    center = (high + low) / 2
    fig.suptitle(f"Sucrose {center:.2f} ppm diagnostic super-zoom — all four fields",
                 fontsize=21, weight="bold", y=1.0)
    fig.tight_layout(rect=(0, 0.02, 1, 0.94))
    png = OUT / f"sucrose_all_fields_superzoom_{name}.png"
    pdf = OUT / f"sucrose_all_fields_superzoom_{name}.pdf"
    fig.savefig(png, dpi=350, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return png, pdf


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, limits in WINDOWS.items():
        png, pdf = panel(name, limits)
        print(png)
        print(pdf)


if __name__ == "__main__":
    main()
