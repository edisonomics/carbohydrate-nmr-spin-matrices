"""Model helpers for single-component and anomer-mixture carbohydrates."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src" / "sucrose" / "bayes_astro"))
from sucrose_sim import load_gissmo_matrix, sucrose_sticks  # noqa: E402
from sucrose_sim import lorentzian_spectrum  # noqa: E402


def _blocks(component: dict[str, Any], matrix) -> list[list[int]]:
    blocks = component.get("blocks", [])
    if not blocks:
        return [list(range(matrix.shape[0]))]
    return [[int(index) - 1 for index in block] for block in blocks]


def component_specs(config: dict[str, Any], repo_root: Path) -> list[dict[str, Any]]:
    """Resolve configured physical components.

    A normal molecule has one component and may contain uncoupled blocks. A
    reducing sugar mixture has multiple independent component matrices (for
    example alpha and beta anomers), each with a population fraction.
    """
    configured = config.get("components")
    if configured:
        specs = []
        for component in configured:
            matrix_path = repo_root / str(component["matrix_file"])
            matrix = load_gissmo_matrix(matrix_path)
            specs.append({
                "name": str(component.get("name", f"component_{len(specs)+1}")),
                "fraction": float(component.get("fraction", 1.0)),
                "matrix_file": matrix_path,
                "matrix": matrix,
                "blocks": _blocks(component, matrix),
                "linewidth_hz": component.get("linewidth_hz"),
                "linewidth_provenance": component.get("linewidth_provenance"),
            })
        total = sum(item["fraction"] for item in specs)
        if total <= 0:
            raise ValueError("Mixture component fractions must sum to a positive value")
        for item in specs:
            item["fraction"] /= total
        return specs

    matrix_path = repo_root / str(config["matrix_file"])
    matrix = load_gissmo_matrix(matrix_path)
    return [{
        "name": str(config.get("name", "single")),
        "fraction": 1.0,
        "matrix_file": matrix_path,
        "matrix": matrix,
        "blocks": _blocks(config, matrix),
    }]


def mixture_sticks(config: dict[str, Any], repo_root: Path,
                   sfo1_mhz: float, carrier_ppm: float):
    """Return population-weighted sticks for a single molecule or mixture."""
    frequencies, intensities = [], []
    for component in component_specs(config, repo_root):
        f, a = sucrose_sticks(component["matrix"], sfo1_mhz, carrier_ppm,
                               blocks=component["blocks"])
        frequencies.append(f)
        intensities.append(component["fraction"] * a)
    import numpy as np
    return np.concatenate(frequencies), np.concatenate(intensities)


def mixture_spectrum(config: dict[str, Any], repo_root: Path, ppm_grid,
                     sfo1_mhz: float, carrier_ppm: float, default_lb_hz: float,
                     offset_ppm: float = 0.0):
    """Return a spectrum with optional component-specific linewidths.

    A single-molecule model uses ``default_lb_hz``.  For an anomer mixture,
    each component may document its own effective linewidth; this preserves
    the validated legacy behavior without hiding a magic number in code.
    """
    import numpy as np
    total = np.zeros_like(ppm_grid, dtype=float)
    for component in component_specs(config, repo_root):
        f, a = sucrose_sticks(component["matrix"], sfo1_mhz, carrier_ppm,
                               blocks=component["blocks"])
        linewidth = float(component.get("linewidth_hz", default_lb_hz))
        total += component["fraction"] * lorentzian_spectrum(
            ppm_grid, f, a, linewidth, sfo1_mhz, offset_ppm
        )
    return total
