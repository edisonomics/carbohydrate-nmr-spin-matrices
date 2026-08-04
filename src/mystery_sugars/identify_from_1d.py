#!/usr/bin/env python3
"""Rank carbohydrate identity hypotheses from complete multifield 1-D spectra.

This is an identity-free screening stage. It does not edit a spin matrix and
it never claims that 1-D evidence alone proves a molecular identity. The
prepared spectra are compared against complete reference proton-shift
fingerprints.  When a local spin matrix exists, the native Python physics
engine also forward-simulates the full multiplet shape and scalar couplings.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any, Iterable

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_BMRB_CATALOG = Path(__file__).with_name("bmrb_candidate_catalog.json")
FULL_PROTON_REGION = (0.50, 10.00)
ANOMERIC_REGION = (4.30, 5.80)
WATER_REGION = (4.65, 4.90)


def _catalog_candidate(item: dict[str, Any]) -> dict[str, Any]:
    """Adapt one explore_BMRB record to the existing ranking contract."""

    centers = [float(value) for value in item.get("reference_anomeric_centers_ppm", [])]
    proton_shifts = [float(value) for value in item.get("reference_proton_shifts_ppm", [])]
    return {
        "id": item["candidate_id"],
        "name": item["name"],
        "class": "chebi_carbohydrate_unreviewed",
        "bmrb_entry": item.get("selected_bmrb_entry"),
        "reference_shift_file": None,
        "reference_anomeric_centers_ppm": centers,
        "reference_proton_shifts_ppm": proton_shifts,
        "expected_anomeric_clusters": len(centers) if centers else None,
        "expected_components": None,
        "forms": [],
        "topology": "unreviewed BMRB/ChEBI carbohydrate candidate",
        "notes": item.get("review_notes"),
        "review_status": item.get("review_status", "needs_review"),
        "bmrb_evidence": item.get("evidence", {}),
        "identity_warning": item.get("identity_warning", ""),
    }


def load_library(
    path: Path | None = None,
    *,
    include_bmrb_catalog: bool = False,
    bmrb_catalog_path: Path | None = None,
) -> list[dict[str, Any]]:
    path = path or Path(__file__).with_name("candidate_library.json")
    payload = json.loads(path.read_text(encoding="utf-8"))
    candidates = list(payload.get("candidates", []))
    for candidate in candidates:
        candidate.setdefault("review_status", "approved")

    if not include_bmrb_catalog:
        return candidates
    catalog_path = bmrb_catalog_path or DEFAULT_BMRB_CATALOG
    catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
    existing_entries = {
        candidate.get("bmrb_entry")
        for candidate in candidates
        if candidate.get("bmrb_entry")
    }
    for item in catalog.get("candidates", []):
        if item.get("selected_bmrb_entry") in existing_entries:
            continue
        adapted = _catalog_candidate(item)
        if adapted["reference_anomeric_centers_ppm"] or adapted["reference_proton_shifts_ppm"]:
            candidates.append(adapted)
    return candidates


def load_prepared(repo_root: Path, molecule: str) -> list[dict[str, Any]]:
    summary = repo_root / "outputs" / molecule / "prepared" / "preparation_summary.csv"
    if not summary.is_file():
        raise FileNotFoundError(
            f"Missing {summary}; prepare the 1-D spectra first with "
            "prepare_carbohydrate_spectra.py"
        )
    rows: list[dict[str, Any]] = []
    with summary.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            spectrum_name = row.get("full_spectrum") or row["fit_spectrum"]
            path = repo_root / "outputs" / molecule / "prepared" / spectrum_name
            with path.open(newline="", encoding="utf-8") as spectrum:
                data = list(csv.DictReader(spectrum))
            ppm = np.array([float(item["ppm_dss"]) for item in data], dtype=float)
            value_key = "intensity_baseline_corrected"
            y = np.array([float(item[value_key]) for item in data], dtype=float)
            order = np.argsort(ppm)
            rows.append({
                "key": str(row["key"]),
                "field_mhz": float(row["field_mhz"]),
                "ppm": ppm[order],
                "intensity": y[order],
            })
    return rows


def _cluster_centers(values: Iterable[float], gap_ppm: float = 0.04) -> list[float]:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return []
    clusters: list[list[float]] = [[ordered[0]]]
    for value in ordered[1:]:
        if value - clusters[-1][-1] <= gap_ppm:
            clusters[-1].append(value)
        else:
            clusters.append([value])
    return [float(np.mean(cluster)) for cluster in clusters]


def detect_anomeric_clusters(
    ppm: np.ndarray,
    intensity: np.ndarray,
    *,
    region: tuple[float, float] = ANOMERIC_REGION,
    water_region: tuple[float, float] = WATER_REGION,
    min_spacing_ppm: float = 0.008,
    cluster_gap_ppm: float = 0.04,
) -> list[float]:
    """Detect reproducible anomeric multiplet clusters in one prepared 1-D trace."""
    mask = (ppm >= region[0]) & (ppm <= region[1])
    mask &= ~((ppm >= water_region[0]) & (ppm <= water_region[1]))
    x, y = ppm[mask], np.asarray(intensity[mask], dtype=float)
    if x.size < 3 or not np.any(np.isfinite(y)):
        return []
    y = np.nan_to_num(y, nan=0.0)
    y -= np.percentile(y, 5)
    peak_height = float(np.max(y))
    if peak_height <= 0:
        return []
    threshold = max(0.015 * peak_height, float(np.percentile(y, 70)) * 0.35)
    candidates = [
        i for i in range(1, x.size - 1)
        if y[i] >= threshold and y[i] >= y[i - 1] and y[i] >= y[i + 1]
    ]
    candidates.sort(key=lambda i: float(y[i]), reverse=True)
    selected: list[int] = []
    for index in candidates:
        if all(abs(float(x[index] - x[other])) >= min_spacing_ppm for other in selected):
            selected.append(index)
    return _cluster_centers((float(x[index]) for index in selected), cluster_gap_ppm)


def detect_full_proton_clusters(
    ppm: np.ndarray,
    intensity: np.ndarray,
    *,
    region: tuple[float, float] = FULL_PROTON_REGION,
    water_region: tuple[float, float] = WATER_REGION,
) -> list[float]:
    """Detect multiplet centers across the complete prepared proton window."""

    return detect_anomeric_clusters(
        ppm,
        intensity,
        region=region,
        water_region=water_region,
        min_spacing_ppm=0.008,
        cluster_gap_ppm=0.04,
    )


def _reference_proton_values(repo_root: Path, candidate: dict[str, Any]) -> list[float]:
    embedded = candidate.get("reference_proton_shifts_ppm")
    if embedded:
        return [float(value) for value in embedded]
    relative = candidate.get("reference_shift_file")
    if not relative:
        return []
    path = repo_root / str(relative)
    if not path.is_file():
        return []
    values: list[float] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("Atom_type") != "H":
                continue
            try:
                values.append(float(row["Val"]))
            except (KeyError, TypeError, ValueError):
                continue
    return values


def _outside_water(value: float, water_region: tuple[float, float] = WATER_REGION) -> bool:
    return not water_region[0] <= value <= water_region[1]


def _reference_full_proton_centers(
    repo_root: Path,
    candidate: dict[str, Any],
    region: tuple[float, float] = FULL_PROTON_REGION,
) -> list[float]:
    values = [
        value
        for value in _reference_proton_values(repo_root, candidate)
        if region[0] <= value <= region[1] and _outside_water(value)
    ]
    return _cluster_centers(values, gap_ppm=0.04)


def _reference_all_anomeric_centers(
    repo_root: Path, candidate: dict[str, Any]
) -> list[float]:
    embedded = candidate.get("reference_anomeric_centers_ppm")
    if embedded:
        return _cluster_centers((float(value) for value in embedded), gap_ppm=0.04)
    values = [
        value
        for value in _reference_proton_values(repo_root, candidate)
        if ANOMERIC_REGION[0] <= value <= ANOMERIC_REGION[1]
    ]
    return _cluster_centers(values, gap_ppm=0.04)


def _reference_anomeric_centers(repo_root: Path, candidate: dict[str, Any]) -> list[float]:
    return [
        value
        for value in _reference_all_anomeric_centers(repo_root, candidate)
        if _outside_water(value)
    ]


def _match_score(observed: list[float], reference: list[float], tolerance_ppm: float = 0.10) -> float | None:
    if not observed or not reference:
        return None
    distances = []
    for target in reference:
        distances.append(min(abs(target - value) for value in observed))
    scores = [math.exp(-0.5 * (distance / tolerance_ppm) ** 2) for distance in distances]
    return float(np.mean(scores))


def _bidirectional_match_score(
    observed: list[float],
    reference: list[float],
    tolerance_ppm: float = 0.08,
) -> float | None:
    """Score both reference coverage and unexplained observed multiplets."""

    forward = _match_score(observed, reference, tolerance_ppm)
    reverse = _match_score(reference, observed, tolerance_ppm)
    if forward is None or reverse is None:
        return None
    return 0.70 * forward + 0.30 * reverse


def _resolved_lines_near(
    field: dict[str, Any],
    center_ppm: float,
    *,
    window_ppm: float = 0.025,
    minimum_distance_hz: float = 0.8,
) -> dict[str, Any]:
    """Find resolved positive lines near one anomeric multiplet center."""

    from scipy.signal import find_peaks  # noqa: PLC0415

    ppm = np.asarray(field["ppm"], dtype=float)
    intensity = np.asarray(field["intensity"], dtype=float)
    mask = (ppm >= center_ppm - window_ppm) & (ppm <= center_ppm + window_ppm)
    x, y = ppm[mask], intensity[mask]
    order = np.argsort(x)
    x, y = x[order], y[order]
    if x.size < 8:
        return {"center_ppm": center_ppm, "lines_ppm": [], "spacings_hz": []}
    z = y - np.median(y)
    differences = np.diff(z)
    noise = (
        float(np.median(np.abs(differences - np.median(differences))))
        / 0.6745
        / math.sqrt(2.0)
    ) or 1.0
    prominence = max(6.0 * noise, 0.05 * float(np.max(z)))
    step_ppm = abs(float(np.median(np.diff(x))))
    distance = max(
        1,
        int(round(minimum_distance_hz / (step_ppm * float(field["field_mhz"])))),
    )
    peaks, properties = find_peaks(z, prominence=prominence, distance=distance)
    lines = sorted(float(x[index]) for index in peaks)
    spacings = sorted(
        (lines[j] - lines[i]) * float(field["field_mhz"])
        for i in range(len(lines))
        for j in range(i + 1, len(lines))
        if 0.5
        <= (lines[j] - lines[i]) * float(field["field_mhz"])
        <= 15.0
    )
    return {
        "center_ppm": center_ppm,
        "lines_ppm": lines,
        "spacings_hz": spacings,
        "noise_proxy": noise,
        "prominence_threshold": prominence,
        "line_prominences": [float(value) for value in properties["prominences"]],
    }


def build_bubb_observations(fields: list[dict[str, Any]]) -> dict[str, Any]:
    """Measure Bubb-style reporter evidence once for all candidates."""

    field_reports: list[dict[str, Any]] = []
    spacing_rows: list[dict[str, float]] = []
    for field in fields:
        centers = detect_anomeric_clusters(field["ppm"], field["intensity"])
        multiplets = []
        for center in centers:
            report = _resolved_lines_near(field, center)
            multiplets.append(report)
            for spacing in report["spacings_hz"]:
                spacing_rows.append({
                    "field_mhz": float(field["field_mhz"]),
                    "center_ppm": float(center),
                    "spacing_hz": float(spacing),
                })
        full_centers = detect_full_proton_clusters(field["ppm"], field["intensity"])
        field_reports.append({
            "field_mhz": float(field["field_mhz"]),
            "anomeric_centers_ppm": centers,
            "anomeric_multiplets": multiplets,
            "methyl_reporters_ppm": [
                value for value in full_centers if 1.1 <= value <= 1.3
            ],
            "acetyl_reporters_ppm": [
                value for value in full_centers if 2.0 <= value <= 2.1
            ],
        })

    center_groups: list[list[dict[str, float]]] = []
    for row in sorted(spacing_rows, key=lambda item: item["center_ppm"]):
        target = next(
            (
                group
                for group in center_groups
                if abs(float(np.mean([item["center_ppm"] for item in group]))
                       - row["center_ppm"]) <= 0.04
            ),
            None,
        )
        if target is None:
            target = []
            center_groups.append(target)
        target.append(row)

    consensus = []
    for center_group in center_groups:
        spacing_groups: list[list[dict[str, float]]] = []
        for row in sorted(center_group, key=lambda item: item["spacing_hz"]):
            target = next(
                (
                    group
                    for group in spacing_groups
                    if abs(float(np.mean([item["spacing_hz"] for item in group]))
                           - row["spacing_hz"]) <= 0.6
                ),
                None,
            )
            if target is None:
                target = []
                spacing_groups.append(target)
            target.append(row)
        supported = [
            group
            for group in spacing_groups
            if len({item["field_mhz"] for item in group}) >= 2
        ]
        if not supported:
            continue
        best = max(
            supported,
            key=lambda group: (
                len({item["field_mhz"] for item in group}),
                len(group),
            ),
        )
        consensus.append({
            "center_ppm": float(np.median([item["center_ppm"] for item in center_group])),
            "spacing_hz": float(np.median([item["spacing_hz"] for item in best])),
            "field_support": len({item["field_mhz"] for item in best}),
            "field_values": best,
        })
    return {
        "field_reports": field_reports,
        "consensus_anomeric_spacings": sorted(
            consensus, key=lambda item: item["center_ppm"]
        ),
        "interpretation": (
            "Cross-field resolved anomeric line spacings are Bubb-style "
            "screening evidence. Crowded or strongly coupled line separations "
            "are not automatically assigned as scalar couplings."
        ),
    }


def _range_score(value: float, expected_range: list[float]) -> float:
    lower, upper = (float(expected_range[0]), float(expected_range[1]))
    distance = max(lower - value, 0.0, value - upper)
    return math.exp(-0.5 * (distance / 0.8) ** 2)


def score_bubb_guidance(
    candidate: dict[str, Any],
    fields: list[dict[str, Any]],
    repo_root: Path,
    observations: dict[str, Any] | None = None,
) -> dict[str, Any] | None:
    """Apply Bubb reporter-group rules without inventing assignments."""

    all_anomeric = _reference_all_anomeric_centers(repo_root, candidate)
    full_reference = _reference_full_proton_centers(repo_root, candidate)
    if not all_anomeric and not full_reference:
        return None

    common_path = repo_root / "src" / "common"
    if str(common_path) not in sys.path:
        sys.path.insert(0, str(common_path))
    from bubb_rules import BUBB_REFERENCE, profile_for  # noqa: PLC0415

    observations = observations or build_bubb_observations(fields)
    profile_name = candidate.get("bubb_profile")
    profile = profile_for(str(profile_name)) if profile_name else None
    visible_reference = [
        value for value in all_anomeric if _outside_water(value)
    ]
    field_counts = [
        len(item["anomeric_centers_ppm"])
        for item in observations["field_reports"]
    ]
    count_score = None
    if visible_reference and field_counts:
        count_score = float(np.mean([
            math.exp(-0.8 * abs(count - len(visible_reference)))
            for count in field_counts
        ]))

    j_checks = []
    j_scores = []
    patterns = list((profile or {}).get("anomeric_j_patterns_hz", []))
    ordered_reference = sorted(all_anomeric, reverse=True)
    for measured in observations["consensus_anomeric_spacings"]:
        if not ordered_reference or not patterns:
            continue
        nearest_index = min(
            range(len(ordered_reference)),
            key=lambda index: abs(ordered_reference[index] - measured["center_ppm"]),
        )
        if nearest_index >= len(patterns):
            continue
        pattern = patterns[nearest_index]
        score = _range_score(float(measured["spacing_hz"]), pattern["range_hz"])
        j_scores.append(score)
        j_checks.append({
            "form": pattern["form"],
            "reference_center_ppm": ordered_reference[nearest_index],
            "observed_center_ppm": measured["center_ppm"],
            "observed_spacing_hz": measured["spacing_hz"],
            "expected_range_hz": pattern["range_hz"],
            "field_support": measured["field_support"],
            "score": score,
        })
    j_score = float(np.mean(j_scores)) if j_scores else None

    diagnostic_checks = []
    diagnostic_scores = []
    diagnostic_regions = {
        "methyl": (1.1, 1.3),
        "acetyl": (2.0, 2.1),
    }
    for label, region in diagnostic_regions.items():
        expected_values = [
            value for value in full_reference if region[0] <= value <= region[1]
        ]
        if not expected_values:
            continue
        observed_values = [
            value
            for field in observations["field_reports"]
            for value in field[f"{label}_reporters_ppm"]
        ]
        score = _bidirectional_match_score(observed_values, expected_values, 0.08)
        diagnostic_scores.append(float(score or 0.0))
        diagnostic_checks.append({
            "reporter": label,
            "reference_ppm": expected_values,
            "observed_ppm": observed_values,
            "score": score,
        })

    available_scores = [
        value
        for value in [
            count_score,
            j_score,
            float(np.mean(diagnostic_scores)) if diagnostic_scores else None,
        ]
        if value is not None
    ]
    if not available_scores:
        return None
    return {
        "reference": BUBB_REFERENCE,
        "profile": profile_name,
        "expected_model": (profile or {}).get("expected_model"),
        "structural_reporter_groups": (
            (profile or {}).get("structural_reporter_groups", ["anomeric H1"])
        ),
        "visible_anomeric_reference_count": len(visible_reference),
        "observed_anomeric_counts": field_counts,
        "anomeric_count_score": count_score,
        "anomeric_j_score": j_score,
        "anomeric_j_checks": j_checks,
        "consensus_anomeric_spacings": observations[
            "consensus_anomeric_spacings"
        ],
        "diagnostic_reporter_checks": diagnostic_checks,
        "score": float(np.mean(available_scores)),
        "warning": (
            "Bubb guidance constrains candidate interpretation. It does not "
            "supply molecule-specific shifts, prove assignments, or convert "
            "crowded-region line separations into signed J values."
        ),
    }


def _matrix_coupling_summary(components: list[dict[str, Any]]) -> list[dict[str, Any]]:
    summary: list[dict[str, Any]] = []
    for component in components:
        matrix = np.asarray(component["matrix"], dtype=float)
        values = [
            abs(float(matrix[i, j]))
            for i in range(matrix.shape[0])
            for j in range(i + 1, matrix.shape[1])
            if 0.1 <= abs(float(matrix[i, j])) <= 25.0
        ]
        summary.append({
            "component": component["name"],
            "matrix_file": str(component["matrix_file"]),
            "scalar_couplings_hz": sorted(round(value, 6) for value in values),
        })
    return summary


def score_physics_model(
    candidate: dict[str, Any],
    fields: list[dict[str, Any]],
    repo_root: Path,
    *,
    region: tuple[float, float] = FULL_PROTON_REGION,
    water_region: tuple[float, float] = WATER_REGION,
) -> dict[str, Any] | None:
    """Compare observed spectra with exact matrix-derived multiplet shapes.

    Chemical shifts and every nonzero scalar coupling in the candidate matrix
    enter the Hamiltonian. Only linewidth and a small global ppm offset are
    searched. The cosine similarity is screening evidence, not a J assignment.
    """

    relative = candidate.get("model_config_file")
    if not relative:
        return None
    config_path = repo_root / str(relative)
    if not config_path.is_file():
        return None

    common_path = repo_root / "src" / "common"
    physics_path = repo_root / "src" / "sucrose" / "bayes_astro"
    for path in (common_path, physics_path):
        if str(path) not in sys.path:
            sys.path.insert(0, str(path))
    from carbohydrate_model import component_specs  # noqa: PLC0415
    from sucrose_sim import lorentzian_spectrum, sucrose_sticks  # noqa: PLC0415

    config = json.loads(config_path.read_text(encoding="utf-8"))
    components = component_specs(config, repo_root)
    field_results: list[dict[str, Any]] = []
    for field in fields:
        ppm = np.asarray(field["ppm"], dtype=float)
        intensity = np.asarray(field["intensity"], dtype=float)
        mask = (ppm >= region[0]) & (ppm <= region[1])
        mask &= ~((ppm >= water_region[0]) & (ppm <= water_region[1]))
        x, observed = ppm[mask], np.maximum(intensity[mask], 0.0)
        if x.size > 2500:
            stride = int(math.ceil(x.size / 2500))
            x, observed = x[::stride], observed[::stride]
        order = np.argsort(x)
        x, observed = x[order], observed[order]
        observed_norm = float(np.linalg.norm(observed))
        if x.size < 3 or observed_norm <= 0:
            field_results.append({
                "field_mhz": field["field_mhz"],
                "multiplet_shape_score": None,
                "reason": "no positive signal in the selected proton window",
            })
            continue
        observed = observed / observed_norm

        sticks = []
        for component in components:
            frequencies, amplitudes = sucrose_sticks(
                component["matrix"],
                float(field["field_mhz"]),
                4.50,
                blocks=component["blocks"],
            )
            sticks.append((component, frequencies, amplitudes))

        trial_widths = (1.0, 1.5, 2.0, 2.5, 3.0, 4.0)
        width_profiles: set[tuple[float, ...]] = set()
        best = {"score": -1.0, "offset_ppm": 0.0, "linewidths_hz": []}
        for trial_width in trial_widths:
            linewidths = tuple(
                float(component.get("linewidth_hz", trial_width))
                for component, _, _ in sticks
            )
            if linewidths in width_profiles:
                continue
            width_profiles.add(linewidths)
            simulated = np.zeros_like(x)
            for linewidth, (component, frequencies, amplitudes) in zip(linewidths, sticks):
                simulated += float(component["fraction"]) * lorentzian_spectrum(
                    x,
                    frequencies,
                    amplitudes,
                    linewidth,
                    float(field["field_mhz"]),
                )
            for offset in np.linspace(-0.03, 0.03, 25):
                shifted = np.interp(x, x + offset, simulated, left=0.0, right=0.0)
                norm = float(np.linalg.norm(shifted))
                score = float(observed @ (shifted / norm)) if norm > 0 else 0.0
                if score > best["score"]:
                    best = {
                        "score": score,
                        "offset_ppm": float(offset),
                        "linewidths_hz": list(linewidths),
                    }
        field_results.append({
            "field_mhz": field["field_mhz"],
            "multiplet_shape_score": best["score"],
            "best_offset_ppm": best["offset_ppm"],
            "best_component_linewidths_hz": best["linewidths_hz"],
            "points_compared": int(x.size),
        })

    scores = [
        float(item["multiplet_shape_score"])
        for item in field_results
        if item.get("multiplet_shape_score") is not None
    ]
    return {
        "model_config_file": str(relative),
        "spin_matrix_status": candidate.get("spin_matrix_status", "unspecified"),
        "matrix_couplings": _matrix_coupling_summary(components),
        "field_results": field_results,
        "mean_multiplet_shape_score": float(np.mean(scores)) if scores else None,
        "interpretation": (
            "Exact full-window forward simulation; all matrix J values affect "
            "the multiplet shape, but no individual line spacing is promoted "
            "to a signed scalar-coupling assignment."
        ),
    }


def score_candidate(
    candidate: dict[str, Any],
    fields: list[dict[str, Any]],
    repo_root: Path,
    *,
    enable_physics: bool = True,
    enable_bubb: bool = True,
    bubb_observations: dict[str, Any] | None = None,
) -> dict[str, Any]:
    anomeric_reference = _reference_anomeric_centers(repo_root, candidate)
    full_reference = _reference_full_proton_centers(repo_root, candidate)
    physics = score_physics_model(candidate, fields, repo_root) if enable_physics else None
    bubb = (
        score_bubb_guidance(
            candidate,
            fields,
            repo_root,
            observations=bubb_observations,
        )
        if enable_bubb
        else None
    )
    physics_by_field = {
        float(item["field_mhz"]): item
        for item in (physics or {}).get("field_results", [])
    }
    field_results: list[dict[str, Any]] = []
    chemical_scores: list[float] = []
    shape_scores: list[float] = []
    for field in fields:
        observed_anomeric = detect_anomeric_clusters(field["ppm"], field["intensity"])
        observed_full = detect_full_proton_clusters(field["ppm"], field["intensity"])
        expected_count = (
            len(anomeric_reference)
            if anomeric_reference
            else candidate.get("expected_anomeric_clusters")
        )
        count_score = None
        if expected_count is not None:
            count_score = math.exp(-0.8 * abs(len(observed_anomeric) - int(expected_count)))
        anomeric_shift_score = _match_score(observed_anomeric, anomeric_reference)
        if anomeric_shift_score is None and count_score is None:
            anomeric_score = None
        elif anomeric_shift_score is None:
            anomeric_score = float(count_score)
        elif count_score is None:
            anomeric_score = float(anomeric_shift_score)
        else:
            anomeric_score = 0.70 * anomeric_shift_score + 0.30 * count_score

        full_shift_score = _bidirectional_match_score(observed_full, full_reference)
        if anomeric_score is not None and full_shift_score is not None:
            chemical_score = 0.35 * anomeric_score + 0.65 * full_shift_score
        elif full_shift_score is not None:
            chemical_score = float(full_shift_score)
        else:
            chemical_score = anomeric_score

        physics_field = physics_by_field.get(float(field["field_mhz"]), {})
        shape_score = physics_field.get("multiplet_shape_score")
        if chemical_score is not None and shape_score is not None:
            score = 0.60 * float(chemical_score) + 0.40 * float(shape_score)
        elif chemical_score is not None:
            score = float(chemical_score)
        elif shape_score is not None:
            score = float(shape_score)
        else:
            score = 0.0
        if chemical_score is not None:
            chemical_scores.append(float(chemical_score))
        if shape_score is not None:
            shape_scores.append(float(shape_score))
        field_results.append({
            "field_mhz": field["field_mhz"],
            "observed_anomeric_clusters_ppm": observed_anomeric,
            "reference_anomeric_clusters_ppm": anomeric_reference,
            "anomeric_shift_match_score": anomeric_shift_score,
            "anomeric_count_score": count_score,
            "anomeric_score": anomeric_score,
            "observed_full_proton_clusters_ppm": observed_full,
            "reference_full_proton_clusters_ppm": full_reference,
            "full_shift_match_score": full_shift_score,
            "chemical_shift_score": chemical_score,
            "multiplet_shape_score": shape_score,
            "physics_fit": physics_field or None,
            "available_field_evidence_score": score,
        })
    mean_chemical = float(np.mean(chemical_scores)) if chemical_scores else None
    mean_shape = float(np.mean(shape_scores)) if shape_scores else None
    bubb_score = float(bubb["score"]) if bubb is not None else None
    if mean_chemical is None:
        mean_score = 0.0
    else:
        # Missing evidence is neutral rather than being confused with failure.
        # This keeps shift-only BMRB candidates comparable with candidates that
        # have matrix and Bubb evidence without rewarding unavailable channels.
        mean_score = (
            0.70 * mean_chemical
            + 0.20 * (mean_shape if mean_shape is not None else 0.50)
            + 0.10 * (bubb_score if bubb_score is not None else 0.50)
        )
    reference_available = bool(anomeric_reference or full_reference)
    return {
        "candidate_id": candidate["id"],
        "name": candidate["name"],
        "class": candidate["class"],
        "bmrb_entry": candidate.get("bmrb_entry"),
        "forms": candidate.get("forms", []),
        "topology": candidate.get("topology"),
        "notes": candidate.get("notes"),
        "review_status": candidate.get("review_status", "approved"),
        "identity_warning": candidate.get("identity_warning", ""),
        "bmrb_evidence": candidate.get("bmrb_evidence", {}),
        "reference_available": reference_available,
        "full_proton_reference_available": bool(full_reference),
        "physics_model_available": physics is not None,
        "bubb_guidance_available": bubb is not None,
        "mean_chemical_shift_score": mean_chemical,
        "mean_multiplet_shape_score": mean_shape,
        "bubb_guidance_score": bubb_score,
        "bubb_guidance": bubb,
        "physics_model": physics,
        "mean_score": mean_score,
        "field_results": field_results,
    }


def rank_candidates(
    fields: list[dict[str, Any]],
    repo_root: Path,
    library: list[dict[str, Any]] | None = None,
    *,
    enable_physics: bool = True,
    enable_bubb: bool = True,
) -> list[dict[str, Any]]:
    bubb_observations = build_bubb_observations(fields) if enable_bubb else None
    results = [
        score_candidate(
            candidate,
            fields,
            repo_root,
            enable_physics=enable_physics,
            enable_bubb=enable_bubb,
            bubb_observations=bubb_observations,
        )
        for candidate in (library or load_library())
    ]
    return sorted(results, key=lambda item: item["mean_score"], reverse=True)


def _direct_j_spacing_summary(repo_root: Path, molecule: str) -> dict[str, Any] | None:
    path = (
        repo_root
        / "outputs"
        / molecule
        / "j_measurements"
        / "j_candidate_measurements.json"
    )
    if not path.is_file():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    grouped: dict[tuple[str, str, float], list[dict[str, Any]]] = {}
    for item in payload.get("coupling_screening", []):
        key = (
            str(item["i"]),
            str(item["j"]),
            abs(float(item.get("provisional_j_hz", 0.0))),
        )
        options = [
            check
            for check in item.get("screening", [])
            if check.get("nearest_spacing_hz") is not None
        ]
        if not options:
            continue
        best = min(options, key=lambda check: float(check["absolute_difference_hz"]))
        grouped.setdefault(key, []).append({
            "field_mhz": item["field_mhz"],
            "nearest_spacing_hz": float(best["nearest_spacing_hz"]),
            "absolute_difference_hz": float(best["absolute_difference_hz"]),
        })
    comparisons = []
    for (left, right, provisional), checks in grouped.items():
        comparisons.append({
            "i": left,
            "j": right,
            "provisional_j_hz": provisional,
            "fields_with_resolved_spacing": len(checks),
            "fields_within_1_hz": sum(
                check["absolute_difference_hz"] <= 1.0 for check in checks
            ),
            "field_checks": checks,
        })
    return {
        "status": payload.get("status", "REVIEW"),
        "report_file": str(path.relative_to(repo_root)),
        "method": payload.get("method"),
        "provisional_seed_warning": (
            "These windows come from the mystery sample's provisional "
            "xylose-like atom centers. They are useful only after the candidate "
            "ranking and do not independently establish identity."
        ),
        "coupling_comparisons": comparisons,
        "interpretation": payload.get("interpretation"),
        "matrix_updated": bool(payload.get("matrix_updated", False)),
    }


def build_report(repo_root: Path, molecule: str, ranked: list[dict[str, Any]]) -> dict[str, Any]:
    if not ranked:
        raise ValueError("No candidates were supplied")
    top = ranked[0]
    second = ranked[1] if len(ranked) > 1 else None
    margin = float(top["mean_score"] - second["mean_score"]) if second else 1.0
    promotable = bool(
        top["reference_available"]
        and top.get("review_status", "approved") == "approved"
        and top["mean_score"] >= 0.75
        and margin >= 0.10
    )
    status = "CANDIDATE" if promotable else "REVIEW"
    return {
        "molecule": molecule,
        "stage": "complete_proton_multifield_1d_screen",
        "proton_window_ppm": list(FULL_PROTON_REGION),
        "excluded_water_window_ppm": list(WATER_REGION),
        "scoring": {
            "chemical_shift_fingerprint": (
                "35% anomeric evidence plus 65% complete-window proton-shift match"
            ),
            "combined_candidate_score": (
                "70% chemical-shift fingerprint, 20% exact matrix multiplet "
                "shape, and 10% Bubb guidance; unavailable optional channels "
                "receive a neutral 0.5 rather than being treated as failure"
            ),
        },
        "status": status,
        "identity_claim": top["name"] if promotable else "unknown carbohydrate",
        "top_candidate": top["name"],
        "topology_hypothesis": top["topology"],
        "top_review_status": top.get("review_status", "approved"),
        "top_score": top["mean_score"],
        "runner_up": second["name"] if second else None,
        "score_margin": margin,
        "identity_confirmation_required": True,
        "confirmation_options": ["verified GISSMO/BMRB matrix", "COSY or TOCSY", "HSQC/HMBC", "matched authentic standard"],
        "direct_j_spacing_evidence": _direct_j_spacing_summary(repo_root, molecule),
        "ranked_candidates": ranked,
        "warning": (
            "Complete-window 1-D multifield screening ranks hypotheses; matrix "
            "simulation tests multiplet/J consistency but does not prove atom "
            "assignments, signed couplings, or molecular identity."
        ),
    }


def write_outputs(repo_root: Path, molecule: str, report: dict[str, Any]) -> tuple[Path, Path]:
    output_dir = repo_root / "outputs" / molecule / "identification"
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "identification_1d_report.json"
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    md_path = output_dir / "identification_1d_report.md"
    lines = [
        f"# Complete 1-D mystery-sugar screen: {molecule}",
        "",
        f"Status: **{report['status']}**",
        f"Top candidate: **{report['top_candidate']}**",
        f"Top score: {report['top_score']:.3f}",
        f"Runner-up: {report['runner_up'] or 'none'}",
        f"Score margin: {report['score_margin']:.3f}",
        "",
        (
            f"Window scored: {report['proton_window_ppm'][0]:.2f}–"
            f"{report['proton_window_ppm'][1]:.2f} ppm; water excluded at "
            f"{report['excluded_water_window_ppm'][0]:.2f}–"
            f"{report['excluded_water_window_ppm'][1]:.2f} ppm."
        ),
        "",
        (
            "This compares the complete proton-shift fingerprint. Where a spin "
            "matrix is available, exact forward simulation also tests multiplet "
            "shapes and all matrix scalar couplings. Bubb reporter-group and "
            "cross-field anomeric-spacing guidance is reported separately."
        ),
        "",
        "This is a multifield 1-D hypothesis ranking, not an identity confirmation.",
        "",
        "| Rank | Candidate | Combined | Shift fingerprint | Multiplet/J shape | Bubb | Matrix |",
        "|---:|---|---:|---:|---:|---:|---|",
    ]
    for rank, item in enumerate(report["ranked_candidates"], 1):
        shift = item.get("mean_chemical_shift_score")
        shape = item.get("mean_multiplet_shape_score")
        bubb = item.get("bubb_guidance_score")
        matrix_status = (
            item["physics_model"]["spin_matrix_status"]
            if item.get("physics_model")
            else "unavailable"
        )
        lines.append(
            f"| {rank} | {item['name']} | {item['mean_score']:.3f} | "
            f"{shift:.3f} | " if shift is not None else
            f"| {rank} | {item['name']} | {item['mean_score']:.3f} | n/a | "
        )
        lines[-1] += (
            f"{shape:.3f} | " if shape is not None else "n/a | "
        )
        lines[-1] += (
            f"{bubb:.3f} | {matrix_status} |"
            if bubb is not None
            else f"n/a | {matrix_status} |"
        )
    lines.extend([
        "",
        (
            "Scalar-coupling note: a matrix-derived shape score tests all J values "
            "together. It does not convert an unresolved 1-D line spacing into an "
            "automatic signed J assignment."
        ),
        "",
        "Required next evidence: " + ", ".join(report["confirmation_options"]) + ".",
    ])
    top_physics = report["ranked_candidates"][0].get("physics_model")
    if top_physics:
        lines.extend(["", "## Top-candidate matrix couplings", ""])
        for component in top_physics["matrix_couplings"]:
            values = ", ".join(
                f"{value:.3f}" for value in component["scalar_couplings_hz"]
            )
            lines.append(f"- {component['component']}: {values} Hz")
    top_bubb = report["ranked_candidates"][0].get("bubb_guidance")
    if top_bubb:
        lines.extend(["", "## Top-candidate Bubb guidance", ""])
        lines.append(
            f"- Bubb score: {top_bubb['score']:.3f}; profile: "
            f"{top_bubb['profile'] or 'general reporter rules'}."
        )
        lines.append(
            "- Observed anomeric counts by field: "
            + ", ".join(str(value) for value in top_bubb["observed_anomeric_counts"])
            + "."
        )
        for check in top_bubb["anomeric_j_checks"]:
            lines.append(
                f"- {check['form']}: {check['observed_spacing_hz']:.3f} Hz "
                f"across {check['field_support']} fields; expected "
                f"{check['expected_range_hz'][0]:.1f}-"
                f"{check['expected_range_hz'][1]:.1f} Hz."
            )
        lines.append(f"- Caution: {top_bubb['warning']}")
    direct_j = report.get("direct_j_spacing_evidence")
    if direct_j:
        lines.extend([
            "",
            "## Direct resolved-line-spacing screen",
            "",
            f"Status: **{direct_j['status']}**",
            "",
            f"Report: `{direct_j['report_file']}`",
            "",
            direct_j["provisional_seed_warning"],
        ])
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, md_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--molecule", required=True)
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--library", type=Path)
    parser.add_argument(
        "--include-bmrb-catalog",
        action="store_true",
        help="Include review-gated candidates synced from explore_BMRB",
    )
    parser.add_argument(
        "--bmrb-catalog",
        type=Path,
        default=DEFAULT_BMRB_CATALOG,
        help=f"Synced explore_BMRB catalog (default: {DEFAULT_BMRB_CATALOG})",
    )
    parser.add_argument(
        "--no-physics-model",
        action="store_true",
        help="Skip exact matrix-derived multiplet/J forward simulations",
    )
    parser.add_argument(
        "--no-bubb-guidance",
        action="store_true",
        help="Skip Bubb reporter-group and anomeric-spacing guidance",
    )
    args = parser.parse_args()
    root = args.repo_root.resolve()
    fields = load_prepared(root, args.molecule)
    library = load_library(
        args.library,
        include_bmrb_catalog=args.include_bmrb_catalog,
        bmrb_catalog_path=args.bmrb_catalog,
    )
    ranked = rank_candidates(
        fields,
        root,
        library,
        enable_physics=not args.no_physics_model,
        enable_bubb=not args.no_bubb_guidance,
    )
    report = build_report(root, args.molecule, ranked)
    json_path, md_path = write_outputs(root, args.molecule, report)
    print(f"1-D identity screen: {report['status']}")
    print(f"Top candidate: {report['top_candidate']} (score={report['top_score']:.3f})")
    top = report["ranked_candidates"][0]
    print(
        "Top evidence: "
        f"full shifts={top.get('mean_chemical_shift_score')}, "
        f"multiplet/J shape={top.get('mean_multiplet_shape_score')}, "
        f"Bubb={top.get('bubb_guidance_score')}"
    )
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
