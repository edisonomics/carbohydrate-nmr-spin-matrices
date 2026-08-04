#!/usr/bin/env python3
"""Machine-readable carbohydrate NMR guidance distilled from Bubb (2003).

This module is deliberately a *guardrail*, not a structure predictor.  Bubb's
review tells us what to look for (anomeric reporter groups, coupling patterns,
crowded proton regions, and the need for 2-D confirmation); BMRB/GISSMO or
student spectra still supply the molecule-specific numbers.
"""

from __future__ import annotations

from typing import Any


BUBB_REFERENCE = {
    "citation": "W. A. Bubb, NMR Spectroscopy in the Study of Carbohydrates (2003)",
    "doi": "10.1002/cmr.a.10080",
    "role": "chemical-interpretation-guidance",
}

# Companion sources are recorded as roles so the workflow can tell students
# why a source is being used.  They are not treated as interchangeable numeric
# matrix sources.
RELATED_REFERENCES = [
    {
        "key": "duus_2000",
        "citation": "Duus, Gotfredsen & Bock, Carbohydrate structural determination by NMR spectroscopy (2000)",
        "doi": "10.1021/cr990302n",
        "role": "2D assignment strategy and limitations",
    },
    {
        "key": "karplus_carbohydrates_2009",
        "citation": "Developments in the Karplus equation as they relate to carbohydrate coupling constants (2009)",
        "doi": "10.1016/S0065-2318(09)00003-1",
        "role": "stereochemistry-informed coupling priors",
    },
    {
        "key": "casper",
        "citation": "Dorst & Widmalm, NMR chemical shift prediction and structural elucidation using CASPER (2023)",
        "doi": "10.1016/j.carres.2023.108937",
        "role": "optional structure/chemical-shift cross-check",
    },
    {
        "key": "gissmo",
        "citation": "Dashti et al., Spin System Modeling for Metabolomics (2017); Applications of Parametrized Spin Systems (2018)",
        "doi": "10.1021/acs.analchem.7b02884; 10.1021/acs.analchem.8b02660",
        "role": "numeric spin-matrix seed and multifield simulation",
    },
]


EVIDENCE_POLICY = {
    "numeric_authority": [
        "verified GISSMO/BMRB matrix",
        "BMRB assigned chemical shifts",
        "student-collected COSY/TOCSY/HSQC/HMBC evidence",
        "resolved multifield experimental peaks",
    ],
    "interpretive_guidance": [
        "Bubb carbohydrate NMR rules",
        "Duus carbohydrate 2-D assignment guidance",
        "Karplus carbohydrate coupling relationships",
        "CASPER structure/shift cross-check",
    ],
    "context_only": [
        "old exploratory fits",
        "single-field refinements without holdout validation",
        "unverified literature assignments",
    ],
    "rule": "Only numeric-authority evidence may change a final matrix; guidance can constrain or flag it, and context-only work cannot change it automatically.",
}


PROFILES: dict[str, dict[str, Any]] = {
    "glucose": {
        "class": "reducing_aldose",
        "expected_model": "anomer_mixture",
        "minimum_components": 2,
        "structural_reporter_groups": ["H1"],
        "anomeric_j_patterns_hz": [
            {"form": "alpha", "range_hz": [2.0, 4.0]},
            {"form": "beta", "range_hz": [7.0, 9.0]},
        ],
        "possible_forms": ["alpha_pyranose", "beta_pyranose", "minor_furanose", "open_chain"],
        "notes": "Aqueous glucose commonly contains distinguishable alpha and beta forms; minor forms should be considered when residuals require them.",
    },
    "xylose": {
        "class": "reducing_aldopentose",
        "expected_model": "anomer_mixture",
        "minimum_components": 2,
        "structural_reporter_groups": ["H1"],
        "anomeric_j_patterns_hz": [
            {"form": "alpha", "range_hz": [2.0, 4.0]},
            {"form": "beta", "range_hz": [7.0, 9.0]},
        ],
        "possible_forms": ["alpha_pyranose", "beta_pyranose", "furanose", "open_chain"],
        "notes": "Aldopentoses can have appreciable furanose populations; an alpha/beta-only model is a provisional starting point.",
    },
    "mannose": {
        "class": "reducing_aldose",
        "expected_model": "anomer_mixture",
        "minimum_components": 2,
        "structural_reporter_groups": ["H1", "H2", "H3"],
        "anomeric_j_patterns_hz": [
            {"form": "alpha", "range_hz": [1.1, 2.1], "typical_hz": 1.6},
            {"form": "beta", "range_hz": [0.4, 1.2], "typical_hz": 0.8},
        ],
        "possible_forms": ["alpha_pyranose", "beta_pyranose", "minor_furanose", "open_chain"],
        "notes": "Treat alpha and beta forms as separate species unless the data demonstrate otherwise.",
    },
    "sucrose": {
        "class": "nonreducing_disaccharide",
        "expected_model": "single_molecule",
        "minimum_components": 1,
        "structural_reporter_groups": ["glucosyl H1"],
        "anomeric_j_patterns_hz": [
            {"form": "fixed_alpha_glucosyl", "range_hz": [2.0, 4.0]},
        ],
        "possible_forms": ["fixed_glucose_anomer", "fixed_fructose_anomer"],
        "notes": "Sucrose has no free anomeric center for ordinary mutarotation; use one linked disaccharide model, not an alpha/beta mixture.",
    },
}


COMMON_CHECKS = {
    "anomeric_reporter": True,
    "anomeric_proton_region_ppm": [4.4, 5.5],
    "proton_crowding_ppm": [3.4, 4.0],
    "acetyl_proton_region_ppm": [2.0, 2.1],
    "methyl_proton_region_ppm": [1.1, 1.3],
    "use_2d_for_assignment": True,
    "warn_1d_splittings_are_not_automatically_J": True,
    "j12_alpha_hz": [2.0, 4.0],
    "j12_beta_hz": [7.0, 9.0],
    "mannose_j12_alpha_typical_hz": 1.6,
    "mannose_j12_beta_typical_hz": 0.8,
}


def profile_for(molecule: str, config: dict[str, Any] | None = None) -> dict[str, Any]:
    """Return an explicit profile, allowing a config override for new sugars."""
    config = config or {}
    chemistry = config.get("chemistry", {})
    profile_name = str(chemistry.get("bubb_profile", molecule)).lower()
    profile = dict(PROFILES.get(profile_name, {
        "class": "unknown_carbohydrate",
        "expected_model": "assignment_required",
        "minimum_components": 1,
        "possible_forms": ["unknown"],
        "notes": "Use Bubb reporter/coupling guidance and 2-D evidence; do not invent a component model.",
    }))
    profile["molecule"] = molecule
    profile["bubb_profile"] = profile_name
    profile["guidance"] = dict(COMMON_CHECKS)
    profile["reference"] = dict(BUBB_REFERENCE)
    profile["companion_references"] = [dict(item) for item in RELATED_REFERENCES]
    profile["evidence_policy"] = dict(EVIDENCE_POLICY)
    return profile


def assess_config(molecule: str, config: dict[str, Any], *, bmrb: dict[str, Any] | None = None) -> dict[str, Any]:
    """Check whether the configured physical model matches Bubb expectations."""
    profile = profile_for(molecule, config)
    components = config.get("components") or []
    model_type = "mixture" if components else "single"
    warnings: list[str] = []
    checks: dict[str, Any] = {}
    expected = profile["expected_model"]
    if expected == "anomer_mixture":
        ok = model_type == "mixture" and len(components) >= profile["minimum_components"]
        checks["reducing_sugar_component_model"] = ok
        if not ok:
            warnings.append("Bubb expects separate anomer components; a single flattened matrix can hide mixture physics.")
    elif expected == "single_molecule":
        ok = model_type == "single"
        checks["nonreducing_single_model"] = ok
        if not ok:
            warnings.append("Bubb chemistry indicates a fixed nonreducing molecule; verify why multiple components were configured.")
    else:
        checks["component_model"] = bool(config.get("matrix_file") or components)
        if not checks["component_model"]:
            warnings.append("No numeric seed or component matrices are configured yet.")

    bmrb = bmrb or {}
    checks["bmrb_shifts_available"] = int(bmrb.get("proton_shift_count", 0)) > 0
    checks["bmrb_2d_available"] = bool(bmrb.get("two_d_experiments"))
    checks["gissmo_matrix_available"] = bool(bmrb.get("gissmo_matrix_file"))
    checks["numeric_seed_available"] = bool(
        checks["gissmo_matrix_available"] or config.get("matrix_file") or components
    )
    if not checks["bmrb_shifts_available"]:
        warnings.append("No BMRB proton-shift artifact is recorded; use assigned spectra or a documented provisional seed.")
    if profile["expected_model"] == "anomer_mixture" and not checks["bmrb_2d_available"]:
        warnings.append("No BMRB COSY/TOCSY/HSQC/HMBC inventory is recorded for this reducing-sugar mixture.")
    checks["bubb_guidance_recorded"] = True
    status = "PASS" if not warnings else "REVIEW"
    if checks["gissmo_matrix_available"]:
        seed_status = "GISSMO_VERIFIED"
    elif checks["bmrb_shifts_available"] and checks["numeric_seed_available"]:
        seed_status = "BMRB_PROVISIONAL"
    elif checks["numeric_seed_available"]:
        # A locally assigned/constructed seed has a different provenance from
        # a BMRB-derived seed, even when both are provisional.
        seed_status = "LOCAL_PROVISIONAL"
    else:
        seed_status = "NO_NUMERIC_SEED"
    return {
        "molecule": molecule,
        "status": status,
        "model_type": model_type,
        "seed_status": seed_status,
        "profile": profile,
        "checks": checks,
        "warnings": warnings,
        "next_step": (
            "Proceed to multifield fitting and transfer validation."
            if status == "PASS"
            else "Resolve the listed chemistry/model warnings before treating the fit as publishable."
        ),
    }
