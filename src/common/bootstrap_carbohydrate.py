#!/usr/bin/env python3
"""Validate reviewed carbohydrate inputs and generate provisional fit artifacts.

The workflow is intentionally semi-automatic.  A scientist supplies explicit
forms, proton assignments, coupling priors, and a multifield spectrum manifest.
This command checks those inputs, then deterministically writes provisional
matrices and a configuration in ``outputs/<molecule>/bootstrap``.  It never
overwrites a deposited or canonical matrix under ``data/``.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import re
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


SCHEMA_VERSION = 1
REQUIRED_FILES = (
    "molecule.json",
    "proton_assignments.csv",
    "couplings.csv",
    "spectra.csv",
)
PROTON_HEADERS = (
    "component",
    "spin_id",
    "atom_label",
    "attached_atom",
    "bmrb_label",
    "shift_ppm",
    "assignment_status",
    "source",
    "uncertainty_ppm",
    "fit",
)
COUPLING_HEADERS = (
    "component",
    "spin_i",
    "spin_j",
    "J_hz",
    "evidence_status",
    "source",
    "uncertainty_hz",
    "fit",
)
SPECTRUM_HEADERS = (
    "key",
    "field_mhz",
    "path",
    "procno",
    "role",
    "sample_id",
    "tube_id",
    "nucleus",
    "pulse_program",
)
EVIDENCE_STATUSES = {
    "DEPOSITED",
    "ASSIGNED",
    "LITERATURE_PRIOR",
    "ANALOG_PRIOR",
    "UNKNOWN",
}
FRACTION_STATUSES = {"MEASURED", "ASSIGNED", "PROVISIONAL", "UNKNOWN"}
ROLES = {"training", "validation"}
TRUE_VALUES = {"true", "yes", "1"}
FALSE_VALUES = {"false", "no", "0"}
SLUG = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


def _issue(items: list[dict[str, Any]], code: str, message: str, **context: Any) -> None:
    items.append({"code": code, "message": message, **context})


def _jcamp_value(path: Path, key: str) -> str | None:
    pattern = re.compile(
        rf"^##\${re.escape(key)}=\s*<?([^>\r\n]+)", re.MULTILINE
    )
    match = pattern.search(path.read_text(encoding="utf-8", errors="ignore"))
    return match.group(1).strip() if match else None


def _read_csv(path: Path, headers: tuple[str, ...]) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        actual = tuple(reader.fieldnames or ())
        missing = [header for header in headers if header not in actual]
        if missing:
            raise ValueError(f"{path.name} is missing columns: {', '.join(missing)}")
        rows = []
        for row_number, row in enumerate(reader, 2):
            if not any(str(value or "").strip() for value in row.values()):
                continue
            cleaned = {key: str(value or "").strip() for key, value in row.items()}
            cleaned["_row"] = row_number
            rows.append(cleaned)
        return rows


def load_bundle(input_dir: Path) -> dict[str, Any]:
    """Load the four-file reviewed-input bundle."""

    input_dir = input_dir.resolve()
    missing = [name for name in REQUIRED_FILES if not (input_dir / name).is_file()]
    if missing:
        raise FileNotFoundError(
            f"Missing input files in {input_dir}: {', '.join(missing)}"
        )
    molecule = json.loads((input_dir / "molecule.json").read_text(encoding="utf-8"))
    return {
        "input_dir": input_dir,
        "molecule": molecule,
        "protons": _read_csv(input_dir / "proton_assignments.csv", PROTON_HEADERS),
        "couplings": _read_csv(input_dir / "couplings.csv", COUPLING_HEADERS),
        "spectra": _read_csv(input_dir / "spectra.csv", SPECTRUM_HEADERS),
    }


def _as_float(
    value: Any,
    *,
    errors: list[dict[str, Any]],
    code: str,
    label: str,
    context: dict[str, Any],
) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        _issue(errors, code, f"{label} must be numeric, got {value!r}", **context)
        return None
    if not math.isfinite(number):
        _issue(errors, code, f"{label} must be finite, got {value!r}", **context)
        return None
    return number


def _as_bool(
    value: Any,
    *,
    errors: list[dict[str, Any]],
    code: str,
    label: str,
    context: dict[str, Any],
) -> bool | None:
    normalized = str(value).strip().lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False
    _issue(errors, code, f"{label} must be true or false, got {value!r}", **context)
    return None


def validate_bundle(bundle: dict[str, Any], repo_root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    """Validate and normalize a bundle without guessing missing science."""

    normalized = copy.deepcopy(bundle)
    molecule = normalized["molecule"]
    protons = normalized["protons"]
    couplings = normalized["couplings"]
    spectra = normalized["spectra"]
    errors: list[dict[str, Any]] = []
    warnings: list[dict[str, Any]] = []

    if molecule.get("schema_version") != SCHEMA_VERSION:
        _issue(
            errors,
            "molecule.schema_version",
            f"molecule.json must declare schema_version={SCHEMA_VERSION}",
        )
    molecule_id = str(molecule.get("molecule_id", "")).strip()
    if not molecule_id or not SLUG.fullmatch(molecule_id):
        _issue(
            errors,
            "molecule.id",
            "molecule_id must contain only letters, numbers, underscores, or hyphens",
        )
    if not str(molecule.get("name", "")).strip():
        _issue(errors, "molecule.name", "molecule.json must provide a display name")

    conditions = molecule.get("conditions")
    if not isinstance(conditions, dict):
        conditions = {}
        _issue(errors, "conditions.missing", "molecule.json must provide conditions")
    for key in ("solvent", "chemical_shift_reference"):
        if not str(conditions.get(key, "")).strip():
            _issue(errors, f"conditions.{key}", f"conditions.{key} is required")
    temperature = _as_float(
        conditions.get("temperature_K"),
        errors=errors,
        code="conditions.temperature",
        label="conditions.temperature_K",
        context={},
    )
    if temperature is not None and not 150.0 <= temperature <= 400.0:
        _issue(errors, "conditions.temperature.range", "temperature_K must be between 150 and 400 K")
    ph = _as_float(
        conditions.get("pH"),
        errors=errors,
        code="conditions.pH",
        label="conditions.pH",
        context={},
    )
    if ph is not None and not 0.0 <= ph <= 14.0:
        _issue(errors, "conditions.pH.range", "pH must be between 0 and 14")
    conditions["temperature_K"] = temperature
    conditions["pH"] = ph
    molecule["conditions"] = conditions

    forms = molecule.get("forms")
    if not isinstance(forms, list) or not forms:
        forms = []
        _issue(errors, "forms.empty", "At least one explicit solution form is required")
    form_ids: list[str] = []
    fractions: list[float] = []
    for index, form in enumerate(forms, 1):
        context = {"form_index": index}
        if not isinstance(form, dict):
            _issue(errors, "forms.type", "Each form must be an object", **context)
            continue
        form_id = str(form.get("id", "")).strip()
        if not form_id or not SLUG.fullmatch(form_id):
            _issue(errors, "forms.id", "Form id must be a safe nonempty identifier", **context)
        form_ids.append(form_id)
        fraction = _as_float(
            form.get("fraction"),
            errors=errors,
            code="forms.fraction",
            label="form fraction",
            context=context,
        )
        if fraction is not None:
            if not 0.0 <= fraction <= 1.0:
                _issue(errors, "forms.fraction.range", "Form fraction must be between 0 and 1", **context)
            fractions.append(fraction)
            form["fraction"] = fraction
        status = str(form.get("fraction_status", "")).strip().upper()
        if status not in FRACTION_STATUSES:
            _issue(
                errors,
                "forms.fraction_status",
                f"fraction_status must be one of {sorted(FRACTION_STATUSES)}",
                **context,
            )
        form["fraction_status"] = status
    if len(set(form_ids)) != len(form_ids):
        _issue(errors, "forms.duplicate", "Form ids must be unique")
    if fractions and len(fractions) == len(forms) and not math.isclose(sum(fractions), 1.0, abs_tol=1e-3):
        _issue(errors, "forms.fraction.sum", "Form fractions must sum to 1.0")
    form_set = set(form_ids)

    if not protons:
        _issue(errors, "protons.empty", "proton_assignments.csv has no assignments")
    spin_ids: set[str] = set()
    component_spins: dict[str, set[str]] = {form: set() for form in form_ids}
    component_labels: set[tuple[str, str]] = set()
    bmrb_labels: Counter[tuple[str, str]] = Counter()
    for row in protons:
        row_number = int(row["_row"])
        context = {"file": "proton_assignments.csv", "row": row_number}
        component = row["component"]
        spin_id = row["spin_id"]
        atom_label = row["atom_label"]
        if component not in form_set:
            _issue(errors, "protons.component", f"Unknown component {component!r}", **context)
        if not spin_id or not SLUG.fullmatch(spin_id):
            _issue(errors, "protons.spin_id", "spin_id must be a safe nonempty identifier", **context)
        elif spin_id in spin_ids:
            _issue(errors, "protons.spin_id.duplicate", f"Duplicate spin_id {spin_id!r}", **context)
        else:
            spin_ids.add(spin_id)
            component_spins.setdefault(component, set()).add(spin_id)
        if not atom_label:
            _issue(errors, "protons.atom_label", "atom_label is required", **context)
        elif (component, atom_label) in component_labels:
            _issue(errors, "protons.atom_label.duplicate", f"Duplicate atom label {atom_label!r} in {component}", **context)
        else:
            component_labels.add((component, atom_label))
        if not row["attached_atom"]:
            _issue(errors, "protons.attached_atom", "attached_atom is required", **context)
        if row["bmrb_label"]:
            bmrb_labels[(component, row["bmrb_label"])] += 1
        shift = _as_float(
            row["shift_ppm"],
            errors=errors,
            code="protons.shift",
            label="shift_ppm",
            context=context,
        )
        if shift is not None and not -1.0 <= shift <= 15.0:
            _issue(errors, "protons.shift.range", "Proton shift must be between -1 and 15 ppm", **context)
        row["shift_ppm"] = shift
        uncertainty = _as_float(
            row["uncertainty_ppm"],
            errors=errors,
            code="protons.uncertainty",
            label="uncertainty_ppm",
            context=context,
        )
        if uncertainty is not None and uncertainty < 0:
            _issue(errors, "protons.uncertainty.range", "uncertainty_ppm cannot be negative", **context)
        row["uncertainty_ppm"] = uncertainty
        status = row["assignment_status"].upper()
        if status not in EVIDENCE_STATUSES:
            _issue(
                errors,
                "protons.status",
                f"assignment_status must be one of {sorted(EVIDENCE_STATUSES)}",
                **context,
            )
        row["assignment_status"] = status
        if status in {"UNKNOWN", "ANALOG_PRIOR"}:
            _issue(
                warnings,
                "protons.review",
                f"{spin_id} uses {status} shift evidence and requires review",
                **context,
            )
        if not row["source"]:
            _issue(errors, "protons.source", "Every shift requires a source", **context)
        row["fit"] = _as_bool(
            row["fit"],
            errors=errors,
            code="protons.fit",
            label="fit",
            context=context,
        )
    for (component, label), count in bmrb_labels.items():
        if count > 1:
            _issue(
                warnings,
                "protons.bmrb_label.duplicate",
                f"BMRB label {label!r} appears {count} times in {component}; preserve the ambiguity explicitly",
            )
    for form_id in form_ids:
        if not component_spins.get(form_id):
            _issue(errors, "protons.component.empty", f"No spins were supplied for component {form_id}")

    coupling_pairs: set[tuple[str, str, str]] = set()
    coupled_spins: dict[str, set[str]] = {form: set() for form in form_ids}
    for row in couplings:
        row_number = int(row["_row"])
        context = {"file": "couplings.csv", "row": row_number}
        component = row["component"]
        left, right = row["spin_i"], row["spin_j"]
        if component not in form_set:
            _issue(errors, "couplings.component", f"Unknown component {component!r}", **context)
        known = component_spins.get(component, set())
        if left not in known or right not in known:
            _issue(
                errors,
                "couplings.spin",
                f"Coupling {left!r}-{right!r} must reference spins in component {component!r}",
                **context,
            )
        if left == right:
            _issue(errors, "couplings.self", "A spin cannot be coupled to itself", **context)
        pair = (component, *sorted((left, right)))
        if pair in coupling_pairs:
            _issue(errors, "couplings.duplicate", f"Duplicate coupling pair {left!r}-{right!r}", **context)
        coupling_pairs.add(pair)
        coupled_spins.setdefault(component, set()).update((left, right))
        value = _as_float(
            row["J_hz"],
            errors=errors,
            code="couplings.value",
            label="J_hz",
            context=context,
        )
        if value is not None and not -30.0 <= value <= 30.0:
            _issue(errors, "couplings.value.range", "J_hz must be between -30 and 30 Hz", **context)
        row["J_hz"] = value
        uncertainty = _as_float(
            row["uncertainty_hz"],
            errors=errors,
            code="couplings.uncertainty",
            label="uncertainty_hz",
            context=context,
        )
        if uncertainty is not None and uncertainty < 0:
            _issue(errors, "couplings.uncertainty.range", "uncertainty_hz cannot be negative", **context)
        row["uncertainty_hz"] = uncertainty
        status = row["evidence_status"].upper()
        if status not in EVIDENCE_STATUSES:
            _issue(
                errors,
                "couplings.status",
                f"evidence_status must be one of {sorted(EVIDENCE_STATUSES)}",
                **context,
            )
        row["evidence_status"] = status
        if status in {"UNKNOWN", "ANALOG_PRIOR", "LITERATURE_PRIOR"}:
            _issue(
                warnings,
                "couplings.review",
                f"{left}-{right} uses {status} coupling evidence and remains provisional",
                **context,
            )
        if not row["source"]:
            _issue(errors, "couplings.source", "Every coupling requires a source", **context)
        row["fit"] = _as_bool(
            row["fit"],
            errors=errors,
            code="couplings.fit",
            label="fit",
            context=context,
        )
    if not couplings:
        _issue(warnings, "couplings.empty", "No scalar couplings were supplied; matrices will contain uncoupled spins")
    for form_id, spins in component_spins.items():
        isolated = sorted(spins - coupled_spins.get(form_id, set()))
        if isolated:
            _issue(
                warnings,
                "couplings.isolated_spins",
                f"Component {form_id} contains spins with no modeled couplings: {', '.join(isolated)}",
            )

    if not spectra:
        _issue(errors, "spectra.empty", "spectra.csv has no datasets")
    spectrum_keys: set[str] = set()
    training_fields: set[float] = set()
    validation_fields: set[float] = set()
    sample_ids: set[str] = set()
    tube_ids: set[str] = set()
    repo_root = repo_root.resolve()
    for row in spectra:
        row_number = int(row["_row"])
        context = {"file": "spectra.csv", "row": row_number}
        key = row["key"]
        if not key:
            _issue(errors, "spectra.key", "Dataset key is required", **context)
        elif key in spectrum_keys:
            _issue(errors, "spectra.key.duplicate", f"Duplicate dataset key {key!r}", **context)
        spectrum_keys.add(key)
        field = _as_float(
            row["field_mhz"],
            errors=errors,
            code="spectra.field",
            label="field_mhz",
            context=context,
        )
        if field is not None and not 40.0 <= field <= 1500.0:
            _issue(errors, "spectra.field.range", "field_mhz must be between 40 and 1500 MHz", **context)
        row["field_mhz"] = field
        role = row["role"].lower()
        if role not in ROLES:
            _issue(errors, "spectra.role", f"role must be one of {sorted(ROLES)}", **context)
        row["role"] = role
        if field is not None and role == "training":
            training_fields.add(field)
        if field is not None and role == "validation":
            validation_fields.add(field)
        path_text = row["path"]
        if not path_text:
            _issue(errors, "spectra.path", "Spectrum path is required", **context)
        else:
            path = Path(path_text)
            resolved = path if path.is_absolute() else repo_root / path
            if not resolved.is_dir():
                _issue(errors, "spectra.path.missing", f"Spectrum path does not exist: {resolved}", **context)
            else:
                procno = row["procno"]
                required_bruker = (
                    resolved / "acqus",
                    resolved / "pdata" / procno / "procs",
                    resolved / "pdata" / procno / "1r",
                )
                for required_path in required_bruker:
                    if not required_path.is_file():
                        _issue(
                            errors,
                            "spectra.bruker_file.missing",
                            f"Required processed Bruker file is missing: {required_path}",
                            **context,
                        )
                acqus_path = resolved / "acqus"
                if acqus_path.is_file():
                    actual_nucleus = (_jcamp_value(acqus_path, "NUC1") or "").upper()
                    if actual_nucleus and actual_nucleus not in {"1H", "H1", "PROTON"}:
                        _issue(
                            errors,
                            "spectra.metadata.nucleus",
                            f"Bruker acqus declares NUC1={actual_nucleus!r}, not 1H",
                            **context,
                        )
                    sfo1_text = _jcamp_value(acqus_path, "SFO1")
                    if sfo1_text and field is not None:
                        actual_field = _as_float(
                            sfo1_text,
                            errors=errors,
                            code="spectra.metadata.sfo1",
                            label="Bruker SFO1",
                            context=context,
                        )
                        tolerance = max(2.0, 0.01 * field)
                        if actual_field is not None and abs(actual_field - field) > tolerance:
                            _issue(
                                errors,
                                "spectra.metadata.field_mismatch",
                                (
                                    f"Manifest field {field:g} MHz does not match "
                                    f"Bruker SFO1 {actual_field:g} MHz"
                                ),
                                **context,
                            )
                    actual_program = _jcamp_value(acqus_path, "PULPROG")
                    if (
                        actual_program
                        and row["pulse_program"]
                        and actual_program.lower() != row["pulse_program"].lower()
                    ):
                        _issue(
                            errors,
                            "spectra.metadata.pulse_program_mismatch",
                            (
                                f"Manifest pulse program {row['pulse_program']!r} does not "
                                f"match Bruker PULPROG {actual_program!r}"
                            ),
                            **context,
                        )
        if not row["procno"]:
            _issue(errors, "spectra.procno", "procno is required", **context)
        if not row["sample_id"]:
            _issue(errors, "spectra.sample_id", "sample_id is required", **context)
        else:
            sample_ids.add(row["sample_id"])
        if not row["tube_id"]:
            _issue(errors, "spectra.tube_id", "tube_id is required", **context)
        else:
            tube_ids.add(row["tube_id"])
        if row["nucleus"] not in {"1H", "H1", "proton"}:
            _issue(errors, "spectra.nucleus", "This workflow currently requires 1H spectra", **context)
        row["nucleus"] = "1H"
        if not row["pulse_program"]:
            _issue(errors, "spectra.pulse_program", "pulse_program is required", **context)
    if len(training_fields) < 2:
        _issue(errors, "spectra.training", "At least two distinct training fields are required")
    if len(validation_fields) < 1:
        _issue(errors, "spectra.validation", "At least one independent validation field is required")
    if training_fields & validation_fields:
        _issue(errors, "spectra.field_overlap", "A field cannot be both training and validation")
    if len(sample_ids) > 1:
        _issue(errors, "spectra.sample_mismatch", "All fields must use the same sample_id")
    if len(tube_ids) > 1:
        _issue(errors, "spectra.tube_mismatch", "All fields must use the same tube_id")

    report = {
        "schema_version": SCHEMA_VERSION,
        "molecule_id": molecule_id or None,
        "status": "PASS" if not errors else "FAIL",
        "errors": errors,
        "warnings": warnings,
        "summary": {
            "forms": len(forms),
            "protons": len(protons),
            "couplings": len(couplings),
            "spectra": len(spectra),
            "training_fields_mhz": sorted(training_fields),
            "validation_fields_mhz": sorted(validation_fields),
            "sample_ids": sorted(sample_ids),
            "tube_ids": sorted(tube_ids),
        },
    }
    return normalized, report


def _connected_blocks(spin_ids: list[str], coupling_rows: list[dict[str, Any]]) -> list[list[int]]:
    adjacency = {spin_id: set() for spin_id in spin_ids}
    for row in coupling_rows:
        left, right = row["spin_i"], row["spin_j"]
        adjacency[left].add(right)
        adjacency[right].add(left)
    index = {spin_id: position + 1 for position, spin_id in enumerate(spin_ids)}
    unseen = set(spin_ids)
    blocks: list[list[int]] = []
    for start in spin_ids:
        if start not in unseen:
            continue
        stack = [start]
        unseen.remove(start)
        component = []
        while stack:
            current = stack.pop()
            component.append(index[current])
            for neighbor in sorted(adjacency[current]):
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        blocks.append(sorted(component))
    return blocks


def build_component_matrix(
    component: str,
    protons: list[dict[str, Any]],
    couplings: list[dict[str, Any]],
) -> tuple[list[str], list[list[float]], list[list[int]]]:
    component_protons = [row for row in protons if row["component"] == component]
    component_couplings = [row for row in couplings if row["component"] == component]
    spin_ids = [row["spin_id"] for row in component_protons]
    index = {spin_id: position for position, spin_id in enumerate(spin_ids)}
    matrix = [[0.0 for _ in spin_ids] for _ in spin_ids]
    for row in component_protons:
        matrix[index[row["spin_id"]]][index[row["spin_id"]]] = float(row["shift_ppm"])
    for row in component_couplings:
        left, right = index[row["spin_i"]], index[row["spin_j"]]
        value = float(row["J_hz"])
        matrix[left][right] = value
        matrix[right][left] = value
    return spin_ids, matrix, _connected_blocks(spin_ids, component_couplings)


def _write_matrix(path: Path, matrix: list[list[float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(" ".join(f"{value:.9f}" for value in row) for row in matrix) + "\n",
        encoding="utf-8",
    )


def _relative_or_absolute(path: Path, repo_root: Path) -> str:
    try:
        return str(path.resolve().relative_to(repo_root.resolve()))
    except ValueError:
        return str(path.resolve())


def _dataset_path_for_config(path_text: str, repo_root: Path, molecule_id: str) -> str:
    """Return a path compatible with the repository's molecule data loaders."""

    path = Path(path_text)
    resolved = path if path.is_absolute() else repo_root / path
    molecule_data = (repo_root / "data" / molecule_id).resolve()
    try:
        return str(resolved.resolve().relative_to(molecule_data))
    except ValueError:
        # Path joins used by the existing loaders preserve an absolute operand.
        return str(resolved.resolve())


def _input_hashes(input_dir: Path) -> dict[str, str]:
    return {
        name: hashlib.sha256((input_dir / name).read_bytes()).hexdigest()
        for name in REQUIRED_FILES
    }


def _default_config(molecule: dict[str, Any]) -> dict[str, Any]:
    processing = {
        "fit_region_ppm": [2.50, 6.20],
        "water_region_ppm": [4.65, 4.90],
        "artifact_region_ppm": [],
        "crowded_region_ppm": [3.40, 4.00],
        "anomeric_region_ppm": [4.30, 5.80],
        "dss_search_region_ppm": [-0.20, 0.20],
        "anomeric_reference_ppm": 5.40,
        "baseline": "median_fit_region",
        "normalization": "anomeric_peak",
    }
    processing.update(molecule.get("processing", {}))
    optimization = molecule.get("optimization", {})
    model = {
        "grid_points": 20000,
        "lb_hz": 1.0,
        "fit_stride": 3,
        "noise_sigma": 0.05,
        "shift_bound_ppm": 0.05,
        "coupling_bound_hz": 3.0,
        "linewidth_bounds_hz": [0.3, 5.0],
        "offset_bound_ppm": 0.03,
        "initial_linewidth_hz": 1.5,
        "max_nfev": 2000,
    }
    model.update(optimization)
    return {"processing": processing, "independent_model": model}


def generate_artifacts(
    bundle: dict[str, Any],
    report: dict[str, Any],
    repo_root: Path,
    output_dir: Path,
    *,
    force: bool = False,
) -> dict[str, Path]:
    """Generate review-only matrices, configuration, and provenance manifest."""

    if report["status"] != "PASS":
        raise ValueError("Cannot generate artifacts from a failing validation report")
    repo_root = repo_root.resolve()
    output_dir = output_dir.resolve()
    molecule = bundle["molecule"]
    molecule_id = molecule["molecule_id"]
    config_path = output_dir / f"{molecule_id}_config.provisional.json"
    manifest_path = output_dir / "seed_manifest.json"
    fit_path = output_dir / "fit_parameter_policy.json"
    expected_matrix_paths = [
        output_dir / "matrices" / f"{form['id']}_provisional_matrix.txt"
        for form in molecule["forms"]
    ]
    protected = [config_path, manifest_path, fit_path, *expected_matrix_paths]
    if not force and any(path.exists() for path in protected):
        raise FileExistsError(
            f"Bootstrap artifacts already exist in {output_dir}; use --force only after review"
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    matrix_dir = output_dir / "matrices"

    forms = molecule["forms"]
    component_records = []
    all_atom_ids: list[str] = []
    matrix_paths: list[Path] = []
    for form in forms:
        component = form["id"]
        spin_ids, matrix, blocks = build_component_matrix(
            component, bundle["protons"], bundle["couplings"]
        )
        matrix_path = matrix_dir / f"{component}_provisional_matrix.txt"
        if matrix_path.exists() and not force:
            raise FileExistsError(matrix_path)
        _write_matrix(matrix_path, matrix)
        matrix_paths.append(matrix_path)
        all_atom_ids.extend(spin_ids)
        component_records.append({
            "name": component,
            "fraction": float(form["fraction"]),
            "fraction_status": form["fraction_status"],
            "linewidth_hz": float(molecule.get("initial_linewidth_hz", 1.5)),
            "linewidth_provenance": "bootstrap nuisance prior; refine per field",
            "matrix_file": _relative_or_absolute(matrix_path, repo_root),
            "matrix_status": "PROVISIONAL_SEED",
            "atom_ids": spin_ids,
            "blocks": blocks,
        })

    datasets = []
    training_keys = []
    validation_keys = []
    for row in bundle["spectra"]:
        dataset = {
            "key": row["key"],
            "field_mhz": float(row["field_mhz"]),
            "relative_dir": _dataset_path_for_config(
                row["path"], repo_root, molecule_id
            ),
            "procno": row["procno"],
            "sample_id": row["sample_id"],
            "tube_id": row["tube_id"],
            "role": row["role"],
            "nucleus": "1H",
            "pulse_program": row["pulse_program"],
        }
        datasets.append(dataset)
        (training_keys if row["role"] == "training" else validation_keys).append(row["key"])

    defaults = _default_config(molecule)
    config: dict[str, Any] = {
        "name": molecule_id,
        "display_name": molecule["name"],
        "status": "PROVISIONAL_INPUT_REVIEW_REQUIRED",
        "conditions": molecule["conditions"],
        "chemistry": {
            **molecule.get("chemistry", {}),
            "forms": [form["id"] for form in forms],
        },
        "seed_selection": {
            "default_source": "reviewed_input_bundle",
            "require_provenance": True,
            "allow_provisional": True,
        },
        "j_measurement": {
            "required_for_new_matrix": True,
            "allowed_methods": ["resolved_1h_multiplet", "j_resolved_2d"],
            "matrix_update_requires_manual_review": True,
        },
        "model_type": "mixture" if len(forms) > 1 else "single",
        "atom_ids": all_atom_ids,
        "datasets": datasets,
        "processing": defaults["processing"],
        "independent_model": {
            **defaults["independent_model"],
            "fields_for_joint_fit": training_keys,
            "validation_fields": validation_keys,
        },
        "quality_gate": {
            "training_min_r": 0.90,
            "training_max_rmse": 0.10,
            "validation_min_r": 0.90,
            "validation_max_rmse": 0.10,
            "minimum_delta_r": 0.0,
            "minimum_delta_rmse": 0.0,
            "borderline_r": 0.85,
            "borderline_rmse": 0.15,
            "max_abs_offset_ppm": 0.03,
            "min_validation_fields": 1,
        },
        "publication": {
            "target": "GISSMO",
            "required_seed_status": "READY",
            "required_quality_status": "PASS",
            "bootstrap_artifacts_are_publishable": False,
        },
        "input_workflow": {
            "schema_version": SCHEMA_VERSION,
            "input_dir": str(bundle["input_dir"]),
            "validation_report": _relative_or_absolute(output_dir / "input_validation_report.json", repo_root),
            "seed_manifest": _relative_or_absolute(manifest_path, repo_root),
        },
    }
    if len(component_records) == 1:
        component = component_records[0]
        config.update({
            "matrix_file": component["matrix_file"],
            "blocks": component["blocks"],
            "component": component,
        })
    else:
        config.update({"matrix_file": None, "blocks": [], "components": component_records})

    fit_parameters = {
        "shifts": [
            {
                "component": row["component"],
                "spin_id": row["spin_id"],
                "fit": row["fit"],
                "uncertainty_ppm": row["uncertainty_ppm"],
                "status": row["assignment_status"],
            }
            for row in bundle["protons"]
        ],
        "couplings": [
            {
                "component": row["component"],
                "spin_i": row["spin_i"],
                "spin_j": row["spin_j"],
                "fit": row["fit"],
                "uncertainty_hz": row["uncertainty_hz"],
                "status": row["evidence_status"],
            }
            for row in bundle["couplings"]
        ],
    }
    fit_path.write_text(json.dumps(fit_parameters, indent=2) + "\n", encoding="utf-8")
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")

    status_counts = Counter(
        row["assignment_status"] for row in bundle["protons"]
    ) + Counter(row["evidence_status"] for row in bundle["couplings"])
    manifest = {
        "schema_version": SCHEMA_VERSION,
        "molecule": molecule_id,
        "status": "PROVISIONAL_SEED",
        "confidence": "provisional",
        "config_file": _relative_or_absolute(config_path, repo_root),
        "matrix_files": [_relative_or_absolute(path, repo_root) for path in matrix_paths],
        "fit_parameter_policy": _relative_or_absolute(fit_path, repo_root),
        "atom_ids": all_atom_ids,
        "input_hashes_sha256": _input_hashes(bundle["input_dir"]),
        "evidence_status_counts": dict(sorted(status_counts.items())),
        "validation_summary": report["summary"],
        "provenance_required": True,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "next_step": (
            "Review atom mappings and priors, measure independent J evidence, "
            "then refine on training fields and pass held-out multifield validation."
        ),
        "publication_block": (
            "Bootstrap matrices are not verified and must not replace a deposited "
            "GISSMO matrix or be published without review and validation."
        ),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return {
        "config": config_path,
        "manifest": manifest_path,
        "fit_policy": fit_path,
        **{f"matrix_{index + 1}": path for index, path in enumerate(matrix_paths)},
    }


def write_validation_report(report: dict[str, Any], output_dir: Path) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "input_validation_report.json"
    json_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    md_path = output_dir / "input_validation_report.md"
    lines = [
        f"# Carbohydrate input validation: {report.get('molecule_id') or 'unknown'}",
        "",
        f"Status: **{report['status']}**",
        "",
        f"Errors: {len(report['errors'])}; warnings: {len(report['warnings'])}.",
    ]
    for heading, key in (("Errors", "errors"), ("Warnings", "warnings")):
        lines.extend(["", f"## {heading}", ""])
        if not report[key]:
            lines.append("None.")
        else:
            for item in report[key]:
                location = ""
                if item.get("file"):
                    location = f" ({item['file']} row {item.get('row', '?')})"
                lines.append(f"- `{item['code']}`{location}: {item['message']}")
    lines.extend([
        "",
        "## Interpretation",
        "",
        (
            "PASS means the inputs are structurally consistent enough to generate a "
            "provisional seed. It does not verify atom assignments, coupling signs, "
            "solution forms, or molecular identity."
        ),
    ])
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, md_path


def initialize_input_directory(
    destination: Path,
    molecule: str,
    forms: list[str],
    *,
    force: bool = False,
) -> list[Path]:
    """Create empty reviewed-input templates without inventing assignments."""

    if not SLUG.fullmatch(molecule):
        raise ValueError("molecule must be a safe identifier")
    if not forms or any(not SLUG.fullmatch(form) for form in forms):
        raise ValueError("At least one safe form identifier is required")
    destination = destination.resolve()
    paths = [destination / name for name in REQUIRED_FILES]
    if not force and any(path.exists() for path in paths):
        raise FileExistsError(f"Input templates already exist in {destination}")
    destination.mkdir(parents=True, exist_ok=True)
    fraction = 1.0 / len(forms)
    molecule_payload = {
        "schema_version": SCHEMA_VERSION,
        "molecule_id": molecule,
        "name": molecule.replace("_", " "),
        "conditions": {
            "solvent": "FILL_ME",
            "temperature_K": None,
            "pH": None,
            "chemical_shift_reference": "FILL_ME",
        },
        "forms": [
            {"id": form, "fraction": fraction, "fraction_status": "PROVISIONAL"}
            for form in forms
        ],
        "chemistry": {"reducing": len(forms) > 1},
        "optimization": {"shift_bound_ppm": 0.05, "coupling_bound_hz": 3.0},
    }
    paths[0].write_text(json.dumps(molecule_payload, indent=2) + "\n", encoding="utf-8")
    for path, headers in zip(paths[1:], (PROTON_HEADERS, COUPLING_HEADERS, SPECTRUM_HEADERS)):
        with path.open("w", newline="", encoding="utf-8") as handle:
            csv.writer(handle).writerow(headers)
    return paths


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repo-root",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    init_parser = subparsers.add_parser("init", help="Create the four input templates")
    init_parser.add_argument("--molecule", required=True)
    init_parser.add_argument("--forms", nargs="+", required=True)
    init_parser.add_argument("--input-dir", type=Path)
    init_parser.add_argument("--force", action="store_true")

    build_parser = subparsers.add_parser("build", help="Validate inputs and generate provisional artifacts")
    build_parser.add_argument("--input-dir", type=Path, required=True)
    build_parser.add_argument("--output-dir", type=Path)
    build_parser.add_argument("--validate-only", action="store_true")
    build_parser.add_argument("--force", action="store_true")

    args = parser.parse_args()
    repo_root = args.repo_root.resolve()
    if args.command == "init":
        input_dir = args.input_dir or repo_root / "data" / args.molecule / "input"
        paths = initialize_input_directory(
            input_dir, args.molecule, args.forms, force=args.force
        )
        print(f"Created reviewed-input templates in {input_dir}")
        for path in paths:
            print(f"  {path.name}")
        print("Fill every field, then run the build command. No matrix was generated.")
        return 0

    try:
        bundle = load_bundle(args.input_dir)
    except (FileNotFoundError, ValueError, json.JSONDecodeError) as error:
        raise SystemExit(f"Input bundle could not be loaded: {error}") from error
    molecule_id = str(bundle["molecule"].get("molecule_id") or args.input_dir.name)
    output_dir = args.output_dir or repo_root / "outputs" / molecule_id / "bootstrap"
    normalized, report = validate_bundle(bundle, repo_root)
    json_report, md_report = write_validation_report(report, output_dir)
    print(f"Input validation: {report['status']}")
    print(f"Wrote {json_report}")
    print(f"Wrote {md_report}")
    if report["status"] != "PASS":
        for item in report["errors"]:
            print(f"ERROR {item['code']}: {item['message']}")
        return 2
    if args.validate_only:
        print("Validation only; no matrices were generated.")
        return 0
    try:
        artifacts = generate_artifacts(
            normalized, report, repo_root, output_dir, force=args.force
        )
    except FileExistsError as error:
        raise SystemExit(str(error)) from error
    for label, path in artifacts.items():
        print(f"Wrote {label}: {path}")
    print("Status: PROVISIONAL_SEED — manual review and held-out validation required.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
