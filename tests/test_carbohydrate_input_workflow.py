#!/usr/bin/env python3
"""Regression tests for the reviewed carbohydrate input workflow."""

from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "common" / "bootstrap_carbohydrate.py"


def load_module():
    spec = importlib.util.spec_from_file_location("bootstrap_carbohydrate", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_csv(path: Path, headers: tuple[str, ...], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(rows)


class CarbohydrateInputWorkflowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.workflow = load_module()

    def make_bundle(self, root: Path) -> Path:
        input_dir = root / "inputs"
        spectrum_dir = root / "spectra"
        input_dir.mkdir(parents=True)
        spectrum_dir.mkdir(parents=True)
        for key in ("f500", "f600", "f700"):
            experiment = spectrum_dir / key
            processed = experiment / "pdata" / "1"
            processed.mkdir(parents=True)
            field = key.removeprefix("f")
            (experiment / "acqus").write_text(
                f"##$NUC1= <1H>\n##$SFO1= {field}.1\n##$PULPROG= <zg30>\n",
                encoding="utf-8",
            )
            (processed / "procs").write_text("##$SI= 1\n", encoding="utf-8")
            (processed / "1r").write_bytes(b"\x00\x00\x00\x00")

        molecule = {
            "schema_version": 1,
            "molecule_id": "test_sugar",
            "name": "Test sugar",
            "conditions": {
                "solvent": "D2O",
                "temperature_K": 298.15,
                "pH": 7.0,
                "chemical_shift_reference": "DSS, 0.00 ppm",
            },
            "forms": [
                {
                    "id": "alpha",
                    "fraction": 0.4,
                    "fraction_status": "MEASURED",
                },
                {
                    "id": "beta",
                    "fraction": 0.6,
                    "fraction_status": "MEASURED",
                },
            ],
            "chemistry": {"reducing": True},
        }
        (input_dir / "molecule.json").write_text(
            json.dumps(molecule, indent=2) + "\n", encoding="utf-8"
        )

        protons = [
            {
                "component": "alpha",
                "spin_id": "alpha_H1",
                "atom_label": "H1",
                "attached_atom": "C1",
                "bmrb_label": "H1",
                "shift_ppm": 5.20,
                "assignment_status": "ASSIGNED",
                "source": "lab HSQC peak A1",
                "uncertainty_ppm": 0.01,
                "fit": "true",
            },
            {
                "component": "alpha",
                "spin_id": "alpha_H2",
                "atom_label": "H2",
                "attached_atom": "C2",
                "bmrb_label": "H2",
                "shift_ppm": 3.50,
                "assignment_status": "ASSIGNED",
                "source": "lab HSQC peak A2",
                "uncertainty_ppm": 0.01,
                "fit": "true",
            },
            {
                "component": "beta",
                "spin_id": "beta_H1",
                "atom_label": "H1",
                "attached_atom": "C1",
                "bmrb_label": "H1",
                "shift_ppm": 4.60,
                "assignment_status": "DEPOSITED",
                "source": "BMRB test entry",
                "uncertainty_ppm": 0.02,
                "fit": "true",
            },
            {
                "component": "beta",
                "spin_id": "beta_H2",
                "atom_label": "H2",
                "attached_atom": "C2",
                "bmrb_label": "H2",
                "shift_ppm": 3.20,
                "assignment_status": "DEPOSITED",
                "source": "BMRB test entry",
                "uncertainty_ppm": 0.02,
                "fit": "true",
            },
        ]
        write_csv(
            input_dir / "proton_assignments.csv",
            self.workflow.PROTON_HEADERS,
            protons,
        )

        couplings = [
            {
                "component": "alpha",
                "spin_i": "alpha_H1",
                "spin_j": "alpha_H2",
                "J_hz": 3.5,
                "evidence_status": "ASSIGNED",
                "source": "resolved 600/700 MHz multiplet A1",
                "uncertainty_hz": 0.2,
                "fit": "true",
            },
            {
                "component": "beta",
                "spin_i": "beta_H1",
                "spin_j": "beta_H2",
                "J_hz": 8.0,
                "evidence_status": "ASSIGNED",
                "source": "resolved 600/700 MHz multiplet B1",
                "uncertainty_hz": 0.2,
                "fit": "true",
            },
        ]
        write_csv(
            input_dir / "couplings.csv",
            self.workflow.COUPLING_HEADERS,
            couplings,
        )

        spectra = [
            {
                "key": "f500",
                "field_mhz": 500,
                "path": str(spectrum_dir / "f500"),
                "procno": "1",
                "role": "training",
                "sample_id": "sample_1",
                "tube_id": "tube_1",
                "nucleus": "1H",
                "pulse_program": "zg30",
            },
            {
                "key": "f600",
                "field_mhz": 600,
                "path": str(spectrum_dir / "f600"),
                "procno": "1",
                "role": "training",
                "sample_id": "sample_1",
                "tube_id": "tube_1",
                "nucleus": "1H",
                "pulse_program": "zg30",
            },
            {
                "key": "f700",
                "field_mhz": 700,
                "path": str(spectrum_dir / "f700"),
                "procno": "1",
                "role": "validation",
                "sample_id": "sample_1",
                "tube_id": "tube_1",
                "nucleus": "1H",
                "pulse_program": "zg30",
            },
        ]
        write_csv(input_dir / "spectra.csv", self.workflow.SPECTRUM_HEADERS, spectra)
        return input_dir

    def validated_bundle(self, root: Path):
        input_dir = self.make_bundle(root)
        loaded = self.workflow.load_bundle(input_dir)
        normalized, report = self.workflow.validate_bundle(loaded, root)
        return input_dir, normalized, report

    def test_valid_bundle_generates_symmetric_component_matrices_and_config(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_dir, bundle, report = self.validated_bundle(root)
            self.assertEqual(report["status"], "PASS")
            self.assertEqual(report["summary"]["training_fields_mhz"], [500.0, 600.0])
            self.assertEqual(report["summary"]["validation_fields_mhz"], [700.0])

            output_dir = root / "outputs" / "test_sugar" / "bootstrap"
            self.workflow.write_validation_report(report, output_dir)
            artifacts = self.workflow.generate_artifacts(
                bundle, report, root, output_dir
            )

            alpha = [
                [float(value) for value in line.split()]
                for line in artifacts["matrix_1"].read_text(encoding="utf-8").splitlines()
            ]
            self.assertEqual(alpha, [[5.2, 3.5], [3.5, 3.5]])
            config = json.loads(artifacts["config"].read_text(encoding="utf-8"))
            self.assertEqual(config["model_type"], "mixture")
            self.assertEqual(config["components"][0]["blocks"], [[1, 2]])
            self.assertEqual(
                config["independent_model"]["fields_for_joint_fit"],
                ["f500", "f600"],
            )
            self.assertEqual(config["independent_model"]["validation_fields"], ["f700"])
            self.assertFalse(config["publication"]["bootstrap_artifacts_are_publishable"])
            manifest = json.loads(artifacts["manifest"].read_text(encoding="utf-8"))
            self.assertEqual(set(manifest["input_hashes_sha256"]), set(self.workflow.REQUIRED_FILES))
            self.assertIn("must not replace", manifest["publication_block"])
            self.assertEqual(Path(bundle["input_dir"]), input_dir.resolve())

    def test_existing_bootstrap_is_not_overwritten_without_force(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, bundle, report = self.validated_bundle(root)
            output_dir = root / "outputs" / "test_sugar" / "bootstrap"
            self.workflow.generate_artifacts(bundle, report, root, output_dir)
            with self.assertRaises(FileExistsError):
                self.workflow.generate_artifacts(bundle, report, root, output_dir)

    def test_duplicate_spin_and_unknown_coupling_fail_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_dir = self.make_bundle(root)
            bundle = self.workflow.load_bundle(input_dir)
            bundle["protons"][1]["spin_id"] = "alpha_H1"
            bundle["couplings"][0]["spin_j"] = "not_a_spin"
            _, report = self.workflow.validate_bundle(bundle, root)
            codes = {item["code"] for item in report["errors"]}
            self.assertEqual(report["status"], "FAIL")
            self.assertIn("protons.spin_id.duplicate", codes)
            self.assertIn("couplings.spin", codes)

    def test_sample_mismatch_and_missing_held_out_field_fail_validation(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_dir = self.make_bundle(root)
            bundle = self.workflow.load_bundle(input_dir)
            bundle["spectra"][2]["role"] = "training"
            bundle["spectra"][2]["sample_id"] = "different_sample"
            _, report = self.workflow.validate_bundle(bundle, root)
            codes = {item["code"] for item in report["errors"]}
            self.assertEqual(report["status"], "FAIL")
            self.assertIn("spectra.validation", codes)
            self.assertIn("spectra.sample_mismatch", codes)

    def test_manifest_field_must_match_bruker_metadata(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_dir = self.make_bundle(root)
            bundle = self.workflow.load_bundle(input_dir)
            bundle["spectra"][0]["field_mhz"] = "900"
            _, report = self.workflow.validate_bundle(bundle, root)
            self.assertIn(
                "spectra.metadata.field_mismatch",
                {item["code"] for item in report["errors"]},
            )

    def test_form_fractions_must_sum_to_one(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_dir = self.make_bundle(root)
            bundle = self.workflow.load_bundle(input_dir)
            bundle["molecule"]["forms"][0]["fraction"] = 0.2
            _, report = self.workflow.validate_bundle(bundle, root)
            self.assertIn(
                "forms.fraction.sum",
                {item["code"] for item in report["errors"]},
            )

    def test_analog_prior_passes_with_explicit_review_warning(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_dir = self.make_bundle(root)
            bundle = self.workflow.load_bundle(input_dir)
            bundle["protons"][0]["assignment_status"] = "ANALOG_PRIOR"
            _, report = self.workflow.validate_bundle(bundle, root)
            self.assertEqual(report["status"], "PASS")
            self.assertIn(
                "protons.review",
                {item["code"] for item in report["warnings"]},
            )

    def test_initializer_creates_four_templates_and_refuses_silent_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "input"
            paths = self.workflow.initialize_input_directory(
                destination, "new_sugar", ["alpha", "beta"]
            )
            self.assertEqual({path.name for path in paths}, set(self.workflow.REQUIRED_FILES))
            molecule = json.loads((destination / "molecule.json").read_text(encoding="utf-8"))
            self.assertAlmostEqual(sum(form["fraction"] for form in molecule["forms"]), 1.0)
            with self.assertRaises(FileExistsError):
                self.workflow.initialize_input_directory(
                    destination, "new_sugar", ["alpha", "beta"]
                )


if __name__ == "__main__":
    unittest.main(verbosity=2)
