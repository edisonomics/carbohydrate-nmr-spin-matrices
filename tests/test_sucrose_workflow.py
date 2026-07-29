#!/usr/bin/env python3
"""Fast, dependency-light regression tests for the sucrose repository."""

from __future__ import annotations

import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
COMMON = ROOT / "src" / "common"
SUCROSE = ROOT / "src" / "sucrose" / "bayes_astro"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class SucroseConfigurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads((ROOT / "data/sucrose/sucrose_config.json").read_text())

    def test_config_has_four_official_fields(self):
        self.assertEqual({str(d["key"]) for d in self.config["datasets"]}, {"600", "800", "900", "1100"})

    def test_training_and_validation_groups_are_explicit(self):
        model = self.config["independent_model"]
        self.assertEqual(model["fields_for_joint_fit"], ["600", "900"])
        self.assertEqual(model["validation_fields"], ["800", "1100"])

    def test_publication_target_is_gissmo(self):
        self.assertEqual(self.config["publication"]["target"], "GISSMO")

    def test_dataset_metadata_and_roles_are_present(self):
        self.assertTrue(all("sample_id" in d and "tube_id" in d and "role" in d for d in self.config["datasets"]))


class XyloseConfigurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.config = json.loads((ROOT / "data/xylose/xylose_config.json").read_text())

    def test_same_5mm_tube_is_recorded_for_every_field(self):
        self.assertEqual({item["tube_id"] for item in self.config["datasets"]}, {"5mm_same_tube"})

    def test_600_900_train_and_1100_validates(self):
        model = self.config["independent_model"]
        self.assertEqual(model["fields_for_joint_fit"], ["600", "900"])
        self.assertEqual(model["validation_fields"], ["1100"])
        roles = {str(item["key"]): item["role"] for item in self.config["datasets"]}
        self.assertEqual(roles, {"600": "training", "900": "training", "1100": "validation"})


class SucroseMatrixTests(unittest.TestCase):
    def test_downloaded_matrix_is_14_by_14_and_matches_canonical(self):
        downloaded = ROOT / "data/sucrose/bmrb/bmse000119/gissmo_spin_matrix.txt"
        canonical = ROOT / "data/sucrose/matrix/sucrose_spin_matrix_GISSMO_14x14.txt"
        self.assertTrue(downloaded.is_file(), downloaded)
        self.assertTrue(canonical.is_file(), canonical)
        a = [[float(x) for x in line.split()] for line in downloaded.read_text().splitlines() if line.strip()]
        b = [[float(x) for x in line.split()] for line in canonical.read_text().splitlines() if line.strip()]
        self.assertEqual((len(a), len(a[0])), (14, 14))
        self.assertEqual((len(b), len(b[0])), (14, 14))
        for row_a, row_b in zip(a, b):
            self.assertEqual(len(row_a), len(row_b))
            for x, y in zip(row_a, row_b):
                self.assertAlmostEqual(x, y, places=8)

    def test_atom_order_has_fourteen_labels(self):
        labels = json.loads((ROOT / "data/sucrose/bmrb/bmse000119/gissmo_atom_ids.json").read_text())
        self.assertEqual(len(labels), 14)
        self.assertEqual(labels, ["24", "25", "30", "32", "34", "35", "37", "36", "33", "31", "26", "27", "28", "29"])


class ProvenanceTests(unittest.TestCase):
    def test_seed_manifest_resolves_downloaded_matrix(self):
        manifest = json.loads((ROOT / "outputs/sucrose/seed_selection.json").read_text())
        self.assertEqual(manifest["status"], "READY")
        self.assertEqual(manifest["selected_source"], "gissmo_bmrb_endpoint")
        self.assertTrue((ROOT / manifest["matrix_file"]).is_file())

    def test_bmrb_provenance_artifacts_exist(self):
        base = ROOT / "data/sucrose/bmrb/bmse000119"
        for name in ("bmse000119.str", "chemical_shifts.csv", "spectral_observations.json", "provenance.json"):
            self.assertTrue((base / name).is_file(), name)
        provenance = json.loads((base / "provenance.json").read_text())
        self.assertEqual(provenance["source"], "BMRB")
        self.assertGreater(provenance["proton_shift_count"], 0)


class PhysicsAndParserTests(unittest.TestCase):
    def test_bubb_profile_distinguishes_reducing_and_nonreducing_sugars(self):
        bubb = load_module("bubb_rules", COMMON / "bubb_rules.py")
        sucrose = bubb.profile_for("sucrose")
        xylose = bubb.profile_for("xylose")
        self.assertEqual(sucrose["expected_model"], "single_molecule")
        self.assertEqual(xylose["expected_model"], "anomer_mixture")
        self.assertTrue(xylose["guidance"]["use_2d_for_assignment"])

    def test_bubb_assessment_flags_flattened_reducing_sugar(self):
        bubb = load_module("bubb_rules_assessment", COMMON / "bubb_rules.py")
        result = bubb.assess_config(
            "xylose",
            {"name": "xylose", "matrix_file": "matrix.txt", "chemistry": {"bubb_profile": "xylose"}},
            bmrb={"proton_shift_count": 12, "gissmo_matrix_file": "matrix.txt"},
        )
        self.assertEqual(result["status"], "REVIEW")
        self.assertTrue(any("separate anomer" in warning for warning in result["warnings"]))

    def test_bubb_assessment_accepts_sucrose_gissmo_seed(self):
        bubb = load_module("bubb_rules_sucrose", COMMON / "bubb_rules.py")
        result = bubb.assess_config(
            "sucrose",
            {"name": "sucrose", "matrix_file": "matrix.txt", "chemistry": {"bubb_profile": "sucrose"}},
            bmrb={"proton_shift_count": 14, "gissmo_matrix_file": "matrix.txt"},
        )
        self.assertEqual(result["status"], "PASS")

    def test_bubb_assessment_does_not_call_component_seed_gissmo(self):
        bubb = load_module("bubb_rules_seed_status", COMMON / "bubb_rules.py")
        result = bubb.assess_config(
            "xylose",
            {"name": "xylose", "components": [{"fraction": 1}], "chemistry": {"bubb_profile": "xylose"}},
            bmrb={"proton_shift_count": 17},
        )
        self.assertFalse(result["checks"]["gissmo_matrix_available"])
        self.assertEqual(result["seed_status"], "BMRB_PROVISIONAL")

    def test_evidence_policy_does_not_promote_old_exploratory_fits(self):
        bubb = load_module("bubb_rules_policy", COMMON / "bubb_rules.py")
        policy = bubb.profile_for("sucrose")["evidence_policy"]
        self.assertIn("old exploratory fits", policy["context_only"])
        self.assertIn("verified GISSMO/BMRB matrix", policy["numeric_authority"])

    def test_empty_artifact_window_is_supported_for_anomeric_peaks(self):
        prep = load_module("prepare_carbohydrate_spectra", COMMON / "prepare_carbohydrate_spectra.py")
        self.assertEqual(prep.reason_for(5.18, {"fit_region_ppm": [3.0, 5.5], "water_region_ppm": [4.65, 4.90], "artifact_region_ppm": []}), "included")

    def test_mixture_component_linewidths_are_preserved(self):
        model = load_module("carbohydrate_model_linewidths", COMMON / "carbohydrate_model.py")
        config = json.loads((ROOT / "data/xylose/xylose_config.json").read_text())
        specs = model.component_specs(config, ROOT)
        self.assertEqual([round(float(spec["linewidth_hz"]), 4) for spec in specs], [1.8357, 1.6329])

    def test_gissmo_html_matrix_parser(self):
        query = load_module("query_bmrb_entry", COMMON / "query_bmrb_entry.py")
        page = """<h2>Spin System Matrix</h2><table><tr><th></th><th>1</th><th>2</th></tr><tr><td>1</td><td>3.5</td><td>7</td></tr><tr><td>2</td><td>0</td><td>4.2</td></tr></table>"""
        self.assertEqual(query.extract_gissmo_matrix(page), (["1", "2"], [[3.5, 7.0], [0.0, 4.2]]))

    def test_bmrb_molecule_name_search(self):
        query = load_module("query_bmrb_name_search", COMMON / "query_bmrb_entry.py")
        page = b'<a href="mol_summary/index.php?id=bmse000119&whichTab=0">Sucrose</a>'

        class Response:
            def __enter__(self):
                return self

            def __exit__(self, *args):
                return False

            def read(self):
                return page

        original = query.urlopen
        query.urlopen = lambda *args, **kwargs: Response()
        try:
            self.assertEqual(query.discover_bmrb_entries_by_molecule("sucrose"), ["bmse000119"])
        finally:
            query.urlopen = original

    def test_bmrb_xylose_experiment_inventory_records_2d_evidence(self):
        query = load_module("query_bmrb_xylose_inventory", COMMON / "query_bmrb_entry.py")
        star = ROOT / "data/xylose/bmrb/bmse000026/bmse000026.str"
        if not star.is_file():
            self.skipTest("xylose BMRB entry is not present")
        text = star.read_text(encoding="utf-8", errors="replace")
        experiments = query.extract_bmrb_experiment_inventory(text)
        names = {item["name"] for item in experiments}
        self.assertIn("2D [1H,1H]-TOCSY", names)
        self.assertIn("2D [1H,1H]-COSY", names)
        peaks = query.extract_bmrb_peak_inventory(text)
        self.assertTrue(any(item["experiment_name"] == "2D [1H,13C]-HSQC" and item["peak_count"] > 0 for item in peaks))

    def test_bmrb_mannose_numeric_hsqc_coordinates_are_preserved(self):
        query = load_module("query_bmrb_mannose_peak_data", COMMON / "query_bmrb_entry.py")
        star = ROOT / "data/mannose/bmrb/bmse000018/bmse000018.str"
        if not star.is_file():
            self.skipTest("mannose BMRB entry is not present")
        data = query.extract_bmrb_peak_data(star.read_text(encoding="utf-8", errors="replace"))
        hsqc = next(item for item in data if item["experiment_name"] == "2D [1H,13C]-HSQC")
        self.assertEqual(hsqc["dimensions"], "2")
        self.assertEqual(hsqc["peak_count"], 11)
        self.assertEqual(hsqc["peaks"][0]["coordinate_values"], [5.17, 96.701])
        assigned = hsqc["peaks"][0]["assigned_atoms"]
        self.assertIn("H19", assigned[0]["atom_ids"])
        self.assertIn("C6", assigned[1]["atom_ids"])

    def test_hsqc_refinement_reads_mannose_proton_targets(self):
        refine = load_module("refine_hsqc_constraints_test", COMMON / "refine_hsqc_constraints.py")
        observations = ROOT / "data/mannose/bmrb/bmse000018/spectral_observations.json"
        if not observations.is_file():
            self.skipTest("mannose BMRB observations are not present")
        targets = refine._read_hsqc_proton_targets(observations)
        self.assertEqual(len(targets), 11)
        self.assertIn(5.17, targets)
        self.assertIn(4.887, targets)

    def test_provisional_seed_builder_preserves_observations(self):
        builder = load_module("build_provisional_seed", COMMON / "build_provisional_seed.py")
        ids, matrix = builder.build_matrix({
            "atoms": [{"id": "H1", "shift_ppm": 5.4}, {"id": "H2", "shift_ppm": 3.8}],
            "couplings": [{"i": "H1", "j": "H2", "j_hz": 7.0}],
        })
        self.assertEqual(ids, ["H1", "H2"])
        self.assertEqual(matrix[0][1], 7.0)
        self.assertEqual(matrix[1][0], 7.0)

    def test_metadata_split_planner_uses_sample_and_tube_groups(self):
        planner = load_module("plan_multifield_split", COMMON / "plan_multifield_split.py")
        result = planner.plan(self_config())
        self.assertEqual(result["status"], "READY_FOR_VALIDATION")
        self.assertEqual(set(result["training_fields"]), {"600", "900"})
        self.assertEqual(set(result["validation_fields"]), {"800", "1100"})

    def test_importer_copies_bruker_files_and_updates_config(self):
        importer = load_module("import_bruker_dataset", COMMON / "import_bruker_dataset.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmp_root = Path(tmp)
            source = tmp_root / "source" / "12"; processed = source / "pdata" / "1"
            processed.mkdir(parents=True)
            (source / "acqus").write_text("##$SFO1= 1000\n")
            (processed / "1r").write_bytes(b"1234")
            (processed / "procs").write_text("##$SI= 1\n")
            molecule = tmp_root / "data" / "demo"
            molecule.mkdir(parents=True)
            (molecule / "demo_config.json").write_text(json.dumps({"name":"demo","datasets":[]}))
            args = type("Args", (), {"repo_root": tmp_root, "source": source, "procno":"1", "molecule":"demo", "field_mhz":1000.0, "experiment":"12", "sample_id":"A", "tube_id":"5mm", "concentration_mM":100.0, "role":"training", "replace":False})()
            result = importer.import_dataset(args)
            self.assertTrue((tmp_root / "data/demo/1000_MHz/12/pdata/1/1r").is_file())
            self.assertEqual(result["dataset"]["sample_id"], "A")

    def test_initializer_uses_seed_manifest(self):
        initializer = load_module("init_carbohydrate", COMMON / "init_carbohydrate.py")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "outputs/demo").mkdir(parents=True)
            (root / "outputs/demo/seed_selection.json").write_text(json.dumps({"matrix_file":"data/demo/matrix.txt","atom_ids":["1","2"]}))
            old_argv = sys.argv
            try:
                sys.argv = ["init_carbohydrate", "--molecule", "demo", "--repo-root", str(root)]
                self.assertEqual(initializer.main(), 0)
            finally:
                sys.argv = old_argv
            config = json.loads((root / "data/demo/demo_config.json").read_text())
            self.assertEqual(config["matrix_file"], "data/demo/matrix.txt")
            self.assertEqual(config["atom_ids"], ["1", "2"])

    def test_new_carbohydrate_does_not_inherit_sucrose_artifact_mask(self):
        initializer = load_module("init_carbohydrate_defaults", COMMON / "init_carbohydrate.py")
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            initializer.create_config(root, "new_sugar")
            config = json.loads((root / "data/new_sugar/new_sugar_config.json").read_text())
            self.assertEqual(config["processing"]["artifact_region_ppm"], [])
            self.assertEqual(config["chemistry"]["bubb_profile"], "new_sugar")

    def test_new_carbohydrate_declares_independent_j_measurement_stage(self):
        initializer = load_module("init_carbohydrate_j_policy", COMMON / "init_carbohydrate.py")
        with tempfile.TemporaryDirectory() as tmp:
            config = json.loads(initializer.create_config(Path(tmp), "new_sugar").read_text())
            policy = config["j_measurement"]
            self.assertTrue(policy["required_for_new_matrix"])
            self.assertTrue(policy["matrix_update_requires_manual_review"])
            self.assertIn("resolved_1h_multiplet", policy["allowed_methods"])

    def test_independent_simulator_selftests(self):
        try:
            simulator = load_module("sucrose_sim", SUCROSE / "sucrose_sim.py")
        except ModuleNotFoundError as error:
            if error.name == "numpy":
                self.skipTest("numpy is not installed in this Python environment")
            raise
        simulator._selftests()


if __name__ == "__main__":
    unittest.main(verbosity=2)


def self_config():
    return json.loads((ROOT / "data/sucrose/sucrose_config.json").read_text())
