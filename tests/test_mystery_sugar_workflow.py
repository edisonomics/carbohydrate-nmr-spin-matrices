#!/usr/bin/env python3
"""Regression tests for identity-free multifield 1-D screening."""

from __future__ import annotations

import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "src" / "mystery_sugars" / "identify_from_1d.py"


def load_module():
    spec = importlib.util.spec_from_file_location("mystery_sugar_identity_screen", MODULE_PATH)
    if spec is None or spec.loader is None:
        raise ImportError(MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MysterySugarScreenTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        try:
            cls.screen = load_module()
        except ModuleNotFoundError as error:
            if error.name == "numpy":
                raise unittest.SkipTest("numpy is not installed in this Python environment")
            raise
        cls.np = cls.screen.np

    def synthetic_fields(self):
        ppm = self.np.linspace(3.0, 5.7, 4000)
        centers = [3.216, 3.312, 3.420, 3.512, 3.625, 3.916, 4.564, 5.182]
        y = sum(
            (1.0 - 0.04 * index) * self.np.exp(-((ppm - center) / 0.0025) ** 2)
            for index, center in enumerate(centers)
        )
        return [
            {"field_mhz": 600.0, "ppm": ppm, "intensity": y},
            {"field_mhz": 900.0, "ppm": ppm + 0.0003, "intensity": y},
            {"field_mhz": 1100.0, "ppm": ppm - 0.0002, "intensity": y},
        ]

    def test_two_anomeric_clusters_are_detected(self):
        fields = self.synthetic_fields()
        clusters = self.screen.detect_anomeric_clusters(fields[0]["ppm"], fields[0]["intensity"])
        self.assertEqual(len(clusters), 2)
        self.assertAlmostEqual(clusters[0], 4.564, places=2)
        self.assertAlmostEqual(clusters[1], 5.182, places=2)

    def test_xylose_candidate_ranks_above_other_reference_candidates(self):
        ranked = self.screen.rank_candidates(
            self.synthetic_fields(), ROOT, enable_physics=False
        )
        self.assertEqual(ranked[0]["candidate_id"], "d_xylose")
        self.assertTrue(ranked[0]["reference_available"])
        self.assertGreater(ranked[0]["mean_score"], 0.8)

    def test_screen_never_calls_one_d_identity_confirmed(self):
        ranked = self.screen.rank_candidates(
            self.synthetic_fields(), ROOT, enable_physics=False
        )
        report = self.screen.build_report(ROOT, "mystery_sugar", ranked)
        self.assertIn(report["status"], {"CANDIDATE", "REVIEW"})
        self.assertNotEqual(report["status"], "CONFIRMED")
        self.assertTrue(report["identity_confirmation_required"])

    def test_reference_free_fructose_stays_review(self):
        fields = self.synthetic_fields()
        fructose = next(item for item in self.screen.load_library() if item["id"] == "d_fructose")
        result = self.screen.score_candidate(fructose, fields, ROOT)
        report = self.screen.build_report(ROOT, "mystery_sugar", [result])
        self.assertFalse(result["reference_available"])
        self.assertEqual(report["status"], "REVIEW")

    def test_embedded_bmrb_centers_can_score_without_a_local_shift_file(self):
        candidate = {
            "id": "bmrb_example",
            "name": "Example sugar",
            "class": "chebi_carbohydrate_unreviewed",
            "bmrb_entry": "bmse999999",
            "reference_anomeric_centers_ppm": [4.564, 5.182],
            "expected_anomeric_clusters": 2,
            "forms": [],
            "topology": "unreviewed",
            "review_status": "needs_review",
        }
        result = self.screen.score_candidate(candidate, self.synthetic_fields(), ROOT)

        self.assertTrue(result["reference_available"])
        self.assertGreater(result["mean_score"], 0.8)
        self.assertEqual(result["review_status"], "needs_review")

    def test_unreviewed_bmrb_candidate_cannot_be_promoted(self):
        candidate = {
            "id": "bmrb_example",
            "name": "Unreviewed sugar",
            "class": "chebi_carbohydrate_unreviewed",
            "bmrb_entry": "bmse999999",
            "reference_anomeric_centers_ppm": [4.564, 5.182],
            "expected_anomeric_clusters": 2,
            "forms": [],
            "topology": "unreviewed",
            "review_status": "needs_review",
        }
        result = self.screen.score_candidate(candidate, self.synthetic_fields(), ROOT)
        report = self.screen.build_report(ROOT, "mystery_sugar", [result])

        self.assertGreater(result["mean_score"], 0.8)
        self.assertEqual(report["status"], "REVIEW")
        self.assertEqual(report["top_review_status"], "needs_review")

    def test_catalog_mode_adds_only_screenable_nonduplicate_candidates(self):
        catalog = {
            "schema_version": 1,
            "candidates": [
                {
                    "candidate_id": "bmrb_bmse000015",
                    "name": "Duplicate glucose",
                    "review_status": "needs_review",
                    "selected_bmrb_entry": "bmse000015",
                    "reference_anomeric_centers_ppm": [4.6, 5.2],
                },
                {
                    "candidate_id": "bmrb_bmse999998",
                    "name": "Full-spectrum-only reference",
                    "review_status": "needs_review",
                    "selected_bmrb_entry": "bmse999998",
                    "reference_anomeric_centers_ppm": [],
                    "reference_proton_shifts_ppm": [3.2, 3.6],
                },
                {
                    "candidate_id": "bmrb_bmse999999",
                    "name": "New screenable sugar",
                    "review_status": "needs_review",
                    "selected_bmrb_entry": "bmse999999",
                    "reference_anomeric_centers_ppm": [5.1],
                    "evidence": {"has_2d": True},
                },
            ],
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "catalog.json"
            path.write_text(json.dumps(catalog), encoding="utf-8")
            candidates = self.screen.load_library(
                include_bmrb_catalog=True, bmrb_catalog_path=path
            )

        ids = {candidate["id"] for candidate in candidates}
        self.assertIn("bmrb_bmse999999", ids)
        self.assertIn("bmrb_bmse999998", ids)
        self.assertNotIn("bmrb_bmse000015", ids)

    def test_full_fingerprint_distinguishes_same_anomeric_centers(self):
        base = {
            "class": "test",
            "reference_anomeric_centers_ppm": [4.564, 5.182],
            "expected_anomeric_clusters": 2,
            "forms": [],
            "topology": "test",
            "review_status": "needs_review",
        }
        matching = {
            **base,
            "id": "matching",
            "name": "Matching",
            "reference_proton_shifts_ppm": [
                3.216, 3.312, 3.420, 3.512, 3.625, 3.916, 4.564, 5.182
            ],
        }
        wrong_body = {
            **base,
            "id": "wrong",
            "name": "Wrong body",
            "reference_proton_shifts_ppm": [
                3.05, 3.12, 3.28, 3.75, 4.05, 4.20, 4.564, 5.182
            ],
        }
        good = self.screen.score_candidate(
            matching, self.synthetic_fields(), ROOT, enable_physics=False
        )
        bad = self.screen.score_candidate(
            wrong_body, self.synthetic_fields(), ROOT, enable_physics=False
        )
        self.assertGreater(good["mean_score"], bad["mean_score"] + 0.2)

    def test_exact_matrix_shape_favors_xylose_for_real_mystery_spectra(self):
        fields = self.screen.load_prepared(ROOT, "mystery_sugar")
        library = self.screen.load_library()
        xylose = next(item for item in library if item["id"] == "d_xylose")
        glucose = next(item for item in library if item["id"] == "d_glucose")
        xylose_result = self.screen.score_physics_model(xylose, fields, ROOT)
        glucose_result = self.screen.score_physics_model(glucose, fields, ROOT)

        self.assertIsNotNone(xylose_result)
        self.assertIsNotNone(glucose_result)
        self.assertGreater(
            xylose_result["mean_multiplet_shape_score"],
            glucose_result["mean_multiplet_shape_score"] + 0.5,
        )
        self.assertTrue(xylose_result["matrix_couplings"])

    def test_bubb_guidance_recovers_xylose_anomeric_j_classes(self):
        fields = self.screen.load_prepared(ROOT, "mystery_sugar")
        xylose = next(
            item for item in self.screen.load_library() if item["id"] == "d_xylose"
        )
        result = self.screen.score_bubb_guidance(xylose, fields, ROOT)

        self.assertIsNotNone(result)
        checks = {item["form"]: item for item in result["anomeric_j_checks"]}
        self.assertAlmostEqual(checks["alpha"]["observed_spacing_hz"], 3.7, delta=0.3)
        self.assertAlmostEqual(checks["beta"]["observed_spacing_hz"], 8.0, delta=0.3)
        self.assertEqual(checks["alpha"]["field_support"], 3)
        self.assertEqual(checks["beta"]["field_support"], 3)
        self.assertEqual(result["score"], 1.0)

    def test_bubb_guidance_recovers_mannose_small_anomeric_j(self):
        fields = self.screen.load_prepared(ROOT, "mannose")
        mannose = next(
            item for item in self.screen.load_library() if item["id"] == "d_mannose"
        )
        result = self.screen.score_bubb_guidance(mannose, fields, ROOT)

        self.assertIsNotNone(result)
        alpha = next(
            item for item in result["anomeric_j_checks"] if item["form"] == "alpha"
        )
        self.assertAlmostEqual(alpha["observed_spacing_hz"], 1.6, delta=0.2)
        self.assertGreaterEqual(alpha["field_support"], 2)
        self.assertEqual(alpha["expected_range_hz"], [1.1, 2.1])

    def test_expanded_known_mannose_screen_beats_shift_only_trehalose(self):
        fields = self.screen.load_prepared(ROOT, "mannose")
        expanded = self.screen.load_library(include_bmrb_catalog=True)
        library = [
            item
            for item in expanded
            if item["id"] == "d_mannose" or item["name"] == "D-Trehalose"
        ]
        ranked = self.screen.rank_candidates(fields, ROOT, library)

        self.assertEqual(ranked[0]["candidate_id"], "d_mannose")
        self.assertGreater(
            ranked[0]["mean_score"], ranked[1]["mean_score"] + 0.05
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
