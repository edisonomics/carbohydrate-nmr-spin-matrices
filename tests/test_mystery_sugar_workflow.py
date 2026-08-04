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
        y = (
            self.np.exp(-((ppm - 5.182) / 0.0025) ** 2)
            + 0.9 * self.np.exp(-((ppm - 4.564) / 0.0025) ** 2)
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
        ranked = self.screen.rank_candidates(self.synthetic_fields(), ROOT)
        self.assertEqual(ranked[0]["candidate_id"], "d_xylose")
        self.assertTrue(ranked[0]["reference_available"])
        self.assertGreater(ranked[0]["mean_score"], 0.8)

    def test_screen_never_calls_one_d_identity_confirmed(self):
        ranked = self.screen.rank_candidates(self.synthetic_fields(), ROOT)
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
                    "name": "No anomeric reference",
                    "review_status": "needs_review",
                    "selected_bmrb_entry": "bmse999998",
                    "reference_anomeric_centers_ppm": [],
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
        self.assertNotIn("bmrb_bmse000015", ids)
        self.assertNotIn("bmrb_bmse999998", ids)


if __name__ == "__main__":
    unittest.main(verbosity=2)
