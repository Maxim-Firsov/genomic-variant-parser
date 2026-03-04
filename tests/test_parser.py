from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import gzip

from src.parser import (
    build_run_report,
    parse_vcf,
    summarize_by_gene,
    summarize_variant_types,
    write_dataframe_output,
)


class ParserTests(unittest.TestCase):
    def setUp(self) -> None:
        self.example_path = Path(__file__).resolve().parents[1] / "data" / "example.vcf"

    def test_parse_vcf_extracts_expected_columns(self) -> None:
        variants = parse_vcf(self.example_path)
        self.assertEqual(len(variants.index), 5)
        self.assertIn("impact", variants.columns)
        self.assertIn("transcript", variants.columns)

    def test_gene_summary_uses_highest_severity_impact(self) -> None:
        variants = parse_vcf(self.example_path)
        summary = summarize_by_gene(variants)
        brca1_row = summary.loc[summary["gene"] == "BRCA1"].iloc[0]
        self.assertEqual(brca1_row["top_impact"], "MODERATE")

    def test_run_report_is_json_ready(self) -> None:
        variants = parse_vcf(self.example_path)
        report = build_run_report(variants)
        self.assertEqual(report["variant_count"], 5)
        self.assertEqual(report["gene_count"], 4)
        self.assertEqual(report["variant_type_counts"], summarize_variant_types(variants).to_dict(orient="records"))

    def test_parse_vcf_supports_gzipped_input(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            gz_path = Path(temp_dir) / "example.vcf.gz"
            gz_path.write_bytes(b"")
            with gzip.open(gz_path, "wt", encoding="utf-8") as handle:
                handle.write(self.example_path.read_text(encoding="utf-8"))
            variants = parse_vcf(gz_path)
            self.assertEqual(len(variants.index), 5)

    def test_parse_vcf_requires_header(self) -> None:
        with tempfile.TemporaryDirectory() as temp_dir:
            invalid_path = Path(temp_dir) / "invalid.vcf"
            invalid_path.write_text("1 10 . A T . PASS GENE=TP53\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                parse_vcf(invalid_path)

    def test_write_dataframe_output_supports_csv_and_json(self) -> None:
        variants = parse_vcf(self.example_path)
        with tempfile.TemporaryDirectory() as temp_dir:
            csv_path = Path(temp_dir) / "variants.csv"
            json_path = Path(temp_dir) / "variants.json"
            write_dataframe_output(variants, csv_path)
            write_dataframe_output(variants.head(1), json_path)
            self.assertTrue(csv_path.exists())
            self.assertTrue(json_path.exists())


if __name__ == "__main__":
    unittest.main()
