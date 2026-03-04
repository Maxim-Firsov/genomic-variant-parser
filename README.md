# Genomic Variant Parser

VCF parsing and summarization project for small-variant review workflows. The repository converts raw variant rows into gene-level and variant-type summaries, and it can emit stable tabular outputs for downstream inspection.

## Problem Statement

Variant call files are compact but not analyst-friendly. Before deeper interpretation, it is useful to normalize records into a tabular representation and produce quick summaries such as gene hit counts, impact levels, and mutation-type counts.

## Technical Approach

- Reads plain-text or gzipped VCF input.
- Requires a valid `#CHROM` header and rejects malformed short records.
- Parses `INFO` tags into structured fields such as `GENE`, `IMPACT`, and `TRANSCRIPT`.
- Expands multi-allelic rows so each alternate allele becomes its own record.
- Classifies each variant as `SNP`, `insertion`, `deletion`, or `complex`.
- Exports normalized variants and summaries as `.csv` or `.json`.

## Architecture

```text
VCF / VCF.GZ input
      |
      v
Header validation + record parsing
      |
      v
INFO normalization
      |
      v
Pandas DataFrame
      |
      +--> normalized variant table
      |
      +--> gene-level summary
      |
      +--> variant-type summary
      |
      +--> machine-readable run report
```

## Repository Layout

- `src/parser.py`: parsing logic, summarization, validation, and CLI
- `data/example.vcf`: sample input with gene and impact annotations
- `tests/test_parser.py`: tests for parsing, gzip support, output writing, and header validation
- `.github/workflows/ci.yml`: automated test workflow

## Example Usage

```powershell
python src\parser.py --input data\example.vcf --summary-json outputs\summary.json --variants-out outputs\variants.csv --gene-summary-out outputs\gene_summary.csv --type-summary-out outputs\type_summary.json
```

Observed CLI summary:

```text
Summary by gene:
 gene  variant_count top_impact
BRCA1              2   MODERATE
 EGFR              1   MODERATE
 KRAS              1       HIGH
 TP53              1       HIGH

Summary by variant type:
variant_type  count
         SNP      3
   insertion      2
```

## Design Decisions

- The parser fails fast on invalid structure because silent record skipping undermines trust in genomic pipelines.
- Outputs are suffix-driven (`.csv` or `.json`) so generated artifacts are predictable and automation-friendly.
- Severity summary uses the most severe observed impact rather than the mode, which is a better default for triage workflows.

## Limitations

- This is a targeted parser for simple small-variant summaries, not a full VCF normalization or annotation engine.
- INFO extraction is key-based and currently focuses on a small set of tags.
- The sample dataset is intentionally compact for portability and test speed.

## Running The Project

```powershell
python -m pip install -r requirements.txt
python src\parser.py --input data\example.vcf --summary-json outputs\summary.json
python -m unittest discover -s tests
```
