import argparse
import gzip
import json
from pathlib import Path

import pandas as pd

IMPACT_RANK = {
    "HIGH": 3,
    "MODERATE": 2,
    "LOW": 1,
    "MODIFIER": 0,
    "UNSPECIFIED": -1,
}


def parse_info_field(info_field: str) -> dict[str, str]:
    """Extract key-value entries from a VCF INFO field."""
    parsed: dict[str, str] = {}
    for entry in info_field.split(";"):
        if "=" not in entry:
            continue
        key, value = entry.split("=", 1)
        parsed[key] = value
    return parsed


def open_variant_stream(input_path: Path):
    """Open plain-text or gzipped VCF input."""
    if input_path.suffix == ".gz":
        return gzip.open(input_path, "rt", encoding="utf-8")
    return input_path.open("r", encoding="utf-8")


def classify_variant(ref: str, alt: str) -> str:
    """Classify a variant based on reference and alternate allele lengths."""
    # This is a lightweight structural classification rather than a full normalization pass.
    if len(ref) == len(alt) == 1:
        return "SNP"
    if len(alt) > len(ref):
        return "insertion"
    if len(ref) > len(alt):
        return "deletion"
    return "complex"


def normalize_filter_status(filter_value: str) -> str:
    """Convert VCF FILTER placeholders into stable analyst-facing labels."""
    normalized = filter_value.strip()
    if not normalized or normalized == ".":
        return "UNSPECIFIED"
    return normalized


def _parse_position(raw_position: str, line: str) -> int:
    """Parse and validate a 1-based VCF position."""
    try:
        position = int(raw_position)
    except ValueError as exc:
        raise ValueError(f"VCF record has a non-integer POS value: {line}") from exc

    if position <= 0:
        raise ValueError(f"VCF record has a non-positive POS value: {line}")
    return position


def _validate_alleles(ref: str, alt: str, line: str) -> list[str]:
    """Validate reference and alternate alleles and expand ALT values."""
    if not ref or ref == ".":
        raise ValueError(f"VCF record is missing a reference allele: {line}")

    alternate_alleles = [allele.strip() for allele in alt.split(",")]
    if not alternate_alleles or any(not allele or allele == "." for allele in alternate_alleles):
        raise ValueError(f"VCF record contains a missing alternate allele: {line}")
    return alternate_alleles


def parse_vcf(input_path: Path) -> pd.DataFrame:
    """Read a VCF file and return parsed variant records."""
    records = []
    saw_header = False

    with open_variant_stream(input_path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                if line.startswith("#CHROM"):
                    saw_header = True
                continue

            fields = line.split()
            if len(fields) < 8:
                raise ValueError(f"VCF record has fewer than 8 columns: {line}")

            chrom, pos, _, ref, alt, qual, filter_value, info = fields[:8]
            info_map = parse_info_field(info)
            gene = info_map.get("GENE", "UNKNOWN")
            impact = info_map.get("IMPACT", "UNSPECIFIED")
            transcript = info_map.get("TRANSCRIPT", "NA")
            position = _parse_position(pos, line)
            alternate_alleles = _validate_alleles(ref, alt, line)

            for alt_allele in alternate_alleles:
                # Multi-allelic rows are expanded so downstream summaries operate on
                # one alternate allele per record instead of a mixed ALT field.
                records.append(
                    {
                        "chromosome": chrom,
                        "position": position,
                        "reference_allele": ref,
                        "alternate_allele": alt_allele,
                        "gene": gene,
                        "impact": impact,
                        "transcript": transcript,
                        "quality": qual,
                        "filter": normalize_filter_status(filter_value),
                        "variant_type": classify_variant(ref, alt_allele),
                    }
                )

    if not saw_header:
        raise ValueError("Input does not contain a #CHROM header line.")

    return pd.DataFrame(records)


def most_severe_impact(series: pd.Series) -> str:
    """Return the highest-severity impact label present in a series."""
    if series.empty:
        return "UNSPECIFIED"
    return max(series.astype(str), key=lambda value: IMPACT_RANK.get(value, -1))


def summarize_by_gene(variants: pd.DataFrame) -> pd.DataFrame:
    """Count variants per gene and return a summary table."""
    if variants.empty:
        return pd.DataFrame(columns=["gene", "variant_count", "top_impact"])

    # Gene-level summary keeps the count and the most severe observed impact so
    # the output remains compact while still prioritizing potentially important calls.
    summary = (
        variants.groupby("gene", as_index=False)
        .agg(
            variant_count=("gene", "size"),
            top_impact=("impact", most_severe_impact),
        )
        .sort_values(["variant_count", "gene"], ascending=[False, True])
    )
    return summary


def summarize_variant_types(variants: pd.DataFrame) -> pd.DataFrame:
    """Count variants by structural class."""
    if variants.empty:
        return pd.DataFrame(columns=["variant_type", "count"])

    return (
        variants.groupby("variant_type", as_index=False)
        .size()
        .rename(columns={"size": "count"})
        .sort_values(["count", "variant_type"], ascending=[False, True])
    )


def summarize_impacts(variants: pd.DataFrame) -> pd.DataFrame:
    """Count variants by impact label using severity-aware ordering."""
    if variants.empty:
        return pd.DataFrame(columns=["impact", "count"])

    summary = (
        variants.groupby("impact", as_index=False)
        .size()
        .rename(columns={"size": "count"})
    )
    summary["impact_rank"] = summary["impact"].map(lambda value: IMPACT_RANK.get(str(value), -1))
    summary = summary.sort_values(["count", "impact_rank", "impact"], ascending=[False, False, True])
    return summary.drop(columns=["impact_rank"]).reset_index(drop=True)


def summarize_filter_statuses(variants: pd.DataFrame) -> pd.DataFrame:
    """Count variants by normalized FILTER status."""
    if variants.empty:
        return pd.DataFrame(columns=["filter", "count"])

    summary = (
        variants.groupby("filter", as_index=False)
        .size()
        .rename(columns={"size": "count"})
    )
    summary["pass_rank"] = summary["filter"].map(lambda value: 1 if str(value) == "PASS" else 0)
    summary = summary.sort_values(["count", "pass_rank", "filter"], ascending=[False, False, True])
    return summary.drop(columns=["pass_rank"]).reset_index(drop=True)


def build_run_report(variants: pd.DataFrame) -> dict:
    """Generate a compact machine-readable run summary."""
    return {
        "variant_count": int(len(variants.index)),
        "gene_count": int(variants["gene"].nunique()) if not variants.empty else 0,
        "variant_type_counts": summarize_variant_types(variants).to_dict(orient="records"),
        "impact_counts": summarize_impacts(variants).to_dict(orient="records"),
        "filter_counts": summarize_filter_statuses(variants).to_dict(orient="records"),
        "top_genes": summarize_by_gene(variants).head(5).to_dict(orient="records"),
    }


def write_dataframe_output(dataframe: pd.DataFrame, output_path: Path) -> None:
    """Write a tabular output using suffix-driven format selection."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = output_path.suffix.lower()
    if suffix == ".csv":
        dataframe.to_csv(output_path, index=False)
        return
    if suffix == ".json":
        output_path.write_text(dataframe.to_json(orient="records", indent=2), encoding="utf-8")
        return
    raise ValueError(f"Unsupported output format for {output_path}. Use .csv or .json.")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Parse a VCF file and summarize genomic variants.")
    parser.add_argument("--input", required=True, help="Path to the input VCF file.")
    parser.add_argument(
        "--summary-json",
        default=None,
        help="Optional path to write a machine-readable summary JSON file.",
    )
    parser.add_argument("--variants-out", default=None, help="Optional path to write normalized variants (.csv or .json).")
    parser.add_argument("--gene-summary-out", default=None, help="Optional path to write gene summary (.csv or .json).")
    parser.add_argument("--type-summary-out", default=None, help="Optional path to write variant-type summary (.csv or .json).")
    parser.add_argument("--impact-summary-out", default=None, help="Optional path to write impact summary (.csv or .json).")
    parser.add_argument("--filter-summary-out", default=None, help="Optional path to write FILTER summary (.csv or .json).")
    return parser


def main() -> None:
    """Parse the requested VCF, print summaries, and optionally write export artifacts."""
    args = build_parser().parse_args()
    input_path = Path(args.input)

    variants = parse_vcf(input_path)
    gene_summary = summarize_by_gene(variants)
    type_summary = summarize_variant_types(variants)
    impact_summary = summarize_impacts(variants)
    filter_summary = summarize_filter_statuses(variants)

    if variants.empty:
        print("No variant records found.")
        return

    print("Parsed variants:")
    print(variants.to_string(index=False))
    print()
    print("Summary by gene:")
    print(gene_summary.to_string(index=False))
    print()
    print("Summary by variant type:")
    print(type_summary.to_string(index=False))
    print()
    print("Summary by impact:")
    print(impact_summary.to_string(index=False))
    print()
    print("Summary by FILTER status:")
    print(filter_summary.to_string(index=False))

    if args.summary_json:
        summary_path = Path(args.summary_json)
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        summary_path.write_text(json.dumps(build_run_report(variants), indent=2), encoding="utf-8")
        print()
        print(f"Summary JSON written to: {summary_path.resolve()}")

    if args.variants_out:
        variants_path = Path(args.variants_out)
        write_dataframe_output(variants, variants_path)
        print(f"Variants table written to: {variants_path.resolve()}")

    if args.gene_summary_out:
        gene_summary_path = Path(args.gene_summary_out)
        write_dataframe_output(gene_summary, gene_summary_path)
        print(f"Gene summary written to: {gene_summary_path.resolve()}")

    if args.type_summary_out:
        type_summary_path = Path(args.type_summary_out)
        write_dataframe_output(type_summary, type_summary_path)
        print(f"Variant-type summary written to: {type_summary_path.resolve()}")

    if args.impact_summary_out:
        impact_summary_path = Path(args.impact_summary_out)
        write_dataframe_output(impact_summary, impact_summary_path)
        print(f"Impact summary written to: {impact_summary_path.resolve()}")

    if args.filter_summary_out:
        filter_summary_path = Path(args.filter_summary_out)
        write_dataframe_output(filter_summary, filter_summary_path)
        print(f"FILTER summary written to: {filter_summary_path.resolve()}")


if __name__ == "__main__":
    main()
