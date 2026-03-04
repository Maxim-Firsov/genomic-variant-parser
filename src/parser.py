import argparse
from pathlib import Path

import pandas as pd


def parse_info_field(info_field: str) -> str:
    """Extract the gene name from a VCF INFO field."""
    for entry in info_field.split(";"):
        if entry.startswith("GENE="):
            return entry.split("=", 1)[1] or "UNKNOWN"
    return "UNKNOWN"


def classify_variant(ref: str, alt: str) -> str:
    """Classify a variant based on reference and alternate allele lengths."""
    if len(ref) == len(alt) == 1:
        return "SNP"
    if len(alt) > len(ref):
        return "insertion"
    if len(ref) > len(alt):
        return "deletion"
    return "complex"


def parse_vcf(input_path: Path) -> pd.DataFrame:
    """Read a VCF file and return parsed variant records."""
    records = []

    with input_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            fields = line.split()
            if len(fields) < 8:
                continue

            chrom, pos, _, ref, alt, _, _, info = fields[:8]
            gene = parse_info_field(info)

            for alt_allele in alt.split(","):
                records.append(
                    {
                        "chromosome": chrom,
                        "position": int(pos),
                        "reference_allele": ref,
                        "alternate_allele": alt_allele,
                        "gene": gene,
                        "variant_type": classify_variant(ref, alt_allele),
                    }
                )

    return pd.DataFrame(records)


def summarize_variants(variants: pd.DataFrame) -> pd.DataFrame:
    """Count variants per gene and return a summary table."""
    if variants.empty:
        return pd.DataFrame(columns=["gene", "variant_count"])

    summary = (
        variants.groupby("gene", as_index=False)
        .size()
        .rename(columns={"size": "variant_count"})
        .sort_values(["variant_count", "gene"], ascending=[False, True])
    )
    return summary


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Parse a VCF file and summarize variants per gene.")
    parser.add_argument("--input", required=True, help="Path to the input VCF file.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    input_path = Path(args.input)

    variants = parse_vcf(input_path)
    summary = summarize_variants(variants)

    if variants.empty:
        print("No variant records found.")
        return

    print("Parsed variants:")
    print(variants.to_string(index=False))
    print()
    print("Summary by gene:")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
