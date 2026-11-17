#!/usr/bin/env python3
"""
CFTR variant analysis from ClinVar.

Downloads the ClinVar variant_summary file (or uses a local copy), filters it
to CFTR variants on the GRCh38 assembly, and produces summary tables and figures.

Usage:
    python scripts/analyze_cftr.py                  # download + analyse
    python scripts/analyze_cftr.py --input FILE.gz  # use a local copy
"""
import argparse
import re
import sys
import urllib.request
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "variant_summary.txt.gz"
PROCESSED = ROOT / "data" / "cftr_variants.csv"
RESULTS = ROOT / "results"

GENE_RE = re.compile(r"(?:^|;)CFTR(?:;|$)")
COLORS = {
    "Pathogenic/Likely pathogenic": "#d62728",
    "Benign/Likely benign": "#2ca02c",
    "Uncertain significance (VUS)": "#ff7f0e",
    "Conflicting": "#9467bd",
    "Other": "#7f7f7f",
}


def download(dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading ClinVar variant_summary to {dest} (several hundred MB)...")

    def hook(blocks, block_size, total):
        if total > 0:
            pct = min(100, blocks * block_size * 100 / total)
            print(f"\r  {pct:5.1f}%", end="", flush=True)

    urllib.request.urlretrieve(URL, dest, hook)
    print("\nDone.")


def simplify(sig) -> str:
    """Collapse ClinVar's many significance labels into five groups."""
    s = str(sig).lower()
    if "conflicting" in s:
        return "Conflicting"
    if "uncertain" in s:
        return "Uncertain significance (VUS)"
    if "pathogenic" in s and "benign" not in s:
        return "Pathogenic/Likely pathogenic"
    if "benign" in s and "pathogenic" not in s:
        return "Benign/Likely benign"
    return "Other"


def load_cftr(path: Path) -> pd.DataFrame:
    header = pd.read_csv(path, sep="\t", nrows=0).columns
    sig_col = next(
        (c for c in ("ClinicalSignificance", "GermlineClassification") if c in header),
        None,
    )
    if sig_col is None:
        sys.exit(f"Could not find a clinical significance column. Columns: {list(header)}")
    wanted = ["Type", "Name", "GeneSymbol", sig_col, "ReviewStatus",
              "Assembly", "Chromosome", "Start", "Stop", "PhenotypeList"]
    usecols = [c for c in wanted if c in header]

    chunks = []
    for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=200_000, low_memory=False):
        mask = chunk["GeneSymbol"].fillna("").str.contains(GENE_RE) & (chunk["Assembly"] == "GRCh38")
        chunks.append(chunk[mask])
    df = pd.concat(chunks, ignore_index=True).rename(columns={sig_col: "ClinicalSignificance"})
    df["Category"] = df["ClinicalSignificance"].map(simplify)
    return df


def make_figures(df: pd.DataFrame) -> None:
    RESULTS.mkdir(exist_ok=True)

    # Figure 1: clinical significance
    counts = df["Category"].value_counts()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.barh(counts.index[::-1], counts.values[::-1],
            color=[COLORS.get(c, "#7f7f7f") for c in counts.index[::-1]])
    for i, v in enumerate(counts.values[::-1]):
        ax.text(v, i, f" {v}", va="center")
    ax.set_xlabel("Number of ClinVar records")
    ax.set_title("CFTR variants in ClinVar by clinical significance (GRCh38)")
    fig.tight_layout()
    fig.savefig(RESULTS / "fig1_clinical_significance.png", dpi=200)
    plt.close(fig)

    # Figure 2: variant types
    types = df["Type"].value_counts()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.barh(types.index[::-1], types.values[::-1], color="#1f77b4")
    ax.set_xlabel("Number of ClinVar records")
    ax.set_title("CFTR variants by molecular type")
    fig.tight_layout()
    fig.savefig(RESULTS / "fig2_variant_types.png", dpi=200)
    plt.close(fig)

    # Figure 3: position along the gene
    pos = df.dropna(subset=["Start"]).copy()
    pos = pos[pos["Start"] > 0]
    fig, ax = plt.subplots(figsize=(10, 4.5))
    for cat in ("Pathogenic/Likely pathogenic", "Benign/Likely benign"):
        sub = pos[pos["Category"] == cat]
        ax.hist(sub["Start"] / 1e6, bins=80, alpha=0.6, label=f"{cat} (n={len(sub)})",
                color=COLORS[cat])
    ax.set_xlabel("Position on chromosome 7 (Mb, GRCh38)")
    ax.set_ylabel("Number of variants")
    ax.set_title("Distribution of CFTR variants along the gene")
    ax.legend()
    fig.tight_layout()
    fig.savefig(RESULTS / "fig3_position_along_gene.png", dpi=200)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", type=Path, help="local variant_summary(.txt/.txt.gz); skips download")
    args = ap.parse_args()

    path = args.input
    if path is None:
        path = RAW
        if not path.exists():
            download(path)

    df = load_cftr(path)
    if df.empty:
        sys.exit("No CFTR GRCh38 records found. Check the input file.")

    PROCESSED.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(PROCESSED, index=False)
    RESULTS.mkdir(exist_ok=True)
    summary = df["Category"].value_counts().rename_axis("Category").reset_index(name="Count")
    summary.to_csv(RESULTS / "summary_by_category.csv", index=False)
    make_figures(df)

    print(f"\nCFTR GRCh38 records: {len(df)}")
    print(summary.to_string(index=False))
    print(f"\nSaved table to {PROCESSED.relative_to(ROOT)} and figures to results/")


if __name__ == "__main__":
    main()
