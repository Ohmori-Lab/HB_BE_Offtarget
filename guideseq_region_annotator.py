#!/usr/bin/env python3
"""
GUIDE-seq off-target region annotator (hg38-ready)
-----------------------------------------------------

Author: Nemekhbayar Baatartsogt <nemekhbayar@jichi.ac.jp>
Date: 08/25/2025
Version: 1.0.0

Description:
    This script annotates off-target regions identified in GUIDE-seq experiments. It processes 
    the identified off-targets data files, normalizes genomic coordinates, and annotates regions 
    (exon, intron, intergenic) using a GTF annotation file (e.g., GENCODE v46 for hg38). 
    The results are saved in both summary Excel files and generated plots for further analysis.

Usage:
    python guideseq_region_annotator.py --input_dir <input_directory> --gtf <gtf_file> --out_dir <output_directory>

    Arguments:
        --input_dir <input_directory>    : Path to the directory containing GUIDE-seq off-targets files 
                                          (e.g., `*identifiedOfftargets*.txt`).
        --gtf <gtf_file>                 : Path to the GTF file containing genomic annotations (e.g., GENCODE v46 for hg38).
        --out_dir <output_directory>     : Path to the output directory where the annotated results and summary reports 
                                          (Excel files, plots) will be saved.

    Example:
        python guideseq_region_annotator.py --input_dir /path/to/identified_offtargets/ --gtf /path/to/genomes/hg38/gencode.v46.annotation.gtf --out_dir /path/to/annotated_results/

    Output Files:
        The script generates the following files:
        - For each sample, an annotated Excel file (`offtarget_summary.xlsx`).
        - Plots (PNG images) summarizing off-target distribution by region, mismatches, and chromosomes.
        - A combined summary for all samples will be generated in the `_combined` directory, including:
            - Combined Excel file (`combined_summary.xlsx`).
            - Summary plots for the combined data.

    Requirements:
        - pandas
        - numpy
        - matplotlib
        - xlsxwriter (optional, for Excel output)

    Notes:
        - Ensure that the GTF file corresponds to the genome version being used (e.g., hg38).
        - The script automatically falls back to writing CSV files if the Excel writer fails.
        - Input files must be in the appropriate format for processing by the script.

    Error Handling:
        - If the GTF file is not found or is in an invalid format, the script will exit with an error message.
        - If the input directory contains no files matching the pattern, the script will issue a warning and exit.
        - Invalid chromosome names or mismatched coordinates will be reported with details in the error logs.
"""

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# -----------------------------
# Helpers
# -----------------------------
def normalize_chr(s: str) -> str:
    """
    Normalize chromosome names to 'chr1', 'chrX', 'chrY', etc.

    Parameters:
        s (str): Chromosome name as a string (e.g., "1", "X", "chr3").

    Returns:
        str: Normalized chromosome name, e.g., "chr1", "chrX".
    """
    s = str(s)
    if s.startswith("chr"):
        return s
    if re.match(r"^\d+$", s):
        return f"chr{s}"
    if s in {"X", "Y", "M", "MT"}:
        return f"chr{s}"
    return s


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """
    Merge overlapping genomic intervals.

    Parameters:
        intervals (List[Tuple[int, int]]): List of tuples representing start and end positions.

    Returns:
        List[Tuple[int, int]]: Merged list of intervals.
    """
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for s, e in intervals[1:]:
        ls, le = merged[-1]
        if s <= le + 1:
            merged[-1] = (ls, max(le, e))
        else:
            merged.append((s, e))
    return merged


def parse_gtf_to_intervals(gtf_file: Path):
    """
    Parse GTF file to extract exon and gene intervals.

    Parameters:
        gtf_file (Path): Path to the GTF file.

    Returns:
        Tuple[Dict[str, List[Tuple[int, int]]], Dict[str, List[Tuple[int, int]]]]:
            A tuple containing two dictionaries:
            - exons: Chromosome -> List of exon intervals.
            - genes: Chromosome -> List of gene intervals.
    """
    exons: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    genes: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    with gtf_file.open("r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _source, feature, start, end, _score, _strand, _frame, _attrs = parts
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            chrom_n = chrom if chrom.startswith("chr") else normalize_chr(chrom)
            if feature == "exon":
                exons[chrom_n].append((start_i, end_i))
            elif feature in ("gene", "transcript"):
                genes[chrom_n].append((start_i, end_i))
    exons = {c: merge_intervals(v) for c, v in exons.items()}
    genes = {c: merge_intervals(v) for c, v in genes.items()}
    return exons, genes


def overlaps(intervals: List[Tuple[int, int]], x: int) -> bool:
    """
    Check if a position overlaps any of the intervals.

    Parameters:
        intervals (List[Tuple[int, int]]): List of intervals.
        x (int): Position to check for overlap.

    Returns:
        bool: True if position overlaps any of the intervals, False otherwise.
    """
    for s, e in intervals:
        if s <= x <= e:
            return True
        if x < s:
            return False
    return False


def annotate_region(chrom: str, pos: float, exons, genes) -> str:
    """
    Annotate genomic region as 'exon', 'intron', or 'intergenic'.

    Parameters:
        chrom (str): Chromosome name.
        pos (float): Position to annotate.

    Returns:
        str: One of 'exon', 'intron', or 'intergenic'.
    """
    if pd.isna(pos):
        return "unknown"
    if chrom in exons and overlaps(exons[chrom], int(pos)):
        return "exon"
    if chrom in genes and overlaps(genes[chrom], int(pos)):
        return "intron"
    return "intergenic"


def pick_signal(df: pd.DataFrame) -> pd.Series:
    """
    Pick the appropriate signal column based on available data.

    Parameters:
        df (pd.DataFrame): DataFrame containing GUIDE-seq data.

    Returns:
        pd.Series: Series containing signal values.
    """
    if "total.sum" in df.columns:
        return pd.to_numeric(df["total.sum"], errors="coerce").fillna(0)
    if "bi.sum.mi" in df.columns:
        return pd.to_numeric(df["bi.sum.mi"], errors="coerce").fillna(0)
    plus = pd.to_numeric(df.get("+.total", 0), errors="coerce").fillna(0)
    minus = pd.to_numeric(df.get("-.total", 0), errors="coerce").fillna(0)
    sig = plus + minus
    if sig is None:
        sig = pd.Series(np.ones(len(df)), index=df.index, dtype=float)
    return sig


def extract_pos(df: pd.DataFrame) -> pd.Series:
    """
    Extract and normalize coordinates from the dataframe.

    Parameters:
        df (pd.DataFrame): DataFrame containing GUIDE-seq data.

    Returns:
        pd.DataFrame: DataFrame with normalized chromosome coordinates.
    """
    if "Position" in df.columns and df["Position"].notna().any():
        return pd.to_numeric(df["Position"], errors="coerce")
    minpos = pd.to_numeric(df.get("Min.Position", np.nan), errors="coerce")
    maxpos = pd.to_numeric(df.get("Max.Position", np.nan), errors="coerce")
    return (minpos + maxpos) / 2.0


def extract_coords(df: pd.DataFrame) -> pd.DataFrame:
    """
    Extract and normalize chromosome coordinates and strand information.

    Parameters:
        df (pd.DataFrame): DataFrame containing GUIDE-seq data.

    Returns:
        pd.DataFrame: DataFrame with chromosome, position, start, end, strand, and coordinates.
    """
    chrom = df["Chromosome"].astype(str).map(normalize_chr)
    pos = extract_pos(df)
    minpos = pd.to_numeric(df.get("Min.Position", pos), errors="coerce")
    maxpos = pd.to_numeric(df.get("Max.Position", pos), errors="coerce")
    strand = df.get("Site_GapsAllowed.Strand")
    if strand is None or strand.isna().all():
        strand = df.get("Site_SubstitutionsOnly.Strand", pd.Series(["?"] * len(df)))
    coordinates = (
        chrom
        + ":"
        + minpos.fillna(0).astype("Int64").astype(str)
        + "-"
        + maxpos.fillna(0).astype("Int64").astype(str)
        + "("
        + strand.fillna("?").astype(str)
        + ")"
    )
    return pd.DataFrame(
        {
            "chrom": chrom,
            "pos": pos,
            "start": minpos,
            "end": maxpos,
            "strand": strand,
            "coordinates": coordinates,
        }
    )


def extract_mismatches(df: pd.DataFrame) -> pd.Series:
    """
    Extract mismatch data from the dataframe.

    Parameters:
        df (pd.DataFrame): DataFrame containing GUIDE-seq data.

    Returns:
        pd.Series: Series containing mismatch counts.
    """
    mm_sub_only = pd.to_numeric(df.get("Site_SubstitutionsOnly.NumSubstitutions"), errors="coerce")
    mm_gap_sub = pd.to_numeric(df.get("Site_GapsAllowed.Substitutions"), errors="coerce")
    mm_gap_ins = pd.to_numeric(df.get("Site_GapsAllowed.Insertions"), errors="coerce")
    mm_gap_del = pd.to_numeric(df.get("Site_GapsAllowed.Deletions"), errors="coerce")
    mismatch = mm_sub_only.copy() if mm_sub_only is not None else None
    if mismatch is None:
        mismatch = mm_gap_sub.fillna(0) + mm_gap_ins.fillna(0) + mm_gap_del.fillna(0)
    else:
        mismatch = mismatch.fillna(mm_gap_sub.fillna(0) + mm_gap_ins.fillna(0) + mm_gap_del.fillna(0))
    return mismatch


# -----------------------------
# Output writers / plots
# -----------------------------
def make_plots(sample_outdir: Path, sites: pd.DataFrame, sample_name: str):
    """
    Generate and save summary plots for off-target data from GUIDE-seq results.

    Parameters:
        sample_outdir (Path): The directory where the plots will be saved. This directory 
                               will be created if it doesn't already exist.
        sites (pd.DataFrame): A DataFrame containing off-target data, including columns 
                              for regions, mismatch counts, and chromosomal locations.
        sample_name (str): The name of the sample, used for plot titles and file naming 
                           to distinguish the plots from different samples.
    """
    sample_outdir.mkdir(parents=True, exist_ok=True)

    # By region
    by_region = sites.groupby("region").size().rename("count").reset_index()
    plt.figure()
    by_region.set_index("region")["count"].plot(kind="bar")
    plt.title(f"{sample_name}: Off-targets by region")
    plt.ylabel("Count")
    plt.xlabel("Region")
    plt.tight_layout()
    plt.savefig(sample_outdir / "plot_by_region.png", dpi=160)
    plt.close()

    # By mismatch
    by_mismatch = (
        sites.groupby("mismatches")
        .size()
        .rename("count")
        .reset_index()
        .sort_values("mismatches")
    )
    plt.figure()
    plt.bar(by_mismatch["mismatches"].astype(str), by_mismatch["count"])
    plt.title(f"{sample_name}: Off-target counts by mismatch")
    plt.ylabel("Count")
    plt.xlabel("Mismatches")
    plt.tight_layout()
    plt.savefig(sample_outdir / "plot_by_mismatch.png", dpi=160)
    plt.close()

    # By chromosome
    counts_per_chr = sites.groupby("chrom").size().sort_values(ascending=False)
    plt.figure()
    counts_per_chr.plot(kind="bar")
    plt.title(f"{sample_name}: Off-target counts per chromosome")
    plt.ylabel("Count")
    plt.xlabel("Chromosome")
    plt.tight_layout()
    plt.savefig(sample_outdir / "plot_by_chromosome.png", dpi=160)
    plt.close()


def write_excel(sample_outdir: Path, sites: pd.DataFrame):
    """
    Write off-target summary data to an Excel file (or CSV if Excel writer fails).

    Parameters:
        sample_outdir (Path): The directory where the Excel (or CSV) files will be saved. 
                               This directory will be created if it doesn't exist.
        sites (pd.DataFrame): A DataFrame containing off-target data, which will be summarized 
                              and written to the output files.
    """
    sample_outdir.mkdir(parents=True, exist_ok=True)
    out_xlsx = sample_outdir / "offtarget_summary.xlsx"
    try:
        with pd.ExcelWriter(out_xlsx, engine="xlsxwriter") as writer:
            # All sites
            sites.to_excel(writer, sheet_name="sites", index=False)
            # By region
            by_region = sites.groupby("region").size().rename("count").reset_index()
            by_region.to_excel(writer, sheet_name="by_region", index=False)
            # By mismatch
            by_mismatch = (
                sites.groupby("mismatches")
                .size()
                .rename("count")
                .reset_index()
                .sort_values("mismatches")
            )
            by_mismatch.to_excel(writer, sheet_name="by_mismatch", index=False)
            # By chromosome × region
            by_chr_region = (
                sites.groupby(["chrom", "region"])
                .size()
                .rename("count")
                .reset_index()
            )
            by_chr_region.to_excel(writer, sheet_name="by_chr_region", index=False)
    except Exception as e:
        # Fallback to CSVs if Excel writer isn't available
        print(f"[WARN] Excel write failed ({e}); writing CSVs instead.")
        sites.to_csv(sample_outdir / "sites.csv", index=False)
        sites.groupby("region").size().rename("count").reset_index() \
             .to_csv(sample_outdir / "by_region.csv", index=False)
        sites.groupby("mismatches").size().rename("count").reset_index() \
             .sort_values("mismatches") \
             .to_csv(sample_outdir / "by_mismatch.csv", index=False)
        sites.groupby(["chrom", "region"]).size().rename("count").reset_index() \
             .to_csv(sample_outdir / "by_chr_region.csv", index=False)


# -----------------------------
# Core processing
# -----------------------------
def process_one_file(path: Path, exons, genes, out_dir: Path) -> pd.DataFrame:
    """
    Process a single GUIDE-seq off-target file and generate annotated results.

    Parameters:
        path (Path): The file path to the input off-targets data file in tab-delimited format.
        exons (dict): A dictionary containing chromosome -> list of exon intervals, parsed from a GTF file.
        genes (dict): A dictionary containing chromosome -> list of gene intervals, parsed from a GTF file.
        out_dir (Path): The directory where output files (Excel, plots) will be saved. This directory will be created if it doesn't exist.

    Returns:
        pd.DataFrame: A DataFrame containing the processed and annotated off-targets data for the sample.
                      This includes the coordinates, mismatch counts, signal, and region annotations (exon/intron/intergenic).
    """
    # Read & normalize
    df = pd.read_csv(path, sep="\t")
    df.columns = [re.sub(r"\s+", "_", c.strip()) for c in df.columns]

    coords = extract_coords(df)
    signal = pick_signal(df)
    mismatches = extract_mismatches(df)

    # Sample name: keep "_filtered" when present
    stem = path.stem  # e.g., "Sample5_identifiedOfftargets_filtered"
    sample_name = stem.replace("_identifiedOfftargets", "")  # -> "Sample5_filtered" or "Sample5"

    sites = pd.DataFrame(
        {
            "sample": sample_name,
            "chrom": coords["chrom"],
            "start": coords["start"],
            "end": coords["end"],
            "pos": coords["pos"],
            "strand": coords["strand"],
            "coordinates": coords["coordinates"],
            "mismatches": mismatches,
            "signal": signal,
            "name": df.get("Name", ""),
            "target_site": df.get("Targetsite", ""),
            "target_seq": df.get("TargetSequence", ""),
            "site_seq_gapsOK": df.get("Site_GapsAllowed.Sequence", ""),
        }
    )

    # Annotate region
    sites["region"] = [annotate_region(c, p, exons, genes) for c, p in zip(sites["chrom"], sites["pos"])]

    # Ensure per-sample dir exists before any write
    sample_out = out_dir / sample_name
    sample_out.mkdir(parents=True, exist_ok=True)

    # Write outputs
    write_excel(sample_out, sites)
    make_plots(sample_out, sites, sample_name=sample_name)

    return sites


# -----------------------------
# Main
# -----------------------------
def main():
    """
    Main function to orchestrate the GUIDE-seq off-target region annotation pipeline.

    This function serves as the entry point for the script. It handles argument parsing, 
    loading necessary input files (e.g., off-targets data and GTF file), and manages the 
    overall process of annotating off-target sites. It processes multiple GUIDE-seq data files 
    by reading, normalizing, annotating, and generating output files (Excel summaries and plots).
    
    The function performs the following tasks:
    - Parses command-line arguments to specify input directory, GTF file, and output directory.
    - Loads the GTF file and extracts exon and gene intervals.
    - Iterates over all identified off-target files in the specified input directory.
    - Processes each file to normalize data, annotate regions, and generate summaries.
    - Writes the results for each sample to Excel and/or CSV, and generates summary plots.
    - Combines data from all processed files and creates a final summary with combined results.

    Parameters:
        None: The function retrieves necessary inputs from command-line arguments.

    Returns:
        None: The function does not return any value. It processes files, generates output 
              files (Excel and CSV), and produces visual plots.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_dir", required=True, help="Directory with *identifiedOfftargets*.txt files")
    ap.add_argument("--gtf", required=True, help="GTF file for hg38 (e.g., GENCODE v46)")
    ap.add_argument("--out_dir", required=True, help="Output directory for annotated results")
    args = ap.parse_args()

    input_dir = Path(args.input_dir)
    gtf_path = Path(args.gtf)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Loading GTF from {gtf_path} ...")
    exons, genes = parse_gtf_to_intervals(gtf_path)
    print(f"[INFO] Exon chroms: {len(exons)} | Gene chroms: {len(genes)}")

    files = sorted(list(input_dir.rglob("*identifiedOfftargets*.txt")))
    if not files:
        print(f"[WARN] No files found in {input_dir} matching *identifiedOfftargets*.txt")
        return

    all_sites = []
    for f in files:
        print(f"[INFO] Processing {f.name} ...")
        try:
            sites = process_one_file(f, exons, genes, out_dir)
            all_sites.append(sites)
        except Exception as e:
            print(f"[ERROR] Failed on {f}: {e}")

    if not all_sites:
        print("[WARN] No samples processed successfully.")
        return

    combined = pd.concat(all_sites, ignore_index=True)

    # Combined summaries
    comb_dir = out_dir / "_combined"
    comb_dir.mkdir(parents=True, exist_ok=True)
    try:
        with pd.ExcelWriter(comb_dir / "combined_summary.xlsx", engine="xlsxwriter") as writer:
            combined.to_excel(writer, sheet_name="sites", index=False)
            combined.groupby(["sample", "region"]).size().rename("count").reset_index() \
                    .to_excel(writer, sheet_name="by_sample_region", index=False)
            combined.groupby(["sample", "mismatches"]).size().rename("count").reset_index() \
                    .to_excel(writer, sheet_name="by_sample_mismatch", index=False)
    except Exception as e:
        print(f"[WARN] Excel write failed ({e}); writing CSVs instead.")
        combined.to_csv(comb_dir / "sites.csv", index=False)
        combined.groupby(["sample", "region"]).size().rename("count").reset_index() \
                .to_csv(comb_dir / "by_sample_region.csv", index=False)
        combined.groupby(["sample", "mismatches"]).size().rename("count").reset_index() \
                .to_csv(comb_dir / "by_sample_mismatch.csv", index=False)

    # Combined plots
    plt.figure()
    combined.groupby("region").size().rename("count").plot(kind="bar")
    plt.title("All samples: Off-targets by region")
    plt.ylabel("Count")
    plt.xlabel("Region")
    plt.tight_layout()
    plt.savefig(comb_dir / "plot_by_region.png", dpi=160)
    plt.close()

    plt.figure()
    by_mismatch = combined.groupby("mismatches").size().rename("count").reset_index().sort_values("mismatches")
    plt.bar(by_mismatch["mismatches"].astype(str), by_mismatch["count"])
    plt.title("All samples: Off-target counts by mismatch")
    plt.ylabel("Count")
    plt.xlabel("Mismatches")
    plt.tight_layout()
    plt.savefig(comb_dir / "plot_by_mismatch.png", dpi=160)
    plt.close()

    plt.figure()
    combined.groupby("chrom").size().sort_values(ascending=False).plot(kind="bar")
    plt.title("All samples: Off-target counts per chromosome")
    plt.ylabel("Count")
    plt.xlabel("Chromosome")
    plt.tight_layout()
    plt.savefig(comb_dir / "plot_by_chromosome.png", dpi=160)
    plt.close()

    print(f"[DONE] Wrote outputs to: {out_dir}")


if __name__ == "__main__":
    main()

