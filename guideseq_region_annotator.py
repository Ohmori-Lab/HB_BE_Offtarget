#!/usr/bin/env python3
"""
GUIDE-seq off-target region annotator (hg38-ready)_Created by Nemekhbayar B. 

- Scans an input directory for "*identifiedOfftargets*.txt" (GUIDE-seq outputs)
- Normalizes coordinates and signal
- Annotates each site as exon / intron / intergenic using a GTF (e.g., GENCODE v46 for hg38)
- Writes per-sample Excel summaries and PNG plots
- Writes a combined summary across all samples

Usage:
  python guideseq_region_annotator.py \
    --input_dir /home/user/Desktop/guideseq/output/identified/Data_PAM_NNN \
    --gtf /home/user/Desktop/Guideseq_Tsai/Human_genome/gencode.v46.annotation.gtf \
    --out_dir /home/user/Desktop/guideseq/output/identified/annotated_hg38

Requires: pandas, numpy, matplotlib, (optional) xlsxwriter for Excel output
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
    s = str(s)
    if s.startswith("chr"):
        return s
    if re.match(r"^\d+$", s):
        return f"chr{s}"
    if s in {"X", "Y", "M", "MT"}:
        return f"chr{s}"
    return s


def merge_intervals(intervals: List[Tuple[int, int]]) -> List[Tuple[int, int]]:
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
    # Linear scan; OK once merged. (Could be bisect for huge sets.)
    for s, e in intervals:
        if s <= x <= e:
            return True
        if x < s:
            return False
    return False


def annotate_region(chrom: str, pos: float, exons, genes) -> str:
    if pd.isna(pos):
        return "unknown"
    if chrom in exons and overlaps(exons[chrom], int(pos)):
        return "exon"
    if chrom in genes and overlaps(genes[chrom], int(pos)):
        return "intron"
    return "intergenic"


def pick_signal(df: pd.DataFrame) -> pd.Series:
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
    if "Position" in df.columns and df["Position"].notna().any():
        return pd.to_numeric(df["Position"], errors="coerce")
    minpos = pd.to_numeric(df.get("Min.Position", np.nan), errors="coerce")
    maxpos = pd.to_numeric(df.get("Max.Position", np.nan), errors="coerce")
    return (minpos + maxpos) / 2.0


def extract_coords(df: pd.DataFrame) -> pd.DataFrame:
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

    # By chromosome
    counts_per_chr = sites.groupby("chrom").size().sort_values(ascending=False)
    plt.figure()
    counts_per_chr.plot(kind="bar")
    plt.title(f"{sample_name}: Off-target counts per chromosome")
    plt.ylabel("Count")
    plt.xlabel("Chromosome")
    plt.tight_layout()
    plt.savefig(sample_outdir / "plot_by_chromosome.png", dpi=160)


def write_excel(sample_outdir: Path, sites: pd.DataFrame):
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

    plt.figure()
    by_mismatch = combined.groupby("mismatches").size().rename("count").reset_index().sort_values("mismatches")
    plt.bar(by_mismatch["mismatches"].astype(str), by_mismatch["count"])
    plt.title("All samples: Off-target counts by mismatch")
    plt.ylabel("Count")
    plt.xlabel("Mismatches")
    plt.tight_layout()
    plt.savefig(comb_dir / "plot_by_mismatch.png", dpi=160)

    plt.figure()
    combined.groupby("chrom").size().sort_values(ascending=False).plot(kind="bar")
    plt.title("All samples: Off-target counts per chromosome")
    plt.ylabel("Count")
    plt.xlabel("Chromosome")
    plt.tight_layout()
    plt.savefig(comb_dir / "plot_by_chromosome.png", dpi=160)

    print(f"[DONE] Wrote outputs to: {out_dir}")


if __name__ == "__main__":
    main()

