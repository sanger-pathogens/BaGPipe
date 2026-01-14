#!/usr/bin/env python3
"""
validation_summaries.py - Pure Python (no pandas)
"""

import argparse
import sys
from pathlib import Path
import csv
import re
from collections import defaultdict


def load_unitig_gene_map(annotated_kmers):
    u2g = {}
    with open(annotated_kmers) as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            if not row:
                continue
            unitig = row[0].strip()
            if not unitig:
                continue
            ann = row[-1]
            genes = re.findall(r'"([^"]+)"', ann)
            if not genes:
                continue
            gene = genes[0]
            u2g[unitig] = gene
    return u2g if u2g else None


def load_rtab(rtab_path):
    samples = []
    unitigs = {}
    with open(rtab_path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if header:
            samples = [s.strip() for s in header[1:]]
        for row in reader:
            if not row or not row[0].strip():
                continue
            unitig = row[0].strip()
            counts = {}
            for i, val in enumerate(row[1:]):
                if i < len(samples):
                    try:
                        counts[samples[i]] = int(val) if int(val) > 0 else 0
                    except (ValueError, IndexError):
                        counts[samples[i]] = 0
            unitigs[unitig] = counts
    return samples, unitigs


def summarise(rtab, annotated_kmers, outdir):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    samples, unitigs = load_rtab(rtab)
    
    if not samples or not unitigs:
        print("Warning: empty Rtab matrix", file=sys.stderr)
        return

    # Per-sample
    per_sample = {}
    for sample in samples:
        count = sum(1 for unitig_counts in unitigs.values() if unitig_counts.get(sample, 0) > 0)
        per_sample[sample] = count

    with open(outdir / "validation_per_sample.tsv", "w") as f:
        f.write("sample\tn_unitigs_present\n")
        for sample in samples:
            f.write(f"{sample}\t{per_sample[sample]}\n")

    # Per-unitig
    per_unitig = {}
    for unitig, counts in unitigs.items():
        prevalence = sum(1 for c in counts.values() if c > 0)
        per_unitig[unitig] = prevalence

    with open(outdir / "validation_per_unitig.tsv", "w") as f:
        f.write("unitig\tn_samples_present\n")
        for unitig, prev in per_unitig.items():
            f.write(f"{unitig}\t{prev}\n")

    # Per-gene
    if annotated_kmers is None or not Path(annotated_kmers).exists():
        return

    u2g = load_unitig_gene_map(annotated_kmers)
    if u2g is None:
        return

    per_gene = defaultdict(set)
    for unitig, gene in u2g.items():
        if unitig in unitigs:
            for sample in samples:
                if unitigs[unitig].get(sample, 0) > 0:
                    per_gene[gene].add(sample)

    with open(outdir / "validation_per_gene.tsv", "w") as f:
        f.write("gene\tn_samples_present\n")
        for gene in sorted(per_gene.keys()):
            f.write(f"{gene}\t{len(per_gene[gene])}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rtab", required=True)
    ap.add_argument("--annotated_kmers", help="annotated_kmers.tsv from discovery")
    ap.add_argument("--outdir", default="validation")
    args = ap.parse_args()
    summarise(args.rtab, args.annotated_kmers, args.outdir)


if __name__ == "__main__":
    main()
