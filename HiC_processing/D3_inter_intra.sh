#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/A2_hicExplorer
OUTDIR=$MAINDIR/HiC/D3_inter_intra
cd $OUTDIR

# Remove sex chrom
# Need normalized and mapq30 filtered ginteraction file

cat > "inter_intra_stats.py" <<'EOF'
#!/usr/bin/env python3
import argparse
import re
import sys
from collections import defaultdict


def parse_sex_chroms(path, species):
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if parts[0] == species:
                scaffolds = {p for p in parts[1:] if p != "NONE"}
                print(f"  Sex chromosomes excluded for {species}: {scaffolds}", file=sys.stderr)
                return scaffolds
    return set()


def parse_chrom_sizes(path, exclude):
    sizes = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            chrom, size = parts[0], int(parts[1])
            if chrom in exclude:
                continue
            sizes[chrom] = size
    return sizes


def parse_ginteractions_binned(path, chrom_sizes, binsize):
    intra     = defaultdict(int)
    inter     = defaultdict(int)
    bin_intra = defaultdict(int)
    bin_inter = defaultdict(int)

    known_chroms = set(chrom_sizes.keys())
    skipped = 0

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                skipped += 1
                continue

            c1, s1, c2, s2 = parts[0], int(parts[1]), parts[3], int(parts[4])

            if c1 not in known_chroms or c2 not in known_chroms:
                skipped += 1
                continue

            bin1 = (s1 // binsize) * binsize
            bin2 = (s2 // binsize) * binsize

            if c1 == c2:
                intra[c1] += 1
                bin_intra[(c1, bin1)] += 1
            else:
                inter[c1] += 1
                inter[c2] += 1
                bin_inter[(c1, bin1)] += 1
                bin_inter[(c2, bin2)] += 1

    if skipped:
        print(f"  [info] {skipped} lines skipped", file=sys.stderr)

    return intra, inter, bin_intra, bin_inter


def write_stats(path, all_chroms, chrom_sizes, intra, inter):
    with open(path, "w") as out:
        out.write("chrom\tsize\tintra\tinter\tintra_per_Mb\tinter_per_Mb\tratio_inter_intra\n")
        for chrom in all_chroms:
            size = chrom_sizes[chrom]
            n_intra = intra.get(chrom, 0)
            n_inter = inter.get(chrom, 0)
            mb = size / 1e6

            intra_per_mb = n_intra / mb if mb > 0 else 0
            inter_per_mb = n_inter / mb if mb > 0 else 0
            ratio = f"{n_inter / n_intra:.4f}" if n_intra > 0 else "NA"

            out.write(f"{chrom}\t{size}\t{n_intra}\t{n_inter}\t"
                      f"{intra_per_mb:.4f}\t{inter_per_mb:.4f}\t{ratio}\n")
    print(f"Stats written to {path}", file=sys.stderr)


def write_binned_stats(path, chrom_sizes, bin_intra, bin_inter, binsize):
    all_chroms = sorted(chrom_sizes.keys(),
                         key=lambda x: [int(c) if c.isdigit() else c for c in re.split(r'(\d+)', x)])
    with open(path, "w") as out:
        out.write("chrom\tstart\tend\tintra\tinter\tratio_inter_intra\n")
        for chrom in all_chroms:
            chrom_len = chrom_sizes[chrom]
            bin_start = 0
            while bin_start < chrom_len:
                bin_end = min(bin_start + binsize, chrom_len)
                n_intra = bin_intra.get((chrom, bin_start), 0)
                n_inter = bin_inter.get((chrom, bin_start), 0)
                ratio = f"{n_inter / n_intra:.4f}" if n_intra > 0 else "NA"
                out.write(f"{chrom}\t{bin_start}\t{bin_end}\t{n_intra}\t{n_inter}\t{ratio}\n")
                bin_start += binsize
    print(f"Binned stats written to {path}", file=sys.stderr)


parser = argparse.ArgumentParser(description="Per-chromosome inter/intra contact stats (sex chroms excluded)")
parser.add_argument("--gint", required=True)
parser.add_argument("--chrom", required=True)
parser.add_argument("--sexchrom", required=True)
parser.add_argument("--species", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--out-bins", required=True)
parser.add_argument("--binsize", type=int, default=100_000)
args = parser.parse_args()

sex_chroms = parse_sex_chroms(args.sexchrom, args.species)
chrom_sizes = parse_chrom_sizes(args.chrom, exclude=sex_chroms)
intra, inter, bin_intra, bin_inter = parse_ginteractions_binned(args.gint, chrom_sizes, args.binsize)
all_chroms = sorted(chrom_sizes.keys())

write_stats(args.out, all_chroms, chrom_sizes, intra, inter)
write_binned_stats(args.out_bins, chrom_sizes, bin_intra, bin_inter, args.binsize)
EOF

interintra() {
    SAMPLE=$1
    mkdir -p $OUTDIR/${SAMPLE}
    SIZES=$MAINDIR/HiC/A1_map_HiC/${SAMPLE}/ref.sizes
    GINT=$INDIR/${SAMPLE}/${SAMPLE}_mapq30_norm_corr_100kb.ginteractions.tsv

    python3 inter_intra_stats.py \
        --gint $GINT \
        --chrom $SIZES \
        --sexchrom $OUTDIR/sexchromlist.txt \
        --species $SAMPLE \
        --out $OUTDIR/${SAMPLE}/${SAMPLE}_inter_intra.tsv \
        --out-bins $OUTDIR/${SAMPLE}/${SAMPLE}_inter_intra_perbin.tsv
}

### Parallelize over individuals
export INDIR OUTDIR
export -f interintra
parallel --colsep '\t' 'interintra {1}' :::: $MAINDIR/HiC/A2_hicExplorer/species_list.txt
