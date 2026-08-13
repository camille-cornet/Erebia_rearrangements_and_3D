#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/A2_hicExplorer
OUTDIR=$MAINDIR/HiC/D3_inter_intra_Z
cd $OUTDIR

# Inter/intra calculations per chromosome but include the Z this time

cat > "inter_intra_stats.py" <<'EOF'
#!/usr/bin/env python3
import argparse
import sys
from collections import defaultdict

def parse_chrom_sizes(path):
    sizes = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            sizes[parts[0]] = int(parts[1])
    return sizes

def parse_ginteractions(path, chrom_sizes):
    intra = defaultdict(int)
    inter = defaultdict(int)
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

            c1, c2 = parts[0], parts[3]
            if c1 not in known_chroms or c2 not in known_chroms:
                skipped += 1
                continue

            if c1 == c2:
                intra[c1] += 1
            else:
                inter[c1] += 1
                inter[c2] += 1

    if skipped:
        print(f"  [info] {skipped} lines skipped", file=sys.stderr)

    return intra, inter

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

parser = argparse.ArgumentParser(description="Per-chromosome inter/intra contact stats (all chroms)")
parser.add_argument("--gint", required=True)
parser.add_argument("--chrom", required=True)
parser.add_argument("--out", required=True)
args = parser.parse_args()

chrom_sizes = parse_chrom_sizes(args.chrom)
intra, inter = parse_ginteractions(args.gint, chrom_sizes)
all_chroms = sorted(chrom_sizes.keys())

write_stats(args.out, all_chroms, chrom_sizes, intra, inter)
EOF

interintra() {
    SAMPLE=$1
    mkdir -p $OUTDIR/${SAMPLE}
    SIZES=$MAINDIR/HiC/A1_map_HiC/${SAMPLE}/ref.sizes
    GINT=$INDIR/${SAMPLE}/${SAMPLE}_mapq30_norm_corr_100kb.ginteractions.tsv

    python3 inter_intra_stats.py --gint $GINT --chrom $SIZES --out $OUTDIR/${SAMPLE}/${SAMPLE}_inter_intra_allchroms.tsv
}

### Parallelize over species
export INDIR OUTDIR
export -f interintra
parallel --colsep '\t' 'interintra {1}' :::: $MAINDIR/HiC/A2_hicExplorer/species_list.txt

### Combine per-chrom all-chrom stats across species
header=1
while read SAMPLE; do
    f="$OUTDIR/${SAMPLE}/${SAMPLE}_inter_intra_allchroms.tsv"
    if [ "$header" -eq 1 ]; then
        echo -e "species\t$(head -1 "$f")"
        header=0
    fi
    tail -n +2 "$f" | awk -v sp="$SAMPLE" '{print sp"\t"$0}'
done < $MAINDIR/HiC/A2_hicExplorer/species_list.txt \
> $OUTDIR/inter_intra_allchroms_allspecies.tsv
