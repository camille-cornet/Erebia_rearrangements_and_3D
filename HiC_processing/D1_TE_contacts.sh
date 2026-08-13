#!/bin/bash

INDIR=$MAINDIR/HiC/A2_hicExplorer
OUTDIR=$MAINDIR/HiC/D1_TE_contacts
TEDIR=$MAINDIR/Annotations/C8_repeat_beds
SEXCHROM=$OUTDIR/sexchromlist.txt
cd $OUTDIR

### Count number of interchromosomal contacts between repeats of the same family
# With and without rDNA associated repeats

source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
conda activate cooler_env

mkdir $OUTDIR/repeats

TEcontacts() {
    SP=$1
    sed 's/\bscaffold_01\b/scaffold_1/g' $TEDIR/${SP}/${SP}_full.bed > $OUTDIR/repeats/${SP}_full.bed
    sed -i 's/\bscaffold_02\b/scaffold_2/g' $OUTDIR/repeats/${SP}_full.bed
    sed -i 's/\bscaffold_03\b/scaffold_3/g' $OUTDIR/repeats/${SP}_full.bed
    sed -i 's/\bscaffold_04\b/scaffold_4/g' $OUTDIR/repeats/${SP}_full.bed
    sed -i 's/\bscaffold_05\b/scaffold_5/g' $OUTDIR/repeats/${SP}_full.bed
    sed -i 's/\bscaffold_06\b/scaffold_6/g' $OUTDIR/repeats/${SP}_full.bed
    sed -i 's/\bscaffold_07\b/scaffold_7/g' $OUTDIR/repeats/${SP}_full.bed
    sed -i 's/\bscaffold_08\b/scaffold_8/g' $OUTDIR/repeats/${SP}_full.bed
    sed -i 's/\bscaffold_09\b/scaffold_9/g' $OUTDIR/repeats/${SP}_full.bed

    mkdir $OUTDIR/${SP}
    cd $OUTDIR/${SP}

    # Get sex chromosomes for this species (skip NONE entries)
    SEX_CHROMS=$(grep "^${SP}" $SEXCHROM | tr '\t' '\n' | tail -n +2 | grep -v "^NONE$")

    # Build awk filter expression for sex chromosomes
    # Matches column 1 in bed files, columns 1 and 4 in HiC pixels
    if [ -n "$SEX_CHROMS" ]; then
        AWK_BED=$(echo "$SEX_CHROMS" | awk '{printf "($1 != \"%s\") && ", $1}' | sed 's/ && $//')
        AWK_HIC=$(echo "$SEX_CHROMS" | awk '{printf "($1 != \"%s\" && $4 != \"%s\") && ", $1, $1}' | sed 's/ && $//')
    else
        AWK_BED="1"
        AWK_HIC="1"
    fi

    HICNORM=$INDIR/${SP}/${SP}_norm_corr_10kb.cool

    # Convert HiC to useable format and remove sex chromosomes
    cooler dump --join $HICNORM | awk "$AWK_HIC" > $OUTDIR/${SP}/${SP}_hic_pixels_norm.tsv

    # Remove sex chromosomes from TE bed
    awk "$AWK_BED" $OUTDIR/repeats/${SP}_full.bed > $OUTDIR/repeats/${SP}_full_autosomes.bed

    # Removing rDNA-associated TEs
    cat $OUTDIR/repeats/${SP}_full_autosomes.bed | grep rRNA | sort -k1,1 -k2,2n > $OUTDIR/repeats/${SP}_rRNA_sorted.bed
    bedtools slop -i $OUTDIR/repeats/${SP}_rRNA_sorted.bed -g $MAINDIR/assemblies/chromosomes_only/${SP}.fasta.fai -b 5000 \
        > $OUTDIR/repeats/${SP}_rRNA_10kb_windows.bed
    bedtools subtract -a $OUTDIR/repeats/${SP}_full_autosomes.bed -b $OUTDIR/repeats/${SP}_rRNA_10kb_windows.bed \
        > $OUTDIR/repeats/${SP}_no_rRNA_clusters.bed

    # Run with rDNA
    cat > "run_${SP}_norm.py" <<EOF
from collections import defaultdict

BIN = 10_000

te_bins = defaultdict(lambda: defaultdict(int))

with open("$OUTDIR/repeats/${SP}_full_autosomes.bed") as f:
    for line in f:
        chrom, start, end, te_type = line.strip().split()
        start = int(start)
        end = int(end)
        bin_start = (start // BIN) * BIN
        if start >= bin_start and end <= bin_start + BIN:
            te_bins[(chrom, bin_start)][te_type] += 1

type_sum   = defaultdict(float)
type_count = defaultdict(int)

with open("$OUTDIR/${SP}/${SP}_hic_pixels_norm.tsv") as f:
    for line in f:
        fields = line.strip().split()
        chrA = fields[0]
        binA = int(fields[1])
        chrB = fields[3]
        binB = int(fields[4])
        contacts = float(fields[6])

        if chrA == chrB:
            continue

        keyA = (chrA, binA)
        keyB = (chrB, binB)

        if keyA not in te_bins or keyB not in te_bins:
            continue

        binsA = te_bins[keyA]
        binsB = te_bins[keyB]

        for te_type in binsA.keys() & binsB.keys():
            n_pairs = binsA[te_type] * binsB[te_type]
            type_sum[te_type]   += contacts * n_pairs
            type_count[te_type] += n_pairs

with open(f"$OUTDIR/${SP}/${SP}_TEtype_interchrom_avg_norm.tsv", "w") as out:
    for te_type in sorted(type_sum):
        avg = type_sum[te_type] / type_count[te_type]
        out.write(f"{te_type}\t{avg:.6f}\n")

EOF

    python run_${SP}_norm.py

    # Run without rRNA clusters
    cat > "run_${SP}_no_rRNA_clusters_norm.py" <<EOF
from collections import defaultdict

BIN = 10_000

te_bins = defaultdict(lambda: defaultdict(int))

with open("$OUTDIR/repeats/${SP}_no_rRNA_clusters.bed") as f:
    for line in f:
        chrom, start, end, te_type = line.strip().split()
        start = int(start)
        end = int(end)
        bin_start = (start // BIN) * BIN
        if start >= bin_start and end <= bin_start + BIN:
            te_bins[(chrom, bin_start)][te_type] += 1

type_sum   = defaultdict(float)
type_count = defaultdict(int)

with open("$OUTDIR/${SP}/${SP}_hic_pixels_norm.tsv") as f:
    for line in f:
        fields = line.strip().split()
        chrA = fields[0]
        binA = int(fields[1])
        chrB = fields[3]
        binB = int(fields[4])
        contacts = float(fields[6])

        if chrA == chrB:
            continue

        keyA = (chrA, binA)
        keyB = (chrB, binB)

        if keyA not in te_bins or keyB not in te_bins:
            continue

        binsA = te_bins[keyA]
        binsB = te_bins[keyB]

        for te_type in binsA.keys() & binsB.keys():
            n_pairs = binsA[te_type] * binsB[te_type]
            type_sum[te_type]   += contacts * n_pairs
            type_count[te_type] += n_pairs

with open(f"$OUTDIR/${SP}/${SP}_TEtype_interchrom_avg_no_rRNA_clusters_norm.tsv", "w") as out:
    for te_type in sorted(type_sum):
        avg = type_sum[te_type] / type_count[te_type]
        out.write(f"{te_type}\t{avg:.6f}\n")

EOF

    python run_${SP}_no_rRNA_clusters_norm.py
}

export INDIR OUTDIR TEDIR SEXCHROM SP
export -f TEcontacts
parallel --colsep '\t' 'TEcontacts {1}' :::: $MAINDIR/HiC/A2_hicExplorer/species_list.txt

conda deactivate


### SOFTWARE VERSIONS
# cooler v0.10.4
# bedtools v2.30.0
