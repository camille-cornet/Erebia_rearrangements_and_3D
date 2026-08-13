#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/A3_callTADs
OUTDIR=$MAINDIR/HiC/B5_plot_TAD_boundaries
cd $OUTDIR

# Script modified from Shi et al. (2026) Dynamic reorganization of three-dimensional genome architecture during Populus diversification. Nature Ecology & Evolution
# https://github.com/jingwanglab/Populus_3D_genome_evolution/

for SP in $(cat $MAINDIR/sample_list_full.txt); do
    mkdir $OUTDIR/$SP
    TADS=$INDIR/$SP/${SP}_10kb_TADs_domains.bed
    GENES=$MAINDIR/HiC/B6_genes_TAD_regioneR/genes_beds/${SP}.bed
    REPEATS=$MAINDIR/HiC/B6_genes_TAD_regioneR/repeats_beds/${SP}_full_unknowns.bed

    # Flank is how far left and right of the TAD I want to plot
    # nf is how many bins each flanking region should be averaged over
    # nb same but for inside the TAD
    awk -v flank=50000 -v nf=20 -v nb=20 'BEGIN{OFS="\t"; tad=0}
    {
    tad++
    chrom=$1; s=$2; e=$3

    up_w = flank / nf
    for (i = 0; i < nf; i++) {
        ws = s - flank + i * up_w
        we = ws + up_w
        if (ws < 0) ws = 0
        if (we > ws) print chrom, int(ws), int(we), "tad"tad, i + 1
    }

    body_w = (e - s) / nb
    for (i = 0; i < nb; i++) {
        ws = s + i * body_w
        we = ws + body_w
        print chrom, int(ws), int(we), "tad"tad, nf + i + 1
    }

    down_w = flank / nf
    for (i = 0; i < nf; i++) {
        ws = e + i * down_w
        we = ws + down_w
        print chrom, int(ws), int(we), "tad"tad, nf + nb + i + 1
    }
    }' $TADS > $OUTDIR/$SP/windows.bed

    sort -k1,1 -k2,2n $OUTDIR/$SP/windows.bed > $OUTDIR/$SP/windows_sorted.bed

    # ---- gene density ----
    # nbins is 30, total number of bins
    sort -k1,1 -k2,2n $GENES > $OUTDIR/$SP/genes_sorted.bed
    bedtools coverage -a $OUTDIR/$SP/windows_sorted.bed -b $OUTDIR/$SP/genes_sorted.bed \
    > $OUTDIR/$SP/gene_coverage.bed

    echo "bin,gene_density" > $OUTDIR/$SP/${SP}_gene_density.csv
    awk 'BEGIN{OFS=","} {bin=$5; sum[bin]+=$NF; n[bin]++}
        END{for (b=1; b<=nbins; b++) if (n[b]>0) print b, sum[b]/n[b]}' \
    nbins=60 $OUTDIR/$SP/gene_coverage.bed \
    | sort -t, -k1,1n >> $OUTDIR/$SP/${SP}_gene_density.csv

    # ---- repeat density ----
    sort -k1,1 -k2,2n $REPEATS > $OUTDIR/$SP/repeats_sorted.bed
    bedtools coverage -a $OUTDIR/$SP/windows_sorted.bed -b $OUTDIR/$SP/repeats_sorted.bed \
    > $OUTDIR/$SP/repeats_coverage.bed

    echo "bin,repeat_density" > $OUTDIR/$SP/${SP}_repeats_density.csv
    awk 'BEGIN{OFS=","} {bin=$5; sum[bin]+=$NF; n[bin]++}
        END{for (b=1; b<=nbins; b++) if (n[b]>0) print b, sum[b]/n[b]}' \
    nbins=60 $OUTDIR/$SP/repeats_coverage.bed \
    | sort -t, -k1,1n >> $OUTDIR/$SP/${SP}_repeats_density.csv

done

### SOFTWARE VERSIONS
# bedtools v2.30.0