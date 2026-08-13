#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/A2_hicExplorer
OUTDIR=$MAINDIR/HiC/A3_callTADs

tads() {
    SAMPLE=$1
    mkdir $OUTDIR/${SAMPLE}

    source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
    conda activate hicexplorer_env

    ### Matrix corrected but not normalized, with mapq30 and 10kb
    hicFindTADs --matrix $INDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_10kb.cool --correctForMultipleTesting fdr \
        --outPrefix $OUTDIR/${SAMPLE}/${SAMPLE}_10kb_TADs --thresholdComparisons 0.01 --delta 0.01

    conda deactivate

    conda activate clodius_env

    ### Aggregate bed for visualisation in HiGlass
    clodius aggregate bedgraph          \
        --output-file $OUTDIR/${SAMPLE}/${SAMPLE}_10kb_TADs_score.hitile \
        --chromosome-col 1              \
        --from-pos-col 2                \
        --to-pos-col 3                  \
        --value-col 4                   \
        --chromsizes-filename $MAINDIR/HiC/A1_map_HiC/${SAMPLE}/ref.sizes \
        --nan-value NA                  \
        --no-header $OUTDIR/${SAMPLE}/${SAMPLE}_10kb_TADs_score.bedgraph
    clodius aggregate bedgraph          \
        --output-file $OUTDIR/${SAMPLE}/${SAMPLE}_20kb_TADs_score.hitile \
        --chromosome-col 1              \
        --from-pos-col 2                \
        --to-pos-col 3                  \
        --value-col 4                   \
        --chromsizes-filename $MAINDIR/HiC/A1_map_HiC/${SAMPLE}/ref.sizes \
        --nan-value NA                  \
        --no-header $OUTDIR/${SAMPLE}/${SAMPLE}_20kb_TADs_score.bedgraph
    clodius aggregate bedpe \
        --output-file $OUTDIR/${SAMPLE}/${SAMPLE}_10kb_TADs_domains_agreg.bed2db \
        --chromsizes-filename $MAINDIR/HiC/A1_map_HiC/${SAMPLE}/ref.sizes \
        --chr1-col 1 --chr2-col 1 \
        --from1-col 2 --to1-col 3 \
        --from2-col 2 --to2-col 3 \
        --no-header $OUTDIR/${SAMPLE}/${SAMPLE}_10kb_TADs_domains.bed
    clodius aggregate bedpe \
        --output-file $OUTDIR/${SAMPLE}/${SAMPLE}_20kb_TADs_domains_agreg.bed2db \
        --chromsizes-filename $MAINDIR/HiC/A1_map_HiC/${SAMPLE}/ref.sizes \
        --chr1-col 1 --chr2-col 1 \
        --from1-col 2 --to1-col 3 \
        --from2-col 2 --to2-col 3 \
        --no-header $OUTDIR/${SAMPLE}/${SAMPLE}_20kb_TADs_domains.bed
    
    conda deactivate
}

### Parallelize over individuals
export INDIR OUTDIR SAMPLE
export -f tads
parallel --colsep '\t' 'tads {1}' :::: $MAINDIR/HiC/A2_hicExplorer/species_list.txt

### Get a file with all species average TAD size per species
while read SAMPLE; do
    bedfile="$OUTDIR/${SAMPLE}/${SAMPLE}_10kb_TADs_domains.bed"

    avg=$(awk '{sum+=($3-$2); n++} END{print sum/n}' "$bedfile")

    printf "%s\t%s\n" "$SAMPLE" "$avg"

done < $MAINDIR/HiC/A2_hicExplorer/species_list.txt \
> $OUTDIR/average_TAD_size.tsv

### Get average TAD size per species per chromosome
while read SAMPLE; do
    bedfile="$OUTDIR/${SAMPLE}/${SAMPLE}_10kb_TADs_domains.bed"

    awk -v sample="$SAMPLE" '
        {
            chrom = $1
            size  = $3 - $2
            sum[chrom] += size
            n[chrom]++
        }
        END {
            for (chrom in sum)
                printf "%s\t%s\t%.2f\n", sample, chrom, sum[chrom] / n[chrom]
        }
    ' "$bedfile"

done < $MAINDIR/HiC/A2_hicExplorer/species_list.txt \
> $OUTDIR/average_TAD_size_per_chrom.tsv

### SOFTWARE VERSIONS
# hicExplorer v3.7.5
# clodius v0.14.3
