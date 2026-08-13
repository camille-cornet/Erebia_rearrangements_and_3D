#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/
OUTDIR=$MAINDIR/HiC/C1_plot_HiC
cd $OUTDIR

conda activate pyGenomeTracks_env

#### HiC maps for fig 1 with compartments and TADS
# Use E. cassioides (X3531) scaffold 4 which has
# One fusion breakpoint at TAD boundary
# One fission breakpoint in B comparment

make_tracks_file --trackFiles $INDIR/A2_hicExplorer/X3531/X3531_corr_100kb.cool \
    -o X3531_tracks_full.ini

make_tracks_file --trackFiles $INDIR/A2_hicExplorer/X3531/X3531_corr_40kb.cool \
    $INDIR/A4_callAB_cool/X3531/X3531_40kb_AB.bedgraph \
    $INDIR/A4_callAB_cool/X3531/X3531_40kb_A.bed \
    $INDIR/A4_callAB_cool/X3531/X3531_40kb_B.bed \
    $INDIR/B1_AB_regioneR/bp_persp/X3531_Erondoui_bp.bed \
    $INDIR/A4_callAB_cool/X3531/X3531_gene_density_10kb.bedgraph \
    -o X3531_tracks_AB.ini

make_tracks_file --trackFiles $INDIR/A2_hicExplorer/X3531/X3531_corr_10kb.cool \
    $INDIR/A3_callTADs/X3531/X3531_10kb_TADs_domains.bed \
    $INDIR/A3_callTADs/X3531/X3531_10kb_TADs_score.bedgraph \
    $INDIR/B1_AB_regioneR/bp_persp/X3531_bp.bed \
    -o X3531_tracks_TAD.ini

pyGenomeTracks --tracks X3531_tracks_full.ini --region scaffold_4:1-58421957 -o X3531_chr4_full.png
pyGenomeTracks --tracks X3531_tracks_AB.ini --region scaffold_4:20000000-42000000 -o X3531_chr4_AB.png
pyGenomeTracks --tracks X3531_tracks_TAD.ini --region scaffold_4:21000000-24000000 -o X3531_chr4_TAD.png

### Supp fig: examples of a good (X3531 = E. cassioides) and a bad (Eepiphron) AB calling
make_tracks_file --trackFiles $INDIR/A2_hicExplorer/X3531/X3531_corr_40kb.cool \
    $INDIR/A4_callAB_cool/X3531/X3531_40kb_AB.bedgraph \
    $INDIR/A4_callAB_cool/X3531/X3531_40kb_A.bed \
    $INDIR/A4_callAB_cool/X3531/X3531_40kb_B.bed \
    $INDIR/A4_callAB_cool/X3531/X3531_gene_density_10kb.bedgraph \
    -o goodAB.ini

make_tracks_file --trackFiles $INDIR/A2_hicExplorer/Eepiphron/Eepiphron_corr_40kb.cool \
    $INDIR/A4_callAB_cool/Eepiphron/Eepiphron_40kb_AB.bedgraph \
    $INDIR/A4_callAB_cool/Eepiphron/Eepiphron_40kb_A.bed \
    $INDIR/A4_callAB_cool/Eepiphron/Eepiphron_40kb_B.bed \
    $INDIR/A4_callAB_cool/Eepiphron/Eepiphron_gene_density_10kb.bedgraph \
    -o badAB.ini

pyGenomeTracks --tracks goodAB.ini --region scaffold_4:1-58421957 -o goodAB_X3531_chr4.png
pyGenomeTracks --tracks badAB.ini --region scaffold_5:1-35217598 -o badAB_Eepiphron_chr5.png

conda deactivate

### Full HiC maps for all species
source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
conda activate hicexplorer_env

cd $OUTDIR/full_matrix

for SAMPLE in $(cat $MAINDIR/HiC/A2_hicExplorer/species_list.txt | cut -f 1); do

    # Do no plot sex chromosomes
    SEX_CHRS=$(grep "^${SAMPLE}" $OUTDIR/sexchromlist.txt | tr '\t' '\n' | grep "scaffold_")
    ALL_CHRS=$(hicInfo -m $INDIR/A2_hicExplorer/${SAMPLE}/${SAMPLE}_corr_40kb.cool 2>&1 | \
        grep "Chromosomes:length:" | \
        grep -o "scaffold_[0-9]*")
    AUTOSOMES=$(echo "$ALL_CHRS" | tr ' ' '\n' | while read chr; do
        echo "$SEX_CHRS" | grep -qx "$chr" || echo "$chr"
    done | tr '\n' ' ')

    hicPlotMatrix -m $INDIR/A2_hicExplorer/${SAMPLE}/${SAMPLE}_corr_100kb.cool \
        --outFileName ${SAMPLE}_full_matrix.pdf \
        --chromosomeOrder $AUTOSOMES \
        --colorMap YlOrRd \
        --log1p \
        --dpi 300

done

conda deactivate

### SOFTWARE VERSIONS
# pyGenomeTracks 3.9
# hicExplorer v3.7.5
