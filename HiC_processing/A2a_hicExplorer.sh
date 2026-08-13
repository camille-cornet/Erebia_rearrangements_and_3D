#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/A1_map_HiC/
OUTDIR=$MAINDIR/HiC/A2_hicExplorer/
cd $OUTDIR

source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
conda activate hicexplorer_env

### Function to parallelize over samples
hic_processing1() {
    SAMPLE=$1
    mkdir $OUTDIR/${SAMPLE}

    # Lower resolution to 10kb (TADs), 40kb (compartments), 100kb (plotting)
    # Original .cool matrix is 1kb resolution
    # Do all for the mapq 0 and mapq30
    hicMergeMatrixBins --matrix $INDIR/${SAMPLE}/${SAMPLE}.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_10kb.cool --numBins 10
    hicMergeMatrixBins --matrix $INDIR/${SAMPLE}/${SAMPLE}.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_40kb.cool --numBins 40
    hicMergeMatrixBins --matrix $INDIR/${SAMPLE}/${SAMPLE}.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_100kb.cool --numBins 100
    hicMergeMatrixBins --matrix $INDIR/${SAMPLE}/${SAMPLE}_mapq30.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_10kb.cool --numBins 10
    hicMergeMatrixBins --matrix $INDIR/${SAMPLE}/${SAMPLE}_mapq30.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_40kb.cool --numBins 40
    hicMergeMatrixBins --matrix $INDIR/${SAMPLE}/${SAMPLE}_mapq30.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_100kb.cool --numBins 100

    # Do a diagnostic plot to set the thresholds for the correction
    hicCorrectMatrix diagnostic_plot --matrix $OUTDIR/${SAMPLE}/${SAMPLE}_10kb.cool -o $OUTDIR/${SAMPLE}/${SAMPLE}_10kb_diagnostic_plot.png
    hicCorrectMatrix diagnostic_plot --matrix $OUTDIR/${SAMPLE}/${SAMPLE}_40kb.cool -o $OUTDIR/${SAMPLE}/${SAMPLE}_40kb_diagnostic_plot.png
    hicCorrectMatrix diagnostic_plot --matrix $OUTDIR/${SAMPLE}/${SAMPLE}_100kb.cool -o $OUTDIR/${SAMPLE}/${SAMPLE}_100kb_diagnostic_plot.png
    hicCorrectMatrix diagnostic_plot --matrix $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_10kb.cool -o $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_10kb_diagnostic_plot.png
    hicCorrectMatrix diagnostic_plot --matrix $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_40kb.cool -o $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_40kb_diagnostic_plot.png
    hicCorrectMatrix diagnostic_plot --matrix $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_100kb.cool -o $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_100kb_diagnostic_plot.png
    # Visually check correction thresholds before next step!

}

### Parallelize over individuals
export INDIR OUTDIR SAMPLE
export -f hic_processing1
parallel --colsep '\t' 'hic_processing1 {1}' :::: $OUTDIR/species_list.txt

conda deactivate

### SOFTWARE VERSIONS
# hicExplorer v3.7.5
