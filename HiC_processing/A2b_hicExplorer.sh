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
    RES=$2
    FILTER=$3
    mkdir $OUTDIR/${SAMPLE}

    # Correct matrix with ICE method (based on thresholds in plots generated in A2a)
    hicCorrectMatrix correct --matrix $OUTDIR/${SAMPLE}/${SAMPLE}_${RES}.cool --correctionMethod ICE --filterThreshold ${FILTER} 5 --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_corr_${RES}.cool
    hicCorrectMatrix correct --matrix $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_${RES}.cool --correctionMethod ICE --filterThreshold ${FILTER} 5 --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_${RES}.cool
}

### Parallelize over individuals
export INDIR OUTDIR
export -f hic_processing1
parallel --colsep '\t' 'hic_processing1 {1} {2} {3}' :::: $OUTDIR/species_list_filters.txt

### Normalize to a 0-1 range for comparisons between species
MATRICES=($OUTDIR/*/*_corr_10kb.cool)
OUTFILES=()
for FILE in "${MATRICES[@]}"; do
    OUTFILES+=("${FILE/_corr_10kb.cool/_norm_corr_10kb.cool}")
done
hicNormalize --matrices ${MATRICES[@]} --normalize norm_range --outFileName ${OUTFILES[@]}
MATRICES=($OUTDIR/*/*_corr_40kb.cool)
OUTFILES=()
for FILE in "${MATRICES[@]}"; do
    OUTFILES+=("${FILE/_corr_40kb.cool/_norm_corr_40kb.cool}")
done
hicNormalize --matrices ${MATRICES[@]} --normalize norm_range --outFileName ${OUTFILES[@]}
MATRICES=($OUTDIR/*/*_corr_100kb.cool)
OUTFILES=()
for FILE in "${MATRICES[@]}"; do
    OUTFILES+=("${FILE/_corr_100kb.cool/_norm_corr_100kb.cool}")
done
hicNormalize --matrices ${MATRICES[@]} --normalize norm_range --outFileName ${OUTFILES[@]}
MATRICES=($OUTDIR/*/*_mapq30_corr_10kb.cool)
OUTFILES=()
for FILE in "${MATRICES[@]}"; do
    OUTFILES+=("${FILE/_mapq30_corr_10kb.cool/_mapq30_norm_corr_10kb.cool}")
done
hicNormalize --matrices ${MATRICES[@]} --normalize norm_range --outFileName ${OUTFILES[@]}
MATRICES=($OUTDIR/*/*_mapq30_corr_40kb.cool)
OUTFILES=()
for FILE in "${MATRICES[@]}"; do
    OUTFILES+=("${FILE/_mapq30_corr_40kb.cool/_mapq30_norm_corr_40kb.cool}")
done
hicNormalize --matrices ${MATRICES[@]} --normalize norm_range --outFileName ${OUTFILES[@]}
MATRICES=($OUTDIR/*/*_mapq30_corr_100kb.cool)
OUTFILES=()
for FILE in "${MATRICES[@]}"; do
    OUTFILES+=("${FILE/_mapq30_corr_100kb.cool/_mapq30_norm_corr_100kb.cool}")
done
hicNormalize --matrices ${MATRICES[@]} --normalize norm_range --outFileName ${OUTFILES[@]}

### Function to parallelize over samples
hic_processing2() {
    SAMPLE=$1
    # Convert to ginteraction to have readable matrix (only for 40kb and 100kb)
    hicConvertFormat --matrices $OUTDIR/${SAMPLE}/${SAMPLE}_norm_corr_40kb.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_norm_corr_40kb.ginteractions --inputFormat cool --outputFormat ginteractions
    hicConvertFormat --matrices $OUTDIR/${SAMPLE}/${SAMPLE}_norm_corr_100kb.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_norm_corr_100kb.ginteractions --inputFormat cool --outputFormat ginteractions
    hicConvertFormat --matrices $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_norm_corr_40kb.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_norm_corr_40kb.ginteractions --inputFormat cool --outputFormat ginteractions
    hicConvertFormat --matrices $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_norm_corr_100kb.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_norm_corr_100kb.ginteractions --inputFormat cool --outputFormat ginteractions
    hicConvertFormat --matrices $OUTDIR/${SAMPLE}/${SAMPLE}_corr_40kb.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_corr_40kb.ginteractions --inputFormat cool --outputFormat ginteractions
    hicConvertFormat --matrices $OUTDIR/${SAMPLE}/${SAMPLE}_corr_100kb.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_corr_100kb.ginteractions --inputFormat cool --outputFormat ginteractions
    hicConvertFormat --matrices $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_40kb.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_40kb.ginteractions --inputFormat cool --outputFormat ginteractions
    hicConvertFormat --matrices $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_100kb.cool --outFileName $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_corr_100kb.ginteractions --inputFormat cool --outputFormat ginteractions

}

### Parallelize over individuals
export INDIR OUTDIR SAMPLE
export -f hic_processing2
parallel --colsep '\t' 'hic_processing2 {1}' :::: $OUTDIR/species_list.txt

conda deactivate

### SOFTWARE VERSIONS
# hicExplorer version 3.7.5
