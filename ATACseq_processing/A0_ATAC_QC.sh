#!/bin/bash

### Set directories
INDIR=$MAINDIR/raw_data/ATACseq/
OUTDIR=$MAINDIR/ATACseq/A0_ATAC_QC

### Run FASTQC
atacqc() {
	SAMPLE=$1
    echo "### Processing sample ${SAMPLE}"
    mkdir $OUTDIR/${SAMPLE}
    conda run -n fastqc_env fastqc -o $OUTDIR/${SAMPLE} $INDIR/${SAMPLE}*R1*.fastq.gz -t 16
    conda run -n fastqc_env fastqc -o $OUTDIR/${SAMPLE} $INDIR/${SAMPLE}*R2*.fastq.gz -t 16
}

### Parallelize over individuals
export INDIR OUTDIR SAMPLE
export -f atacqc
parallel --colsep '\t' 'atacqc {1}' :::: $MAINDIR/ATACseq/samples_atac.txt

### SOFTWARE VERSIONS
# FastQC v0.12.1
