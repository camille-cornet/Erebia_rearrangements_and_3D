#!/bin/bash

### Set directories
INDIR=$MAINDIR/ATACseq/A1_map_ATAC
OUTDIR=$MAINDIR/ATACseq/A3_peak_calling

source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
conda activate macs3_env

while IFS=$'\t' read -r SAMPLE SAMPLEREF COL3 COL4 COL5 SIZE; do
    mkdir -p $OUTDIR/${SAMPLE}
    cd $OUTDIR/${SAMPLE}

    ### MACS3 needs the bam files sorted by name
    samtools sort -n -o $INDIR/${SAMPLE}/${SAMPLE}_nodup_sorted.bam $INDIR/${SAMPLE}/${SAMPLE}_nodup.bam

    macs3 callpeak --keep-dup all --nomodel --broad --broad-cutoff 0.1 -q 0.1 \
        -f BAMPE -t $INDIR/${SAMPLE}/${SAMPLE}_nodup_sorted.bam \
        -g $SIZE --name $SAMPLE --outdir $OUTDIR/${SAMPLE} --cutoff-analysis

    ### Calculate the FRiP score (fraction of reads in peaks, default 1 nucleotide overlap is enough)
    ### Need a SAF file for featureCounts
    echo "GeneID    Chr Start   End Strand" > $OUTDIR/${SAMPLE}/${SAMPLE}_peaks.saf
    awk '{OFS = "\t"} {print "Interval_"NR,$1,$2,$3,"."}' \
        $OUTDIR/${SAMPLE}/${SAMPLE}_peaks.broadPeak >> $OUTDIR/${SAMPLE}/${SAMPLE}_peaks.saf

    featureCounts -F SAF -p \
        -a $OUTDIR/${SAMPLE}/${SAMPLE}_peaks.saf \
        -o $OUTDIR/${SAMPLE}/${SAMPLE}_feature_Counts.txt \
        $INDIR/${SAMPLE}/${SAMPLE}_nodup_sorted.bam \
        2> $OUTDIR/${SAMPLE}/${SAMPLE}_featureCounts.log
    cat $OUTDIR/${SAMPLE}/${SAMPLE}_featureCounts.log

done < $MAINDIR/ATACseq/samples_atac.txt

conda deactivate

### SOFTWARE VERSIONS
# macs3 v.3.0.3
# featureCounts v2.1.1