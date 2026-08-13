#!/bin/bash

### Set directories
INDIR=$MAINDIR/raw_data/ATACseq/
OUTDIR=$MAINDIR/ATACseq/A1_map_ATAC

mapatac() {
	SAMPLE=$1
    SAMPLEREF=$2
    REF=$MAINDIR/assemblies/chromosomes_only/${SAMPLEREF}.fasta
    mkdir $OUTDIR/${SAMPLE}
    
    ### trim reads before mapping
    # fastp also removes Nextera transposase adapters automatically
    fastp --trim_poly_g -l 45 -q 30 -w 16 \
    -i $INDIR/${SAMPLE}_*_R1_*.fastq.gz -o $OUTDIR/${SAMPLE}/${SAMPLE}_R1_trimmed.fastq.gz \
    -I $INDIR/${SAMPLE}_*_R2_*.fastq.gz -O $OUTDIR/${SAMPLE}/${SAMPLE}_R2_trimmed.fastq.gz

    ### Index ref
    bwa index $REF

    ### Mapping to ref (without mitochondria), removing low mapping quality, and keeping only properly paired reads and converting to bam
    bwa mem -t 8 $REF $OUTDIR/${SAMPLE}/${SAMPLE}_R1_trimmed.fastq.gz $OUTDIR/${SAMPLE}/${SAMPLE}_R2_trimmed.fastq.gz \
        | samtools view -f 0x02 -@ 8 -Sbq 20 - | samtools sort -@ 8 -o $OUTDIR/${SAMPLE}/${SAMPLE}.bam
    # Samtools view flags: https://broadinstitute.github.io/picard/explain-flags.html

    ### Removing PCR and optical duplicates
    mkdir $OUTDIR/${SAMPLE}/tmpdir
    java -Xmx30G -XX:-UseGCOverheadLimit -jar $SOFTWAREDIR/picard.jar MarkDuplicates \
    INPUT=$OUTDIR/${SAMPLE}/${SAMPLE}.bam OUTPUT=$OUTDIR/${SAMPLE}/${SAMPLE}_nodup.bam \
    METRICS_FILE=$OUTDIR/${SAMPLE}/metrics.replicates.txt TMP_DIR=$OUTDIR/${SAMPLE}/tmpdir/ \
    ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=TRUE

    ### Get stats
    samtools flagstat $OUTDIR/${SAMPLE}/${SAMPLE}_nodup.bam > $OUTDIR/${SAMPLE}/${SAMPLE}_flagstat.txt

    ### Convert to BigWig for visualisation into IGV (bamCoverage from package deeptools)
    source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
    conda activate deeptools_env
    bamCoverage -b $OUTDIR/${SAMPLE}/${SAMPLE}_nodup.bam -o $OUTDIR/${SAMPLE}/${SAMPLE}_nodup.bw --normalizeUsing CPM
    conda deactivate

}

### Parallelize over individuals
export INDIR OUTDIR SAMPLE SAMPLEREF
export -f mapatac
parallel --colsep '\t' 'mapatac {1} {2}' :::: $MAINDIR/ATACseq/samples_atac.txt

### SOFTWARE VERSIONS
# fastp v0.23.4
# bwa 0.7.17
# samtools 1.16
# Picard 3.1.1
# deeptools 3.5.6
