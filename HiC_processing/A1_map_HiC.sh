#!/bin/bash

### Set directories
INDIR=$MAINDIR/assemblies/chromosomes_only/
OUTDIR=$MAINDIR/HiC/A1_map_HiC/

### SAMPLE and SAMPLESANGER are known from the job submission script

mkdir -p $OUTDIR/${SAMPLE}
cd $OUTDIR/${SAMPLE}
mkdir -p $OUTDIR/${SAMPLE}/tmp

### Index ref, both with samtools and bwa
REF=$INDIR/${SAMPLE}.fasta
bwa index -a bwtsw $REF
samtools faidx $REF

### Get chrom size file for pairtools and cooler
cut -f1,2 $INDIR/${SAMPLE}.fasta.fai | LC_ALL=C sort -k2rn,2 > $OUTDIR/${SAMPLE}/ref.sizes

### Get HiC fastq files (to delete after)
if [[ "$SAMPLE" == "Epalarica" ]]; then
    R1=$MAINDIR/HiC/Epalarica_R1.fastq.gz
    R2=$MAINDIR/HiC/Epalarica_R2.fastq.gz
elif grep -qF "$SAMPLE" $MAINDIR/HiC/sample_list_psyche.txt; then
    R1=$OUTDIR/${SAMPLE}/${SAMPLE}_R1.fastq.gz
    R2=$OUTDIR/${SAMPLE}/${SAMPLE}_R2.fastq.gz
    for FILE in $(ls $TOLDIR/$SAMPLESANGER/genomic_data/il*/hic-arima*/*.cram)
    do 
        samtools sort -n -@6 $FILE | \
            samtools fastq -@2 -s /dev/null -0 /dev/null - \
            -1 >(gzip >> $R1) \
            -2 >(gzip >> $R2)
    done
else
    RAWDIR=$MAINDIR/raw_data/HiC
    zcat $RAWDIR/iSeq/${SAMPLE}_*_R1_*.fastq.gz $RAWDIR/${SAMPLE}_L1_R1_*.fastq.gz \
    $RAWDIR/${SAMPLE}_L2_R1_*.fastq.gz | gzip > $OUTDIR/${SAMPLE}/${SAMPLE}_R1.fastq.gz
    zcat $RAWDIR/iSeq/${SAMPLE}_*_R2_*.fastq.gz $RAWDIR/${SAMPLE}_L1_R2_*.fastq.gz \
    $RAWDIR/${SAMPLE}_L2_R2_*.fastq.gz | gzip > $OUTDIR/${SAMPLE}/${SAMPLE}_R2.fastq.gz
    R1=$OUTDIR/${SAMPLE}/${SAMPLE}_R1.fastq.gz
    R2=$OUTDIR/${SAMPLE}/${SAMPLE}_R2.fastq.gz
fi

### Map to reference
bwa mem -t 8 -SP $REF $R1 $R2 | samtools view -b - > $OUTDIR/${SAMPLE}/tmp/mapped.bam

### Pairtools to remove non-valid pairs and duplicates
source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
conda activate pairtools_env
pairtools parse -c $OUTDIR/${SAMPLE}/ref.sizes --nproc-in 4 --nproc-out 4 \
    $OUTDIR/${SAMPLE}/tmp/mapped.bam --walks-policy 5unique --add-columns mapq | \
pairtools sort --nproc 8 --nproc-in 4 --nproc-out 4 | \
pairtools dedup --nproc-in 4 --nproc-out 4 \
    --output \
        >( pairtools split \
            --output-pairs $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.pairs.gz \
            --output-sam $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.bam \
        ) \
    --output-stats $OUTDIR/${SAMPLE}/${SAMPLE}.dedup.stats.txt
pairtools select "(mapq1>=30) and (mapq2>=30)" $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.pairs.gz -o $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.pairs.mapq30.gz
pairtools stats --nproc-in 8 --nproc-out 8 \
    $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.pairs.gz -o $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.stats.txt
pairtools stats --nproc-in 8 --nproc-out 8 \
    $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.pairs.mapq30.gz -o $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.mapq30.stats.txt
conda deactivate

### Mapping stats after deduplication, and also after quality filter
samtools view -@8 -q 30 -b $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.bam -o $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.mapq30.bam
samtools flagstat -@8 $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.bam >> $OUTDIR/${SAMPLE}/${SAMPLE}_flagstat.txt
samtools flagstat -@8 $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.mapq30.bam >> $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30_flagstat.txt

### Generate .cool with cooler
conda activate cooler_env
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $OUTDIR/${SAMPLE}/ref.sizes:1000 \
    $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.pairs.gz $OUTDIR/${SAMPLE}/${SAMPLE}.cool
cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 $OUTDIR/${SAMPLE}/ref.sizes:1000 \
    $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.pairs.mapq30.gz $OUTDIR/${SAMPLE}/${SAMPLE}_mapq30.cool
# Generate a multiresolution .mcool for visualisation in HiGlass
cooler zoomify $OUTDIR/${SAMPLE}/${SAMPLE}.cool
conda deactivate

### Clean
if [[ "$SAMPLE" != "Epalarica" ]]; then
    rm -f $R1 $R2
fi
rm -rf $OUTDIR/${SAMPLE}/tmp/
rm -f $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.pairs.gz
rm -f $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.pairs.mapq30.gz
rm -f $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.bam
rm -f $OUTDIR/${SAMPLE}/${SAMPLE}.nodups.mapq30.bam

### SOFTWARE VERSIONS
# bwa mem v0.7.18
# pairtools v1.1.3
# cooler v0.10.4
# samtools v1.18
