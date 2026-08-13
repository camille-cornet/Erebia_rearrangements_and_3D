#!/bin/bash

### Set directories
INDIR=$MAINDIR/ATACseq/A1_map_ATAC
OUTDIR=$MAINDIR/ATACseq/A2_ATAC_QC
cd $OUTDIR

atacqc() {
	SAMPLE=$1
  SAMPLEREF=$2
  SAMPLEANNOT=$3
  CHROM=$5
  REF=$MAINDIR/assemblies/chromosomes_only/${SAMPLEREF}.fasta
  echo "### Processing sample ${SAMPLE}"
  mkdir $OUTDIR/${SAMPLE}

  ### Index bam
  samtools index $INDIR/${SAMPLE}/${SAMPLE}_nodup.bam

  ### Fragment length distribution
  java -Xmx32G -jar $SOFTWAREDIR/picard.jar CollectInsertSizeMetrics -I $INDIR/${SAMPLE}/${SAMPLE}_nodup.bam \
  -O $OUTDIR/${SAMPLE}/${SAMPLE}_fraglen.stats -H $OUTDIR/${SAMPLE}/${SAMPLE}_fraglen.pdf -M 0.5

  ### Generate R scripts for next step
  cat <<EOF > run_ATAC_QC_${SAMPLE}.R
library(ATACseqQC)
library(Rsamtools)
library(BSgenome)
library(BSgenomeForge)
library(Biostrings)
library(GenomeInfoDb)
library(rtracklayer)
library(GenomicRanges)
library(ChIPpeakAnno)

bamFile <- "/lustre/scratch122/tol/teams/blaxter/projects/lepidoptera/erebia/ATACseq/A1_map_ATAC/${SAMPLE}/${SAMPLE}_nodup.bam"

possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))
bamTop100 <- scanBam(BamFile(bamFile, yieldSize = 100), param = ScanBamParam(tag = possibleTag))[[1]]\$tag
tags <- names(bamTop100)[lengths(bamTop100) > 0]
bam <- readBamFile(bamFile, tag = tags, asMates = TRUE, bigFile = TRUE)

outPath <- "${SAMPLE}/splitBam"
dir.create(outPath)
saveRDS(bam, file = "${SAMPLE}/splitBam/bam.rds", ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL)

gtf_file <- "${SAMPLEANNOT}/braker_transcripts.gtf"
txs <- import(gtf_file, format = "gtf")

fasta <- "${SAMPLEREF}.fasta"
ref <- readDNAStringSet(fasta, format = "fasta")
output_dir <- "${SAMPLEREF}"
dir.create(output_dir)
for (seqname in names(ref)) {
  writeXStringSet(ref[seqname], filepath = file.path(output_dir, paste0(seqname, ".fa")))
}
seed_content <- paste(
  "Package: BSgenome.Erebia.${SAMPLEREF}\n",
  "Title: BSgenome for Erebia ${SAMPLEREF}\n",
  "Description: Full genome sequences for ${SAMPLEREF} as a BSgenome data package\n",
  "Version: 1.0.0\n",
  "organism: Erebia\n",
  "common_name: Erebia\n",
  "provider: Unine\n",
  "provider_version: 1.0\n",
  "release_date: 2025-03-20\n",
  "genome: ${SAMPLEREF}\n",
  "seqnames: c(", paste(sprintf('\"%s\"', names(ref)), collapse = ", "), ")\n",
  "seqs_srcdir: ${SAMPLEREF}\n",
  "circ_seqs: character(0)\n",
  "BSgenomeObjname: Erebia${SAMPLEREF}\n",
  "organism_biocview: Erebia\n",
  sep = ""
)
writeLines(seed_content, "${SAMPLE}/seed_file.txt")
seed_file <- "${SAMPLE}/seed_file.txt"
forgeBSgenomeDataPkg(seed_file, destdir = "${SAMPLE}", replace = TRUE)
install.packages("${SAMPLE}/BSgenome.Erebia.${SAMPLEREF}", repos = NULL, type = "source")
library(BSgenome.Erebia.${SAMPLEREF})
genome <- getBSgenome("BSgenome.Erebia.${SAMPLEREF}")

split_bam <- splitGAlignmentsByCut(bam, txs=txs, genome=genome, outPath = outPath)
saveRDS(split_bam, file = "${SAMPLE}/splitBam/split_bam.rds", ascii = FALSE, version = NULL,compress = TRUE, refhook = NULL)
#split_bam <- readRDS(file = "${SAMPLE}/splitBam/split_bam.rds")

bamFiles <- file.path(outPath, c("NucleosomeFree.bam", "mononucleosome.bam", "dinucleosome.bam", "trinucleosome.bam"))
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)
librarySize <- estLibSize(bamFiles)
NTILE <- 101
dws <- ups <- 2000
sigs <- enrichedFragments(gal=split_bam[c("NucleosomeFree", "mononucleosome", "dinucleosome", "trinucleosome")],
                          TSS=TSS,librarySize=librarySize, TSS.filter=0.5, n.tile = NTILE, upstream = ups, downstream = dws, seqlev = paste0("scaffold_", c(1:${CHROM})))
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
pdf("${SAMPLE}/${SAMPLE}_Heatmap_splitbam.pdf")
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws), zeroAt=.5, n.tile=NTILE)
dev.off()

out <- featureAlignedDistribution(sigs, reCenterPeaks(TSS, width=ups+dws), zeroAt=.5, n.tile=NTILE, type="l", ylab="Averaged coverage")
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)
pdf("${SAMPLE}/${SAMPLE}_TSSprofile_splitbam.pdf")
        matplot(out, type="l", xaxt="n", xlab="Position (bp)", ylab="Fraction of signal")
        axis(1, at=seq(0, 100, by=10)+1, labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
        abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()

EOF

  source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
  conda activate R_env

  Rscript ./run_ATAC_QC_${SAMPLE}.R | tee ATAC_QC_${SAMPLE}.log

  conda deactivate

  ### Also use ataqv to generate a common html report for all samples
  cut $MAINDIR/assemblies/chromosomes_only/${SAMPLEREF}.fasta.fai -f 1 > $OUTDIR/${SAMPLE}/${SAMPLE}_chrom_list.txt
  CHRLIST=$OUTDIR/${SAMPLE}/${SAMPLE}_chrom_list.txt
  mkdir $OUTDIR/tss
  cp $MAINDIR/Annotations/C3_braker_reads/*/*tss.bed $OUTDIR/tss
  cp $MAINDIR/Annotations/C3b_braker_psyche/*/*tss.bed $OUTDIR/tss
  TSS=$OUTDIR/tss/${SAMPLEANNOT}_tss.bed
  BAM=$INDIR/${SAMPLE}/${SAMPLE}_nodup.bam
  conda activate ataqv_env
  ataqv --name $SAMPLE --metrics-file $OUTDIR/${SAMPLE}/${SAMPLE}_ataqv.json --tss-file $TSS --autosomal-reference-file $CHRLIST erebia $BAM > $OUTDIR/${SAMPLE}/${SAMPLE}_ataqv.out
  conda deactivate

}

### Parallelize over individuals
export INDIR OUTDIR SAMPLE
export -f atacqc
parallel --colsep '\t' 'atacqc {1} {2} {3} {4} {5} {6}' :::: $MAINDIR/ATACseq/samples_atac.txt

### SOFTWARE VERSIONS
# bwa 0.7.17
# samtools 1.16
# Picard 3.1.1
# AtacQC 1.30.0
# Subread v2.1.1
