#!/bin/bash

### Set directories
INDIR=$MAINDIR/ATACseq/A1_map_ATAC
OUTDIR=$MAINDIR/ATACseq/B2_reads_regioneR
cd $OUTDIR
mkdir $OUTDIR/read_counts
mkdir $OUTDIR/genmap_index
mkdir $OUTDIR/genmap_results

# Test if ATACseq reads are enriched at breakpoints
# Mask genes because otherwise breakpoints might be depleted in ATACseq signal just because they are gene poor by definition

atacreads() {
  CODE=$1
  SP=$2
  ATACCODE=$3

  # First need a bedgraph file with the number of ATACseq reads per 40kb genomic window
  # Remove reads overlapping genes, then count in fixed windows
  READS=$INDIR/$ATACCODE/${ATACCODE}_nodup_sorted.bam
  MASK=$OUTDIR/mask_beds/${SP}_mask.bed
  FAI=$OUTDIR/fai_files/${SP}.fasta.fai
  bedtools makewindows -g $FAI -w 40000 > $OUTDIR/read_counts/${SP}_windows_40kb.bed
  bedtools intersect -v -abam $READS -b $MASK > $OUTDIR/read_counts/${SP}_reads_intergenic.bam
  WINDOWS=$OUTDIR/read_counts/${SP}_windows_40kb.bed
  READSNOGENES=$OUTDIR/read_counts/${SP}_reads_intergenic.bam
  bedtools coverage -a $WINDOWS -b $READSNOGENES -counts > $OUTDIR/read_counts/${SP}_windows_counts.bedgraph

  # Even after removing genes, maybe low read counts in breakpoints reflect low mappability instead of closed chromatin
  # Compute mappability of windows to normalize read counts per mappable bases in windows
  source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
  conda activate genmap_env

  REF=$MAINDIR/assemblies/chromosomes_only/${SP}.fasta
  genmap index -F $REF -I $OUTDIR/genmap_index/${SP}
  genmap map -K 50 -E 2 -I $OUTDIR/genmap_index/${SP} -O $OUTDIR/genmap_results/${SP} -t -w -bg

  conda deactivate

  # Find uniquely mappable regions
  awk '$4==1' $OUTDIR/genmap_results/${SP}.bedgraph > $OUTDIR/genmap_results/${SP}_uniquely_mappable.bed

  # remove genes from the uniquely-mappable regions so the denominator is intergenic only
  bedtools subtract \
    -a $OUTDIR/genmap_results/${SP}_uniquely_mappable.bed -b $MASK \
    > $OUTDIR/genmap_results/${SP}_uniquely_mappable_intergenic.bed

  # calculate the proportion of each 40kb window that is uniquely mapping AND intergenic 
  bedtools coverage \
    -a $OUTDIR/read_counts/${SP}_windows_40kb.bed \
    -b $OUTDIR/genmap_results/${SP}_uniquely_mappable_intergenic.bed \
    > $OUTDIR/read_counts/${SP}_windows_mappability_intergenic.txt    

  # Normalize intergenic ATAC counts per intergenic mappable bp   
  paste \
      $OUTDIR/read_counts/${SP}_windows_counts.bedgraph \
      $OUTDIR/read_counts/${SP}_windows_mappability_intergenic.txt | awk 'BEGIN{OFS="\t"}{
      counts=$4;
      mappable_bp=$9;
      norm=(mappable_bp>0)?counts/mappable_bp:0;
      print $1,$2,$3,norm
  }' \
  > $OUTDIR/read_counts/${SP}_windows_counts_normalized.bedgraph

  # Keep windows whose intergenic-mappable fraction >= 0.3 
  awk '$7 >= 0.3 {print $1"\t"$2"\t"$3}' \
    $OUTDIR/read_counts/${SP}_windows_mappability_intergenic.txt > $OUTDIR/read_counts/${SP}_windows_intergenic_ge03.bed

  bedtools intersect \
    -a $OUTDIR/read_counts/${SP}_windows_counts_normalized.bedgraph \
    -b $OUTDIR/read_counts/${SP}_windows_intergenic_ge03.bed -u > $OUTDIR/read_counts/${SP}_windows_counts_normalized_ge03.bedgraph

  bedtools intersect \
    -a bp_beds/${CODE}_bp_nogenes.bed \
    -b $OUTDIR/read_counts/${SP}_windows_intergenic_ge03.bed \
    -u > bp_beds/${CODE}_bp_nogenes_ge03.bed

  # Write the R script to run regioneR
  cat > "run_${CODE}.R" <<EOF

library(regioneR)
library(rtracklayer)
library(dplyr)
library(plyranges)

mc.set.seed=FALSE
set.seed(123)

### Testing if the breakpoints have different nb of ATACseq reads mapped
### A = breakpoints (that is the region to randomize)
### universe = genomic windows

# Identify species
code <- "$CODE"
sp <- "$SP"

# Directory paths
bp_dir <- "bp_beds"
count_dir <- "read_counts"
out_dir <- "results_bp_reads"

# --- Load data ---
breaks <- toGRanges(file.path(bp_dir, paste0(code, "_bp_nogenes_ge03.bed")))
count <- import(file.path(count_dir, paste0(sp, "_windows_counts_normalized_ge03.bedgraph")))

# --- Run test for difference in number of reads between bp and in whole genome ---

  pt <- permTest(
    A = breaks,
    universe = count,
    x = count,
    randomize.function = resampleRegions,
    evaluate.function = meanInRegions,
    ntimes = 10000,
    per.chromosome = TRUE
  )

  results <- data.frame(
    pval = pt\$meanInRegions\$pval,
    zscore = pt\$meanInRegions\$zscore,
    observed = pt\$meanInRegions\$observed,
    expected = mean(pt\$meanInRegions\$permuted),
    stringsAsFactors = FALSE
  )

# --- Save results as TSV ---
out_path <- file.path(out_dir, paste0(code, "_results_read_count_regioneR.tsv"))
write.table(results, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

# Save full permutation distribution for plot later
perm_df <- data.frame(
  permuted_bp_overlap = as.numeric(pt\$meanInRegions\$permuted)
)
perm_path <- file.path(out_dir, paste0(code, "_read_count_permuted_distribution_regioneR.tsv"))
write.table(perm_df, perm_path, sep = "\t", quote = FALSE, row.names = FALSE)

EOF

  conda activate R_env

  Rscript ./run_${CODE}.R | tee logs/regioneR_${CODE}.log

  conda deactivate

}

### Parallelize over individuals
export INDIR OUTDIR 
export -f atacreads
parallel --colsep '\t' 'atacreads {1} {2} {3}' :::: $OUTDIR/bp_atac_list.txt

### Make one result file per test with all species
# write header
echo -e "species\tpval\tzscore\tobserved\texpected" > results_atacseq_read_count_regioneR.tsv
# loop over species
cut -f1 bp_atac_list.txt | while read CODE; do
    FILE="results_bp_reads/${CODE}_results_read_count_regioneR.tsv"
    # skip header and prepend species name
    tail -n +2 "$FILE" | awk -v s="$CODE" 'BEGIN{OFS="\t"}{print s,$0}'
done >> results_atacseq_read_count_regioneR.tsv

### SOFTWARE VERSIONS
# genmap v1.3.0
# regioneR v1.38.0
# GenomicRanges v1.56.2
