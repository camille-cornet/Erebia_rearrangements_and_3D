#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/A3_callTADs
OUTDIR=$MAINDIR/HiC/B4_TAD_features_regioneR
cd $OUTDIR

# Test if genes/all repeats/R1 are enriched in TAD boundaries

source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
conda activate R_env

# Function to parallelize over species
regioneR() {
  SP=$1
  
  ### Create the R script for each species pair (including fusions, fissions)
  # For genes and repeats
cat > "run_${SP}.R" <<EOF

library(regioneR)
library(rtracklayer)
library(dplyr)

mc.set.seed=FALSE
set.seed(123)

### Testing if genes are enriched at tad boundaries
### A = TAD boundaries (that is the region to randomize)
### B = genes

### Load genome (chr lengths) as GR to randomize over
fai_files <- list.files("fai_files", pattern = "\\\\.fai\$", full.names = TRUE)
for (file in fai_files) {
  species_name <- sub("\\\\.fasta\\\\.fai\$", "", basename(file))
  fai <- read.table(file, header = FALSE, sep = "\t",
                    col.names = c("chr", "length", "offset", "linebases", "linewidth"))
  genome_gr <- GRanges(seqnames = fai\$chr, ranges = IRanges(1, fai\$length))
  genome_var <- paste0(species_name, "_genome")
  assign(genome_var, genome_gr, envir = .GlobalEnv)
}

# Identify species
sp <- "$SP"

# Directory paths
tad_dir <- "tad_beds"
gene_dir <- "genes_beds"
repeat_dir <- "repeats_beds"
out_dir <- "results_features_tads"

# --- Load data ---
tads <- toGRanges(file.path(tad_dir, paste0(sp, "_10kb_TADs_boundaries.bed")))
repeats <- toGRanges(file.path(repeat_dir, paste0(sp, "_full_unknowns.bed")))
genes <- toGRanges(file.path(gene_dir, paste0(sp, ".bed")))
genome <- get(paste0(sp, "_genome"))

# --- Make sure that sex chromosomes are removed from the universe (they are absent from fai files) ---
# And because we test per chromosome, we should restrict the analysis to the chromosomes that do rearrange 
common <- intersect(
  seqlevels(R1),
  intersect(seqlevels(tads), seqlevels(genome)))
repeats  <- keepSeqlevels(repeats,  common, pruning.mode="coarse")
genes  <- keepSeqlevels(genes,  common, pruning.mode="coarse")
tads   <- keepSeqlevels(tads,   common, pruning.mode="coarse")
genome  <- keepSeqlevels(genome,  common, pruning.mode="coarse")

# --- Run enrichment of overlaps between repeats and TAD boundaries ---
  pt <- permTest(
    A = tads,
    B = repeats,
    randomize.function = randomizeRegions,
    evaluate.function = numOverlaps,
    ntimes = 10000,
    genome = genome,
    verbose = FALSE,
    per.chromosome = TRUE
  )
  results <- data.frame(
    pval = pt\$numOverlaps\$pval,
    zscore = pt\$numOverlaps\$zscore,
    observed = pt\$numOverlaps\$observed,
    expected = mean(pt\$numOverlaps\$permuted),
    stringsAsFactors = FALSE
  )
# --- Save results as TSV ---
out_path <- file.path(out_dir, paste0(sp, "_results_10kb_repeats_tads_regioneR.tsv"))
write.table(results, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

# --- Run enrichment of overlaps between genes and TAD boundaries ---
  pt <- permTest(
    A = tads,
    B = genes,
    randomize.function = randomizeRegions,
    evaluate.function = numOverlaps,
    ntimes = 10000,
    genome = genome,
    verbose = FALSE,
    per.chromosome = TRUE
  )
  results <- data.frame(
    pval = pt\$numOverlaps\$pval,
    zscore = pt\$numOverlaps\$zscore,
    observed = pt\$numOverlaps\$observed,
    expected = mean(pt\$numOverlaps\$permuted),
    stringsAsFactors = FALSE
  )
# --- Save results as TSV ---
out_path <- file.path(out_dir, paste0(sp, "_results_10kb_genes_tads_regioneR.tsv"))
write.table(results, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

EOF

  Rscript ./run_${SP}.R | tee logs/regioneR_${SP}.log

}

### Parallelize over individuals
export SP
export -f regioneR
parallel --colsep '\t' 'regioneR {1}' :::: $OUTDIR/species_list.txt

### Make one result file per test with all species
# write header
echo -e "species\tpval\tzscore\tobserved\texpected" > results_10kb_repeats_tads_regioneR.tsv
# loop over species
cut -f1 $MAINDIR/HiC/A2_hicExplorer/species_list.txt | while read SP; do
    FILE="results_features_tads/${SP}_results_10kb_repeats_tads_regioneR.tsv"
    # skip header and prepend species name
    tail -n +2 "$FILE" | awk -v s="$SP" 'BEGIN{OFS="\t"}{print s,$0}'
done >> results_10kb_repeats_tads_regioneR.tsv

# write header
echo -e "species\tpval\tzscore\tobserved\texpected" > results_10kb_genes_tads_regioneR.tsv
# loop over species
cut -f1 $MAINDIR/HiC/A2_hicExplorer/species_list.txt | while read SP; do
    FILE="results_features_tads/${SP}_results_10kb_genes_tads_regioneR.tsv"
    # skip header and prepend species name
    tail -n +2 "$FILE" | awk -v s="$SP" 'BEGIN{OFS="\t"}{print s,$0}'
done >> results_10kb_genes_tads_regioneR.tsv

conda deactivate

### SOFTWARE VERSIONS
# regioneR v1.38.0
# GenomicRanges v1.56.2
