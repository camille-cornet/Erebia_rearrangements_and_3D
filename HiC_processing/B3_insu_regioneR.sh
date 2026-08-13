#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/
OUTDIR=$MAINDIR/HiC/B3_insu_regioneR
cd $OUTDIR

source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
conda activate R_env

# Remove the sex chromosomes
while IFS=$'\t' read -r SP CHROM1 CHROM2 CHROM3; do
  cat tads_bedgraphs/${SP}_10kb_TADs_score.bedgraph | grep -v "${CHROM1}\b" | grep -v "${CHROM2}\b" | grep -v "${CHROM3}\b" \
  > tads_bedgraphs/${SP}_10kb_TADs_score_nosex.bedgraph
done < sexchromlist.txt

# Function to parallelize over species
regioneR() {
  CODE=$1
  SP=$2

  ### Create the R script for each species pair for insulation score (including fusions, fissions, heteroz Egorge)
cat > "run_${CODE}_bp_insu.R" <<EOF

library(regioneR)
library(rtracklayer)
library(dplyr)
library(plyranges)

mc.set.seed=FALSE
set.seed(123)

### Testing if the breakpoints have different insulation scores
### A = breakpoints (that is the region to randomize)
### universe = genomic windows

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
code <- "$CODE"
sp <- "$SP"

# Directory paths
bp_dir <- "bp_persp_full"
insu_dir <- "tads_bedgraphs"
out_dir <- "results_bp_insu"

# --- Load data ---
breaks <- toGRanges(file.path(bp_dir, paste0(code, "_bp.bed")))
tads_insu <- import(file.path(insu_dir, paste0(code, "_10kb_TADs_score_nosex.bedgraph")))
genome <- get(paste0(sp, "_genome"))

# --- Make sure that sex chromosomes are removed from the universe (they are absent from fai files) ---
# And because we test per chromosome, we should restrict the analysis to the chromosomes that do rearrange 
common <- intersect(seqlevels(breaks), seqlevels(genome))
breaks  <- keepSeqlevels(breaks,  common, pruning.mode="coarse")
genome  <- keepSeqlevels(genome,  common, pruning.mode="coarse")

# --- Run test for difference in PC1 Score between bp and in whole genome ---

  pt <- permTest(
    A = breaks,
    x = tads_insu,
    randomize.function = randomizeRegions,
    evaluate.function = meanInRegions,
    genome = genome,
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
out_path <- file.path(out_dir, paste0(code, "_results_10kb_bp_insu_regioneR.tsv"))
write.table(results, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save full permutation distribution for plot later ---
perm_df <- data.frame(
  permuted_bp_overlap = as.numeric(pt\$meanInRegions\$permuted)
)
perm_path <- file.path(out_dir, paste0(code, "_insu_permuted_distribution_regioneR.tsv"))
write.table(perm_df, perm_path, sep = "\t", quote = FALSE, row.names = FALSE)

EOF

  Rscript ./run_${CODE}_bp_insu.R | tee logs/regioneR_${CODE}_bp_insu.log
}

### Parallelize over individuals
export CODE SP
export -f regioneR
parallel --colsep '\t' 'regioneR {1} {2}' :::: $MAINDIR/HiC/B4_scores_regioneR/species_pairs_list.txt

### Make one result file per test with all species
# write header
echo -e "species\tpval\tzscore\tobserved\texpected" > results_10kb_bp_insu_regioneR.tsv
# loop over species
cut -f1 species_pairs_list.txt | while read CODE; do
    FILE="results_bp_insu/${CODE}_results_10kb_bp_insu_regioneR.tsv"
    # skip header and prepend species name
    tail -n +2 "$FILE" | awk -v s="$CODE" 'BEGIN{OFS="\t"}{print s,$0}'
done >> results_10kb_bp_insu_regioneR.tsv

conda deactivate

### SOFTWARE VERSIONS
# regioneR v1.38.0
# GenomicRanges v1.56.2
