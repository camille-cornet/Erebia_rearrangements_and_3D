#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/A3_callTADs
OUTDIR=$MAINDIR/HiC/B2_TAD_regioneR
cd $OUTDIR

source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
conda activate R_env

# Function to parallelize over species
regioneR() {
  CODE=$1
  SP=$2
  
  ### Create the R script for each species pair (including fusions, fissions, heteroz Egorge)
cat > "run_${CODE}.R" <<EOF

library(regioneR)
library(rtracklayer)
library(dplyr)

mc.set.seed=FALSE
set.seed(123)

### Testing if tad boundaries are enriched at breakpoints
### A = breakpoints (that is the region to randomize)
### B = TAD boundaries

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
tad_dir <- "tad_beds"
bp_dir <- "bp_persp"
out_dir <- "results_bp_tads"

# --- Load data ---
tads <- toGRanges(file.path(tad_dir, paste0(sp, "_10kb_TADs_boundaries.bed")))
breaks  <- toGRanges(file.path(bp_dir, paste0(code, "_bp.bed")))
genome <- get(paste0(sp, "_genome"))

# --- Make sure that sex chromosomes are removed from the universe (they are absent from fai files) ---
# And because we test per chromosome, we should restrict the analysis to the chromosomes that do rearrange 
common <- intersect(
  seqlevels(breaks),
  intersect(seqlevels(tads), seqlevels(genome)))
breaks  <- keepSeqlevels(breaks,  common, pruning.mode="coarse")
tads   <- keepSeqlevels(tads,   common, pruning.mode="coarse")
genome  <- keepSeqlevels(genome,  common, pruning.mode="coarse")

# --- Custom evaluate function: test the number of bp that overlap ---
basepairOverlap <- function(A, B, ...) {
  if (length(A) == 0 || length(B) == 0) return(0)
  ov <- GenomicRanges::intersect(A, B, ignore.strand = TRUE)
  if (length(ov) == 0) return(0)
  sum(width(ov))
}

# --- Run enrichment of overlaps between TAD boundaries and breakpoints ---

  pt <- permTest(
    A = breaks,
    B = tads,
    randomize.function = randomizeRegions,
    evaluate.function = basepairOverlap,
    ntimes = 10000,
    genome = genome,
    verbose = FALSE,
    per.chromosome = TRUE
  )

  results <- data.frame(
    pval     = pt\$basepairOverlap\$pval,
    zscore   = pt\$basepairOverlap\$zscore,
    observed = pt\$basepairOverlap\$observed,
    expected = mean(as.numeric(pt\$basepairOverlap\$permuted)),
    stringsAsFactors = FALSE
  )

# --- Save results as TSV ---
out_path <- file.path(out_dir, paste0(code, "_results_10kb_tads_bpoverlap_regioneR.tsv"))
write.table(results, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save full permutation distribution for plot later ---
perm_df <- data.frame(
  permuted_bp_overlap = as.numeric(pt\$basepairOverlap\$permuted)
)
perm_path <- file.path(out_dir, paste0(code, "_tads_permuted_distribution_regioneR.tsv"))
write.table(perm_df, perm_path, sep = "\t", quote = FALSE, row.names = FALSE)

EOF

  Rscript ./run_${CODE}.R | tee logs/regioneR_${CODE}.log
}

### Parallelize over individuals
export CODE SP
export -f regioneR
parallel --colsep '\t' 'regioneR {1} {2}' :::: $MAINDIR/HiC/B2_TAD_regioneR/species_pairs_list.txt

### Make one result file per test with all species
# write header
echo -e "species\tpval\tzscore\tobserved\texpected" > results_10kb_tads_regioneR_bpoverlap.tsv
# loop over species
cut -f1 species_pairs_list.txt | while read CODE; do
    FILE="results_bp_tads/${CODE}_results_10kb_tads_bpoverlap_regioneR.tsv"
    # skip header and prepend species name
    tail -n +2 "$FILE" | awk -v s="$CODE" 'BEGIN{OFS="\t"}{print s,$0}'
done >> results_10kb_tads_regioneR_bpoverlap.tsv

conda deactivate

### SOFTWARE VERSIONS
# regioneR v1.38.0
# GenomicRanges v1.56.2
