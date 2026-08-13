#!/bin/bash

### Set directories
INDIR=$MAINDIR/HiC/A4_callAB
OUTDIR=$MAINDIR/HiC/B1_AB_regioneR
cd $OUTDIR

source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
conda activate R_env

# Function to parallelize over species
regioneR() {
  CODE=$1
  SP=$2
  
  ### Create the R script for each species pair, A compartment (including fusions, fissions, heteroz Egorge)
cat > "run_${CODE}_A.R" <<EOF

library(regioneR)
library(rtracklayer)
library(dplyr)

mc.set.seed=FALSE
set.seed(123)

### Testing if compA are enriched at breakpoints
### A = breakpoints (that is the region to randomize)
### B = A compartment

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
comp_dir <- "comp_beds"
bp_dir <- "bp_persp"
out_dir <- "results_bp_Acomp"

# --- Load data ---
Acomp <- toGRanges(file.path(comp_dir, paste0(sp, "_40kb_A.bed")))
breaks  <- toGRanges(file.path(bp_dir, paste0(code, "_bp.bed")))
genome <- get(paste0(sp, "_genome"))

# --- Make sure that sex chromosomes are removed from the universe (they are absent from fai files) ---
# And because we test per chromosome, we should restrict the analysis to the chromosomes that do rearrange 
common <- intersect(
  seqlevels(breaks),
  intersect(seqlevels(Acomp), seqlevels(genome)))
breaks  <- keepSeqlevels(breaks,  common, pruning.mode="coarse")
Acomp   <- keepSeqlevels(Acomp,   common, pruning.mode="coarse")
genome  <- keepSeqlevels(genome,  common, pruning.mode="coarse")

# --- Custom evaluate function: test the number of bp that overlap ---
basepairOverlap <- function(A, B, ...) {
  if (length(A) == 0 || length(B) == 0) return(0)
  ov <- GenomicRanges::intersect(A, B, ignore.strand = TRUE)
  if (length(ov) == 0) return(0)
  sum(width(ov))
}

# --- Run enrichment of overlaps between A compartment and breakpoints ---

  pt <- permTest(
    A = breaks,
    B = Acomp,
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
out_path <- file.path(out_dir, paste0(code, "_results_40kb_Acomp_bpoverlap_regioneR.tsv"))
write.table(results, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save full permutation distribution for plot later ---
perm_df <- data.frame(
  permuted_bp_overlap = as.numeric(pt\$basepairOverlap\$permuted)
)
perm_path <- file.path(out_dir, paste0(code, "_Acomp_permuted_distribution_regioneR.tsv"))
write.table(perm_df, perm_path, sep = "\t", quote = FALSE, row.names = FALSE)

EOF

  ### Create the R script for each species pair, B compartment (including fusions, fissions, heteroz Egorge)
cat > "run_${CODE}_B.R" <<EOF

library(regioneR)
library(rtracklayer)
library(dplyr)

mc.set.seed=FALSE
set.seed(123)

### Testing if compA are enriched at breakpoints
### A = breakpoints (that is the region to randomize)
### B = B compartment

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
comp_dir <- "comp_beds"
bp_dir <- "bp_persp"
out_dir <- "results_bp_Bcomp"

# --- Load data ---
Bcomp <- toGRanges(file.path(comp_dir, paste0(sp, "_40kb_B.bed")))
breaks  <- toGRanges(file.path(bp_dir, paste0(code, "_bp.bed")))
genome <- get(paste0(sp, "_genome"))

# --- Make sure that sex chromosomes are removed from the universe (they are absent from fai files) ---
# And because we test per chromosome, we should restrict the analysis to the chromosomes that do rearrange 
common <- intersect(
  seqlevels(breaks),
  intersect(seqlevels(Bcomp), seqlevels(genome)))
breaks  <- keepSeqlevels(breaks,  common, pruning.mode="coarse")
Bcomp   <- keepSeqlevels(Bcomp,   common, pruning.mode="coarse")
genome  <- keepSeqlevels(genome,  common, pruning.mode="coarse")

# --- Custom evaluate function: test the number of bp that overlap ---
basepairOverlap <- function(A, B, ...) {
  if (length(A) == 0 || length(B) == 0) return(0)
  ov <- GenomicRanges::intersect(A, B, ignore.strand = TRUE)
  if (length(ov) == 0) return(0)
  sum(width(ov))
}

# --- Run enrichment of overlaps between B compartment and breakpoints ---

  pt <- permTest(
    A = breaks,
    B = Bcomp,
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
out_path <- file.path(out_dir, paste0(code, "_results_40kb_Bcomp_bpoverlap_regioneR.tsv"))
write.table(results, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save full permutation distribution for plot later ---
perm_df <- data.frame(
  permuted_bp_overlap = as.numeric(pt\$basepairOverlap\$permuted)
)
perm_path <- file.path(out_dir, paste0(code, "_Bcomp_permuted_distribution_regioneR.tsv"))
write.table(perm_df, perm_path, sep = "\t", quote = FALSE, row.names = FALSE)

EOF

  Rscript ./run_${CODE}_A.R | tee logs/regioneR_${CODE}_A.log
  Rscript ./run_${CODE}_B.R | tee logs/regioneR_${CODE}_B.log

}

### Parallelize over individuals
export CODE SP
export -f regioneR
parallel --colsep '\t' 'regioneR {1} {2}' :::: $MAINDIR/HiC/B1_AB_regioneR/species_pairs_list.txt

### Make one result file per test with all species
# write header
echo -e "species\tpval\tzscore\tobserved\texpected" > results_40kb_ABbound_regioneR_bpoverlap.tsv
# loop over species
cut -f1 species_pairs_list.txt | while read CODE; do
    FILE="results_bp_ABbound/${CODE}_results_40kb_ABbound_bpoverlap_regioneR.tsv"
    # skip header and prepend species name
    tail -n +2 "$FILE" | awk -v s="$CODE" 'BEGIN{OFS="\t"}{print s,$0}'
done >> results_40kb_ABbound_regioneR_bpoverlap.tsv

# write header
echo -e "species\tpval\tzscore\tobserved\texpected" > results_40kb_Acomp_regioneR_bpoverlap.tsv
# loop over species
cut -f1 species_pairs_list.txt | while read CODE; do
    FILE="results_bp_Acomp/${CODE}_results_40kb_Acomp_bpoverlap_regioneR.tsv"
    # skip header and prepend species name
    tail -n +2 "$FILE" | awk -v s="$CODE" 'BEGIN{OFS="\t"}{print s,$0}'
done >> results_40kb_Acomp_regioneR_bpoverlap.tsv

conda deactivate

### SOFTWARE VERSIONS
# regioneR v1.38.0
# GenomicRanges v1.56.2
