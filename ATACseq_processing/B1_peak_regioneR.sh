#!/bin/bash

### Set directories
INDIR=$MAINDIR/ATACseq/A3_peak_calling
OUTDIR=$MAINDIR/ATACseq/B1_peak_regioneR
cd $OUTDIR

# Test if ATACseq peaks are enriched at bp
# Mask genes because otherwise breakpoints might be depleted in ATACseq signal just because they are gene poor by definition

atacreads() {
  CODE=$1
  SP=$2
  ATACCODE=$3
  
  cat > "run_${CODE}.R" <<EOF

library(regioneR)
library(rtracklayer)
library(dplyr)

### Testing if ATACseq peaks are enriched at breakpoints
### A = breakpoints (that is the region to randomize)
### B = peaks

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
bp_code <- "$CODE"
sp <- "$SP"
ATAC_code <- "$ATACCODE"

# Directory paths
peak_dir <- "peak_beds"
bp_dir <- "bp_beds"
mask_dir <- "mask_beds"
out_dir <- "results_bp_peaks"
rds_saved <- "rds_saved"

# --- Load data ---
peaks <- toGRanges(file.path(peak_dir, paste0(ATAC_code, "_peaks.broadPeak")))
breaks  <- toGRanges(file.path(bp_dir, paste0(bp_code, "_bp_nogenes.bed")))
mask <- toGRanges(file.path(mask_dir, paste0(sp, "_mask.bed")))
genome <- get(paste0(sp, "_genome"))

# --- Run enrichment based on number of basepairs that overlap ---

  basepairOverlap <- function(A, B, ...) {
    if (length(A) == 0 || length(B) == 0) return(0)
    ov <- GenomicRanges::intersect(A, B, ignore.strand = TRUE)
    if (length(ov) == 0) return(0)
    sum(width(ov))
  }

  pt <- permTest(
    A = breaks,
    B = peaks,
    randomize.function = randomizeRegions,
    evaluate.function = basepairOverlap,
    ntimes = 10000,
    genome = genome,
    mask = mask,
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
out_path <- file.path(out_dir, paste0(bp_code, "_results_peak_nuc.tsv"))
saveRDS(pt, file = file.path("rds_saved", paste0(bp_code, "_peaks.rds")))
write.table(results, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save full permutation distribution for plot later ---
perm_df <- data.frame(
  permuted_bp_overlap = as.numeric(pt\$basepairOverlap\$permuted)
)
perm_path <- file.path(out_dir, paste0(bp_code, "_nuc_atac_distribution_regioneR.tsv"))
write.table(perm_df, perm_path, sep = "\t", quote = FALSE, row.names = FALSE)

EOF

  source $SOFTWAREDIR/miniconda3/etc/profile.d/conda.sh
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
echo -e "species\tpval\tzscore\tobserved\texpected" > results_atacseq_regioneR_nuc.tsv
# loop over species
cut -f1 bp_atac_list.txt | while read CODE; do
    FILE="results_bp_peaks/${CODE}_results_peak_nuc.tsv"
    # skip header and prepend species name
    tail -n +2 "$FILE" | awk -v s="$CODE" 'BEGIN{OFS="\t"}{print s,$0}'
done >> results_atacseq_regioneR_nuc.tsv

### SOFTWARE VERSIONS
# regioneR v1.38.0
# GenomicRanges v1.56.2
