library(rtracklayer)
library(karyoploteR)
library(tidyverse)

species <- c("Eligea", "X3531", "Eepiphron", "Egorge", "Emedusa", "Epharte",
             "X3258", "X3737")

# Display names
display_names <- ifelse(
  grepl("^E", species),
  paste0("E. ", sub("^E", "", species)),
  species
)
display_names[species == "X3531"] <- "E. cassioides"
display_names[species == "X3737"] <- "E. tyndarus"
display_names[species == "X3258"] <- "E. nivalis"

# Sex chromosomes to exclude
sex_chr <- read.table("sexchrom.txt", header = FALSE,
                      col.names = c("species", "chr1", "chr2", "chr3")) |>
  pivot_longer(-species, values_to = "chr") |>
  filter(chr != "NONE") |>
  mutate(chr = gsub("^scaffold_", "", chr)) |>
  group_by(species) |>
  summarise(sex_chrs = list(chr))

# Load genomes
genome_files <- paste0("genomes/", species, "_genome.txt")
genomes <- lapply(genome_files, function(f) {
  gr <- toGRanges(f)
  seqlevels(gr) <- gsub("^scaffold_", "", seqlevels(gr))
  gr
})

# Load ATAC tracks
atac_files <- paste0("ATACseq_reads/", species, "_windows_counts_normalized_ge03.bedgraph")
atac_tracks <- lapply(atac_files, function(f) {
  bed <- import(f, format = "BEDGraph")
  seqlevels(bed) <- gsub("^scaffold_", "", seqlevels(bed))
  bed
})

# Plot
for (i in seq_along(species)) {

  # Get sex chromosomes for this species
  sex_idx <- which(sex_chr$species == species[i])
  sex_to_exclude <- if (length(sex_idx) > 0) sex_chr$sex_chrs[[sex_idx]] else character(0)

  # Filter genome to autosomes
  chr_df <- as.data.frame(genomes[[i]])
  chr_df <- chr_df[!chr_df$seqnames %in% sex_to_exclude, ]
  autosomes <- as.character(chr_df$seqnames)
  chr_seqlengths <- setNames(chr_df$end, chr_df$seqnames)

  # Filter ATAC to autosomes
  atac_filt <- atac_tracks[[i]][!as.character(seqnames(atac_tracks[[i]])) %in% sex_to_exclude]

  # Average to 500kb windows
  tiles <- tileGenome(chr_seqlengths,
                      tilewidth = 1000000,
                      cut.last.tile.in.chrom = TRUE)
  win_gr <- unlist(GRangesList(tiles))

  overlaps <- findOverlaps(win_gr, atac_filt)
  avg_counts <- tapply(atac_filt$score[subjectHits(overlaps)],
                       queryHits(overlaps), mean)

  win_gr$score <- 0
  win_gr$score[as.integer(names(avg_counts))] <- avg_counts

  # Plot
  pdf(paste0("ATACseq_reads_plots/ATAC_seq_reads_", species[i], ".pdf"), width = 15, height = 2)
  pp <- getDefaultPlotParams(plot.type = 4)
  pp$data1inmargin <- 0
  pp$bottommargin <- 50
  pp$leftmargin <- 0.12

  kp <- plotKaryotype(genome = genomes[[i]], chromosomes = autosomes, plot.type = 4,
                      plot.params = pp, cex = 0.6)
  kpAddLabels(kp, labels = bquote(italic(.(display_names[i]))),
              srt = 0, pos = 4, cex = 1, r0 = 1.15)
  kpAddCytobandsAsLine(kp)
  kpBars(kp, win_gr, y1 = win_gr$score / max(win_gr$score), col = "#114B5F", border = NA)
  kpAxis(kp, ymax = max(win_gr$score), cex = 0.7)
  kpAddLabels(kp, labels = "ATAC-seq reads",
              side = "left", cex = 1.2, srt = 90,
              r0 = 0.9, r1 = 1,
              label.margin = 0.05)
  dev.off()
}
