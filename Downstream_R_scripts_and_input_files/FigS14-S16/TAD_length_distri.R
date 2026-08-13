library(rtracklayer)
library(karyoploteR)
library(tidyverse)

species <- c("Emeolans", "Epalarica", "Emedusa", "Etriaria", "Ealbergana",
             "Emnestra", "Egorge", "C0001", "Emelancholica", "Estirius", "Emontana",
             "Epronoe", "X3506", "X3311", "Eepiphron", "Epharte", "Echristi", "Eflavofasciata",
             "Eaethiops", "Esudetica", "Emelampus", "Erondoui", "X3531", "X3737", "X3258", "C0080",
             "C0100", "C0055", "Eligea", "Eeriphyle", "X2575", "X2576", "Ebubastis", "Emanto", "X3252",
             "Edisa", "Eembla")

# Display names
display_names <- ifelse(
  grepl("^E", species),
  paste0("E. ", sub("^E", "", species)),
  species
)
display_names[species == "C0001"] <- "E. pluto"
display_names[species == "X2575"] <- "E. euryale adyte"
display_names[species == "X2576"] <- "E. euryale isarica"
display_names[species == "X3252"] <- "E. pandrose"
display_names[species == "X3258"] <- "E. nivalis"
display_names[species == "X3531"] <- "E. cassioides"
display_names[species == "X3311"] <- "E. oeme"
display_names[species == "X3506"] <- "E. styx"
display_names[species == "X3737"] <- "E. tyndarus"
display_names[species == "C0055"] <- "E. ottomana"
display_names[species == "C0100"] <- "E. graucasica"
display_names[species == "C0080"] <- "E. calcaria"

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

# Load TAD tracks
tad_files <- paste0("tad_beds/", species, "_10kb_TADs_domains.bed")
tad_tracks <- lapply(tad_files, function(f) {
  bed <- import(f, format = "BED")
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

  # Filter TADs to autosomes
  tad_filt <- tad_tracks[[i]][!as.character(seqnames(tad_tracks[[i]])) %in% sex_to_exclude]

  # Compute windowed average TAD size
  window_size <- 1e6
  tiles <- tileGenome(chr_seqlengths,
                      tilewidth = window_size,
                      cut.last.tile.in.chrom = TRUE)
  win_gr <- unlist(GRangesList(tiles))

  overlaps <- findOverlaps(win_gr, tad_filt)
  avg_sizes <- tapply(width(tad_filt)[subjectHits(overlaps)],
                      queryHits(overlaps), mean)

  win_gr$score <- 0
  win_gr$score[as.integer(names(avg_sizes))] <- avg_sizes / 1000  # convert to kb

  # Plot
  pdf(paste0("TAD_sizes_plots/TAD_size_", species[i], ".pdf"), width = 15, height = 2)
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
  kpAddLabels(kp, labels = "TAD size (kb)",
              side = "left", cex = 1.2, srt = 90,
              r0 = 0.9, r1 = 1,
              label.margin = 0.05)
  dev.off()
}
