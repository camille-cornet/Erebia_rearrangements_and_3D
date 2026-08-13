library(tidyverse)
library(ggtext)

# Group definitions (fusions and fissions)
group1 <- c("Eligea_C0111", "Eligea_C0056", "X3529_Erondoui")

group_labels <- c(group1 = "Fissions", group2 = "Fusions")

species_labels <- c(
  "Eligea_C0111" = "*E. graucasica* (fissions)",
  "Eepiphron" = "*E. epiphron*",
  "Egorge" = "*E. gorge*",
  "Egorge_hap1" = "*E. gorge* (haplotypes)",
  "Eligea_C0056" = "*E. ottomana* (fissions)",
  "Emedusa" = "*E. medusa*",
  "Epharte" = "*E. pharte*",
  "X3336" = "*E. nivalis*",
  "X3529" = "*E. cassioides*",
  "X3529_Erondoui" = "*E. rondoui* (fissions)",
  "X3738" = "*E. tyndarus*",
  "Eaethiops" = "*E. aethiops*",
  "X3349" = "*E. oeme*",
  "C0088" = "*E. euryale*"
)

# Analysis definitions
analyses <- list(
  list(
    name     = "ATAC-seq peaks overlap",
    pattern  = "_nuc_atac_distribution_regioneR\\.tsv$",
    exclude  = "_R\\d+_nuc_atac_",
    strip    = "_nuc_atac_distribution_regioneR.tsv",
    obs_file = "results_atacseq_regioneR_nuc.tsv",
    xlab     = "ATAC-seq peaks overlap (kb)"
  ),
  list(
    name     = "ATAC-seq reads",
    pattern  = "_read_count_permuted_distribution_regioneR\\.tsv$",
    exclude  = "_R\\d+_read_count_",
    strip    = "_read_count_permuted_distribution_regioneR.tsv",
    obs_file = "results_atacseq_read_count_regioneR.tsv",
    xlab     = "Number of uniquely mapping ATAC-seq reads"
  )
)

# Main loop
for (a in analyses) {
  cat("Processing:", a$name, "\n")

  # Load permuted distributions
  perm_files <- list.files("perm_distri", pattern = a$pattern, full.names = TRUE)
  perm_files <- perm_files[!grepl(a$exclude, perm_files)]

  if (length(perm_files) == 0) {
    cat("  No files found for", a$name, "- skipping\n")
    next
  }

  perm_data <- map_dfr(perm_files, function(f) {
    sp <- sub(a$strip, "", basename(f))
    read_tsv(f, show_col_types = FALSE) %>% mutate(species = sp)
  })

  obs_data <- read_tsv(a$obs_file, show_col_types = FALSE) %>%
    filter(species %in% unique(perm_data$species))

  # Scale if needed
  if (a$name %in% c("ATAC-seq peaks overlap")) {
    perm_data <- perm_data %>% mutate(permuted_bp_overlap = permuted_bp_overlap / 1000)
    obs_data  <- obs_data  %>% mutate(observed = observed / 1000)
  }

  # Assign groups and labels
  perm_data <- perm_data %>%
    mutate(
      group = case_when(species %in% group1 ~ "group1", TRUE ~ "group2"),
      group_label = factor(group_labels[group], levels = c("Fissions", "Fusions")),
      species = factor(species, levels = c(setdiff(unique(species), group1), group1))
    )

  obs_data <- obs_data %>%
    mutate(
      group = case_when(species %in% group1 ~ "group1", TRUE ~ "group2"),
      group_label = factor(group_labels[group], levels = c("Fissions", "Fusions")),
      species = factor(species, levels = levels(perm_data$species))
    )


  # one facet per species pair
  p_species <- ggplot() +
    geom_density(
      data = perm_data,
      aes(x = permuted_bp_overlap),
      fill = "grey70", colour = "grey50", alpha = 0.4
    ) +
    geom_vline(
      data = obs_data,
      aes(xintercept = observed, colour = group_label),
      linewidth = 1
    ) +
    scale_color_manual(
      values = c("Fissions" = "#3D7A96", "Fusions" = "#D07D59"),
      name = NULL
    ) +
    geom_richtext(
      data = obs_data,
      aes(x = Inf, y = Inf,
          label = paste0("*p*-value = ", round(pval, 3))),
      hjust = 1.1, vjust = 1.5,
      size = 3,
      fill = NA, label.color = NA
    ) +
    facet_wrap(~ species, scales = "free", labeller = as_labeller(species_labels)) +
    labs(
      x = a$xlab, y = "Density",
      title = paste0(a$name, " — Distribution across random windows")
    ) +
    theme_classic(base_size = 10) +
    theme(
      plot.title       = element_text(hjust = 0, face = "bold"),
      strip.background = element_blank(),
      strip.text       = element_markdown(face = "bold", size = 8, hjust = 0),
      legend.position  = "bottom",
      legend.text      = element_text(size = 10),
    )

  n_sp  <- n_distinct(perm_data$species)
  ncols <- 3
  nrows <- ceiling(n_sp / ncols)

  p_species <- p_species + facet_wrap(~ species, scales = "free",
                                      labeller = as_labeller(species_labels),
                                      ncol = ncols)

  ggsave(paste0(a$name, "_by_species.pdf"), p_species,
         width = 8.27, height = 11.69,
         device = cairo_pdf)

}
