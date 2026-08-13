library(tidyverse)
library(ggtext)

# Define groups (fusions and fissions)
group1 <- c("C0055_C0100", "Eligea_C0055", "X3531_Erondoui")

group_labels <- c(group1 = "Fissions", group2 = "Fusions")

species_labels <- c(
  "C0001" = "*E. pluto*",
  "C0055_C0100" = "*E. graucasica* (fissions)",
  "C0080" = "*E. calcaria*",
  "Eepiphron" = "*E. epiphron*",
  "Egorge" = "*E. gorge*",
  "Egorge_hap1" = "*E. gorge* (haplotypes)",
  "Eligea_C0055" = "*E. ottomana* (fissions)",
  "Emedusa" = "*E. medusa*",
  "Epharte" = "*E. pharte*",
  "Estirius" = "*E. stirius*",
  "Estirius_hap1" = "*E. stirius* (haplotypes)",
  "Etriaria" = "*E. triaria*",
  "X3258" = "*E. nivalis*",
  "X3531" = "*E. cassioides*",
  "X3531_Erondoui" = "*E. rondoui* (fissions)",
  "X3737" = "*E. tyndarus*"
)

# Which analyses to plot for
analyses <- list(
  list(
    name     = "Insulation score",
    pattern  = "_insu_permuted_distribution_regioneR\\.tsv$",
    exclude  = "_R\\d+_insu_",
    strip    = "_insu_permuted_distribution_regioneR.tsv",
    obs_file = "results_10kb_bp_insu_regioneR.tsv",
    xlab     = "Insulation score"
  ),
  list(
    name     = "TAD boundaries",
    pattern  = "_tads_permuted_distribution_regioneR\\.tsv$",
    exclude  = "_R\\d+_tads_",
    strip    = "_tads_permuted_distribution_regioneR.tsv",
    obs_file = "results_10kb_tads_regioneR_bpoverlap.tsv",
    xlab     = "TAD boundaries overlap (kb)"
  ),
  list(
    name     = "A compartment",
    pattern  = "_Acomp_permuted_distribution_regioneR\\.tsv$",
    exclude  = "_R\\d+_Acomp_",
    strip    = "_Acomp_permuted_distribution_regioneR.tsv",
    obs_file = "results_40kb_Acomp_regioneR_bpoverlap.tsv",
    xlab     = "A compartment overlap (kb)"
  ),
  list(
    name     = "B compartment",
    pattern  = "_Bcomp_permuted_distribution_regioneR\\.tsv$",
    exclude  = "_R\\d+_Bcomp_",
    strip    = "_Bcomp_permuted_distribution_regioneR.tsv",
    obs_file = "results_40kb_Bcomp_regioneR_bpoverlap.tsv",
    xlab     = "B compartment overlap (kb)"
  )
)

# Main loop
for (a in analyses) {

  # Load permuted distributions
  perm_files <- list.files("perm_distri", pattern = a$pattern, full.names = TRUE)
  perm_files <- perm_files[!grepl(a$exclude, perm_files)]

  perm_data <- map_dfr(perm_files, function(f) {
    sp <- sub(a$strip, "", basename(f))
    read_tsv(f, show_col_types = FALSE) %>% mutate(species = sp)
  })

  obs_data <- read_tsv(a$obs_file, show_col_types = FALSE) %>%
    filter(species %in% unique(perm_data$species))

  # Scale for bp overlap analyses
  if (a$name %in% c("TAD boundaries", "A compartment", "B compartment")) {
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

  # plot one facet per species pair
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
