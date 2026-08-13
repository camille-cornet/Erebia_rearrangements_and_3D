library(tidyverse)
library(readxl)
library(patchwork)
library(ggtext)
library(pals)

regioneR_results <- read.delim("regioneR_results_hic.tsv")
regioneR_results$TAD_boundaries_kb_observed <- regioneR_results$TAD_boundaries_observed/1000
regioneR_results$A_compartment_kb_observed <- regioneR_results$A_compartment_observed/1000
regioneR_results$B_compartment_kb_observed <- regioneR_results$B_compartment_observed/1000
regioneR_results$TAD_boundaries_kb_expected <- regioneR_results$TAD_boundaries_expected/1000
regioneR_results$A_compartment_kb_expected <- regioneR_results$A_compartment_expected/1000
regioneR_results$B_compartment_kb_expected <- regioneR_results$B_compartment_expected/1000
atac_results <- read_excel("regioneR_results_atac_seq.xlsx")
atac_results$peak_nuc_kb_observed <- atac_results$peak_nuc_observed/1000
atac_results$peak_nuc_kb_expected <- atac_results$peak_nuc_expected/1000

#### Plots ####
plot_df <- regioneR_results %>%
  pivot_longer(
    cols = ends_with(c("_observed", "_expected")),
    names_to = "name",
    values_to = "value") %>%
  mutate(
    obs_type = ifelse(grepl("_observed$", name), "observed", "expected"),
    compartment = sub("_(observed|expected)$", "", name))

plot_df <- plot_df %>%
  mutate(compartment = factor(compartment, levels = unique(compartment))) %>%
  arrange(desc(Type))

plot_atac <- atac_results %>%
  pivot_longer(
    cols = ends_with(c("_observed", "_expected")),
    names_to = "name",
    values_to = "value") %>%
  mutate(
    obs_type = ifelse(grepl("_observed$", name), "observed", "expected"),
    compartment = sub("_(observed|expected)$", "", name))

p1 <- ggplot(
  plot_df |> filter(compartment == "TAD_boundaries_kb"),
  aes(x = obs_type, y = value, group = Species_pair, color = Type)
) +
  geom_line(alpha = 1) +
  geom_point(size = 1.5) +
  scale_x_discrete(
    limits = c("expected", "observed"),
    labels = c("expected" = "Random\nWindows", "observed" = "At\nBreakpoints"),
    expand = expansion(add = 0.2)
  ) +
  scale_color_manual(values = c("Fission" = "#3D7A96", "Fusion" = "#D07D59"),
                     labels = c("Fission" = "Fissions", "Fusion" = "Fusions")) +
  labs(x = NULL, y = "Length of overlap\nwith TAD boundaries (kb)",
       color = "Type", title = "A) TAD Boundaries",
       subtitle = "*p* < 0.001") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 2), size = 8),
        axis.text.y = element_markdown(size = 6),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = -4)),
        plot.subtitle = element_markdown(size = 7, hjust = 1, margin = margin(t = 4, b = -4)),
        plot.title.position = "plot")

p2 <- ggplot(
  plot_df |> filter(compartment == "Insulation_score"),
  aes(x = obs_type, y = value, group = Species_pair, color = Type)
) +
  geom_line(alpha = 1) +
  geom_point(size = 1.5) +
  scale_x_discrete(
    limits = c("expected", "observed"),
    labels = c("expected" = "Random\nWindows", "observed" = "At\nBreakpoints"),
    expand = expansion(add = 0.2)
  ) +
  scale_color_manual(values = c("Fission" = "#3D7A96", "Fusion" = "#D07D59"),
                     labels = c("Fission" = "Fissions", "Fusion" = "Fusions")) +
  labs(x = NULL, y = "Insulation score", color = "Type",
       title = "B) Insulation score",
       subtitle = "*p* = 0.002") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 2), size = 8),
        axis.text.y = element_markdown(size = 6),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = -4)),
        plot.subtitle = element_markdown(size = 7, hjust = 1, margin = margin(t = 4, b = -4)),
        plot.title.position = "plot")

p3 <- ggplot(
  plot_df |> filter(compartment == "A_compartment_kb"),
  aes(x = obs_type, y = value, group = Species_pair, color = Type)
) +
  geom_line(alpha = 1) +
  geom_point(size = 1.5) +
  scale_x_discrete(
    limits = c("expected", "observed"),
    labels = c("expected" = "Random\nWindows", "observed" = "At\nBreakpoints"),
    expand = expansion(add = 0.2)
  ) +
  ylim(0, 2100) +
  scale_color_manual(values = c("Fission" = "#3D7A96", "Fusion" = "#D07D59"),
                     labels = c("Fission" = "Fissions", "Fusion" = "Fusions")) +
  labs(x = NULL, y = "Length of overlap\nwith A compartment (kb)",
       color = "Type", title = "C) A compartment",
       subtitle = "*p* = 0.037") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 2), size = 8),
        axis.text.y = element_markdown(size = 6),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = -4)),
        plot.subtitle = element_markdown(size = 7, hjust = 1, margin = margin(t = 4, b = -4)),
        plot.title.position = "plot")

p4 <- ggplot(
  plot_df |> filter(compartment == "B_compartment_kb"),
  aes(x = obs_type, y = value, group = Species_pair, color = Type)
) +
  geom_line(alpha = 1) +
  geom_point(size = 1.5) +
  scale_x_discrete(
    limits = c("expected", "observed"),
    labels = c("expected" = "Random\nWindows", "observed" = "At\nBreakpoints"),
    expand = expansion(add = 0.2)
  ) +
  ylim(0, 2100) +
  scale_color_manual(values = c("Fission" = "#3D7A96", "Fusion" = "#D07D59"),
                     labels = c("Fission" = "Fissions", "Fusion" = "Fusions")) +
  labs(x = NULL, y = "Length of overlap\nwith B compartment (kb)",
       color = "Type", title = "D) B compartment",
       subtitle = "*p* = 0.305") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 2), size = 8),
        axis.text.y = element_markdown(size = 6),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = -4)),
        plot.subtitle = element_markdown(size = 7, hjust = 1, margin = margin(t = 4, b = -4)),
        plot.title.position = "plot")

p6 <- ggplot(
  plot_atac |> filter(compartment == "reads"),
  aes(x = obs_type, y = value, group = species, color = Type)
) +
  geom_line(alpha = 1) +
  geom_point(size = 1.5) +
  scale_x_discrete(
    limits = c("expected", "observed"),
    labels = c("expected" = "Random\nWindows", "observed" = "At\nBreakpoints"),
    expand = expansion(add = 0.2)
  ) +
  scale_color_manual(values = c("Fission" = "#3D7A96", "Fusion" = "#D07D59"),
                     labels = c("Fission" = "Fissions", "Fusion" = "Fusions")) +
  labs(x = NULL, y = "Number of uniquely\nmappping ATAC-seq reads", color = "Type",
       title = "E) ATAC-seq reads",
       subtitle = "*p* = 0.022") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title.y = element_text(margin = margin(r = 2), size = 8),
        axis.text.y = element_markdown(size = 6),
        plot.title = element_text(size = 8, face = "bold", margin = margin(b = -4)),
        plot.subtitle = element_markdown(size = 7, hjust = 1, margin = margin(t = 4, b = -4)),
        plot.title.position = "plot")

full_p <- (p1 | p2 | p3 | p4 | p6) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom",
                                          legend.margin = margin(t = -10),
                                          legend.text = element_text(size = 8, face = "bold"),
                                          legend.title = element_blank())
full_p
ggsave("ttest_full_p.pdf", full_p, device = cairo_pdf, width = 8, height = 2.5)

#### Colour according to species pair for supp. fig. ####
species_levels <- plot_df %>% pull(Species_pair) %>% as.character() %>% unique()
atac_levels <- plot_atac %>% pull(species) %>% unique()
all_levels  <- union(species_levels, atac_levels)

n_lev <- length(all_levels)
pal   <- if (n_lev <= 25) pals::cols25(n_lev) else pals::glasbey(n_lev)
named_colors <- setNames(pal, all_levels)

sp_labels <- c(
  "Ecalcaria" = "*E. calcaria* (4 fusions)",
  "Ecassioides" = "*E. cassioides* (6 fusions)",
  "Ecassioides_Erondoui" = "*E. rondoui* (7 fissions)",
  "Eepiphron" = "*E. epiphron* (10 fusions)",
  "Egorge" = "*E. gorge* (2 fusions)",
  "Egorge_hap1" = "*E. gorge* (2 heterozygous fusions)",
  "Eligea_Eottomana" = "*E. ottomana* (20 fissions)",
  "Emedusa" = "*E. medusa* (5 fusions)",
  "Enivalis" = "*E. nivalis* (1 fusion)",
  "Eottomana_Egraucasica" = "*E. graucasica* (37 fissions)",
  "Eligea_Egraucasica" = "*E. graucasica* (40 fissions)",
  "Epharte" = "*E. pharte* (8 fusions)",
  "Epluto" = "*E. pluto* (4 fusions)",
  "Estirius" = "*E. stirius* (1 fusion)",
  "Etriaria" = "*E. triaria* (1 fusion)",
  "Estirius_hap1" = "*E. stirius* (1 heterozygous fusion)",
  "Etyndarus" = "*E. tyndarus* (2 fusions)",
  "Eeuryale" = "*E. euryale* (1 fusion)",
  "Eaethiops" = "*E. aethiops* (9 fusions)",
  "Eoeme" = "*E. oeme* (15 fusions)"
)
all_labels <- setNames(all_levels, all_levels)
all_labels[names(sp_labels)] <- sp_labels

plot_df$Species_pair   <- factor(as.character(plot_df$Species_pair), levels = all_levels)
plot_atac$Species_pair <- factor(plot_atac$species,                  levels = all_levels)

color_scale <- function(legend_levels) {
  scale_color_manual(values = named_colors,
                     limits = legend_levels,
                     labels = all_labels[legend_levels],
                     drop   = FALSE,
                     name   = "Species pair")
}

x_scale <- scale_x_discrete(
  limits = c("expected", "observed"),
  labels = c("expected" = "Random\nWindows", "observed" = "At\nBreakpoints"),
  expand = expansion(add = 0.2)
)

base_theme <- function(markdown_y = TRUE) {
  theme_minimal(base_size = 8) +
    theme(panel.grid.major.x = element_blank(),
          legend.position = "none",
          axis.title.y = element_text(margin = margin(r = 2), size = 8),
          axis.text.y  = if (markdown_y) element_markdown(size = 6) else element_text(size = 6),
          plot.title = element_text(size = 8, face = "bold", margin = margin(b = 1)),
          plot.subtitle = element_markdown(size = 7, hjust = 1, margin = margin(t = 0, b = -4)),
          plot.title.position = "plot",
          legend.text = element_markdown())
}

make_panel <- function(data, comp, ylab, title, legend_levels,
                       ylim = NULL, markdown_y = TRUE, pval = NULL) {
  p <- ggplot(filter(data, compartment == comp),
              aes(x = obs_type, y = value, group = Species_pair, color = Species_pair)) +
    geom_line(alpha = 1) +
    geom_point(size = 1.5) +
    x_scale + color_scale(legend_levels) +
    labs(x = NULL, y = ylab, title = title,
         subtitle = if (!is.null(pval)) paste0("*p* = ", pval) else NULL) +
    base_theme(markdown_y)
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}

## HiC supp figure
comp_levels <- species_levels

c1 <- make_panel(plot_df, "TAD_boundaries_kb", "Length of overlap\nwith TAD boundaries (kb)", "A) TAD Boundaries", comp_levels, pval = "0.0002")
c2 <- make_panel(plot_df, "Insulation_score",  "Insulation score",                            "B) Insulation score", comp_levels, markdown_y = FALSE, pval = "0.002")
c3 <- make_panel(plot_df, "A_compartment_kb",  "Length of overlap\nwith A compartment (kb)",  "C) A compartment", comp_levels, ylim = c(0, 2100), pval = "0.037")
c4 <- make_panel(plot_df, "B_compartment_kb",  "Length of overlap\nwith B compartment (kb)",  "D) B compartment", comp_levels, ylim = c(0, 2100), pval = "0.305")
c1 <- c1 + guides(color = guide_legend(ncol = 4))
c2 <- c2 + guides(color = "none")
c3 <- c3 + guides(color = "none")
c4 <- c4 + guides(color = "none")

fig_compartments <- (c1 | c2 | c3 | c4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.spacing.y = unit(-2, "pt"),
        legend.key.height = unit(12, "pt"))

fig_compartments
ggsave("hic_supp.pdf", fig_compartments, device = cairo_pdf, width = 8, height = 4.5)

## ATACseq supp figure
atac_legend <- atac_levels

a1 <- make_panel(plot_atac, "peak_nuc_kb", "Length of overlap\nwith ATAC-seq peaks (kb)", "A) ATAC-seq peaks", atac_legend, pval = "0.00003") +
  theme(plot.title.position = "panel")
a2 <- make_panel(plot_atac, "reads", "Number of uniquely\nmappping ATAC-seq reads", "B) ATAC-seq reads", atac_legend, pval = "0.022") +
  theme(plot.title.position = "panel")
a1 <- a1 + guides(color = guide_legend(ncol = 1))
a2 <- a2 + guides(color = "none")

fig_atac <- (a1 | a2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right",
        legend.key.spacing.y = unit(-2, "pt"),
        legend.key.height = unit(12, "pt"))

fig_atac
ggsave("atac_supp.pdf", fig_atac, device = cairo_pdf, width = 6, height = 4)


#### Stats ####
# Stouffer's Z across species, using signed regioneR z-scores.
# Combined sign = overall direction; opposite-direction species cancel.
stouffer_combined <- function(zscores) {
  z <- zscores[is.finite(zscores)]
  z_comb <- sum(z) / sqrt(length(z))
  data.frame(
    z_combined = z_comb,
    p.value    = 2 * pnorm(-abs(z_comb)),   # two-sided: enrichment OR depletion
    direction  = ifelse(z_comb > 0, "enriched", "depleted"),
    n          = length(z),
    n_enriched = sum(z > 0),
    n_depleted = sum(z < 0)
  )
}

zscore_cols <- c("A_compartment_zscore", "B_compartment_zscore",
                 "TAD_boundaries_zscore", "Insulation_score_zscore")

run_stouffer <- function(data, group_name) {
  map_dfr(zscore_cols, ~ {
    stouffer_combined(data[[.x]]) %>%
      mutate(compartment = gsub("_zscore", "", .x), group = group_name)
  })
}

stouffer_all <- bind_rows(
  run_stouffer(regioneR_results, "all"),
  run_stouffer(filter(regioneR_results, Type == "Fusion"),  "fusions"),
  run_stouffer(filter(regioneR_results, Type == "Fission"), "fissions")
) %>%
  select(group, compartment, z_combined, p.value,
         direction, n, n_enriched, n_depleted)

stouffer_all

write.table(stouffer_all, file = "stouffer_results.tsv", sep = "\t")
