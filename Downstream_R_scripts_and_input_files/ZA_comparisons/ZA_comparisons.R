library(tidyverse)
library(phytools)
library(phylolm)
library(ggtext)

AB <- read.delim("AB_per_chrom.tsv")
TAD <- read.delim("average_TAD_size_per_chrom.tsv")
slope <- read.delim("avg_slope_loglog_perchrom_allchroms_allspecies.tsv")
ratio <- read.delim("inter_intra_allchroms_allspecies.tsv")
strength <- read.delim("compartment_strength_per_chrom.tsv")

merged <- merge(AB, TAD, by = c("species","chrom"))
merged <- merge(merged, slope, by = c("species","chrom"))
merged <- merge(merged, ratio, by = c("species","chrom"))
merged <- merge(merged, strength, by = c("species","chrom"))

sexchromlist <- read.delim("sexchromlist.txt")

sex_long <- sexchromlist %>%
  pivot_longer(cols = c(Z1, W, Z2, W2), names_to = "type", values_to = "chrom") %>%
  filter(chrom != "NONE") %>%
  mutate(chrom_type = case_when(
    grepl("^Z", type) ~ "Z",
    grepl("^W", type) ~ "W"
  )) %>%
  select(species, chrom, chrom_type)

merged <- merged %>%
  left_join(sex_long, by = c("species", "chrom")) %>%
  mutate(chrom_type = replace_na(chrom_type, "autosome")) %>%
  filter(!chrom_type %in% c("W")) %>% filter(!species %in% c("Mcinxia", "Mjurtina"))

tree <- read.tree("erebia_dated_phylo_clean_full.newick")
tree <- drop.tip(tree, c("Mjurtina", "Mcinxia"))

# paired Z - autosome difference
run_paired_phylo_cov <- function(df, tree, metric, size_col = "size") {
  df_paired <- df %>%
    group_by(species) %>%
    summarise(
      mean_auto   = mean(.data[[metric]][chrom_type == "autosome"], na.rm = TRUE),
      Z_val       = mean(.data[[metric]][chrom_type == "Z"], na.rm = TRUE),
      auto_size   = mean(.data[[size_col]][chrom_type == "autosome"], na.rm = TRUE),
      Z_size      = mean(.data[[size_col]][chrom_type == "Z"], na.rm = TRUE),
      diff_Z_auto = Z_val  - mean_auto,
      diff_size   = Z_size - auto_size, # to correct for chrom size
      .groups = "drop"
    ) %>%
    filter(is.finite(diff_Z_auto), is.finite(diff_size)) %>%
    mutate(diff_size_c = diff_size - mean(diff_size))  # interpretable intercept

  df_paired <- as.data.frame(df_paired)
  rownames(df_paired) <- df_paired$species
  keep <- intersect(tree$tip.label, df_paired$species)
  this_tree <- drop.tip(tree, setdiff(tree$tip.label, keep))
  df_paired <- df_paired[this_tree$tip.label, ]

  phylolm(diff_Z_auto ~ diff_size_c,
          data = df_paired,
          phy = this_tree,
          model = "lambda")
}

metrics <- c("A_perc", "B_perc", "TAD_size_kb",
             "avg_slope_loglog", "strength", "ratio_inter_intra")

models_cov <- lapply(metrics, function(m) run_paired_phylo_cov(merged, tree, m))
names(models_cov) <- metrics

for (m in metrics) {
  cat("\n==================", m, "==================\n")
  print(summary(models_cov[[m]]))
}

# Result table: intercept = size-adjusted mean Z-autosome diff, slope = size effect
results_cov <- map_dfr(metrics, function(m) {
  co <- summary(models_cov[[m]])$coefficients
  tibble(
    metric        = m,
    n_species     = models_cov[[m]]$n,
    adj_diff      = co["(Intercept)", "Estimate"], # Z-auto diff, size held at mean
    pval          = co["(Intercept)", "p.value"],
    size_slope    = co["diff_size_c", "Estimate"],
    size_slope_p  = co["diff_size_c", "p.value"],
    lambda        = models_cov[[m]]$optpar
  )
})

results_cov

#### Plot ####
plot_df <- merged %>%
  filter(chrom_type %in% c("autosome", "Z")) %>%
  select(species, chrom_type,
         A_perc, B_perc, TAD_size_kb,
         avg_slope_loglog, strength, ratio_inter_intra) %>%
  pivot_longer(
    cols = c(A_perc, B_perc, TAD_size_kb,
             avg_slope_loglog, strength, ratio_inter_intra),
    names_to = "metric", values_to = "value"
  ) %>%
  # average per species: one autosome mean and one Z mean per species
  group_by(species, chrom_type, metric) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    chrom_type = factor(chrom_type, levels = c("Z", "autosome")),
    metric     = factor(metric, levels = c("A_perc", "B_perc", "TAD_size_kb",
                                           "avg_slope_loglog", "strength",
                                           "ratio_inter_intra"))
  )

metric_labels <- c(
  "A_perc" = "%A compartment",
  "B_perc" = "%B compartment",
  "avg_slope_loglog" = "P(s) slope (0 - 10 Mb range)",
  "TAD_size_kb" = "TAD size (kb)",
  "strength" = "Compartment strength",
  "ratio_inter_intra" = "Inter-/intra chromosomal contact ratio"
)

highlight <- c("C0055","C0080","C0100","Erondoui","X3258","X3531","X3737")
plot_df$col_group <- ifelse(plot_df$species %in% highlight, plot_df$species, "Other")
plot_df$col_group <- factor(plot_df$col_group, levels = c(highlight, "Other"))

# Plot "Other" first so highlighted species draw on top
plot_df <- plot_df %>% arrange(desc(col_group == "Other"))

colors <- c("Other" = "gray60",
            "C0055"="chocolate4",
            "C0080"="skyblue",
            "C0100"="#F0E442",
            "Erondoui"="#E69F00",
            "X3258"="#009E73",
            "X3531"="#0072B2",
            "X3737"="royalblue4")

labs <- c("Z", "Autosomes")

ggplot(plot_df, aes(x = chrom_type, y = value)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.52, fill = "grey90") +
  geom_point(aes(color = col_group),
             position = position_jitter(width = 0.15, seed = 1),
             size = 1.8) +
  facet_wrap(~ metric, scales = "free_y", nrow = 2,
             strip.position = "left",
             labeller = as_labeller(metric_labels)) +
  scale_color_manual(values = colors,
                     labels = c(
                       "C0055"    = "*E. ottomana*",
                       "C0080"    = "*E. calcaria*",
                       "C0100"    = "*E. graucasica*",
                       "Erondoui" = "*E. rondoui*",
                       "X3258"    = "*E. nivalis*",
                       "X3531"    = "*E. cassioides*",
                       "X3737"    = "*E. tyndarus*",
                       "Other"    = "Other species"
                     ),
                     name = "Species") +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(labels = labs) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    legend.text = element_markdown(),
    strip.background = element_blank(),
    strip.placement  = "outside"
  )

ggsave("boxplot_ZA.pdf",
       width  = 297/1.5,
       height = 210/1.5,
       units  = "mm",
       device = "pdf")
