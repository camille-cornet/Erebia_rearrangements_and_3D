library(tidyverse)
library(ggplot2)
library(ggtext)
library(ape)

set.seed(123)

# Load data
inter_intra <- read.delim("inter_intra_df.tsv")
slopes      <- read.delim("slopes_perchrom.tsv")
comp_strength <- read.delim("compartment_strength_per_chrom.tsv")
tree        <- read.tree("erebia_dated_phylo_clean_full.newick")
tree        <- drop.tip(tree, c("Mjurtina", "Mcinxia"))


# Species name mapping
all_names <- c(
  "C0001" = "*E. pluto*",
  "C0055" = "*E. ottomana*",
  "C0080" = "*E. calcaria*",
  "C0100" = "*E. graucasica*",
  "X2575" = "*E. euryale adyte*",
  "X2576" = "*E. euryale isarica*",
  "X3252" = "*E. pandrose*",
  "X3258" = "*E. nivalis*",
  "X3311" = "*E. oeme*",
  "X3506" = "*E. styx*",
  "X3531" = "*E. cassioides*",
  "X3737" = "*E. tyndarus*",
  "Eaethiops"      = "*E. aethiops*",
  "Ealbergana"     = "*E. albergana*",
  "Ebubastis"      = "*E. bubastis*",
  "Echristi"       = "*E. christi*",
  "Edisa"          = "*E. disa*",
  "Eembla"         = "*E. embla*",
  "Eepiphron"      = "*E. epiphron*",
  "Eeriphyle"      = "*E. eriphyle*",
  "Eflavofasciata" = "*E. flavofasciata*",
  "Egorge"         = "*E. gorge*",
  "Eligea"         = "*E. ligea*",
  "Emanto"         = "*E. manto*",
  "Emedusa"        = "*E. medusa*",
  "Emelampus"      = "*E. melampus*",
  "Emelancholica"  = "*E. melancholica*",
  "Emeolans"       = "*E. meolans*",
  "Emnestra"       = "*E. mnestra*",
  "Emontana"       = "*E. montana*",
  "Epalarica"      = "*E. palarica*",
  "Epharte"        = "*E. pharte*",
  "Epronoe"        = "*E. pronoe*",
  "Erondoui"       = "*E. rondoui*",
  "Estirius"       = "*E. stirius*",
  "Esudetica"      = "*E. sudetica*",
  "Etriaria"       = "*E. triaria*"
)

# Phylogeny tip order
is_tip          <- tree$edge[, 2] <= Ntip(tree)
tip_order_codes <- tree$tip.label[tree$edge[is_tip, 2]]
label_order     <- ifelse(tip_order_codes %in% names(all_names),
                          all_names[tip_order_codes],
                          tip_order_codes)

# Merge and reshape
combined <- inter_intra %>%
  select(species, chrom, ratio_inter_intra) %>%
  full_join(slopes %>% select(species, chrom, mean_slope_0_10Mb),
            by = c("species", "chrom")) %>%
  full_join(comp_strength %>% select(species, chrom, strength),
            by = c("species", "chrom")) %>%
  filter(!species %in% c("Mcinxia", "Mjurtina")) %>%  # remove outgroups
  pivot_longer(cols = c(ratio_inter_intra, mean_slope_0_10Mb, strength),
               names_to  = "metric",
               values_to = "value") %>%
  mutate(
    species_label = ifelse(species %in% names(all_names),
                           all_names[species], species),
    species_label = factor(species_label, levels = label_order),
    metric = factor(metric,
                    levels = c("ratio_inter_intra", "strength", "mean_slope_0_10Mb"),
                    labels = c("Inter/Intra ratio", "strength", "Mean slope (0-10Mb)"))
  )

# Plot
metric_labels <- c(
  "Mean slope (0-10Mb)" = "P(s) slope (0 - 10 Mb range)",
  "Inter/Intra ratio" = "Inter-/intrachromosomal\ncontact ratio",
  "strength" = "Compartment strength"
)

# Reorder so slope is on top
combined <- combined %>%
  mutate(metric = factor(metric,
                         levels = c("Mean slope (0-10Mb)",
                                    "strength",
                                    "Inter/Intra ratio"))) %>%
  filter(!(species == "Egorge" & chrom == "scaffold_19"))

ggplot(combined, aes(x = species_label, y = value)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  facet_wrap(~ metric, ncol = 1, scales = "free_y",
             strip.position = "left",
             labeller = as_labeller(metric_labels)) +
  theme_bw() +
  theme(
    axis.text.x      = element_markdown(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.placement  = "outside",
    strip.text.y.left = element_text(angle = 90),
    panel.spacing    = unit(0.5, "lines"),
    axis.title.y     = element_blank()
  ) +
  labs(x = "Species")

ggsave("boxplot_persp.pdf",
       width  = 297,
       height = 210,
       units  = "mm",
       device = "pdf")

