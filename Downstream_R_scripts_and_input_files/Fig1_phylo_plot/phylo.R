library(ggplot2)
library(dplyr)
library(ggtree)
library(ape)
library(tidyr)
library(aplot)
library(ggtext)
library(ggh4x)

format_species <- function(x) {
  x <- gsub("^Eeuryaleadyte$",   "E. euryale adyte",   x)
  x <- gsub("^Eeuryaleisarica$", "E. euryale isarica", x)
  ifelse(
    grepl("^[A-Z][a-z]+$", x),
    sub("^(.)(.+)$", "\\1. \\2", x),
    x
  )
}

# Load and prepare data
tree <- read.tree("erebia_dated_phylo_clean_full.newick")
tree <- drop.tip(tree, c("Mjurtina", "Mcinxia"))
tree$tip.label <- format_species(tree$tip.label)

RM_summary <- read.table("RM_summary_nfufi.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>%
  filter(!species %in% c("Mcinxia", "Mjurtina")) %>%
  mutate(display_name = format_species(species))

# Plot the phylogeny
p_tree <- ggtree(tree, linewidth = 0.3) %>%
  flip(44, 51) %>%
  flip(11, 57) %>%
  flip(6, 50) %>%
  flip(48, 47) %>%
  flip(1, 2) %>%
  flip(3, 4) %>%
  flip(65, 66) %>%
  flip(31, 32) %>%
  flip(7, 8) %>%
  flip(29, 30)

p_tree$data <- p_tree$data %>%
  left_join(RM_summary %>% select(display_name, in_bp),
            by = c("label" = "display_name"))

p_tree <- p_tree +
  geom_tiplab(aes(fontface = ifelse(in_bp == "yes", "bold.italic", "italic")),
              align = TRUE, size = 3, linetype = 0) +
  geom_treescale(x = 0, y = 0.5, width = 1, label = "1 Myr", fontsize = 3) +
  xlim(0, 12)

# Reshape data for bar plots
data_long_plot <- RM_summary %>%
  select(display_name, karyotype, nb_fusions, nb_fissions) %>%
  pivot_longer(-display_name, names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = c("karyotype", "nb_fusions", "nb_fissions")))

# Bar plots
p_bars <- ggplot(data_long_plot, aes(x = value, y = display_name, color = variable)) +
  geom_segment(aes(x = 0, xend = value, yend = display_name), linewidth = 0.4) +
  geom_point(size = 2.5) +
  facet_wrap(~variable, scales = "free_x", nrow = 1, strip.position = "bottom",
             labeller = as_labeller(c(
               "karyotype"   = "Chromosome<br>number (*n*)",
               "nb_fusions"  = "Number of fusions",
               "nb_fissions" = "Number of fissions"
             ))) +
  facetted_pos_scales(x = list(
    variable == "karyotype"   ~ scale_x_continuous(breaks = scales::pretty_breaks(n = 3),
                                                   limits = function(x) c(0, max(x, na.rm = TRUE))),
    variable == "nb_fusions"  ~ scale_x_continuous(breaks = c(0, 20, 40), limits = c(0, 42)),
    variable == "nb_fissions" ~ scale_x_continuous(breaks = c(0, 20, 40), limits = c(0, 42))
  )) +
  scale_color_manual(values = c("karyotype"   = "grey60",
                                "nb_fusions"  = "#D07D59",
                                "nb_fissions" = "#3D7A96")) +
  theme_minimal() +
  theme(
    axis.text.y        = element_blank(),
    axis.title.y       = element_blank(),
    axis.title.x       = element_blank(),
    axis.ticks.y       = element_blank(),
    legend.position    = "none",
    strip.text         = element_markdown(face = "bold"),
    strip.placement    = "outside",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# A/B compartment stacked bar plot
ab_long <- RM_summary %>%
  select(display_name, A_perc, B_perc) %>%
  pivot_longer(-display_name, names_to = "compartment", values_to = "value") %>%
  mutate(
    compartment = factor(compartment, levels = c("A_perc", "B_perc")),
    panel_label = "A & B compartments <br><span style='color:#FFC107'>&#9632;</span> A% &nbsp;<span style='color:#D81B60'>&#9632;</span> B%"
  )

p_ab <- ggplot(ab_long, aes(x = value, y = display_name, fill = compartment)) +
  geom_bar(stat = "identity", position = position_stack(reverse = TRUE)) +
  facet_wrap(~panel_label, nrow = 1, strip.position = "bottom") +
  scale_x_continuous(
    limits = c(0, 100),
    breaks = c(0, 50, 100),
    labels = c("0", "50", "100")
  ) +
  scale_fill_manual(
    values = c("A_perc" = "#FFC107", "B_perc" = "#D81B60")
  ) +
  theme_minimal() +
  theme(
    axis.text.y        = element_blank(),
    axis.title.y       = element_blank(),
    axis.title.x       = element_blank(),
    axis.ticks.y       = element_blank(),
    legend.position    = "none",
    strip.text         = element_markdown(face = "bold", size = 8),
    strip.placement    = "outside",
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.minor.x = element_blank()
  )

# Align tree + bars + A/B panel
phylo_p <- p_ab |>
  insert_left(p_bars, width = 4) |>
  insert_left(p_tree, width = 2.5)
phylo_p

ggsave("phylo_p.pdf", phylo_p, device = cairo_pdf, width = 9, height = 6)

