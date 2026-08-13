library(ggplot2)
library(dplyr)
library(readr)
library(purrr)
library(patchwork)

# Gene density
gene <- read.csv("all_species_gene_density.csv") # average across species
gene$bin <- as.numeric(gene$bin)

p1 <- ggplot(gene, aes(x = bin, y = gene_density)) +
  geom_line(size = 0.8, color = "#fc8d62") +
  geom_vline(xintercept = 21, linetype = "dashed") +
  geom_vline(xintercept = 40, linetype = "dashed") +
  theme_classic() +
  labs(x = "", y = "Gene density", tag = "A)") +
  scale_x_continuous(breaks = c(1, 21, 40, 60),
                     labels = c("-50k", "start", "end", "+50k"))
p1
ggsave(p1, filename = "all_gene_density.pdf", width = 3.7, height = 2.7)

# Repeat density
repeats <- read.csv("all_species_repeats_density.csv") # average across species
repeats$bin <- as.numeric(repeats$bin)

p2 <- ggplot(repeats, aes(x = bin, y = repeat_density)) +
  geom_line(size = 0.8, color = "#fc8d62") +
  geom_vline(xintercept = 21, linetype = "dashed") +
  geom_vline(xintercept = 40, linetype = "dashed") +
  theme_classic() +
  labs(x = "", y = "Repeat density", tag = "B)") +
  scale_x_continuous(breaks = c(1, 21, 40, 60),
                     labels = c("-50k", "start", "end", "+50k"))
p2
ggsave(p2, filename = "all_repeats_density.pdf", width = 3.7, height = 2.7)

# Plot for supp_fig
full_p <- p1 + p2
full_p

ggsave("TAD_bound.pdf", full_p, device = cairo_pdf, width = 8, height = 3.5)
