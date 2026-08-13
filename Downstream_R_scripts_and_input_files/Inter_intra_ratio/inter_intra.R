library(ggplot2)
library(dplyr)
library(ggtext)
library(phylolm)
library(phytools)
library(patchwork)
library(ape)
library(MCMCglmm)
library(mulTree)

set.seed(123)

# Load all species files
df <- read.table("inter_intra_df.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Colour grouping
highlight <- c("C0055","C0080","C0100","Erondoui","X3258","X3531","X3737")

colors <- c(
  "Other"    = "gray60",
  "C0055"    = "chocolate4",
  "C0080"    = "skyblue",
  "C0100"    = "#F0E442",
  "Erondoui" = "#E69F00",
  "X3258"    = "#009E73",
  "X3531"    = "#0072B2",
  "X3737"    = "royalblue4"
)

df <- df |>
  mutate(
    col_group = ifelse(species %in% highlight, species, "Other"),
    col_group = factor(col_group, levels = c(highlight, "Other"))
  ) |>
  # Plot "Other" first so highlighted species appear on top
  arrange(desc(col_group == "Other"))

ordered_highlight <- df |>
  filter(col_group != "Other") |>
  group_by(col_group) |>
  summarise(med = median(ratio_inter_intra, na.rm = TRUE)) |>
  arrange(desc(med)) |>
  pull(col_group) |>
  as.character()

df$col_group <- factor(df$col_group, levels = c(ordered_highlight, "Other"))

# Plot
p <- ggplot(df, aes(x = size / 1e6, y = ratio_inter_intra,
                    colour = col_group, alpha = col_group)) +
  geom_point(size = 3, stroke = 0, shape = 16) +
  scale_colour_manual(
    values = colors, labels = c(
      "C0055"    = "*E. ottomana* (*n* = 40)",
      "C0080"    = "*E. calcaria* (*n* = 8)",
      "C0100"    = "*E. graucasica* (*n* = 51)",
      "Erondoui" = "*E. rondoui* (*n* = 24)",
      "X3258"    = "*E. nivalis* (*n* = 11)",
      "X3531"    = "*E. cassioides* (*n* = 10)",
      "X3737"    = "*E. tyndarus* (*n* = 10)",
      "Other"    = "Other species"),
    name   = "Species",
    breaks = c(ordered_highlight, "Other")   # legend order: highlighted first
  ) +
  scale_alpha_manual(
    values = c(setNames(rep(0.9, length(ordered_highlight)), ordered_highlight), "Other" = 0.3),
    guide  = "none"
  ) +
  scale_x_log10(labels = scales::comma_format(suffix = " Mb"),
                breaks = c(10, 50, 100)) +
  annotation_logticks(sides = "b") +
  labs(
    x     = "Chromosome size (Mb)",
    y     = "Inter-/intrachromosomal contact ratio"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.key.size  = unit(0.4, "cm"),
    legend.text      = element_markdown(size = 10),
    axis.text        = element_text(size = 10),
    legend.position      = c(0.95, 1),
    legend.key.height = unit(0.6, "cm"),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = "grey90", color = NA)
  ) +
  theme(legend.key = element_rect(fill = "grey90", color = NA)) +
  guides(color = guide_legend(title = "Species", ncol = 1,
                              override.aes = list(size = 3, alpha = 1,
                                                  fill = "grey90")))
p

y_limits <- range(df$ratio_inter_intra, na.rm = TRUE)

p_box <- ggplot(df, aes(x = col_group, y = ratio_inter_intra,
                        colour = col_group, fill = col_group)) +
  geom_boxplot(alpha = 0.6, outlier.size = 1) +
  scale_colour_manual(values = colors, guide = "none") +
  scale_fill_manual(values = colors, guide = "none") +
  scale_x_discrete(
    limits = c(ordered_highlight, "Other"),
    labels = c(
    "C0055"    = "*E. ottomana*",
    "C0080"    = "*E. calcaria*",
    "C0100"    = "*E. graucasica*",
    "Erondoui" = "*E. rondoui*",
    "X3258"    = "*E. nivalis*",
    "X3531"    = "*E. cassioides*",
    "X3737"    = "*E. tyndarus*",
    "Other"    = "Other species"
    )) +
  coord_cartesian(ylim = y_limits) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x  = element_markdown(angle = 45, hjust = 1),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank()
  )
p_box

p_full <- (p | p_box) +
  plot_layout(widths = c(3, 1))
p_full

# Save
ggsave("inter_intra_size.pdf", p_full, device = cairo_pdf, width = 10, height = 6)


#### Test ####
RM_summary <- read.table("RM_summary_nfufi.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
RM_summary$median_autosome_size_Mb <- RM_summary$median_autosome_size/1000000

df_sp <- df |>
  group_by(species) |>
  summarise(
    mean_ratio = mean(ratio_inter_intra)
  ) |>
  left_join(RM_summary, by = "species")

tree <- read.tree("erebia_dated_phylo_clean_full.newick")
tree <- drop.tip(tree, c("Mjurtina", "Mcinxia"))
rownames(df_sp) <- df_sp$species
setdiff(df_sp$species, tree$tip.label)  # in data but not tree
setdiff(tree$tip.label, df_sp$species)  # in tree but not data

#### Simple pgls for size
model_ratio_size <- phylolm(
  mean_ratio ~  median_autosome_size_Mb,
  data = df_sp,
  phy = tree,
  model = "lambda"
)
summary(model_ratio_size)


#### MCMCglmm ####
tree <- force.ultrametric(tree)
invA <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv
my_priors_a <- list(
  R = list(V = 1, fix = TRUE),
  G = list(G1 = list(V = 1, nu = 0.002))
)
df_sp$animal <- as.factor(df_sp$species)
df_sp$animal <- factor(as.character(df_sp$animal))

# Karyotype
cl <- makeCluster(4)
clusterExport(cl, list("df_sp", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(karyotype ~ mean_ratio + median_autosome_size_Mb,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = df_sp,
           family   = "poisson",
           prior    = my_priors_a,
           nitt = 200000, burnin = 50000, thin = 100)
})
stopCluster(cl)

# Check convergence
fit_sol <- lapply(fit, function(m) m$Sol)
fit_sol <- do.call(mcmc.list, fit_sol)
print(gelman.diag(fit_sol))

# Save MCMC chains
save(fit_sol, file = "fit_n.RData")

# Summarise pooled chains
combined_Sol <- do.call(rbind, fit_sol)
pMCMC    <- apply(combined_Sol, 2, function(x) 2 * min(mean(x > 0), mean(x < 0)))
post_mean <- apply(combined_Sol, 2, mean)
ci        <- apply(combined_Sol, 2, quantile, probs = c(0.025, 0.975))
eff_samp  <- effectiveSize(fit_sol)

results <- data.frame(
  param    = names(post_mean),
  mean     = post_mean,
  lower    = ci[1, ],
  upper    = ci[2, ],
  eff_samp = eff_samp,
  pMCMC    = pMCMC
)

print(results)
write.table(results, file = paste0("results_n.tsv"), sep = "\t", row.names = FALSE)

# Fusions
cl <- makeCluster(4)
clusterExport(cl, list("df_sp", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(nb_fusions ~ mean_ratio + median_autosome_size_Mb,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = df_sp,
           family   = "poisson",
           prior    = my_priors_a,
           nitt = 200000, burnin = 50000, thin = 100)
})
stopCluster(cl)

# Check convergence
fit_sol <- lapply(fit, function(m) m$Sol)
fit_sol <- do.call(mcmc.list, fit_sol)
print(gelman.diag(fit_sol))

# Save MCMC chains
save(fit_sol, file = "fit_fu.RData")

# Summarise pooled chains
combined_Sol <- do.call(rbind, fit_sol)
pMCMC    <- apply(combined_Sol, 2, function(x) 2 * min(mean(x > 0), mean(x < 0)))
post_mean <- apply(combined_Sol, 2, mean)
ci        <- apply(combined_Sol, 2, quantile, probs = c(0.025, 0.975))
eff_samp  <- effectiveSize(fit_sol)

results <- data.frame(
  param    = names(post_mean),
  mean     = post_mean,
  lower    = ci[1, ],
  upper    = ci[2, ],
  eff_samp = eff_samp,
  pMCMC    = pMCMC
)

print(results)
write.table(results, file = paste0("results_fu.tsv"), sep = "\t", row.names = FALSE)

# Fissions
cl <- makeCluster(4)
clusterExport(cl, list("df_sp", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(nb_fissions ~ mean_ratio + median_autosome_size_Mb,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = df_sp,
           family   = "poisson",
           prior    = my_priors_a,
           nitt = 200000, burnin = 50000, thin = 100)
})
stopCluster(cl)

# Check convergence
fit_sol <- lapply(fit, function(m) m$Sol)
fit_sol <- do.call(mcmc.list, fit_sol)
print(gelman.diag(fit_sol))

# Save MCMC chains
save(fit_sol, file = "fit_fi.RData")

# Summarise pooled chains
combined_Sol <- do.call(rbind, fit_sol)
pMCMC    <- apply(combined_Sol, 2, function(x) 2 * min(mean(x > 0), mean(x < 0)))
post_mean <- apply(combined_Sol, 2, mean)
ci        <- apply(combined_Sol, 2, quantile, probs = c(0.025, 0.975))
eff_samp  <- effectiveSize(fit_sol)

results <- data.frame(
  param    = names(post_mean),
  mean     = post_mean,
  lower    = ci[1, ],
  upper    = ci[2, ],
  eff_samp = eff_samp,
  pMCMC    = pMCMC
)

print(results)
write.table(results, file = paste0("results_fi.tsv"), sep = "\t", row.names = FALSE)


# Does the tyndarus clade have higher inter/intra ratio?
tyndarus <- tree$tip.label[phytools::getDescendants(tree, 62)[
  phytools::getDescendants(tree, 62) <= length(tree$tip.label)]]
medusa <- tree$tip.label[phytools::getDescendants(tree, 45)[
  phytools::getDescendants(tree, 45) <= length(tree$tip.label)]]
pluto <- tree$tip.label[phytools::getDescendants(tree, 49)[
  phytools::getDescendants(tree, 49) <= length(tree$tip.label)]]
pronoe <- tree$tip.label[phytools::getDescendants(tree, 53)[
  phytools::getDescendants(tree, 53) <= length(tree$tip.label)]]
epiphron <- tree$tip.label[phytools::getDescendants(tree, 58)[
  phytools::getDescendants(tree, 58) <= length(tree$tip.label)]]
ligea <- tree$tip.label[phytools::getDescendants(tree, 68)[
  phytools::getDescendants(tree, 68) <= length(tree$tip.label)]]
embla <- tree$tip.label[phytools::getDescendants(tree, 73)[
  phytools::getDescendants(tree, 73) <= length(tree$tip.label)]]

all_tips <- tree$tip.label
group <- rep("other", length(all_tips))
names(group) <- all_tips
group[all_tips %in% tyndarus] <- "tyndarus"
group[all_tips %in% medusa]   <- "medusa"
group[all_tips %in% pluto]    <- "pluto"
group[all_tips %in% pronoe]   <- "pronoe"
group[all_tips %in% epiphron] <- "epiphron"
group[all_tips %in% ligea]    <- "ligea"
group[all_tips %in% embla]    <- "embla"

mean_ratio <- df_sp$mean_ratio
names(mean_ratio) <- df_sp$species
all(names(mean_ratio) %in% tree$tip.label) # to check if all names correspond
ratio_test <- phylANOVA(tree, group, mean_ratio, nsim = 10000, p.adj = "fdr")
ratio_test

# For testing just tyndarus vs others, do:
group <- ifelse(all_tips %in% tyndarus, "focal", "rest")
names(group) <- tree$tip.label
all(names(mean_ratio) %in% tree$tip.label) # to check if all names correspond
ratio_test <- phylANOVA(tree, group, mean_ratio, nsim = 10000, p.adj = "fdr")
ratio_test


#### Plots ####
highlight <- c("C0055","C0080","C0100","Erondoui","X3258","X3531","X3737")

colors <- c(
  "Other"    = "gray60",
  "C0055"    = "chocolate4",
  "C0080"    = "skyblue",
  "C0100"    = "#F0E442",
  "Erondoui" = "#E69F00",
  "X3258"    = "#009E73",
  "X3531"    = "#0072B2",
  "X3737"    = "royalblue4"
)

df_sp <- df_sp |>
  mutate(
    col_group = ifelse(species %in% highlight, as.character(species), "Other"),
    col_group = factor(col_group, levels = c(highlight, "Other"))
  ) |>
  arrange(desc(col_group == "Other"))

p_n <- ggplot(df_sp, aes(x = mean_ratio, y = karyotype, color = col_group)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors, name = "Species",
                     labels = c(
                       "C0055"    = "*E. ottomana* (*n* = 40)",
                       "C0080"    = "*E. calcaria* (*n* = 8)",
                       "C0100"    = "*E. graucasica* (*n* = 51)",
                       "Erondoui" = "*E. rondoui* (*n* = 24)",
                       "X3258"    = "*E. nivalis* (*n* = 11)",
                       "X3531"    = "*E. cassioides* (*n* = 10)",
                       "X3737"    = "*E. tyndarus* (*n* = 10)",
                       "Other"    = "Other species"
                     )) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) == 0.867",
           parse = TRUE, size = 3) +
  theme_classic(base_size = 10) +
  theme(legend.text = element_markdown(),
        axis.title.y = element_markdown()) +
  labs(x = "Inter-/intrachromosomal\ncontact ratio", y = "Chromosome number (*n*)")
p_n

p_fu <- ggplot(df_sp, aes(x = mean_ratio, y = nb_fusions, color = col_group)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors, name = "Species",
                     labels = c(
                       "C0055"    = "*E. ottomana* (*n* = 40)",
                       "C0080"    = "*E. calcaria* (*n* = 8)",
                       "C0100"    = "*E. graucasica* (*n* = 51)",
                       "Erondoui" = "*E. rondoui* (*n* = 24)",
                       "X3258"    = "*E. nivalis* (*n* = 11)",
                       "X3531"    = "*E. cassioides* (*n* = 10)",
                       "X3737"    = "*E. tyndarus* (*n* = 10)",
                       "Other"    = "Other species"
                     )) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) == 0.782",
           parse = TRUE, size = 3) +
  theme_classic(base_size = 10) +
  theme(legend.text = element_markdown()) +
  labs(x = "Inter-/intrachromosomal\ncontact ratio", y = "Number of fusions")
p_fu

p_fi <- ggplot(df_sp, aes(x = mean_ratio, y = nb_fissions, color = col_group)) +
  geom_point(size = 2) +
  scale_color_manual(values = colors, name = "Species",
                     labels = c(
                       "C0055"    = "*E. ottomana* (*n* = 40)",
                       "C0080"    = "*E. calcaria* (*n* = 8)",
                       "C0100"    = "*E. graucasica* (*n* = 51)",
                       "Erondoui" = "*E. rondoui* (*n* = 24)",
                       "X3258"    = "*E. nivalis* (*n* = 11)",
                       "X3531"    = "*E. cassioides* (*n* = 10)",
                       "X3737"    = "*E. tyndarus* (*n* = 10)",
                       "Other"    = "Other species"
                     )) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) == '0.633'",
           parse = TRUE, size = 3) +
  theme_classic(base_size = 10) +
  theme(legend.text = element_markdown(),
        axis.text.x = element_markdown()) +
  labs(x = "Inter-/intrachromosomal\ncontact ratio", y = "Number of fissions")
p_fi

nfufi <- (p_n + p_fu + p_fi) +
  plot_layout(guides = "collect") &
  theme(legend.position = "none")
nfufi

ggsave("full_plot.pdf", nfufi, device = cairo_pdf, width = 8, height = 2.5)

