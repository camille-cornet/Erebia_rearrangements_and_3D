library(tidyverse)
library(multcomp)
library(geiger)
library(phytools)
library(phylolm)
library(MCMCglmm)
library(parallel)
library(ggtext)
library(readr)
library(dplyr)
library(purrr)
library(stringr)

set.seed(123)

#### All TEs, normalization, rDNA clusters included ####
df <- read.table("allsp_interchrom_avg_norm.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
remove_sp <- c("Mjurtina", "Mcinxia")

df_long <- df %>%
  pivot_longer(
    cols = -TEtype,
    names_to = "Species",
    values_to = "Value"
  )

# remove repeats that have interchromosomal contacts in less than 50% of species
df_long_filtered <- df_long %>%
  filter(!Species %in% remove_sp) %>%
  group_by(TEtype) %>%
  filter(mean(is.na(Value)) < 0.5) %>%
  ungroup()

df_long_filtered <- df_long_filtered %>%
  group_by(TEtype) %>%
  mutate(median_value = median(Value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(TEtype = fct_reorder(TEtype, median_value))

# define which species to color
highlight <- c("C0055","C0080","C0100","Erondoui","X3258","X3531","X3737")
df_long_filtered$col_group <- ifelse(df_long_filtered$Species %in% highlight, df_long_filtered$Species, "Other")
df_long_filtered$col_group <- factor(df_long_filtered$col_group,
                                     levels = c(highlight, "Other"))
# Colour palette
colors <- c("Other" = "gray60",
            "C0055" = "chocolate4",
            "C0080" = "skyblue",
            "C0100" = "#F0E442",
            "Erondoui" = "#E69F00",
            "X3258" = "#009E73",
            "X3531" = "#0072B2",
            "X3737" = "royalblue4")

allTEs <- ggplot(df_long_filtered, aes(x = TEtype, y = Value, color = col_group)) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.7, shape = 16) +
  geom_boxplot(outlier.shape = NA, color = "black",
               alpha = 0, width = 0.4, linewidth = 0.3) +
  scale_y_log10() +
  scale_color_manual(values = colors,
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
  scale_x_discrete(labels = function(x) gsub("\\.", "/", gsub("_", " ", x))) +
  coord_flip() +
  labs(x = "Repeat family", y = "Average interchromosomal contacts") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        axis.title.y = element_markdown(),
        axis.title.x = element_text(size = 10),
        legend.text = element_markdown(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "grey90", color = NA),
        legend.key.spacing.y = unit(-2, "pt"),
        legend.key.height = unit(12, "pt"),
        legend.box.margin = margin(l = -90),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  theme(legend.key = element_rect(fill = "grey90", color = NA)) +
  guides(color = guide_legend(title = "Species", ncol = 2,
                              override.aes = list(size = 2, alpha = 1)))
allTEs
ggsave("TE_full_plot.pdf", allTEs, device = cairo_pdf, width = 4, height = 8)

# Stats with MCMCglmm
tree <- read.tree("erebia_dated_phylo_clean_full.newick")
tree <- drop.tip(tree, c("Mjurtina", "Mcinxia"))
tree <- force.ultrametric(tree)
invA <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv
df_mcmc <- as.data.frame(df_long_filtered) %>%
  mutate(Species = factor(as.character(Species),
                          levels = rownames(invA)))
# Also include genomic proportion as a predictor of contact
te_keep <- unique(df_long_filtered$TEtype)
RM_summary_full <- read.delim("RM_summary_full.tsv", check.names = FALSE)
df_prop_long <- RM_summary_full %>%
  dplyr::rename(Species = species) %>%
  dplyr::select(Species, all_of(te_keep)) %>%
  pivot_longer(-Species, names_to = "TEtype", values_to = "genomic_prop")
# Also include average TE length per species as covariate
repeat_length_matrix <- read.delim("repeat_length_matrix.tsv", check.names = FALSE)
repeat_length_matrix <- repeat_length_matrix |> column_to_rownames("repeat_type")
repeat_length_matrix <- as.data.frame(t(repeat_length_matrix)) |>
  rownames_to_column("species")
df_length_long <- repeat_length_matrix %>%
  dplyr::rename(Species = species) %>%
  dplyr::select(Species, all_of(te_keep)) %>%
  pivot_longer(-Species, names_to = "TEtype", values_to = "TE_length")
# Join contacts and proportions and length
df_mcmc <- as.data.frame(df_long_filtered) %>%
  left_join(df_prop_long, by = c("Species", "TEtype")) %>%
  left_join(df_length_long, by = c("Species", "TEtype")) %>%
  mutate(Species = factor(as.character(Species), levels = rownames(invA)))
# DNA.TcMar-Tigger has the lowest average inter contacts: used as ref level
df_mcmc <- df_mcmc %>%
  mutate(TEtype = relevel(factor(TEtype), ref = "DNA.TcMar-Tigger"))

prior <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002))

# Run 4 chains in parallel to check convergence
cl <- makeCluster(4)
clusterExport(cl, list("df_mcmc", "prior", "invA"))
clusterEvalQ(cl, library(MCMCglmm))
model <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(Value ~ TEtype + genomic_prop + TE_length,
           random = ~ Species,
           ginverse = list(Species = invA),
           data = df_mcmc,
           prior = prior,
           nitt = 200000, burnin = 50000, thin = 100)
})
stopCluster(cl)

model <- lapply(model, function(m) m$Sol)
model <- do.call(mcmc.list, model)
gelman.diag(model)
summary(model)
save(model, file = "model.RData")
load("model.RData")
summary(model)
# Get pMCMC values (and CI) from pooled chains
combined_Sol <- do.call(rbind, model)
pMCMC <- apply(combined_Sol, 2, function(x) 2 * min(mean(x > 0), mean(x < 0)))
post_mean  <- apply(combined_Sol, 2, mean)
ci         <- apply(combined_Sol, 2, quantile, probs = c(0.025, 0.975))
eff_samp   <- effectiveSize(model)
results <- data.frame(
  mean     = post_mean,
  lower    = ci[1, ],
  upper    = ci[2, ],
  eff_samp = eff_samp,
  pMCMC    = pMCMC
)
print(results)
write.table(results, file = "results.tsv", sep = "\t")

# Check specifically R1 vs all or all others
load("model.RData")
combined_Sol <- do.call(rbind, model)
te_cols <- grep("^TEtype", colnames(combined_Sol), value = TRUE)
length(te_cols)          # should be n_families - 1
r1_col <- "TEtypeLINE.R1"

# Family effects on the reference scale; reference family is 0 by construction
fam_effects <- cbind(`DNA.TcMar-Tigger` = 0,
                     combined_Sol[, te_cols, drop = FALSE])

# Contrast: R1 vs mean of all other families
others <- fam_effects[, colnames(fam_effects) != r1_col, drop = FALSE]
contrast_others <- combined_Sol[, r1_col] - rowMeans(others)

p <- 2 * min(mean(contrast_others > 0), mean(contrast_others < 0))
p <- max(p, 2 / length(contrast_others))
c(mean  = mean(contrast_others),
  lower = unname(quantile(contrast_others, 0.025)),
  upper = unname(quantile(contrast_others, 0.975)),
  p_dir = mean(contrast_others > 0),
  pMCMC = p)

# And get the % increase in R1 compared to others
fam_means <- df_mcmc %>%
  group_by(TEtype) %>%
  summarise(m = mean(Value, na.rm = TRUE))
others <- mean(fam_means$m[fam_means$TEtype != "LINE.R1"])
others
0.0007 / others * 100 # R1 has 67% highest contact values than the others
               
# Is there an association between R1 contacts and nb of fu/fi?
RM_summary_nfufi <- read.delim("RM_summary_nfufi.tsv")
df_mcmc_R1 <- df_mcmc %>%
  dplyr::filter(TEtype == "LINE.R1")
R1 <- df_mcmc_R1$Value
names(R1) <- df_mcmc_R1$Species
all(names(R1) %in% tree$tip.label) # to check if all names correspond

R1_df <- data.frame(
  species = names(R1),
  R1_value = as.numeric(R1)
) %>%
  dplyr::left_join(RM_summary_nfufi %>% dplyr::select(species, karyotype, nb_fusions, nb_fissions), by = "species")
rownames(R1_df) <- R1_df$species
tree_pruned <- drop.tip(tree, tree$tip.label[!tree$tip.label %in% R1_df$species])
model_R1_n <- phylolm(karyotype ~ R1_value, data = R1_df, phy = tree_pruned, model = "lambda")
summary(model_R1_n)
model_R1_fu <- phylolm(nb_fusions ~ R1_value, data = R1_df, phy = tree_pruned, model = "lambda")
summary(model_R1_fu)
model_R1_fi <- phylolm(nb_fissions ~ R1_value, data = R1_df, phy = tree_pruned, model = "lambda")
summary(model_R1_fi)

# MCMCglmm
tree <- force.ultrametric(tree)
invA <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv
my_priors_a <- list(
  R = list(V = 1, fix = TRUE),
  G = list(G1 = list(V = 1, nu = 0.002))
)
R1_df$animal <- as.factor(R1_df$species)
R1_df$animal <- factor(as.character(R1_df$animal))

# Karyotype
cl <- makeCluster(4)
clusterExport(cl, list("R1_df", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(karyotype ~ R1_value,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = R1_df,
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
clusterExport(cl, list("R1_df", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(nb_fusions ~ R1_value,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = R1_df,
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
clusterExport(cl, list("R1_df", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(nb_fissions ~ R1_value,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = R1_df,
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


# Does the tyndarus clade have higher interchromosomal contact for R1?
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

R1_test <- phylANOVA(tree, group, R1, nsim = 10000, p.adj = "fdr")
R1_test

# For testing just tyndarus vs others, do:
group <- ifelse(all_tips %in% tyndarus, "focal", "rest")
names(group) <- all_tips
R1_test_group <- phylANOVA(tree, group, R1, nsim = 10000, p.adj = "fdr")
R1_test_group

#### All TEs, normalization, rDNA clusters excluded ####
df <- read.table("allsp_interchrom_avg_no_rRNA_clusters_norm.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
remove_sp <- c("Mjurtina", "Mcinxia")

df_long <- df %>%
  pivot_longer(
    cols = -TEtype,
    names_to = "Species",
    values_to = "Value"
  )

# remove repeats that have interchromosomal contacts in less than 50% of species
df_long_filtered <- df_long %>%
  filter(!Species %in% remove_sp) %>%
  group_by(TEtype) %>%
  filter(mean(is.na(Value)) < 0.5) %>%
  ungroup()

df_long_filtered <- df_long_filtered %>%
  group_by(TEtype) %>%
  mutate(median_value = median(Value, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(TEtype = fct_reorder(TEtype, median_value))

# define which species to color
highlight <- c("C0055","C0080","C0100","Erondoui","X3258","X3531","X3737")
df_long_filtered$col_group <- ifelse(df_long_filtered$Species %in% highlight, df_long_filtered$Species, "Other")
df_long_filtered$col_group <- factor(df_long_filtered$col_group,
                                     levels = c(highlight, "Other"))
# Colour palette
colors <- c("Other" = "gray60",
            "C0055"="chocolate4",
            "C0080"="skyblue",
            "C0100"="#F0E442",
            "Erondoui"="#E69F00",
            "X3258"="#009E73",
            "X3531"="#0072B2",
            "X3737"="royalblue4")

allTEs <- ggplot(df_long_filtered, aes(x = TEtype, y = Value, color = col_group)) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.7, shape = 16) +
  geom_boxplot(outlier.shape = NA, color = "black",
               alpha = 0, width = 0.4, linewidth = 0.3) +
  scale_y_log10() +
  scale_color_manual(values = colors,
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
  scale_x_discrete(labels = function(x) gsub("\\.", "/", gsub("_", " ", x))) +
  coord_flip() +
  labs(x = "Repeat family", y = "Average interchromosomal contacts") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 8),
        axis.title.y = element_markdown(),
        axis.title.x = element_text(size = 10),
        legend.text = element_markdown(),
        legend.position = "bottom",
        legend.background = element_rect(fill = "grey90", color = NA),
        legend.key.spacing.y = unit(-2, "pt"),
        legend.key.height = unit(12, "pt"),
        legend.box.margin = margin(l = -90),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  theme(legend.key = element_rect(fill = "grey90", color = NA)) +
  guides(color = guide_legend(title = "Species", ncol = 2,
                              override.aes = list(size = 2, alpha = 1,
                                                  fill = "grey90")))
allTEs
ggsave("TE_full_plot_no_rDNA.pdf", allTEs, device = cairo_pdf, width = 4, height = 8)

# Stats with MCMCglmm
# Reformat data
tree <- read.tree("erebia_dated_phylo_clean_full.newick")
tree <- drop.tip(tree, c("Mjurtina", "Mcinxia"))
tree <- force.ultrametric(tree)
invA <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv
df_mcmc <- as.data.frame(df_long_filtered) %>%
  mutate(Species = factor(as.character(Species),
                          levels = rownames(invA)))
# Also include genomic proportion as a predictor of contact
te_keep <- unique(df_long_filtered$TEtype)
RM_summary_full <- read.delim("RM_summary_full.tsv", check.names = FALSE)
df_prop_long <- RM_summary_full %>%
  dplyr::rename(Species = species) %>%
  dplyr::select(Species, all_of(te_keep)) %>%
  pivot_longer(-Species, names_to = "TEtype", values_to = "genomic_prop")
# Also include average TE length per species as covariate
repeat_length_matrix <- read.delim("repeat_length_matrix.tsv", check.names = FALSE)
repeat_length_matrix <- repeat_length_matrix |> column_to_rownames("repeat_type")
repeat_length_matrix <- as.data.frame(t(repeat_length_matrix)) |>
  rownames_to_column("species")
df_length_long <- repeat_length_matrix %>%
  dplyr::rename(Species = species) %>%
  dplyr::select(Species, all_of(te_keep)) %>%
  pivot_longer(-Species, names_to = "TEtype", values_to = "TE_length")
# Join contacts and proportions and length
df_mcmc <- as.data.frame(df_long_filtered) %>%
  left_join(df_prop_long, by = c("Species", "TEtype")) %>%
  left_join(df_length_long, by = c("Species", "TEtype")) %>%
  mutate(Species = factor(as.character(Species), levels = rownames(invA)))
# DNA.TcMar-Tigger has the lowest average inter contacts: used as ref level
df_mcmc <- df_mcmc %>%
  mutate(TEtype = relevel(factor(TEtype), ref = "DNA.TcMar-Tigger"))

prior <- list(
  G = list(G1 = list(V = 1, nu = 0.002)),
  R = list(V = 1, nu = 0.002))

# Run 4 chains in parallel to check convergence
cl <- makeCluster(4)
clusterExport(cl, list("df_mcmc", "prior", "invA"))
clusterEvalQ(cl, library(MCMCglmm))
model_norDNA <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(Value ~ TEtype + genomic_prop + TE_length,
           random = ~ Species,
           ginverse = list(Species = invA),
           data = df_mcmc,
           prior = prior,
           nitt = 200000, burnin = 50000, thin = 100)
})
stopCluster(cl)

model_norDNA <- lapply(model_norDNA, function(m) m$Sol)
model_norDNA <- do.call(mcmc.list, model_norDNA)
gelman.diag(model_norDNA)
summary(model_norDNA)
save(model_norDNA, file = "model_norDNA.RData")
load("model_norDNA.RData")
summary(model_norDNA)
# Get pMCMC values (and CI) from pooled chains
combined_Sol <- do.call(rbind, model_norDNA)
pMCMC <- apply(combined_Sol, 2, function(x) 2 * min(mean(x > 0), mean(x < 0)))
post_mean  <- apply(combined_Sol, 2, mean)
ci         <- apply(combined_Sol, 2, quantile, probs = c(0.025, 0.975))
eff_samp   <- effectiveSize(model_norDNA)
results_norDNA <- data.frame(
  mean     = post_mean,
  lower    = ci[1, ],
  upper    = ci[2, ],
  eff_samp = eff_samp,
  pMCMC    = pMCMC
)
print(results_norDNA)
write.table(results_norDNA, file = "results_norDNA.tsv", sep = "\t")

