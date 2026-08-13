library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggtext)
library(multcomp)
library(geiger)
library(phytools)
library(phylolm)
library(MCMCglmm)
library(parallel)

set.seed(123)

# Plot
files <- list.files(
  pattern = "_mapq30_expected_cis_100kb_smooth_slope.tsv$"
)
ps_all <- files %>%
  map_dfr(~ {
    read.table(.x, header = TRUE, sep = "\t") %>%
      mutate(species = sub("_mapq30_expected_cis_100kb_smooth_slope.tsv", "", .x))
  })
ps_all$species <- factor(ps_all$species)

# define which species to color
highlight <- c("C0055","C0080","C0100","Erondoui","X3258","X3531","X3737")
ps_all$col_group <- ifelse(ps_all$species %in% highlight, ps_all$species, "Other")

# Colour palette
colors <- c("Other" = "gray60",
            "C0055" = "chocolate4",
            "C0080" = "skyblue",
            "C0100" = "#F0E442",
            "Erondoui" = "#E69F00",
            "X3258" = "#009E73",
            "X3531" = "#0072B2",
            "X3737" = "royalblue4")

legend_labels <- c(
  "C0055"    = "*E. ottomana* (*n* = 40)",
  "C0080"    = "*E. calcaria* (*n* = 8)",
  "C0100"    = "*E. graucasica* (*n* = 51)",
  "Erondoui" = "*E. rondoui* (*n* = 24)",
  "X3258"    = "*E. nivalis* (*n* = 11)",
  "X3531"    = "*E. cassioides* (*n* = 10)",
  "X3737"    = "*E. tyndarus* (*n* = 10)",
  "Other"    = "Other species"
)

ps_all$col_group <- factor(
  ifelse(ps_all$species %in% highlight, as.character(ps_all$species), "Other"),
  levels = c(highlight, "Other")
)

p1 <- ggplot(ps_all, aes(
  x = dist_bp,
  y = contact_prob_smoothed,
  group = species,
  color = col_group
)) +
  geom_line(aes(alpha = col_group), linewidth = 0.7) +
  scale_x_log10(labels = scales::label_number(scale = 1e-6)) +
  scale_y_log10(labels = scales::label_number()) +
  scale_color_manual(values = colors, labels = legend_labels) +
  scale_alpha_manual(values = c("Other" = 0.25, setNames(rep(1, length(highlight)), highlight)), guide = "none") +
  labs(x = "Genomic distance (Mb)", y = "Smoothed P(s)", color = "Species") +
  theme_minimal() +
  theme(legend.text = element_markdown(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p2 <- ggplot(ps_all, aes(
  x = dist_bp,
  y = slope_loglog,
  group = species,
  color = col_group
)) +
  geom_line(aes(alpha = col_group), linewidth = 0.7) +
  scale_x_log10(labels = scales::label_number(scale = 1e-6)) +
  scale_color_manual(values = colors, labels = legend_labels) +
  scale_alpha_manual(values = c("Other" = 0.25, setNames(rep(1, length(highlight)), highlight)), guide = "none") +
  labs(x = "Genomic distance (Mb)", y = "P(s) slope", color = "Species") +
  theme_minimal() +
  theme(legend.text = element_markdown())

full_plot <- p1 / p2 +
  plot_layout(
    guides = "collect",
    heights = c(2, 1)
  ) &
  theme(legend.position = "none")
full_plot

ggsave("full_plot.pdf", full_plot, device = cairo_pdf, width = 5, height = 7)


#### stats ####
tree <- read.tree("erebia_dated_phylo_clean_full.newick")
tree <- drop.tip(tree, c("Mjurtina", "Mcinxia"))
tree <- force.ultrametric(tree)

slopes_summary <- ps_all %>%
  group_by(species) %>%
  summarise(
    mean_slope = mean(slope_loglog, na.rm = TRUE),
    # or restrict to a distance window first:
    mean_slope_0_10Mb = mean(slope_loglog[dist_bp <= 10e6], na.rm = TRUE),
  )
write.table(slopes_summary, file = "slopes_allsp.tsv", sep = "\t")

# Other dataset
RM_summary_nfufi <- read.delim("RM_summary_nfufi.tsv")
combined <- left_join(RM_summary_nfufi, slopes_summary, by = "species")
combined$median_autosome_size_Mb <- combined$median_autosome_size/1000000

# Check if same species in phylo and in dataset
combined <- combined %>% filter(!species %in% c("Mcinxia", "Mjurtina"))
combined$species <- as.factor(combined$species)
row.names(combined) <- combined$species
ere_check <- name.check(tree, combined)
ere_check

# n
fit_n <- phylolm(karyotype ~ mean_slope + median_autosome_size_Mb,
                  data = combined, phy = tree,
                  model = "lambda")
summary(fit_n)
fit_n_range <- phylolm(karyotype ~ mean_slope_0_10Mb + median_autosome_size_Mb,
                        data = combined, phy = tree,
                        model = "lambda")
summary(fit_n_range)

# Fusions
fit_fu <- phylolm(nb_fusions ~ mean_slope + median_autosome_size_Mb,
                        data = combined, phy = tree,
                        model = "lambda")
summary(fit_fu)
fit_fu_range <- phylolm(nb_fusions ~ mean_slope_0_10Mb + median_autosome_size_Mb,
               data = combined, phy = tree,
               model = "lambda")
summary(fit_fu_range)

# Fissions
fit_fi <- phylolm(nb_fissions ~ mean_slope + median_autosome_size_Mb,
                  data = combined, phy = tree,
                  model = "lambda")
summary(fit_fi)
fit_fi_range <- phylolm(nb_fissions ~ mean_slope_0_10Mb + median_autosome_size_Mb,
                        data = combined, phy = tree,
                        model = "lambda")
summary(fit_fi_range)

# Chrom size alone
fit_size <- phylolm(mean_slope ~ median_autosome_size_Mb,
                    data = combined, phy = tree,
                    model = "lambda")
summary(fit_size)
fit_size_range <- phylolm(mean_slope_0_10Mb ~ median_autosome_size_Mb,
                    data = combined, phy = tree,
                    model = "lambda")
summary(fit_size_range)

#### MCMCglmm ####
### Using MCMCglmm so I can implement a Poisson distribution
invA <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv
my_priors_a <- list(
  R = list(V = 1, fix = TRUE),
  G = list(G1 = list(V = 1, nu = 0.002))
  )
combined$animal <- as.factor(combined$species)
combined$animal <- factor(as.character(combined$animal))

# Karyotype 0-10
cl <- makeCluster(4)
clusterExport(cl, list("combined", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(karyotype ~ mean_slope_0_10Mb + median_autosome_size_Mb,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = combined,
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
save(fit_sol, file = "fit_n_0-10.RData")

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
write.table(results, file = paste0("results_n_0-10.tsv"), sep = "\t", row.names = FALSE)

# Karyotype not restricted range
cl <- makeCluster(4)
clusterExport(cl, list("combined", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(karyotype ~ mean_slope + median_autosome_size_Mb,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = combined,
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


# Fusions 0-10
cl <- makeCluster(4)
clusterExport(cl, list("combined", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(nb_fusions ~ mean_slope_0_10Mb + median_autosome_size_Mb,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = combined,
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
save(fit_sol, file = "fit_fu_0-10.RData")

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
write.table(results, file = paste0("results_fu_0-10.tsv"), sep = "\t", row.names = FALSE)

# Fusions not restricted range
cl <- makeCluster(4)
clusterExport(cl, list("combined", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(nb_fusions ~ mean_slope + median_autosome_size_Mb,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = combined,
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

# Fissions 0-10
cl <- makeCluster(4)
clusterExport(cl, list("combined", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(nb_fissions ~ mean_slope_0_10Mb + median_autosome_size_Mb,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = combined,
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
save(fit_sol, file = "fit_fi_0-10.RData")

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
write.table(results, file = paste0("results_fi_0-10.tsv"), sep = "\t", row.names = FALSE)

# Fissions not restricted range
cl <- makeCluster(4)
clusterExport(cl, list("combined", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(nb_fissions ~ mean_slope + median_autosome_size_Mb,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = combined,
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

# Does the tyndarus clade have steeper slope?
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

# For testing just tyndarus vs others, do:
#group <- ifelse(all_tips %in% tyndarus, "focal", "rest")

mean_slope <- combined$mean_slope
mean_slope_range <- combined$mean_slope_0_10Mb
names(mean_slope) <- combined$species
names(mean_slope_range) <- combined$species
all(names(mean_slope) %in% tree$tip.label) # to check if all names correspond
all(names(mean_slope_range) %in% tree$tip.label) # to check if all names correspond

mean_slope_test <- phylANOVA(tree, group, mean_slope, nsim = 10000, p.adj = "fdr")
mean_slope_test
mean_slope_range_test <- phylANOVA(tree, group, mean_slope_range, nsim = 10000, p.adj = "fdr")
mean_slope_range_test

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

combined <- combined |>
  mutate(
    col_group = ifelse(species %in% highlight, as.character(species), "Other"),
    col_group = factor(col_group, levels = c(highlight, "Other"))
  ) |>
  arrange(desc(col_group == "Other"))

p_size <- ggplot(combined, aes(x = mean_slope_0_10Mb, y = median_autosome_size / 1000000, color = col_group)) +
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
           label = "italic(p) < 0.001",
           parse = TRUE, size = 3.5) +
  theme_classic(base_size = 10) +
  theme(legend.text = element_markdown()) +
  labs(x = "Mean P(s) slope (0 - 10 Mb range)", y = "Median autosome size (Mb)")
p_size

p_n <- ggplot(combined, aes(x = mean_slope_0_10Mb, y = karyotype, color = col_group)) +
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
           label = "italic(p) == 0.522",
           parse = TRUE, size = 3.5) +
  labs(x = "Mean P(s) slope (0 - 10 Mb range)", y = "Chromosome number (*n*)") +
  theme_classic(base_size = 10) +
  theme(legend.text = element_markdown(),
        axis.title.y = element_markdown())
p_n

p_fu <- ggplot(combined, aes(x = mean_slope_0_10Mb, y = nb_fusions, color = col_group)) +
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
           label = "italic(p) == 0.093",
           parse = TRUE, size = 3.5) +
  theme_classic(base_size = 10) +
  theme(legend.text = element_markdown()) +
  labs(x = "Mean P(s) slope (0 - 10 Mb range)", y = "Number of fusions")
p_fu

p_fi <- ggplot(combined, aes(x = mean_slope_0_10Mb, y = nb_fissions, color = col_group)) +
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
           label = "italic(p) == 0.529",
           parse = TRUE, size = 3.5) +
  theme_classic(base_size = 10) +
  theme(legend.text = element_markdown()) +
  labs(x = "Mean P(s) slope (0 - 10 Mb range)", y = "Number of fissions")
p_fi

sizenfufi <- (p_size + p_n) / (p_fu + p_fi) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "grey90", color = NA),
        legend.text = element_markdown(size = 10),
        legend.key = element_rect(fill = "grey90", color = NA)) &
  guides(color = guide_legend(nrow = 3),
         fill  = guide_legend(nrow = 3))
sizenfufi

ggsave("sizenfufi.pdf", sizenfufi, device = cairo_pdf, width = 6, height = 7)

