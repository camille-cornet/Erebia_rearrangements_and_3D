library(ggplot2)
library(dplyr)
library(ggtext)
library(phylolm)
library(phytools)
library(patchwork)
library(MCMCglmm)
library(mulTree)

set.seed(123)

### Load data
df <- read.delim("compartment_strength_per_species.tsv")
df <- df %>% filter(!species %in% c("Mcinxia", "Mjurtina"))
df$median_autosome_size <- df$median_autosome_size/1000000

### Load tree
tree <- read.tree("erebia_dated_phylo_clean_full.newick")
tree <- drop.tip(tree, c("Mjurtina", "Mcinxia"))
rownames(df) <- df$species
setdiff(df$species, tree$tip.label)  # in data but not tree
setdiff(tree$tip.label, df$species)  # in tree but not data

### Run phylolm
# Take median chromosome size into account
model_size <- phylolm(strength ~ median_autosome_size,
                    data = df, phy = tree, model = "lambda")
summary(model_size)

#### MCMCglmm ####
tree <- force.ultrametric(tree)
invA <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv
my_priors_a <- list(
  R = list(V = 1, fix = TRUE),
  G = list(G1 = list(V = 1, nu = 0.002))
)
df$animal <- as.factor(df$species)
df$animal <- factor(as.character(df$animal))

# Karyotype
cl <- makeCluster(4)
clusterExport(cl, list("df", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(karyotype ~ strength,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = df,
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
clusterExport(cl, list("df", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(nb_fusions ~ strength,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = df,
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
clusterExport(cl, list("df", "my_priors_a", "invA"))
clusterEvalQ(cl, library(MCMCglmm))

fit <- parLapply(cl, 1:4, function(i) {
  MCMCglmm(nb_fissions ~ strength,
           random   = ~ animal,
           ginverse = list(animal = invA),
           data     = df,
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


### Comparisons between tyndarus and other clades
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

strength <- df$strength
names(strength) <- df$species
all(names(strength) %in% tree$tip.label) # to check if all names correspond
ratio_test <- phylANOVA(tree, group, strength, nsim = 10000, p.adj = "fdr")
ratio_test

# For testing just tyndarus vs others, do:
group <- ifelse(all_tips %in% tyndarus, "focal", "rest")
names(group) <- tree$tip.label
all(names(strength) %in% tree$tip.label) # to check if all names correspond
ratio_test <- phylANOVA(tree, group, strength, nsim = 10000, p.adj = "fdr")
ratio_test


### Plot
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


p_n <- ggplot(df, aes(x = strength, y = karyotype, color = col_group)) +
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
           label = "italic(p) == 0.696",
           parse = TRUE, size = 3.5) +
  theme_classic(base_size = 10) +
  theme(legend.text = element_markdown(),
        axis.title.y = element_markdown()) +
  labs(x = "Compartment strength", y = "Chromosome number (*n*)")
p_n

p_fu <- ggplot(df, aes(x = strength, y = nb_fusions, color = col_group)) +
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
           label = "italic(p) == 0.642",
           parse = TRUE, size = 3.5) +
  theme_classic(base_size = 10) +
  theme(legend.text = element_markdown()) +
  labs(x = "Compartment strength", y = "Number of fusions")
p_fu

p_fi <- ggplot(df, aes(x = strength, y = nb_fissions, color = col_group)) +
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
           label = "italic(p) == 0.908",
           parse = TRUE, size = 3.5) +
  theme_classic(base_size = 10) +
  theme(legend.text = element_markdown()) +
  labs(x = "Compartment strength", y = "Number of fissions")
p_fi

p_size <- ggplot(df, aes(y = median_autosome_size, x = strength, color = col_group)) +
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
           label = "italic(p) == 0.924",
           parse = TRUE, size = 3.5) +
  theme_classic(base_size = 10) +
  theme(legend.text = element_markdown()) +
  labs(x = "Compartment strength", y = "Median autosome size (Mb)")
p_size

full_p <- (p_size + p_n) / (p_fu + p_fi) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.background = element_rect(fill = "grey90", color = NA),
        legend.text = element_markdown(size = 10),
        legend.key = element_rect(fill = "grey90", color = NA)) &
  guides(color = guide_legend(nrow = 4),
         fill  = guide_legend(nrow = 4))
full_p

ggsave("comp_strength_full_plot.pdf", full_p, device = cairo_pdf, width = 6, height = 7)

