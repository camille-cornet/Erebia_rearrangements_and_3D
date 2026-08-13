library(ape)
library(phytools)
library(caper)
library(MCMCglmm)
library(mulTree)
library(tidyverse)
library(ggeffects)
library(patchwork)
library(parallel)
library(ggtext)

set.seed(123)

### controling for extant chromosome size

# Full dataset
data <- read.delim("data_per_merian_full.tsv")

# Remove rows with NA (corresponding to sex chrom)
data <- data %>% drop_na()
# Factors
data$Merian <- as.factor(data$Merian)
data$species <- as.factor(data$species)
data$animal <- as.factor(data$species)
data$animal <- factor(as.character(data$animal))

# Load tree
tree <- read.tree("erebia_dated_phylo_clean_full.newick")
tree <- drop.tip(tree, c("Mcinxia", "Mjurtina"))
tree <- force.ultrametric(tree)
is.ultrametric(tree)
plot(tree, cex = 0.6)

invA <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv

my_priors_a <- list(
  R = list(V = 1, fix = TRUE),
  G = list(
    G1 = list(V = 1, nu = 0.002),
    G2 = list(V = 1, nu = 0.002)
  )
)

# List of interesting variables to loop over
interesting_vars <- c(
  "avg_ratio_inter_intra",
  "strength",
  "avg_slope_loglog_0_10"
)

#### Loop over each variable for fusions ####
for (var in interesting_vars) {

  formula_fu <- as.formula(paste("nb_fusion ~", var, "+ chrom_size"))

  cl <- makeCluster(4)
  clusterExport(cl, list("data", "my_priors_a", "formula_fu", "invA"))
  clusterEvalQ(cl, library(MCMCglmm))

  fit <- parLapply(cl, 1:4, function(i) {
    MCMCglmm(formula_fu,
             random   = ~ animal + Merian,
             ginverse = list(animal = invA),
             data     = data,
             family   = "poisson",
             prior    = my_priors_a,
             nitt = 200000, burnin = 50000, thin = 100)
  })
  stopCluster(cl)

  # Check convergence
  fit_sol <- lapply(fit, function(m) m$Sol)
  fit_sol <- do.call(mcmc.list, fit_sol)

  cat("Gelman-Rubin diagnostics for", var, ":\n")
  print(gelman.diag(fit_sol))

  # Save MCMC chains
  save(fit_sol, file = paste0("fit_fu_", var, ".RData"))

  # Summarise pooled chains
  combined_Sol <- do.call(rbind, fit_sol)
  pMCMC    <- apply(combined_Sol, 2, function(x) 2 * min(mean(x > 0), mean(x < 0)))
  post_mean <- apply(combined_Sol, 2, mean)
  ci        <- apply(combined_Sol, 2, quantile, probs = c(0.025, 0.975))
  eff_samp  <- effectiveSize(fit_sol)

  results <- data.frame(
    variable = var,
    param    = names(post_mean),
    mean     = post_mean,
    lower    = ci[1, ],
    upper    = ci[2, ],
    eff_samp = eff_samp,
    pMCMC    = pMCMC
  )

  print(results)
  write.table(results, file = paste0("results_fu_", var, ".tsv"), sep = "\t", row.names = FALSE)
}


#### Loop over each variable for fissions ####
for (var in interesting_vars) {

  formula_fi <- as.formula(paste("nb_fission ~", var, "+ chrom_size"))

  cl <- makeCluster(4)
  clusterExport(cl, list("data", "my_priors_a", "formula_fi", "invA"))
  clusterEvalQ(cl, library(MCMCglmm))

  fit <- parLapply(cl, 1:4, function(i) {
    MCMCglmm(formula_fi,
             random   = ~ animal + Merian,
             ginverse = list(animal = invA),
             data     = data,
             family   = "poisson",
             prior    = my_priors_a,
             nitt = 200000, burnin = 50000, thin = 100)
  })
  stopCluster(cl)

  # Check convergence
  fit_sol <- lapply(fit, function(m) m$Sol)
  fit_sol <- do.call(mcmc.list, fit_sol)

  cat("Gelman-Rubin diagnostics for", var, ":\n")
  print(gelman.diag(fit_sol))

  # Save MCMC chains
  save(fit_sol, file = paste0("fit_fi_", var, ".RData"))

  # Summarise pooled chains
  combined_Sol <- do.call(rbind, fit_sol)
  pMCMC    <- apply(combined_Sol, 2, function(x) 2 * min(mean(x > 0), mean(x < 0)))
  post_mean <- apply(combined_Sol, 2, mean)
  ci        <- apply(combined_Sol, 2, quantile, probs = c(0.025, 0.975))
  eff_samp  <- effectiveSize(fit_sol)

  results <- data.frame(
    variable = var,
    param    = names(post_mean),
    mean     = post_mean,
    lower    = ci[1, ],
    upper    = ci[2, ],
    eff_samp = eff_samp,
    pMCMC    = pMCMC
  )

  print(results)
  write.table(results, file = paste0("results_fi_", var, ".tsv"), sep = "\t", row.names = FALSE)
}


#### Loop over each variable for extant chromosome size ####
for (var in interesting_vars) {

  formula_fi <- as.formula(paste("chrom_size ~", var))

  cl <- makeCluster(4)
  clusterExport(cl, list("data", "my_priors_a", "formula_fi", "invA"))
  clusterEvalQ(cl, library(MCMCglmm))

  fit <- parLapply(cl, 1:4, function(i) {
    MCMCglmm(formula_fi,
             random   = ~ animal + Merian,
             ginverse = list(animal = invA),
             data     = data,
             family   = "poisson",
             prior    = my_priors_a,
             nitt = 200000, burnin = 50000, thin = 100)
  })
  stopCluster(cl)

  # Check convergence
  fit_sol <- lapply(fit, function(m) m$Sol)
  fit_sol <- do.call(mcmc.list, fit_sol)

  cat("Gelman-Rubin diagnostics for", var, ":\n")
  print(gelman.diag(fit_sol))

  # Save MCMC chains
  save(fit_sol, file = paste0("fit_size_", var, ".RData"))

  # Summarise pooled chains
  combined_Sol <- do.call(rbind, fit_sol)
  pMCMC    <- apply(combined_Sol, 2, function(x) 2 * min(mean(x > 0), mean(x < 0)))
  post_mean <- apply(combined_Sol, 2, mean)
  ci        <- apply(combined_Sol, 2, quantile, probs = c(0.025, 0.975))
  eff_samp  <- effectiveSize(fit_sol)

  results <- data.frame(
    variable = var,
    param    = names(post_mean),
    mean     = post_mean,
    lower    = ci[1, ],
    upper    = ci[2, ],
    eff_samp = eff_samp,
    pMCMC    = pMCMC
  )

  print(results)
  write.table(results, file = paste0("results_size_", var, ".tsv"), sep = "\t", row.names = FALSE)
}



#### Plot ####
data <- data %>%
  mutate(Merian = factor(gsub("M20M17", "M17M20", as.character(Merian)),
                         levels = gsub("M20M17", "M17M20", levels(Merian))))
unique_vals <- unique(as.character(data$Merian))
is_Mnum <- grepl("^M[0-9]", unique_vals) & unique_vals != "MZ"
first_num <- rep(NA_integer_, length(unique_vals))
first_num[is_Mnum] <- as.integer(sub("^M([0-9]+).*", "\\1", unique_vals[is_Mnum]))
order_M <- order(first_num[is_Mnum],
                 nchar(unique_vals[is_Mnum]),
                 unique_vals[is_Mnum],
                 na.last = TRUE)
sorted_Mvals <- unique_vals[is_Mnum][order_M]
levels_final <- if ("MZ" %in% unique_vals) c(sorted_Mvals, "MZ") else sorted_Mvals
data$Merian <- factor(as.character(data$Merian), levels = levels_final)
levels(data$Merian)

# Merian colour palette
hex_colors = c(#"#666666", MZ not shown
  "#710093", "#BC007B", "#EA005B", "#FC1D1D", "#FD9514", "#E9CB19",
  "#87BF13", "#00AB3E", "#00F2A1", "#005C66", "#00589E", "#006DDB", "#0080FF",
  "#A676FF", "#FF1AEF", "#FF82CD", "#FF6B70", "#EE6A15", "#A16C00", #"#4E6400", (we combine M17 and M20)
  "#005200", "#00E251", "#009286", "#00B0E0", "#00C8FF", "#A1D0FF", "#ABA9E5",
  "#AE83B6", "#A56183", "#8E454F", "#6B3122"
)

levels_final <- if ("MZ" %in% unique_vals) c("MZ", sorted_Mvals) else sorted_Mvals
data$Merian <- factor(as.character(data$Merian), levels = levels_final)
names(hex_colors) <- levels_final
merian_color_scale <- scale_colour_manual(values = hex_colors)

# Row 1: chromosome size
size_slope_p <- ggplot(data, aes(x = avg_slope_loglog_all, y = chrom_size/1000000, colour = Merian)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  merian_color_scale +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) < 0.001",
           parse = TRUE, size = 3) +
  ylab("Extant chromosome size (Mb)") + xlab("P(s) slope (0 - 10 Mb range)") +
  labs(colour = "Merian element") +
  theme_minimal(base_size = 10)

size_strength_p <- ggplot(data, aes(x = strength, y = chrom_size/1000000, colour = Merian)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  merian_color_scale +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) == 0.033",
           parse = TRUE, size = 3) +
  ylab("Extant chromosome size (Mb)") + xlab("Compartment strength") +
  labs(colour = "Merian element") +
  theme_minimal(base_size = 10)

size_ratio_p <- ggplot(data, aes(x = avg_ratio_inter_intra, y = chrom_size/1000000, colour = Merian)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  merian_color_scale +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) < 0.001",
           parse = TRUE, size = 3) +
  ylab("Extant chromosome size (Mb)") + xlab("Inter-/intrachromosomal\ncontact ratio") +
  labs(colour = "Merian element") +
  theme_minimal(base_size = 10)

# Row 2: number of fusions
fu_slope_p <- ggplot(data, aes(x = avg_slope_loglog_all, y = nb_fusion, colour = Merian)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  merian_color_scale +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) == 0.763",
           parse = TRUE, size = 3) +
  ylab("Number of fusions") + xlab("P(s) slope (0 - 10 Mb range)") +
  labs(colour = "Merian element") +
  theme_minimal(base_size = 10)

fu_strength_p <- ggplot(data, aes(x = strength, y = nb_fusion, colour = Merian)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  merian_color_scale +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) == 0.547",
           parse = TRUE, size = 3) +
  ylab("Number of fusions") + xlab("Compartment strength") +
  labs(colour = "Merian element") +
  theme_minimal(base_size = 10)

fu_ratio_p <- ggplot(data, aes(x = avg_ratio_inter_intra, y = nb_fusion, colour = Merian)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  merian_color_scale +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) == 0.282",
           parse = TRUE, size = 3) +
  ylab("Number of fusions") + xlab("Inter-/intrachromosomal\ncontact ratio") +
  labs(colour = "Merian element") +
  theme_minimal(base_size = 10)

# Row 3: number of fissions
fi_slope_p <- ggplot(data, aes(x = avg_slope_loglog_all, y = nb_fission, colour = Merian)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  merian_color_scale +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) == 0.070",
           parse = TRUE, size = 3) +
  ylab("Number of fissions") + xlab("P(s) slope (0 - 10 Mb range)") +
  labs(colour = "Merian element") +
  theme_minimal(base_size = 10)

fi_strength_p <- ggplot(data, aes(x = strength, y = nb_fission, colour = Merian)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  merian_color_scale +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) == 0.426",
           parse = TRUE, size = 3) +
  ylab("Number of fissions") + xlab("Compartment strength") +
  labs(colour = "Merian element") +
  theme_minimal(base_size = 10)

fi_ratio_p <- ggplot(data, aes(x = avg_ratio_inter_intra, y = nb_fission, colour = Merian)) +
  geom_point(alpha = 0.3, shape = 16) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  merian_color_scale +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = "italic(p) == 0.781",
           parse = TRUE, size = 3) +
  ylab("Number of fissions") + xlab("Inter-/intrachromosomal\ncontact ratio") +
  labs(colour = "Merian element") +
  theme_minimal(base_size = 10)


# Remove x titles from rows 1 and 2
size_slope_p <- size_slope_p + theme(axis.title.x = element_blank())
size_ratio_p <- size_ratio_p + theme(axis.title.x = element_blank())
size_strength_p     <- size_strength_p     + theme(axis.title.x = element_blank())
fu_slope_p   <- fu_slope_p   + theme(axis.title.x = element_blank())
fu_ratio_p   <- fu_ratio_p   + theme(axis.title.x = element_blank())
fu_strength_p       <- fu_strength_p       + theme(axis.title.x = element_blank())

# Remove y titles from columns 2 and 3
size_ratio_p <- size_ratio_p + theme(axis.title.y = element_blank())
size_strength_p     <- size_strength_p     + theme(axis.title.y = element_blank())
fu_ratio_p   <- fu_ratio_p   + theme(axis.title.y = element_blank())
fu_strength_p       <- fu_strength_p       + theme(axis.title.y = element_blank())
fi_ratio_p   <- fi_ratio_p   + theme(axis.title.y = element_blank())
fi_strength_p       <- fi_strength_p       + theme(axis.title.y = element_blank())

# Combine
full_p <- (size_slope_p + size_strength_p + size_ratio_p) /
  (fu_slope_p + fu_strength_p + fu_ratio_p)   /
  (fi_slope_p + fi_strength_p + fi_ratio_p)   +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
full_p

ggsave("Merian_full_p.pdf", full_p, device = cairo_pdf, width = 6, height = 9)
