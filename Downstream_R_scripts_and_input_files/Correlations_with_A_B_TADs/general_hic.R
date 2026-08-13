library(readxl)
library(ggplot2)
library(dplyr)
library(ggtext)
library(phylolm)
library(phytools)
library(patchwork)
library(MCMCglmm)
library(mulTree)

set.seed(123)

### Load file
df <- read_excel("general_hic_data.xlsx")
df$TAD_size_Mb <- df$TAD_size/1000000
df$median_autosome_size_Mb <- df$median_autosome_size/1000000

### Load tree
tree <- read.tree("erebia_dated_phylo_clean_full.newick")
tree <- drop.tip(tree, c("Mjurtina", "Mcinxia"))
rownames(df) <- df$species
setdiff(df$species, tree$tip.label)  # in data but not tree
setdiff(tree$tip.label, df$species)  # in tree but not data


### Tests
responses  <- c("karyotype", "nb_fusions", "nb_fissions")
predictors <- c("A_perc", "B_perc", "TAD_size_Mb")

### Test with autosome size
size_TAD <- phylolm(TAD_size_Mb ~ median_autosome_size_Mb, data = df, phy = tree, model = "lambda")
summary(size_TAD)
size_A <- phylolm(A_perc ~ median_autosome_size_Mb, data = df, phy = tree, model = "lambda")
summary(size_A)
size_B <- phylolm(B_perc ~ median_autosome_size_Mb, data = df, phy = tree, model = "lambda")
summary(size_B)

#### MCMCglmm ####
tree <- force.ultrametric(tree)
invA <- inverseA(tree, nodes = "TIPS", scale = TRUE)$Ainv
my_priors_a <- list(
  R = list(V = 1, fix = TRUE),
  G = list(G1 = list(V = 1, nu = 0.002))
)
df$animal <- as.factor(df$species)
df$animal <- factor(as.character(df$animal))

run_mcmcglmm <- function(resp, pred, df, invA, prior,
                         nitt = 200000, burnin = 50000, thin = 100, n_chains = 4) {

  formula <- as.formula(paste(resp, "~", pred))

  cl <- makeCluster(n_chains)
  clusterExport(cl, list("df", "prior", "invA", "formula"), envir = environment())
  clusterEvalQ(cl, library(MCMCglmm))
  fit <- parLapply(cl, 1:n_chains, function(i) {
    MCMCglmm(formula,
             random   = ~ animal,
             ginverse = list(animal = invA),
             data     = df,
             family   = "poisson",
             prior    = prior,
             nitt = nitt, burnin = burnin, thin = thin)
  })
  stopCluster(cl)

  # Check convergence
  fit_sol <- lapply(fit, function(m) m$Sol)
  fit_sol <- do.call(mcmc.list, fit_sol)
  gd <- gelman.diag(fit_sol)
  cat("\n---", resp, "~", pred, "---\n")
  print(gd)

  # Save MCMC chains
  save(fit_sol, file = paste0("fit_", resp, "_", pred, ".RData"))

  # Summarise pooled chains
  combined_Sol <- do.call(rbind, fit_sol)
  pMCMC     <- apply(combined_Sol, 2, function(x) 2 * min(mean(x > 0), mean(x < 0)))
  post_mean <- apply(combined_Sol, 2, mean)
  ci        <- apply(combined_Sol, 2, quantile, probs = c(0.025, 0.975))
  eff_samp  <- effectiveSize(fit_sol)

  res <- data.frame(
    response  = resp,
    predictor = pred,
    param     = names(post_mean),
    mean      = post_mean,
    lower     = ci[1, ],
    upper     = ci[2, ],
    eff_samp  = eff_samp,
    pMCMC     = pMCMC,
    row.names = NULL
  )

  write.table(res, file = paste0("results_", resp, "_", pred, ".tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  res
}

all_results <- do.call(rbind, lapply(responses, function(resp) {
  do.call(rbind, lapply(predictors, function(pred) {
    run_mcmcglmm(resp, pred, df, invA, my_priors_a)
  }))
}))
all_results
write.table(all_results, file = "mcmcglmm_all_results.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)


### Is the tyndarus clade different?
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

results <- lapply(predictors, function(v) {
  x <- df[[v]]
  names(x) <- df$species
  test <- phylANOVA(tree, group, x, nsim = 10000, p.adj = "fdr")
  data.frame(
    variable = v,
    Fval = test$F,
    pval = test$Pf
  )
})
results_df <- do.call(rbind, results)
write.table(results_df,
            file = "phylANOVA_clade_results.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# For testing just tyndarus vs others, do:
group <- ifelse(all_tips %in% tyndarus, "focal", "rest")
names(group) <- tree$tip.label
results <- lapply(predictors, function(v) {
  x <- df[[v]]
  names(x) <- df$species
  test <- phylANOVA(tree, group, x, nsim = 10000, p.adj = "fdr")
  data.frame(
    variable = v,
    Fval = test$F,
    pval = test$Pf
  )
})
results_df <- do.call(rbind, results)
write.table(results_df,
            file = "phylANOVA_tyn_results.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
