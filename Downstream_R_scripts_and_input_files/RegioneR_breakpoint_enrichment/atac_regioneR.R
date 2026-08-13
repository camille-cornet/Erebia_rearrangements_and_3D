library(tidyverse)
library(readxl)
library(patchwork)

regioneR_results <- read_excel("regioneR_results_atac_seq.xlsx")

#### Stats ####
# Stouffer's Z across species, using signed regioneR z-scores.
# Combined sign = overall direction; opposite-direction species cancel.
stouffer_combined <- function(zscores) {
  z <- zscores[is.finite(zscores)]
  z_comb <- sum(z) / sqrt(length(z))
  data.frame(
    z_combined = z_comb,
    p.value    = 2 * pnorm(-abs(z_comb)),   # two-sided: enrichment OR depletion
    direction  = ifelse(z_comb > 0, "enriched", "depleted"),
    n          = length(z),
    n_enriched = sum(z > 0),
    n_depleted = sum(z < 0)
  )
}

zscore_cols <- c("peak_nuc_zscore", "reads_zscore")

run_stouffer <- function(data, group_name) {
  map_dfr(zscore_cols, ~ {
    stouffer_combined(data[[.x]]) %>%
      mutate(compartment = gsub("_zscore", "", .x), group = group_name)
  })
}

stouffer_all <- bind_rows(
  run_stouffer(regioneR_results, "all"),
  run_stouffer(filter(regioneR_results, Type == "Fusion"),  "fusions"),
  run_stouffer(filter(regioneR_results, Type == "Fission"), "fissions")
) %>%
  select(group, compartment, z_combined, p.value,
         direction, n, n_enriched, n_depleted)

write.table(stouffer_all, file = "stouffer_results.tsv", sep = "\t")
