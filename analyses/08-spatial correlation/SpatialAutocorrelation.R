# Spatial autocorrelation by permutation test

rm(list = ls())

library(tidyverse)
library(nmgc)

data <- read.csv("data/reflectance.csv", header = TRUE)

# Spatial autocorrelation test
res_cortest <- nspcortest(
  data, variables = paste0("PC", 1:4), nesting = "island", nperm = 1000, seed = 42,
  to_pcomp = paste0("wl", 300:400), keep = "habitat"
)

# Save table
t1 <- res_cortest$res
t1 <- t1 %>% add_signif()
save_table(t1, "table_spatial", digits = c(0, 3, 3, 0))

# Table with multivariate means per site
res_cortest$sites
write.csv(res_cortest$sites, "table_sites.csv", row.names = FALSE)

# Make a Latex table
t2 <- res_cortest$sites
save_table(t2, "table_sites", digits = c(0, 4, 4, 0, 4, 4, 4, 4))
