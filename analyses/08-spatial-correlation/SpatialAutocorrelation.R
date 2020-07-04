# Spatial autocorrelation by permutation test

rm(list = ls())

library(dewlap)
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
t1_fname <- "analyses/08-spatial correlation/table_spatial"
t1_names <- c("Island", "$r$", "$P$", "$N$", "")
save_table(t1, t1_fname, digits = c(0, 3, 3, 0), col.names = t1_names)

# Table with multivariate means per site
t2 <- res_cortest$sites
t2_fname <- "analyses/08-spatial correlation/table_sites"
t2_names <- c("Island", "Longitude", "Latitude", "Habitat", "PC1", "PC2", "PC3", "PC4")
save_table(t2, t2_fname, digits = c(0, 4, 4, 0, 4, 4, 4, 4), col.names = t2_names)
