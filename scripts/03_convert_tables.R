## ---------------------------
##
## Script name: 03_convert_tables.R
##
## Purpose of script: This script loads relevant tables from the results and
## turns them into nice-looking LaTeX tables that can directly be loaded into
## the manuscript.
##
## Author: Raphael Scherrer
##
## Date Created: 2021-04-05
##
## This script comes with no guarantee whatsoever.
##
## Find me on GitHub at https://github.com/rscherrer
##
## Email:
## r.scherrer@rug.nl
## raphael.scherrer@evobio.eu
## raph.rjfs@hotmail.fr
##
## ---------------------------

rm(list = ls())

library(tidyverse)
library(knitr)
library(kableExtra)

#### 0. Accessory functions ####

# Function to add a significance column with asterisk symbols to a tibble
mutate_signif <- function(
  D, col = "pvalue", levels = c(0.05, 0.01, 0.001), dropname = TRUE
) {

  D <- D %>%
    mutate(
      .signif = ifelse(get(col) < levels[1], "*", ""),
      .signif = ifelse(get(col) < levels[2], "**", .signif),
      .signif = ifelse(get(col) < levels[3], "***", .signif)
    )

  if (dropname) D <- D %>% rename(" " = ".signif")

  return(D)

}

# Round and format P-values
round_pvalue <- function(pvalue, digits = 4) {

  threshold <- 1 / 10^digits
  template <- paste0("%.", sprintf("%sf", digits))
  label <- sprintf(paste("<", template), threshold, digits)
  rounded <- sprintf(template, round(pvalue, digits), digits)
  ifelse(pvalue < threshold, label, rounded)

}

#### 1. Random forest table ####

# Random forest results
read_csv("results/group_comparisons/machine_learning/PC/randomforest/summary.csv") %>%
  mutate_signif() %>%
  mutate(pvalue = round_pvalue(pvalue, digits = 4)) %>%
  kable(
    "latex",
    col.names = c("Island", "$N$", "Score", "$P$", ""),
    digits = c(0, 0, 3, 0, 0),
    align = "lrrrl",
    booktabs = TRUE,
    linesep = "",
    escape = FALSE
  ) %>%
  cat(file = "ms/tables/randomforests.tex")

#### 2. LDA table ####

# LDA results
read_csv("results/group_comparisons/machine_learning/PC/lda/summary.csv") %>%
  mutate_signif() %>%
  mutate(pvalue = round_pvalue(pvalue, digits = 4)) %>%
  kable(
    "latex",
    col.names = c(
      "Island", "$N$", "Score", "$P$", ""
    ),
    digits = c(0, 0, 3, 0, 0),
    align = "lrrrl",
    booktabs = TRUE,
    linesep = "",
    escape = FALSE
  ) %>%
  cat(file = "ms/tables/ldas.tex")

#### 3. SVM table ####

# SVM results
read_csv("results/group_comparisons/machine_learning/PC/ksvm/summary.csv") %>%
  mutate_signif() %>%
  mutate(pvalue = round_pvalue(pvalue, digits = 4)) %>%
  kable(
    "latex",
    col.names = c(
      "Island", "$N$", "Score", "$P$", ""
    ),
    digits = c(0, 0, 3, 0, 0),
    align = "lrrrl",
    booktabs = TRUE,
    linesep = "",
    escape = FALSE
  ) %>%
  cat(file = "ms/tables/ksvms.tex")

#### 4. ANOVA table ####

# ANOVA table
read_csv("results/group_comparisons/anovas.csv") %>%
  mutate_signif() %>%
  mutate(pvalue = round_pvalue(pvalue, digits = 4)) %>%
  kable(
    "latex",
    col.names = c(
      "Island", "Variable", "AICc", "$\\Delta$AICc", "AICw", "Model",
      "Log-lik.", "$\\chi^2$", "df", "$P$", ""
    ),
    digits = c(0, 0, 2, 2, 3, 0, 2, 2, 0, 0, 0),
    align = "llrrrlrrrrl",
    booktabs = TRUE,
    linesep = "",
    escape = FALSE
  ) %>%
  cat(file = "ms/tables/anova.tex")

#### 5. Spatial autocorrelation tables ####

# Mantel's test of spatial autocorrelation
read_csv("results/spatial_correlation/mantel_test.csv") %>%
  mutate_signif() %>%
  mutate(pvalue = round_pvalue(pvalue, digits = 3)) %>%
  kable(
    "latex",
    col.names = c(
      "Island", "$\\rho$", "$P$", ""
    ),
    digits = c(0, 3, 0, 0),
    align = "lrrl",
    booktabs = TRUE,
    linesep = "",
    escape = FALSE
  ) %>%
  cat(file = "ms/tables/mantel.tex")

# Spatial autocorrelation (custom permutation test)
read_csv("results/spatial_correlation/spatial_correlation.csv") %>%
  mutate_signif() %>%
  mutate(pvalue = round_pvalue(pvalue, digits = 3)) %>%
  kable(
    "latex",
    col.names = c(
      "Island", "$\\rho$", "$P$", ""
    ),
    digits = c(0, 3, 0, 0),
    align = "lrrl",
    booktabs = TRUE,
    linesep = "",
    escape = FALSE
  ) %>%
  cat(file = "ms/tables/autocorrelation.tex")

#### 6. Sample size table ####

# Sample sizes
read_csv("data/counts.csv") %>%
  rename(" " = "island") %>%
  kable(
    "latex",
    booktabs = TRUE,
    linesep = ""
  ) %>%
  cat(file = "ms/tables/counts.tex")

#### 7. Site metadata table ####

# Site metadata
read_csv("data/sites.csv") %>%
  kable(
    "latex",
    col.names = c(
      "Island", "Longitude", "Latitude", "Habitat", paste0("PC", 1:4)
    ),
    digits = c(0, 1, 1, 0, 3, 3, 3, 3),
    booktabs = TRUE,
    linesep = ""
  ) %>%
  cat(file = "ms/tables/sites.tex")

#### 8. PCA table ####

# PCA explained variance
read_csv("results/pc_expvars/pc_expvars.csv") %>%
  kable(
    "latex",
    col.names = c(
      "Island", "Total", paste0("PC", 1:4)
    ),
    digits = c(0, 3, 3, 3, 3, 3),
    booktabs = TRUE,
    linesep = ""
  ) %>%
  cat(file = "ms/tables/pcavariances.tex")

#### 9. Multivariate normality table ####

# Tests of multivariate normality
read_csv("results/assumptions/multinorm.csv") %>%
  arrange(island, habitat) %>%
  mutate_signif() %>%
  mutate(pvalue = round_pvalue(pvalue, digits = 4)) %>%
  kable(
    "latex",
    col.names = c(
      "Island", "Habitat", "Outliers", "$HZ$", "$P$", ""
    ),
    digits = c(0, 0, 0, 3, 0, 0),
    align = "llrrrl",
    booktabs = TRUE,
    linesep = "",
    escape = FALSE
  ) %>%
  cat(file = "ms/tables/multinorm.tex")

#### 10. Normality table ####

# Univariate normality
read_csv("results/assumptions/uninorm.csv") %>%
  mutate_signif(col = "padj") %>%
  mutate(padj = round_pvalue(padj, digits = 4)) %>%
  kable(
    "latex",
    col.names = c(
      "Island", "Variable", "$W$", "$P$", "$P_{adj}$", ""
    ),
    digits = c(0, 0, 3, 4, 0, 0),
    align = "llrrrl",
    booktabs = TRUE,
    linesep = "",
    escape = FALSE
  ) %>%
  cat(file = "ms/tables/normality.tex")

#### 11. Kruskal-Wallis table ####

# Kruskal-Wallis tests
read_csv("results/group_comparisons/kruskalwallis.csv") %>%
  mutate_signif() %>%
  mutate(pvalue = round_pvalue(pvalue, digits = 4)) %>%
  kable(
    "latex",
    col.names = c(
      "Island", "Variable", "$\\chi^2$", "df", "$P$", ""
    ),
    digits = c(0, 0, 2, 0, 0, 0),
    align = "llrrrl",
    booktabs = TRUE,
    linesep = "",
    escape = FALSE
  ) %>%
  cat(file = "ms/tables/kruskal.tex")
