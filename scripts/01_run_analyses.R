## ---------------------------
##
## Script name: analyses.R
##
## Purpose of script: This script contains all the analyses of this study.
##
## Author: Raphael Scherrer
##
## Date Created: 2021-02-06
##
## Copyright (c) Raphael Scherrer, 2021
##
## Email:
## r.scherrer@rug.nl
## raphael.scherrer@evobio.eu
## raph.rjfs@hotmail.fr
##
## ---------------------------

# Notation used for naming datasets in this script:
# D = the original dataset
# X = placeholder within a function
# P = temporary data frame for plotting or containing plots

rm(list = ls())

library(tidyverse)
library(ggridges) # ridge plots
library(broom) # tidy statistical tables
library(heplots) # Box M test
library(MVN) # multivariate normality
library(assertthat) # assertions within functions
library(rminer) # machine learning
library(geosphere) # geographical distances
library(MuMIn) # AICc
library(nlme) # univariate ANOVAs
library(PMCMRplus) # posthoc tests
library(RColorBrewer) # color palette
library(pairwiseAdonis) # posthoc PERMANOVA
library(scales) # color palette
library(vegan) # PERMANOVA
library(ape) # Mantel test

#### 1. Set-up ####

# Set plotting theme
theme_set(theme_classic())

# Graphical parameters
island_colors <- brewer.pal(8, "Set2")
island_colors <- c(island_colors, "darkgreen")
habitat_colors <- c("goldenrod", "forestgreen", "mediumseagreen")

# Global settings
pc_names <- paste0("PC", 1:4)
wpc_names <- paste0("w", pc_names)
wl <- 300:700
wl_names <- paste0("wl", wl)

# Whether to load machines already fitted to whole spectrum data
# (instead of rerunning the very lengthy classification)
read_fitted <- TRUE

#### 2. Read reflectance data ####

# Read the data
D <- read_csv("data/reflectance.csv")

# Useful things to keep in mind
island_names <- unique(D$island)
habitat_names <- unique(D$habitat)

# Save sample sizes
with(D, table(island, habitat)) %>%
  as_tibble() %>%
  pivot_wider(names_from = habitat, values_from = n) %>%
  write_csv("metadata/counts.csv")

#### 2.1. Plot reflectance profiles ####

# Prepare to plot reflectance
P <- D %>%
  select(island, habitat, wl_names, id) %>%
  pivot_longer(cols = wl_names) %>%
  mutate(wl = as.numeric(str_remove(name, "wl")))

# Function to plot all reflectance profiles
PLOTFUN <- function(P) {

  P %>%
    ggplot(aes(x = wl, y = value, group = id, color = habitat)) +
    geom_line(alpha = 0.5) +
    scale_color_manual(values = habitat_colors) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(x = "Wavelength (nm)", y = "Reflectance (%)", color = "Habitat")

}

# Plot all reflectance profiles
P %>%
  PLOTFUN() +
  facet_wrap(. ~ island)

figname <- "results/reflectance/refl_profiles.png"
ggsave(figname, width = 6, height = 5, dpi = 300)

# Or each island separately
P <- P %>%
  group_by(island) %>%
  nest()

with(P, walk2(island, data, function(island, P) {

  P %>% PLOTFUN()

  figname <- sprintf("results/reflectance/refl_profiles_%s.png", island)
  ggsave(figname, width = 4, height = 2, dpi = 300)

}))

#### 2.2. Plot reflectance ribbons ####

# Compute reflectance ribbons (interquantile range)
P <- D %>%
  select(island, habitat, wl_names) %>%
  pivot_longer(cols = wl_names) %>%
  mutate(wl = as.numeric(str_remove(name, "wl"))) %>%
  group_by(wl, habitat, island) %>%
  summarize(
    q05 = quantile(value, 0.05),
    q95 = quantile(value, 0.95)
  )

# Function to plot the ribbons
PLOTFUN <- function(P) {

  P %>%
    ggplot(aes(x = wl, ymin = q05, ymax = q95, color = habitat, fill = habitat)) +
    geom_ribbon(alpha = 0.3) +
    scale_color_manual(values = habitat_colors) +
    scale_fill_manual(values = habitat_colors) +
    labs(
      x = "Wavelength (nm)", y = "Reflectance (%)", color = "Habitat",
      fill = "Habitat"
    ) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

}

# Plot the ribbons
P %>%
  PLOTFUN() +
  facet_wrap(. ~ island)

figname <- "results/reflectance/refl_ribbons.png"
ggsave(figname, width = 6, height = 5, dpi = 300)

# Or each island separately
P <- P %>%
  group_by(island) %>%
  nest()

with(P, walk2(island, data, function(island, P) {

  P %>% PLOTFUN()

  figname <- sprintf("results/reflectance/refl_ribbons_%s.png", island)
  ggsave(figname, width = 4, height = 2, dpi = 300)

}))

#### 2.3. Boxplots at various wavelengths ####

# Boxplots at various wavelengths
P <- D %>%
  group_by(island) %>%
  nest()

with(P, walk2(island, data, function(island, X) {

  plot <- X %>%
    select(habitat, paste0("wl", seq(350, 650, 100))) %>%
    pivot_longer(wl350:wl650) %>%
    mutate(name = str_c(str_remove(name, "wl"), "nm")) %>%
    ggplot(aes(x = habitat, y = value, color = habitat)) +
    facet_wrap(. ~ name) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_color_manual(values = habitat_colors) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 60, hjust = 1)
    ) +
    labs(x = NULL, y = "Reflectance (%)")

  figname <- sprintf("results/refl_boxplots/boxplot_refl_%s.png", island)
  ggsave(figname, plot, width = 4, height = 4, dpi = 300)

}))

#### 3. Principal component analysis ####

# Perform a PCA first because wavelengths are highly correlated
PCA <- prcomp(D[, wl_names], center = TRUE, scale = TRUE)

# How many PCs to retain?
cumsum(PCA$sdev / sum(PCA$sdev))

# Retain explained variance
expvars <- PCA$sdev[1:4] / sum(PCA$sdev)

# With four PCs we explain more than 90% of the variance
PCs <- PCA$x[, 1:4]

# Scale the principal components to unit variance (they are centered by constr)
PCs <- scale(PCs, scale = TRUE)

# Attach the matrix of principal components to the dataset
D <- cbind(D, PCs)

#### 3.1. Scatterplots ####

saveto <- "results/pc_global"

# Scatterplots showing the PCs against each other in various ways
for (i in 1:3) {
  for (j in (i + 1):4) {

    pci <- paste0("PC", i)
    pcj <- paste0("PC", j)

    # Prepare the labels
    pcilab <- paste0(pci, " (", round(100 * expvars[i], 1), "%)")
    pcjlab <- paste0(pcj, " (", round(100 * expvars[j], 1), "%)")

    # Points colored by island
    ggplot(D, aes(x = get(pci), y = get(pcj), color = island)) +
      geom_point() +
      coord_fixed() +
      scale_color_manual(values = island_colors) +
      labs(x = pcilab, y = pcjlab, color = "Island")

    figname <- sprintf("%s/plot_%s_%s_by_island.png", saveto, pci, pcj)
    ggsave(figname, width = 4, height = 4, dpi = 300)

    # Points colored by habitat and facetted by island
    ggplot(D, aes(x = get(pci), y = get(pcj), color = habitat)) +
      geom_point() +
      coord_fixed() +
      scale_color_manual(values = habitat_colors) +
      facet_wrap(~ island) +
      labs(x = pcilab, y = pcjlab, color = "Habitat")

    figname <- sprintf(
      "%s/plot_%s_%s_by_island_by_habitat_facetted.png", saveto, pci, pcj
    )
    ggsave(figname, width = 8, height = 5, dpi = 300)

    # Points colored by habitat
    ggplot(D, aes(x = get(pci), y = get(pcj), color = habitat)) +
      geom_point() +
      coord_fixed() +
      scale_color_manual(values = habitat_colors) +
      labs(x = pcilab, y = pcjlab, color = "Habitat")

    sprintf("%s/plot_%s_%s_by_habitat.png", saveto, pci, pcj)
    ggsave(figname, width = 4, height = 4, dpi = 300)

    # Points colored by island and shaped by habitat
    ggplot(
      D,
      aes(x = get(pci), y = get(pcj), color = island, shape = habitat)
    ) +
      geom_point() +
      coord_fixed() +
      scale_color_manual(values = island_colors) +
      labs(x = pcilab, y = pcjlab, color = "Island", shape = "Habitat")

    figname <- sprintf(
      "%s/plot_%s_%s_by_island_by_habitat.png", saveto, pci, pcj
    )
    ggsave(figname, width = 4, height = 4, dpi = 300)

  }
}

#### 3.2. Density plots ####

# Check the distributions of each principal component
for (i in 1:4) {

  pc <- paste0("PC", i)

  # Plot density across islands and habitats
  ggplot(D, aes(x = get(pc), y = island, fill = habitat)) +
    geom_density_ridges(alpha = 0.5) +
    scale_fill_manual(values = habitat_colors) +
    labs(y = NULL, fill = "Habitat", x = pc)

  figname <- sprintf("%s/densities_%s.png", saveto, pc)
  ggsave(figname, width = 4, height = 4, dpi = 300)

}

#### 4. Assumption testing ####

# Test multivariate normality and outliers within each group
D %>%
  group_by(habitat, island) %>%
  nest() %>%
  mutate(

    # Perform tests on each group
    mvn = map(data, function(data) {

      mvn(data[, paste0("PC", 1:4)], mvnTest = "hz", showOutliers = TRUE)

    }),

    # Henze-Zirkler test results
    hztest = map(mvn, ~ .x$multivariateNormality),

    # Multivariate outliers
    outliers = map(mvn, ~ .x$multivariateOutliers),

    # Number of outliers
    noutliers = map_int(outliers, length),

    # Extract statistic and P-value from HZ test
    HZ = map_dbl(hztest, ~ .x[["HZ"]]),
    pvalue = map_dbl(hztest, ~ .x[["p value"]])

  ) %>%
  select(-data, -mvn, -hztest, -outliers) %>%
  unnest(cols = c()) %>%
  ungroup() %>%
  select(island, habitat, noutliers, HZ, pvalue) %>%
  write_csv("results/assumptions/multinorm.csv")

# Second, homogeneity of covariance matrices across groups
boxM(PCs, unite(D[, c("island", "habitat")], "")[[1]])

# We now know that our groups are not normally distributed in multivariate space
# and their covariance matrices are not homogeneous (although Box's
# M test is very sensitive to non-normality and to many groups, Warner 2013)
# https://www.researchgate.net/post/MANOVA_violation_of_assumptions-what_now
# That said, visual inspection suggests that covariance matrices are indeed
# not homogeneous

#### 5. Two-way multivariate analysis of variance ####

# Now we use a two-way MANOVA to ask about differences
summary(manova(PCs ~ island * habitat, data = D))

# All effects are very significant
# But the assumption of normality is violated so we back our results with
# a nonparametric equivalent

set.seed(42)

# PERMANOVA
adonis(dist(PCs) ~ island * habitat, data = D, permutations = 1000)

# The PERMANOVA also gives us estimates of effect sizes R2
# They tell us that there is an interaction, so we go for within-island tests

#### 6. Within-island principal component analysis ####

# Perform PCA within each island
D <- D %>%
  group_by(island) %>%
  nest() %>%
  mutate(

    # Within-island principal component analysis
    PCA = map(data, ~ prcomp(.x[, wl_names], center = TRUE, scale = TRUE)),

    # Extract rotated data
    PCs = map(PCA, function(PCA) {

      PCs <- as_tibble(scale(PCA$x[, 1:4], scale = TRUE))
      colnames(PCs) <- paste0("w", colnames(PCs)) # label to recognize them
      return(PCs)

    })
  )

# Extract the fitted within-island PCAs
wPCA <- D$PCA
names(wPCA) <- island_names

# Unnest the data frame
D <- D %>%
  select(-PCA) %>%
  unnest(cols = c(data, PCs))

# Now the dataset is complete
# saveRDS(D, "data.rds")

# Set specific colors for the retained PCs
wpc_colors <- hue_pal()(length(wpc_names))
names(wpc_colors) <- str_remove(wpc_names, "w")

#### 6.1. Boxplots across islands ####

# Nest by island
P <- D %>%
  group_by(island) %>%
  nest()

saveto <- "results/pc_boxplots"

# Naive boxplots of within-island PC
with(P, walk2(data, island, function(X, island) {

  plot <- X %>%
    pivot_longer(cols = wPC1:wPC4) %>%
    mutate(name = str_remove(name, "w")) %>%
    ggplot(aes(x = habitat, y = value, color = habitat)) +
    facet_wrap(. ~ name) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_classic() +
    scale_color_manual(values = habitat_colors) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      legend.position = "none"
    ) +
    labs(x = NULL, y = "Value")

  figname <- sprintf("%s/pc_boxplot_%s.png", saveto, island)
  ggsave(figname, plot, width = 4, height = 4, dpi = 300)

}))

#### 6.2. Correlation between PCs and wavelengths ####

# Make a data frame with correlation between PCs and wavelengths
P <- D %>%
  group_by(island) %>%
  nest() %>%
  mutate(

    # Correlations
    P = map(data, function(X) {

      P <- map_dfc(wpc_names, function(wpc) {
        map_dbl(wl_names, function(wl) {

          cor(X[[wpc]], X[[wl]])

        })
      })

      colnames(P) <- str_remove(wpc_names, "w")

      P <- P %>%
        mutate(wl = wl) %>%
        pivot_longer(cols = PC1:PC4)

    })

  ) %>%
  select(-data) %>%
  unnest(cols = c(P))

saveto <- "results/pc_correlations"

# Plot the correlations
P %>%
  ggplot(aes(x = wl, y = value, color = name)) +
  facet_wrap(. ~ island) +
  geom_line() +
  scale_color_manual(values = wpc_colors) +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Wavelength (nm)", y = "Correlation", color = NULL) +
  ylim(c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

figname <- sprintf("%s/pc_correlations.png", saveto)
ggsave(figname, width = 5, height = 4, dpi = 300)

# Heatmap version
P %>%
  ggplot(aes(x = wl, y = name, fill = value)) +
  facet_wrap(. ~ island) +
  geom_tile() +
  scale_fill_gradient2(
    low = "darkblue", mid = "white", high = "darkgreen", limits = c(-1, 1)
  ) +
  labs(x = "Wavelength (nm)", fill = "Correlation", y = NULL) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

figname <- sprintf("%s/pc_correlation_heatmaps.png", saveto)
ggsave(figname, width = 5, height = 4, dpi = 300)

# Or islands separately
P <- P %>%
  group_by(island) %>%
  nest()

# Make and save the plots
with(P, walk2(data, island, function(X, island) {

  plot <- X %>%
    ggplot(aes(x = wl, y = value, color = name)) +
    geom_line() +
    scale_color_manual(values = wpc_colors) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "Wavelength (nm)", y = "Correlation", color = NULL) +
    ylim(c(-1, 1))

  figname <- sprintf("%s/pc_correlation_%s.png", saveto, island)
  ggsave(figname, width = 4, height = 2, dpi = 300)

  plot <- X %>%
    ggplot(aes(x = wl, y = name, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "darkblue", mid = "white", high = "darkgreen", limits = c(-1, 1)
    ) +
    labs(x = "Wavelength (nm)", fill = "Correlation", y = NULL) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  figname <- sprintf("%s/pc_correlation_heatmap_%s.png", saveto, island)
  ggsave(figname, width = 4, height = 2, dpi = 300)

}))

#### 6.3. Rotation matrices ####

# Prepare to plot the rotation matrices
P <- map2_dfr(wPCA, island_names, function(PCA, island) {

  PCA$rotation[, 1:4] %>%
    as_tibble() %>%
    mutate(island = island) %>%
    mutate(wl = wl)

}) %>%
  pivot_longer(cols = PC1:PC4)

saveto <- "results/pc_rotations"

# Plot rotation coefficients across islands
P %>%
  ggplot(aes(x = wl, y = value, color = name)) +
  geom_line() +
  scale_color_manual(values = wpc_colors) +
  theme_classic() +
  geom_hline(yintercept = 0, lty = 2) +
  labs(x = "Wavelength (nm)", y = "Coefficient", color = NULL) +
  facet_wrap(. ~ island)

figname <- sprintf("%s/pc_rotations.png", saveto)
ggsave(figname, width = 5, height = 4, dpi = 300)

# Heatmap version
P %>%
  ggplot(aes(x = wl, y = name, fill = value)) +
  facet_wrap(. ~ island) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkgreen") +
  labs(x = "Wavelength (nm)", fill = "Coefficient", y = NULL) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

figname <- sprintf("%s/pc_rotation_heatmaps.png", saveto)
ggsave(figname, width = 5, height = 4, dpi = 300)

# Or for each island separately
P <- P %>%
  group_by(island) %>%
  nest()

# Line plots
with(P, walk2(data, island, function(X, island) {

  plot <- X %>%
    ggplot(aes(x = wl, y = value, color = name)) +
    geom_line() +
    scale_color_manual(values = wpc_colors) +
    theme_classic() +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "Wavelength (nm)", y = "Coefficient", color = NULL)

  figname <- sprintf("%s/pc_rotation_%s.png", saveto, island)
  ggsave(figname, width = 4, height = 2, dpi = 300)

  plot <- X %>%
    ggplot(aes(x = wl, y = name, fill = value)) +
    geom_tile() +
    scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkgreen") +
    labs(x = "Wavelength (nm)", fill = "Coefficient", y = NULL) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

  figname <- sprintf("%s/pc_rotation_heatmap_%s.png", saveto, island)
  ggsave(figname, plot, width = 4, height = 2, dpi = 300)

}))

#### 6.4. Explained variance ####

# Make a table with the variance explained by each PC on each island
X <- map(wPCA, ~ with(.x, sdev[1:4] / sum(sdev))) %>%
  do.call("rbind", .) %>%
  as_tibble(.name_repair = ) %>%
  mutate(expvar = rowSums(.)) %>%
  mutate(island = island_names) %>%
  setNames(c(pc_names, "expvar", "island")) %>%
  select(island, expvar, pc_names)

# Record explained variance for the global PCA on all islands
pc_expvar <- with(PCA, sdev[1:4] / sum(sdev))
names(pc_expvar) <- pc_names

# Prepare extra row for all islands
new_row <- c(
  island = "All islands",
  expvar = sum(pc_expvar),
  as.list(pc_expvar)
)

# Append extra row
X <- X %>% add_row(!!!new_row)

# Save the table
write_csv(X, "results/pc_expvars/pc_expvars.csv")

#### 6.5. Mapping of the PC onto the reflectance curves ####

# Plot the reflectance curves and their PC scores
P <- D %>%
  group_by(island) %>%
  nest()

saveto <- "results/pc_reflectance"

# Multiple plots
with(P, walk2(data, island, function(X, island) {

  plot <- X %>%
    pivot_longer(cols = wl_names, names_to = "wl", values_to = "refl") %>%
    mutate(wl = as.numeric(str_remove(wl, "wl"))) %>%
    pivot_longer(cols = wpc_names) %>%
    mutate(name = str_remove(name, "w")) %>%
    ggplot(aes(x = wl, y = refl, group = id, color = value)) +
    geom_line(alpha = 0.5) +
    scale_color_gradient2(
      low = "darkblue", mid = "lightgray", high = "darkgreen"
    ) +
    facet_wrap(. ~ name) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
    labs(x = "Wavelength (nm)", y = "Reflectance (%)", color = "PC score")

  figname <- sprintf("%s/pc_reflectance_%s.png", saveto, island)
  ggsave(figname, width = 4, height = 3, dpi = 300)

}))

#### 7. Within-island assumption testing ####

# Within-island multivariate normality assumption testing
P <- D %>%
  group_by(island, habitat) %>%
  nest() %>%
  mutate(

    # Perform tests on each group
    mvn = map(data, function(X) {

      MVN::mvn(X[, wpc_names], mvnTest = "hz", showOutliers = TRUE)

    }),

    # Henze-Zirkler test results
    hztest = map(mvn, ~ .x$multivariateNormality),

    # Multivariate outliers
    outliers = map(mvn, ~ .x$multivariateOutliers),

    # Number of outliers
    noutliers = map_int(outliers, length),

    # Extract statistic and P-value from HZ test
    HZ = map_dbl(hztest, ~ .x[["HZ"]]),
    pvalue = map_dbl(hztest, ~ .x[["p value"]])

  ) %>%
  select(-data, -mvn, -hztest, -outliers)

P %>% write_csv("results/assumptions/multinorm_within.csv")

# How many violations of normality per island?
P %>%
  group_by(island) %>%
  summarize(NV = length(which(pvalue < 0.05))) %>%
  write_csv("results/assumptions/multinorm_within_nviolations.csv")

# Within-island test of homogeneity of covariance matrices
D %>%
  group_by(island) %>%
  nest() %>%
  mutate(

    # Test of homogeneity of covariance matrices
    boxM = map(data, function(X) boxM(X[wpc_names], X$habitat)),

    # Extract important statistics
    chisq = map_dbl(boxM, ~ .x[["statistic"]]),
    df = map_int(boxM, ~ as.integer(.x[["parameter"]])),
    pvalue = map_dbl(boxM, ~ .x[["p.value"]])

  ) %>%
  select(-data, -boxM) %>%
  write_csv("results/assumptions/boxM_within.csv")

#### 8. Within-island multivariate group comparison ####

#### 8.1. MANOVA ####

D %>%
  group_by(island) %>%
  nest() %>%
  mutate(

  # Perform test
  fit = map(data, ~ manova(as.matrix(.x[, wpc_names]) ~ habitat, data = .x)),

  # Extract important statistics
  pillai = map_dbl(fit, ~ tidy(.x)[["pillai"]][1]),
  pseudoF = map_dbl(fit, ~ tidy(.x)[["statistic"]][1]),
  df1 = map_int(fit, ~ as.integer(tidy(.x)[["num.df"]][1])),
  df2 = map_int(fit, ~ as.integer(tidy(.x)[["den.df"]][1])),
  pvalue = map_dbl(fit, ~ tidy(.x)[["p.value"]][1])

) %>%
  select(pillai, pseudoF, df1, df2, pvalue) %>%
  write_csv("results/group_comparisons/manova_within.csv")

#### 8.2. PERMANOVA ####

set.seed(24)

D %>%
  group_by(island) %>%
  nest() %>%
  mutate(

  # Perform test
  fit = map(data, function(X) {

      vegan::adonis(
        dist(X[, wpc_names]) ~ habitat,
        data = X, permutations = 1000
      )

    }
  ),

  # Extract important statistics
  pseudoF = map_dbl(fit, ~ tidy(.x$aov.tab)[["F.Model"]][1]),
  df = map_int(fit, ~ as.integer(tidy(.x$aov.tab)[["df"]][1])),
  pvalue = map_dbl(fit, ~ tidy(.x$aov.tab)[["p.value"]][1]),
  r2 = map_dbl(fit, ~ tidy(.x$aov.tab)[["R2"]][1])

) %>%
  select(pseudoF, df, pvalue, r2) %>%
  write_csv("results/group_comparisons/permanova_within.csv")

#### 8.3. Machine learning ####

#### 8.3.1. Set-up ####

# Function to run machine learning classification on a dataset
# Note: we do not rely on rminer's mining function because of its limited
# options
MLFUN <- function(

  data,
  variables,
  groups,
  model,
  search = NULL,
  k = 5,
  n = 5,
  minsize = 5

) {

  formula <- formula(paste(groups, "~", paste(variables, collapse = " + ")))

  # Keep only the predictor and grouping variables, in that order
  # This is for the Importance function not to be annoying later on
  data <- data %>% select(variables, all_of(groups))

  # Create a vector of bin labels based on the number of bins
  bins <- rep(seq(k), each = floor(1 / k * nrow(data)))

  # Add some if needed to make sure all observations are included
  if (length(bins) < nrow(data)) {
    bins <- c(bins, seq(nrow(data) - length(bins)))
  }

  # Check that they are all included
  assert_that(length(bins) == nrow(data))

  # For each replicate run
  runs <- map(seq(n), function(r) {

    is_fine <- FALSE # to know whether to re-sample the partition

    while(!is_fine) {

      # Shuffle the labels
      bins <- sample(bins)

      # Check whether each class is sufficiently represented
      is_fine <- map_lgl(seq(k), function(k) {
        all(table(data[bins != k, groups]) >= minsize)
      })
      is_fine <- all(is_fine)

    }

    # For each bin in the cross-validation procedure...
    kfits <- map(seq(k), function(k) {

      # What is the sample of the least represented group in the training set?
      least <- min(table(data[[groups]][bins != k]))

      # Down-sample each group in the training set to ensure balanced sizes
      training <- unlist(map(
        unique(data[[groups]]),
        ~ sample(which(bins != k & data[[groups]] == .x), least)
      ))

      # Fit a single machine on the training set
      # Note: data must be a data.frame, not a tibble!
      fit <- fit(
        formula,
        data = as.data.frame(data[training,]),
        model = model,
        task = "class",
        search = search,
        importance = TRUE
      )

      # Test the machine on the testing set
      pred <- predict(fit, newdata = data[bins == k,])
      pred <- factor(pred, levels = unique(data[[groups]]))
      test <- factor(data[[groups]][bins == k], levels = unique(data[[groups]]))

      # Compute the confusion matrix
      confmat <- table(pred, test)

      RFI <- NULL

      # Random forest-specific feature importance analysis
      if (model == "randomforest") RFI <- fit@object$importance

      # 1D sensitivity analysis
      SA1 <- Importance(fit, data.frame(data[training,]))
      SA1 <- SA1$imp[seq(variables)]

      # 2D sensitivity analysis
      SA2 <- matrix(0, length(variables), length(variables))
      colnames(SA2) <- rownames(SA2) <- variables

      for (i in seq(variables)) {
        for (j in i:length(variables)) {

          SA2[i,j] <- Importance(
            fit,
            data.frame(data[training,]),
            method = "GSA",
            interactions = c(i, j)
          )$value

        }
      }

      return(list(fit, confmat, RFI, SA1, SA2))

    })

    names(kfits) <- seq(k)

    # Sum the confusion matrices over bins
    confmat <- Reduce("+", map(kfits, ~ .x[[2]]))

    # Importance table over bins
    RFI <- map_dfr(
      seq(kfits),
      ~ mutate(as_tibble(kfits[[.x]][[3]]), CV = .x)
    )

    # One-dimensional sensitivity analysis table
    SA1 <- t(map_dfr(kfits, ~ .x[[4]]))
    SA1 <- data.frame(SA1)
    colnames(SA1) <- variables
    SA1$CV <- seq(nrow(SA1))

    # List of two-dimensional sensitivity analysis matrices
    SA2 <- map(kfits, ~ .x[[5]])
    SA2 <- tibble(SA2 = SA2, CV = seq(k), R = r)

    return(list(confmat, RFI, SA1, SA2, kfits))

  })

  # Average confusion matrix over runs
  confmat <- nmgc::mavg(map(runs, ~ .x[[1]]))

  # Importance data frame
  RFI <- map_dfr(runs, ~ .x[[2]], .id = "R")
  RFI <- RFI[, c(2:ncol(RFI), 1)]
  RFI <- tibble(
    cbind(variable = rep(variables, nrow(RFI) / length(variables)), RFI)
  )

  # One dimensional sensitivity analysis data frame
  SA1 <- map_dfr(seq(runs), ~ mutate(tibble(runs[[.x]][[3]]), R = .x))

  # Two dimensional sensitivity analysis data frame
  SA2 <- map_dfr(runs, ~ .x[[4]])

  # Fitted machines
  machines <- map(runs, ~ .x[[5]])

  return(
    list(
      confmat = confmat, RFI = RFI, SA1 = SA1, SA2 = SA2, machines = machines
    )
  )

}

# Machine learning models available
models <- c("lda", "ksvm", "randomforest")

# Hyperparameter search setup
searches <- list(

  lda = NULL,

  ksvm = list(
    smethod = "grid",
    search = list(C = seq(1, 10, 1)),
    method = c("holdouto", 2/3),
    convex = 0
  ),

  randomforest = list(
    smethod = "grid",
    search = list(ntree = c(500, 1000), mtry = 1:4),
    method = c("holdouto", 2/3),
    convex = 0
  )

)

# Function to save a plot for each island separately
PLOTISLANDS <- function(X, PLOTFUN, figname, width = 3.5, height = 2.5) {

  # Split by island
  X <- X %>%
    group_by(island) %>%
    nest()

  # Make space for the island label
  figname <- str_replace(figname, ".png", "_%s.png")

  # Save each plot
  with(X, walk2(data, island, function(X, island) {

    # Make the plot for that island
    plot <- X %>%
      PLOTFUN()

    figname <- sprintf(figname, island)
    ggsave(figname, plot, width = width, height = height, dpi = 300)

  }))

}

#### 8.3.2. Run on principal components ####

set.seed(53)

# Split the data by island
P <- D %>%
  group_by(island) %>%
  nest() %>%
  ungroup()

# For each type of model...
for (i in seq(models)) {

  model <- models[i]
  search <- searches[[model]]

  # Run the machine learning classification
  ML <- map(P$data, MLFUN, wpc_names, "habitat", model, search)
  names(ML) <- island_names

  # Where to save the results?
  saveto <- paste0("results/group_comparisons/machine_learning/PC/", model)

  # Summarize the results
  S <- P %>%
    mutate(

      n = map_int(data, nrow),
      score = map_dbl(ML, ~ sum(diag(.x$confmat)) / sum(.x$confmat)),
      pvalue = map2_dbl(score, n, ~ 1 - pbinom(.x * .y, .y, 1/3))

    ) %>%
    select(-data)

  # Save the summary
  S %>%
    write_csv(sprintf("%s/summary.csv", saveto))

  # Extract average confusion matrix for each island (normalized by column)
  Q <- map2_dfr(ML, names(ML), function(ML, island) {

    ML$confmat %>% '/'(colSums(.)) %>%
      as_tibble() %>%
      mutate(island = island)

  })

  # Function to plot a confusion matrix
  PLOTFUN <- function(X) {

    X %>%
      ggplot(aes(x = test, y = pred, fill = n)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "darkgreen", limits = c(0, 1)) +
      labs(x = "True", y = "Predicted", fill = "Frequency") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1))

  }

  # Plot across islands
  Q %>% PLOTFUN() + facet_wrap(. ~ island)
  figname <- paste0(saveto, "/confusion.png")
  ggsave(figname, width = 5, height = 4.5, dpi = 300)
  PLOTISLANDS(Q, PLOTFUN, figname, width = 3.5, height = 2.5)

  # Reduce to the significant islands
  ML <- ML[S$pvalue < 0.05]

  # Sensitivity analysis
  Q <- map2_dfr(ML, names(ML), ~ .x$SA1 %>% mutate(island = .y)) %>%
    pivot_longer(cols = wpc_names) %>%
    mutate(name = str_remove(name, "w"))

  # Function for plotting the sensitivity analysis
  PLOTFUN <- function(X) {

    X %>%
      ggplot(aes(x = name, y = value)) +
      geom_point(mapping = aes(group = interaction(CV, R))) +
      geom_path(mapping = aes(group = interaction(CV, R)), alpha = 0.5) +
      labs(x = NULL, y = "Importance")

  }

  Q %>% PLOTFUN() + facet_wrap(. ~ island)
  figname <- paste0(saveto, "/importance1D.png")
  ggsave(figname, width = 5, height = 4, dpi = 300)
  PLOTISLANDS(Q, PLOTFUN, figname, width = 3, height = 2)

  # Two-dimensional sensitivity analysis
  Q <- map2_dfr(ML, names(ML), function(ML, island) {

    mavg(ML$SA2$SA2) %>%
      as.data.frame() %>%
      rownames_to_column() %>%
      mutate(island = island)

  }) %>%
    pivot_longer(cols = wpc_names) %>%
    mutate(name = str_remove(name, "w")) %>%
    mutate(rowname = str_remove(rowname, "w"))

  # Function for plotting the 2D sensitivity analysis
  PLOTFUN <- function(X) {

    X %>%
      ggplot(aes(x = rowname, y = name, size = value, alpha = value)) +
      geom_point() +
      scale_alpha(range = c(0, 1)) +
      labs(x = NULL, y = NULL, size = "Importance", alpha = "Importance")

  }

  Q %>% PLOTFUN() + facet_wrap(. ~ island)
  figname <- paste0(saveto, "/importance2D.png")
  ggsave(figname, width = 5, height = 4, dpi = 300)
  PLOTISLANDS(Q, PLOTFUN, figname, width = 3, height = 2)

  if (model == "randomforest") {

    # Random forest-specific feature importance analysis
    Q <- map2_dfr(ML, names(ML), function(ML, island) {

      ML$RFI %>%
        mutate(island = island)

    }) %>%
      pivot_longer(cols = habitat_names) %>%
      mutate(variable = str_remove(variable, "w"))

    # Function to plot random forest-specific feature importance
    PLOTFUN <- function(X) {

      X %>%
        ggplot(
          aes(x = variable, y = value, color = name, group = interaction(CV, R))
        ) +
        geom_line(alpha = 0.5) +
        theme(
          axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "top"
        ) +
        scale_color_manual(values = habitat_colors) +
        labs(color = "Habitat", y = "Importance", x = NULL)

    }

    # Plot forest-specific feature importance across islands
    Q %>% PLOTFUN() + facet_grid(name ~ island)
    figname <- paste0(saveto, "/importanceRF.png")
    ggsave(figname, width = 7, height = 4, dpi = 300)

    # Or split by island
    Q <- Q %>%
      group_by(island) %>%
      nest()

    # Make place for the island label
    figname <- str_replace(figname, ".png", "_%s.png")

    # Save each plot
    with(Q, walk2(data, island, function(X, island) {

      # Make the plot for that island
      plot <- X %>%
        PLOTFUN() +
        facet_grid(name ~ .) +
        theme(legend.position = "right")

      figname <- sprintf(figname, island)
      ggsave(figname, plot, width = 4, height = 3, dpi = 300)

    }))

  }
}

#### 8.3.3. Significant islands ####

# Significant islands
islands <- c(
  "Abaco", "Bimini", "Cayman Brac", "Eleuthera", "Little Cayman",
  "Long Island", "North Andros", "South Andros"
)

#### 8.3.4. Run on reflectance data ####

set.seed(35)

# Random forest only, takes too long
model <- models[3]
search <- searches[[model]]

# Now the same thing but with reflectance instead of PCs
# Note: this can take several hours

P <- D %>%
  group_by(island) %>%
  nest()

# Run the machine learning classification with reflectance every 5nm
wl_names_reduced <- wl_names[seq(1, length(wl_names), 5)]

# Where to save the results?
saveto <- paste0("results/group_comparisons/machine_learning/WL/", model)

if (read_fitted) {

  ML <- readRDS(paste0(saveto, "/fitted.rds"))

} else {

  ML <- map(P$data, MLFUN, wl_names_reduced, "habitat", model, search)
  names(ML) <- island_names

  saveRDS(ML, paste0(saveto, "/fitted.rds"))

}

# Summarize the results
D %>%
  group_by(island) %>%
  nest() %>%
  ungroup() %>%
  mutate(

    n = map_int(data, nrow),
    score = map_dbl(ML, ~ sum(diag(.x$confmat)) / sum(.x$confmat)),
    pvalue = map2_dbl(score, n, ~ 1 - pbinom(.x * .y, .y, 1/3))

  ) %>%
  select(-data) %>%
  write_csv(sprintf("%s/summary.csv", saveto))

# Extract average confusion matrix for each island (normalized by column)
Q <- map2_dfr(ML, names(ML), function(ML, island) {

  ML$confmat %>% '/'(colSums(.)) %>%
    as_tibble() %>%
    mutate(island = island)

})

# Function to plot a confusion matrix
PLOTFUN <- function(X) {

  X %>%
    ggplot(aes(x = test, y = pred, fill = n)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "darkgreen", limits = c(0, 1)) +
    labs(x = "True", y = "Predicted", fill = "Frequency") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

}

# Plot across islands
Q %>% PLOTFUN() + facet_wrap(. ~ island)
figname <- paste0(saveto, "/confusion.png")
ggsave(figname, width = 5, height = 4.5, dpi = 300)
PLOTISLANDS(Q, PLOTFUN, figname, width = 3.5, height = 2.5)

# Reduce to the significant islands
ML <- ML[islands]

# Random forest-specific feature importance analysis
Q <- map2_dfr(ML, names(ML), function(ML, island) {

  ML$RFI %>%
    mutate(island = island)

}) %>%
  pivot_longer(cols = habitat_names) %>%
  mutate(variable = as.numeric(str_remove(variable, "wl")))

# Function to plot random forest-specific feature importance
PLOTFUN <- function(X) {

  X %>%
    ggplot(
      aes(x = variable, y = value, color = name, group = interaction(CV, R))
    ) +
    geom_line(alpha = 0.5) +
    theme(
      axis.text.x = element_text(angle = 60, hjust = 1),
      legend.position = "top"
    ) +
    scale_color_manual(values = habitat_colors) +
    labs(color = "Habitat", y = "Importance", x = "Wavelength (nm)")

}

# Plot forest-specific feature importance across islands
Q %>% PLOTFUN() + facet_grid(name ~ island)
figname <- paste0(saveto, "/importanceRF.png")
ggsave(figname, width = 7, height = 4, dpi = 300)

# Or split by island
Q <- Q %>%
  group_by(island) %>%
  nest()

# Make place for the island label
figname <- str_replace(figname, ".png", "_%s.png")

# Save each plot
with(Q, walk2(data, island, function(X, island) {

  # Make the plot for that island
  plot <- X %>%
    PLOTFUN() +
    facet_grid(name ~ .) +
    theme(legend.position = "right")

  figname <- sprintf(figname, island)
  ggsave(figname, plot, width = 4, height = 3, dpi = 300)

}))

#### 9. Test of spatial autocorrelation ####

# Save site metadata
D %>%
  group_by(island, longitude, latitude, habitat) %>%
  summarize(across(wpc_names, mean)) %>%
  write_csv("metadata/sites.csv")

#### 9.1. Permutation test ####

set.seed(1994)

# Permutation test
D %>%
  group_by(island) %>%
  nest() %>%
  mutate(

    # Get a vector of observed and permuted correlation coefficients
    R = map(data, function(X) {

      map_dbl(seq(1001), function(i) {

        if (i > 1) {

          # Shuffle the individuals across sites on the island
          jj <- sample(nrow(X))
          X[c("longitude", "latitude")] <- X[jj, c("longitude", "latitude")]

        }

        # Compute mean phenotypes per site
        S <- X %>%
          group_by(longitude, latitude) %>%
          nest() %>%
          mutate(

            data = map(data, function(data) {

              tibble(
                variable = wpc_names,
                mean = colMeans(data[wpc_names])
              )

            })

          ) %>%
          unnest(cols = c(data)) %>%
          pivot_wider(names_from = "variable", values_from = "mean")

        # Matrix of geographic distances between sites
        G <- as.dist(
          distm(
            S[c("longitude", "latitude")], fun = distGeo
          )
        )

        # Matrix of euclidean distances between site means in phenotype space
        M <- dist(S[wpc_names])

        # Correlation between the two matrices
        cor(G, M)

      })

    }),

    # Extract observed correlation
    robs = map_dbl(R, ~ .x[1]),

    # Compute P-value
    pvalue = map_dbl(
      R,
      ~ 1 - length(which(robs > .x[2:length(.x)])) / (length(.x) - 1)
    )

  ) %>%
  select(robs, pvalue) %>%
  write_csv("results/spatial_correlation/spatial_correlation.csv")

#### 9.2. Mantel's test ####

set.seed(99)

D %>%
  group_by(island) %>%
  nest() %>%
  mutate(

    # Get a vector of observed and permuted correlation coefficients
    mantel_test = map(data, function(X) {

      # Compute mean phenotypes per site
      S <- X %>%
        group_by(longitude, latitude) %>%
        nest() %>%
        mutate(

          data = map(data, function(data) {

            tibble(
              variable = wpc_names,
              mean = colMeans(data[wpc_names])
            )

          })

        ) %>%
        unnest(cols = c(data)) %>%
        pivot_wider(names_from = "variable", values_from = "mean")

      # Matrix of geographic distances between sites
      G <- as.dist(
        distm(
          S[c("longitude", "latitude")], fun = distGeo
        )
      )

      # Matrix of euclidean distances between site means in phenotype space
      M <- dist(S[wpc_names])

      # Mantel's test of correlation between these two matrices
      test_result <- mantel(G, M)
      test_stat <- test_result$statistic
      p_value <- test_result$signif

      return(list(test_stat, p_value))

    }),

    # Extract the test statistics
    mantel = map_dbl(mantel_test, ~ .x[[1]]),
    pvalue = map_dbl(mantel_test, ~ .x[[2]])

  ) %>%
  select(mantel, pvalue) %>%
  write_csv("results/spatial_correlation/mantel_test.csv")

#### 10. Univariate signal decomposition ####

#### 10.1. Multiple univariate analyses of variance ####

# Multiple analyses of variance
P <- D %>%
  filter(island %in% islands) %>%
  pivot_longer(cols = wpc_names, names_to = "variable") %>%
  mutate(variable = str_remove(variable, "w")) %>%
  group_by(island, variable) %>%
  nest() %>%
  mutate(

    anova = map(data, function(X) {

      # Fit a GLS model with one residual variance per habitat
      M <- gls(
        value ~ habitat,
        data = X,
        weights = varIdent(form = ~ 1|habitat)
      )

      # Make a second model fitted with OLS
      M2 <- update(M, weights = NULL)

      # Compute AICc
      AIC1 <- AICc(M)
      AIC2 <- AICc(M2)

      # Compute AICc weights
      W <- exp(-0.5 * c(0, AIC2 - AIC1))
      W <- W / sum(W)

      # Compute the difference in AIC
      delta <- AIC1 - AIC2

      # Retain the best-fitting model
      if (delta > 0) {
        M <- M2
        W <- W[2]
        AIC <- AIC2
      } else {
        W <- W[1]
        AIC <- AIC1
      }

      # Re-fit the model with maximum likelihood (instead of REML)
      M <- update(M, method = "ML")

      # Test the habitat effect by using a likelihood ratio test
      M0 <- update(M, ~ 1)
      LRT <- anova(M, M0)

      # Note: we cannot update the model outside of this function, so we
      # compute AIC here because that needs to be done prior to model
      # update

      # Check the assumption of normality of the residuals
      shapirotest <- shapiro.test(residuals(M, type = "pearson"))

      # Note: the normality of the raw data around their group-means
      # does not equate to normality of the residuals within a GLS model,
      # because the (standardized) residuals may differ between the groups

      return(

        # Return test statistics
        list(
          AICc = AIC,
          dAIC = delta,
          AICw = W,
          model = ifelse(delta > 0, "OLS", "GLS"),
          loglik = LRT$logLik[1],
          lratio = LRT$L.Ratio[2],
          df = as.integer(with(LRT, df[1] - df[2])),
          pvalue = LRT$"p-value"[2],
          shapirotest = shapirotest
        )
      )
    })
  ) %>%
  select(-data) %>%
  unnest_wider(anova)

P %>%
  select(-shapirotest) %>%
  write_csv("results/group_comparisons/anovas.csv")

# Check univariate normality
P %>%
  select(shapirotest) %>%
  mutate(
    W = map_dbl(shapirotest, ~ .x$statistic),
    pvalue = map_dbl(shapirotest, ~ .x$p.value)
  ) %>%
  select(-shapirotest) %>%
  ungroup() %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  write_csv("results/assumptions/uninorm.csv")

#### 10.2. Univariate PERMANOVAs ####

P <- read_csv("results/assumptions/uninorm.csv")

set.seed(52)

# Nonparametric follow-up for deviations from normality (PERMANOVA)
D %>%
  pivot_longer(cols = wpc_names, names_to = "variable") %>%
  mutate(variable = str_remove(variable, "w")) %>%
  group_by(island, variable) %>%
  nest() %>%
  right_join(P) %>%
  filter(padj < 0.05) %>%
  mutate(nonparametric = map(data, function(X) {

    # Perform univariate PERMANOVA
    permanova <- vegan::adonis(
      dist(X$value) ~ habitat,
      data = X, permutations = 1000
    )

    return(

      # Extract important statistics
      list(
        pseudoF = permanova$aov.tab$F.Model[1],
        df = tidy(permanova$aov.tab)$df[1],
        pvalue = tidy(permanova$aov.tab)$p.value[1],
        r2 = tidy(permanova$aov.tab)$R2[1]
      )
    )

  })) %>%
  select(-data, -W, -pvalue, -padj) %>%
  unnest_wider(nonparametric) %>%
  write_csv("results/group_comparisons/unipermanovas.csv")

#### 10.3. Kruskal-Wallis tests ####

P <- read_csv("results/assumptions/uninorm.csv")

# Nonparametric follow-up for deviations from normality (Kruskal-Wallis)
D %>%
  pivot_longer(cols = wpc_names, names_to = "variable") %>%
  mutate(variable = str_remove(variable, "w")) %>%
  group_by(island, variable) %>%
  nest() %>%
  right_join(P) %>%
  filter(padj < 0.05) %>%
  mutate(nonparametric = map(data, function(X) {

    # Perform a Kruskal-Wallis test
    kw = with(X, kruskal.test(value, habitat))

    return(

      # Extract important statistics
      list(
        chisq = tidy(kw)$statistic[1],
        df = tidy(kw)$parameter[1],
        pvalue = tidy(kw)$p.value[1]
      )
    )

  })) %>%
  select(-data, -W, -pvalue, -padj) %>%
  unnest_wider(nonparametric) %>%
  write_csv("results/group_comparisons/kruskalwallis.csv")

#### 10.4. Comparing parametric and nonparametric results ####

# Load both results
P1 <- read_csv("results/group_comparisons/anovas.csv")
P2 <- read_csv("results/group_comparisons/kruskalwallis.csv")

# Compare the P-values
P2 %>%
  mutate(

    # Get the P-values from the parametric and nonparametric tests aligned
    pparam = map2_dbl(
      island, variable, ~ P1$pvalue[P1$island == .x & P1$variable == .y]
    ),

    # Is the parametric test significant when the nonparametric test is?
    is_same = (pvalue < 0.05) == (pparam < 0.05) |
      (pvalue >= 0.05) == (pparam >= 0.05)

  ) %>%
  rename(pnonparam = "pvalue") %>%
  select(-chisq, -df) %>%
  write_csv("results/group_comparisons/test_comparison.csv")

#### 11. Post-hoc tests: which groups differ from which? ####

# Here we perform different post-hoc tests depending on which model assumptions
# are met: Tukey, Dunnett, Nemenyi and pairwise PERMANOVAs

P <- D %>%
  pivot_longer(wpc_names, names_to = "variable") %>%
  mutate(variable = str_remove(variable, "w")) %>%
  group_by(island, variable) %>%
  nest()

nhabitats <- length(habitat_names)
ncontrasts <- nhabitats * (nhabitats - 1) / 2

# Prepare a data frame for the contrasts
M <- tibble(
  i = rep(habitat_names[-nhabitats], each = ncontrasts - 1),
  j = rep(habitat_names[-1], ncontrasts - 1)
) %>%
  filter(i != j)

# We want to restrict the post-hoc tests to a set of selected contrasts
# which depend on what univariate group comparisons test was performed
# which itself depends on what assumptions were met

# Normality data set
N <- read_csv("results/assumptions/uninorm.csv") %>%
  select(-W, -pvalue)

#### 11.1. Tukey test ####

# Tukey test for normally distributed variables with equal variances
read_csv("results/group_comparisons/anovas.csv") %>%
  group_by(island, variable, model, pvalue) %>%
  nest() %>%
  left_join(N) %>%
  filter(padj >= 0.05) %>%
  select(-data, -padj) %>%
  right_join(P) %>%
  filter(model == "OLS", pvalue < 0.05) %>%
  mutate(posthoc = map(data, function(X) {

    M %>%
      mutate(
        pvalue = c(tukeyTest(X$value, factor(X$habitat))[[3]]) %>%
          .[!is.na(.)]
      )

  })) %>%
  ungroup() %>%
  select(-data, -model, -pvalue) %>%
  unnest(cols = c(posthoc)) %>%
  write_csv("results/group_comparisons/posthoc/tukey.csv")

#### 11.2. Dunnett test ####

# Dunnett test for heteroskedastic variables
read_csv("results/group_comparisons/anovas.csv") %>%
  group_by(island, variable, model, pvalue) %>%
  nest() %>%
  left_join(N) %>%
  filter(padj >= 0.05) %>% # keep only the normally distributed
  select(-data, -padj) %>%
  right_join(P) %>%
  filter(model == "GLS", pvalue < 0.05) %>%
  mutate(posthoc = map(data, function(X) {

    M %>%
      mutate(
        pvalue = c(dunnettT3Test(X$value, factor(X$habitat))[[3]]) %>%
          .[!is.na(.)]
      )

  })) %>%
  ungroup() %>%
  select(-data, -model, -pvalue) %>%
  unnest(cols = c(posthoc)) %>%
  write_csv("results/group_comparisons/posthoc/dunnett.csv")

#### 11.3. Nemenyi test ####

# Nemenyi test for non-normal variables
read_csv("results/group_comparisons/kruskalwallis.csv") %>%
  group_by(island, variable, pvalue) %>%
  nest() %>%
  select(-data) %>%
  right_join(P) %>%
  filter(pvalue < 0.05) %>%
  mutate(posthoc = map(data, function(X) {

    M %>%
      mutate(
        pvalue = c(kwAllPairsNemenyiTest(X$value, factor(X$habitat))[[3]]) %>%
          .[!is.na(.)]
      )

  })) %>%
  ungroup() %>%
  select(-data, -pvalue) %>%
  unnest(cols = c(posthoc)) %>%
  write_csv("results/group_comparisons/posthoc/nemenyi.csv")

#### 11.4. Multiple pairwise PERMANOVAs ####

set.seed(33)

# Multiple pairwise PERMANOVAs
read_csv("results/group_comparisons/kruskalwallis.csv") %>%
  group_by(island, variable, pvalue) %>%
  nest() %>%
  select(-data) %>%
  right_join(P) %>%
  filter(pvalue < 0.05) %>%
  mutate(posthoc = map(data, function(X) {

    M %>%
      mutate(
        pvalue =
          pairwise.adonis(
            dist(X$value), factors = X$habitat, perm = 1000,
            p.adjust.m = "BH"
          )$p.adjusted
      )

  })) %>%
  ungroup() %>%
  select(-data, -pvalue) %>%
  unnest(cols = c(posthoc)) %>%
  write_csv("results/group_comparisons/posthoc/permanova.csv")

#### 11.5. Post-hoc summary ####

map_dfr(c("tukey", "dunnett", "nemenyi"), function(test) {

  read_csv(sprintf("results/group_comparisons/posthoc/%s.csv", test)) %>%
    mutate(test = test)

}) %>%
  write_csv("results/group_comparisons/posthoc/summary.csv")

#### 12. Brightness ####

saveto <- "results/brightness"

# Compute the brightness of each dewlap as mean reflectance
P <- D %>%
  rowwise() %>%
  mutate(brightness = mean(c_across(wl300:wl700)))

# Check the correlation between brightness and global PC1
with(P, cor.test(PC1, brightness))

# Plot it (the P-value is so small we do not need to show it)
P %>%
  ggplot(aes(x = PC1, y = brightness)) +
  geom_point() +
  labs(x = "PC1", y = "Brightness (%)") +
  annotate(
    geom = "text",
    label = paste0("r = ", round(with(P, cor(PC1, brightness)), 3), "***"),
    x = min(P$wPC1),
    y = max(P$brightness),
    hjust = 0
  )

figname <- sprintf("%s/brightness.png", saveto)
ggsave(figname, width = 3, height = 2, dpi = 300)

# Check correlation between brightness and within-island PC1
labels <- P %>%
  group_by(island) %>%
  summarize(
    r = cor(wPC1, brightness),
    pvalue = cor.test(wPC1, brightness)$p.value
  )

# The P-values are again extremely small
labels

# Arrange the labels for pretty display
labels <- labels %>% mutate(r = paste0("r = ", round(r, 3), "***"))

# Function to plot brightness versus PC1
PLOTFUN <- function(P) {

  P %>%
    ggplot(aes(x = wPC1, y = brightness)) +
    geom_point() +
    labs(x = "PC1", y = "Brightness (%)")

}

# Plot brightness versus PC1 across islands
P %>%
  PLOTFUN() +
  geom_text(
    data = labels,
    mapping = aes(label = r),
    x = 0.7 * min(P$wPC1),
    y = 0.9 * max(P$brightness),
    hjust = 0,
    size = 3
  ) +
  facet_wrap(. ~ island)

figname <- sprintf("%s/brightnesses.png", saveto)
ggsave(figname, width = 4, height = 4, dpi = 300)

P <- P %>%
  group_by(island) %>%
  nest()

# Or each island separately
with(P, walk2(island, data, function(island, P) {

  plot <- P %>%
    PLOTFUN() +
    annotate(
      geom = "text",
      label = paste0("r = ", round(with(P, cor(PC1, brightness)), 3), "***"),
      x = min(P$wPC1),
      y = max(P$brightness),
      hjust = 0
    )

  figname <- sprintf("%s/brightness_%s.png", saveto, island)
  ggsave(figname, plot, width = 3, height = 2, dpi = 300)

}))

#### 13. Multiple site comparisons ####

# This is a posthoc analysis asking between which sites, within each island,
# differences in dewlap color were detected. We perform Wilcoxon tests, only
# for islands and pairs of habitats that were found significant in previous
# steps. We look at the distances at which differences are found.

saveto <- "results/site_comparisons"

# Perform pairwise comparisons between all sites within each island
P <- D %>%
  group_by(island) %>%
  nest() %>%
  mutate(wilcoxon = map(data, function(X) {

    # Group the data by site
    X <- X %>%
      group_by(site, longitude, latitude, habitat) %>%
      nest()

    # For each pair of sites and each variable...
    map_dfr(wpc_names, function(variable) {
      map_dfr(seq(nrow(X) - 1), function(i) {
        map_dfr((i + 1):nrow(X), function(j) {

          # Perform the Wilcoxon test
          test <- with(
            X,
            wilcox.test(data[[i]][[variable]], data[[j]][[variable]])
          )

          # Also mesure the absolute phenotypic distance between the sites
          euclid <- with(X, abs(
            mean(data[[i]][[variable]]) - mean(data[[j]][[variable]])
          ))

          # Extract the relevant information
          tibble(
            variable = variable,
            habitat_i = X$habitat[i],
            habitat_j = X$habitat[j],
            longitude_i = X$longitude[i],
            longitude_j = X$longitude[j],
            latitude_i = X$latitude[i],
            latitude_j = X$latitude[j],
            W = test$statistic,
            pvalue = test$p.value,
            euclid = euclid
          )

        })
      })
    })

  })) %>%
  select(-data) %>%
  unnest(cols = c(wilcoxon)) %>%
  mutate(

    # Add the geographical distances between the sites
    distance = pmap_dbl(
      list(longitude_i, latitude_i, longitude_j, latitude_j),
      function(x1, y1, x2, y2) {

        # Measure distance from coordinates
        distGeo(c(x1, y1), c(x2, y2))

      }
    )
  ) %>%
  mutate(variable = str_remove(variable, "w"))

# Load the results of the post-hoc comparisons
PH <- read_csv("results/group_comparisons/posthoc/summary.csv") %>%
  filter(pvalue < 0.05)

# Significant contrasts
contrasts <- with(PH, paste(island, variable, i, j))

# Restrict this analysis to the relevant comparisons (it is a posthoc test)
P <- P %>%
  filter(paste(island, variable, habitat_i, habitat_j) %in% contrasts) %>%
  mutate(pvalue = p.adjust(pvalue, method = "BH")) %>% # corrected P-values
  filter(pvalue < 0.05)

set.seed(33) # for the jitter

# Function to plot the results of the multiple site comparisons
PLOTFUN <- function(P) {

  P %>%
    ggplot(
      aes(
        x = distance / 1000,
        y = paste(habitat_i, "vs", habitat_j),
        fill = str_remove(variable, "w"),
        shape = str_remove(variable, "w")
      )
    ) +
    scale_shape_manual(values = 21:24) +
    geom_hline(
      mapping = aes(yintercept = paste(habitat_i, "vs", habitat_j)),
      linetype = 4,
      color = "grey"
    ) +
    geom_jitter(height = 0.3, size = 2, alpha = 0.8) +
    labs(
      x = "Distance (km)",
      y = NULL,
      fill = NULL,
      shape = NULL
    ) +
    scale_x_log10()

}

# Now phenotypic distances versus geographic distances
PLOTFUN2 <- function(P) {

  P %>%
    ggplot(
      aes(
        x = distance / 1000,
        y = euclid,
        shape = paste(habitat_i, "vs", habitat_j),
        fill = str_remove(variable, "w")
      )
    ) +
    geom_point(size = 2, alpha = 0.8) +
    scale_shape_manual(values = 21:23) +
    labs(
      x = "Distance (km)",
      y = "Phenotypic distance",
      fill = "Variable",
      shape = "Contrast"
    ) +
    scale_x_log10() +
    guides(
      fill = guide_legend(override.aes = list(shape = 21), order = 1),
      shape = guide_legend(order = 2)
    )

}

# Plot the multiple site comparisons
P %>%
  PLOTFUN() +
  facet_wrap(. ~ island) +
  scale_fill_manual(values = wpc_colors)

figname <- sprintf("%s/site_comparisons.png", saveto)
ggsave(figname, width = 6, height = 4, dpi = 300)

P %>%
  PLOTFUN2() +
  facet_wrap(. ~ island) +
  scale_fill_manual(values = wpc_colors)

figname <- sprintf("%s/site_comparisons2.png", saveto)
ggsave(figname, width = 6, height = 4, dpi = 300)

P <- P %>%
  group_by(island) %>%
  nest()

# Or each island separately
with(P, walk2(island, data, function(island, P) {

  # Pick the colors for the current PCs found significant on this island
  curr_colors <- wpc_colors[unique(P$variable)]

  # Distances at which significant differences are detected
  P %>%
    PLOTFUN() +
    scale_fill_manual(values = curr_colors)

  figname <- sprintf("%s/site_comparisons_%s.png", saveto, island)
  ggsave(figname, width = 4, height = 2, dpi = 300)

  # Phenotypic versus geographical distances
  P %>%
    PLOTFUN2() +
    scale_fill_manual(values = curr_colors)

  figname <- sprintf("%s/site_comparisons2_%s.png", saveto, island)
  ggsave(figname, width = 4, height = 3, dpi = 300)

}))
