rm(list = ls())

# Here we explore dewlap coloration between habitats across the whole archipelago

library(nmgc)
library(GGally)
library(MANOVA.RM)
library(heplots)
library(cowplot)

data <- read.csv("data/reflectance.csv", header = TRUE)
wl <- paste0("wl", 300:700)
variables <- paste0("PC", 1:4)

data <- cbind(data, data.frame(npcomp(data, wl)$x)[, variables])

#### 1. Eyeball ####

# We plot the color variables across habitats at the scale of the whole archipelago

colors <- c(coastal = "darkgoldenrod1", coppice = "forestgreen", mangrove = "mediumseagreen")
alpha <- 0.2

lowerplot <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(..., alpha = alpha) +
    stat_ellipse() +
    scale_color_manual(values = colors)
}

diagplot <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_density(..., alpha = alpha) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors)
}

diagplot2 <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_boxplot(..., alpha = alpha) +
    scale_color_manual(values = colors) +
    coord_flip()
}

p <- ggpairs(
  data, columns = variables, mapping = aes(color = habitat),
  upper = NULL, lower = list(continuous = lowerplot), diag = list(continuous = diagplot2)
) +
  theme_bw()

p

# There seems to be no difference between habitats.
# To confirm this we run a two-way MANOVA with islands and habitats as effects

#### 2. 2-way MANOVA ####

# 2.1. Parametric MANOVA

res_manova <- manova(cbind(PC1, PC2, PC3, PC4) ~ island * habitat, data = data)
res_manova %>% summary(test = "Pillai")

# This test suggests that all effects are significant, including an archipelago-wide
# habitat-effect
# How big is this effect?

# 2.2. Effect sizes

etasq(res_manova)

# Very small archipelago-wide effect of habitat

# 2.3. Semi-parametric MANOVA

# We know, however that there are deviations from multivariate normality...
read.csv("analyses/05-assumptions/table_multinorm_pooled.csv", header = TRUE)

# and homogeneity of covariance matrices
read.csv("analyses/05-assumptions/table_covariance_pooled.csv", header = TRUE)

# So we repeat the analysis with a semi-parametric MANOVA

res_manova2 <- MANOVA.wide(cbind(PC1, PC2, PC3, PC4) ~ island * habitat, data, seed = 42, iter = 1000)
res_manova2

# The results hold

#### 3. Univariate ANOVAs ####

# What explains the archipelago-wide differences between habitats?
# We analyze each variable separately
# But first we assess normality and homogeneity of variances

# There is evidence for non-normality for all variables ...
read.csv("analyses/05-assumptions/table_normality_pooled.csv", header = TRUE)

# and no evidence for heterogeneity of variance
read.csv("analyses/05-assumptions/table_variance_pooled.csv", header = TRUE)

# 3.1 Univariate ANOVAs

# We run multiple ANOVAs despite non-normality. We also allow for
# islands to be treated as a random effect.
res_nanova <- nanova(
  data, variables, grouping = "habitat", random = "island", univariate = TRUE,
  parametric = TRUE
)
t1 <- res_nanova$res[, -1]
t1_fname <- "analyses/06-pooled/table_anova_pooled"
t1_names <- c(
  "Variable", "Best fit", "df", "AICc", "$\\Delta$AICc", "AICcw",
  "df$_{\\mbox{LRT}}$", "Log-lik.", "$\\chi^2$", "$P$", ""
)
save_table(t1, t1_fname, digits = c(0, 0, 0, 1, 1, 3, 0, 1, 2, 4, 0), col.names = t1_names)

# Mixed models were always the best fit. They indicate significant differences,
# after accounting for the effect of islands, along PC1, 2 and 3
# But the differences are so small that posthoc tests could not identify what
# groups differ
t2 <- res_nanova$ph[, -1]
t2_fname <- "analyses/06-pooled/table_posthoc_pooled"
t2_names <- c("Variable", "Test", "Contrast", "", "Statistic", "$P$", "")
save_table(t2, t2_fname, digits = c(0, 0, 0, 0, 3, 4, 0), col.names = t2_names)

# 3.2. Kruskal-Wallis tests

# We check with Kruskal-Wallis to account for non-normality, but then we
# cannot account for the effect of islands
res_kw <- nanova(
  data, variables, grouping = "habitat", univariate = TRUE, parametric = FALSE
)
res_kw

t3 <- res_kw$res[, -1]
t3_fname <- "analyses/06-pooled/table_kw_pooled"
t3_names <- c("Variable", "$\\chi^2$", "df", "$P$", "")
save_table(t3, t3_fname, digits = c(0, 2, 0, 4, 0), col.names = t3_names)

t4 <- res_kw$ph[, -1]
t4_fname <- "analyses/06-pooled/table_postkw_pooled"
t4_names <- c("Variable", "Test", "Contrast", "", "Statistic", "$P$", "")
save_table(t4, t4_fname, digits = c(0, 0, 0, 0, 2, 4, 0), col.names = t4_names)

# There seems indeed to be small differences along PC2, and this seems to
# be due to mangroves having higher scores than coastal

#### 4. Full figure with rotation matrix ####

# What wavelengths is that?
# We plot the rotation matrix of the PCA
rotation <- data.frame(npcomp(data, wl)$rotation)
rotation <- rotation %>% rownames_to_column("variable")

p2 <- rotation %>%
  gather_("PC", "loading", variables) %>%
  mutate(variable = variable %>% str_replace("wl", "") %>% as.numeric) %>%
  ggplot(aes(x = variable, y = PC, fill = loading)) +
  geom_tile() +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkgreen") +
  labs(y = NULL, x = "Wavelength (nm)", fill = "Loading") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.key.height = unit(0.2, "in"))

# Combine plots
blank <- ggplot() + theme_void()
p2. <- plot_grid(blank, p2, blank, ncol = 1)
fig <- plot_grid(ggmatrix_gtable(p), p2., labels = c("A", "B"), rel_widths = c(2, 1))

ggsave("figure_pooled.png", fig, width = 8, height = 5, dpi = 400)
