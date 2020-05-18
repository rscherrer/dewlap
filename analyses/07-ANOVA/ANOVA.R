# This pipeline performs the following tests
# 1. MANOVA and tests of MANOVA assumptions
# 2. Multiple ANOVAs accounting for heteroskedasticity and posthoc tests
# 3. Plotting of the results
# 4. Identification of deviations from normality
# 5. Repeat tests with Kruskal-Wallis when such deviations occur

rm(list = ls())

library(nmgc)
library(tidyverse)
library(cowplot)
library(knitr)
library(dewlap)

data <- read.csv("data/reflectance.csv", header = TRUE)
variables <- paste0("PC", 1:4)
wl <- paste0("wl", 300:700)

#### 1. Multiple one-way MANOVAs ####

# Are there multivariate differences between habitats on each island?

# 1.1. Multiple parametric MANOVAs

res_manova <- nanova(
  data, variables, grouping = "habitat", nesting = "island", to_pcomp = wl,
  univariate = FALSE, test = "Pillai"
)

res_manova$res # MANOVA table

# It looks like on some islands, yes
# But we know that multivariate normality and homogeneity of covariance matrices
# are violated on multiple islands
read.csv("analyses/05-assumptions/table_multinorm.csv", header = TRUE)
read.csv("analyses/05-assumptions/table_covariance.csv", header = TRUE)

# 1.2. Multiple semi-paramatric MANOVAs

res_smanova <- nanova(
  data, variables, grouping = "habitat", nesting = "island", to_pcomp = wl,
  univariate = FALSE, parametric = FALSE, seed = 42, iter = 1000
)

t1 <- res_smanova$res
save_table(t1, "table_smanova", digits = c(0, 1, 4, 0))

#### 2. Multiple ANOVAs on significant islands ####

# Keep the most significant islands only (based on machine learning)
sislands <- c("Abaco", "Bimini", "Cayman Brac", "Little Cayman", "Long Island")
data <- data %>% filter(island %in% sislands) %>% droplevels

# 2.1. Parametric one-way ANOVAs (with various variance structures)

# Fit ANOVAs
res_anova <- nanova(
  data, variables, grouping = "habitat", nesting = "island", to_pcomp = wl
)

# Note: it would not make sense to consider sites as a random variable
# because only one habitat was measured per site, so sites are completely
# colinear with our grouping variable

t2 <- res_anova$res
save_table(t2, "table_anova", digits = c(0, 0, 0, 0, 1, 1, 3, 0, 1, 2, 4, 0))

t2. <- res_anova$ph
save_table(t2., "table_posthoc", digits = c(0, 0, 0, 0, 0, 2, 4, 0))

# 2.2. Kruskal-Wallis tests wherever normality is not met

# We know that univariate normality is violated in some cases
norm <- read.csv("analyses/05-assumptions/table_normality.csv", header = TRUE)
norm %>% filter(pvalue < 0.05)

# Perform KW tests and record tests of univariate normality (for later)
res_kw <- nanova(
  data, variables, grouping = "habitat", nesting = "island", to_pcomp = wl,
  parametric = FALSE
)

norm2 <- norm %>%
  group_by(island, variable) %>%
  summarize(dev = any(pvalue < 0.05))

t3 <- res_kw$res %>%
  rename(island = "nesting") %>%
  right_join(norm2) %>%
  filter(dev) %>%
  select(-dev) %>%
  drop_na()

save_table(t3, "table_kw", digits = c(0, 0, 2, 0, 4, 0))

t3. <- res_kw$ph %>%
  rename(island = "nesting") %>%
  right_join(norm2) %>%
  filter(dev) %>%
  select(-dev) %>%
  drop_na()

save_table(t3., "table_postkw", digits = c(0, 0, 0, 0, 0, 2, 4, 0))

#### 3. Plot the results ####

# Prepare PCA data for plotting
pca <- npcomp(data, wl, nesting = "island", combine = TRUE, reduce = 1:4)
pcdata <- pca$x
rotation <- pca$rotation

data <- cbind(data, pcdata[, variables])
data <- data %>% gather_("variable", "score", variables)

# 3.1. Plot the distribution of the data

p <- data %>%
  ggplot(aes(x = habitat, y = score)) +
  geom_violin(aes(color = habitat)) +
  geom_jitter(aes(color = habitat), alpha = 0.5, width = 0.2) +
  facet_grid(variable ~ island, scales = "free_y") +
  theme_bw() +
  labs(x = NULL, y = "Principal component score", color = "Habitat") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_color_manual(values = c(coastal = "darkgoldenrod1", coppice = "forestgreen", mangrove = "mediumseagreen"))

# 3.2. Add ANOVA P-values

# Prepare P-value dataset
p_anova <- res_anova$res %>% rename(island = "nesting")

# Replace with KW where applicable
p_kw <- t3 %>% select(island, variable, pvalue)

p_anova <- p_anova %>%
  full_join(p_kw %>% group_by(island, variable) %>% nest()) %>%
  mutate(pvalue = map2_dbl(pvalue, data, ~ if (is.null(.y)) .x else .y$pvalue)) %>%
  select(-data)

p_anova <- p_anova %>% mutate(plabel = round(pvalue, 4) %>% paste("P =", .)) %>%
  mutate(plabel = ifelse(pvalue < 0.0001, "P < 0.0001", plabel)) %>%
  mutate(plabel = ifelse(pvalue < 0.05, str_replace(plabel, "$", "*"), plabel))

p_anova <- p_anova %>% select(island, variable, plabel)

p_anova <- p_anova %>%
  right_join(
    data %>%
      group_by(variable) %>%
      mutate(ymin = min(score), ymax = max(score)) %>%
      ungroup() %>%
      group_by(island, variable) %>%
      summarize(ymin = first(ymin), ymax = first(ymax))
  )

# Add them to the plot
p <- p + geom_text(
  data = p_anova, aes(label = plabel, y = ymin), x = 3.5, hjust = 1, size = 3
)

# 3.3. Add post-hoc tests results

# Prepare post-hoc bar data
p_posthoc <- res_anova$ph
p_posthoc <- p_posthoc %>% rename(island = "nesting")
p_posthoc <- p_posthoc %>% filter(pvalue < 0.05)

# Replace with nonparametric post-hoc test results when applicable
t3.. <- t3. %>%
  filter(pvalue < 0.05) %>%
  group_by(island, variable) %>%
  nest() %>%
  rename(dataKW = "data")

p_posthoc <- p_posthoc %>%
  group_by(island, variable) %>%
  nest() %>%
  full_join(t3..) %>%
  mutate(data = map2(data, dataKW, function(x, y) {
    if (is.null(y)) x else y
  })) %>%
  select(-dataKW) %>%
  unnest(cols = c(data))

# Set horizontal coordinates based on groups being compared
levs <- levels(data$habitat)
step <- 0.1

p_posthoc <- p_posthoc %>%
  mutate(
    xmin = map_int(contrast1, ~ which(.x == levs)) + step,
    xmax = map_int(contrast2, ~ which(.x == levs)) - step
  )

# Append limits of the vertical axis for each facet
p_posthoc <- p_posthoc %>%
  right_join(p_anova) %>%
  drop_na() %>%
  select(-signif, -plabel)

# Set vertical coordinates based on facet limits
p_posthoc <- p_posthoc %>% unite("contrast", contrast1, contrast2, sep = ":")
contrasts <- levels(factor(p_posthoc$contrast))
p_posthoc <- p_posthoc %>% mutate(contrast = map_int(contrast, ~ which(.x == contrasts)))
ylags <- c(0.05, 0, 0.05) # deviations from top of the plot
p_posthoc <- p_posthoc %>% mutate(y = ymax - ylags[contrast] * (ymax - ymin))

# Add hooks to the posthoc bars
add_geom_posthoc_hook <- function(p, data, x, hook = 0.02) {
  p + geom_segment(data = data, aes(
    x = x, y = y, xend = x, yend = y - hook * (ymax - ymin), color = NULL
  ))
}

# Add posthoc bars
add_geom_posthoc <- function(p, data, hook = 0.02) {
  p <- p + geom_segment(data = data, aes(
      x = xmin, y = y, xend = xmax, yend = y, color = NULL
    ))
  p <- p %>% add_geom_posthoc_hook(data, data$xmin, hook = hook)
  p <- p %>% add_geom_posthoc_hook(data, data$xmax, hook = hook)
  p
}

# Display them all
p <- p %>% add_geom_posthoc(p_posthoc)

# 3.4. Add PCA rotation matrices to the figure

# Add rotation matrices
rotplot <- rotation %>%
  gather_("PC", "loading", variables) %>%
  mutate(variable = variable %>% str_replace("wl", "") %>% as.numeric) %>%
  ggplot(aes(x = variable, y = PC, fill = loading)) +
  geom_tile() +
  facet_grid(. ~ island) +
  scale_fill_gradient2(low = "darkblue", mid = "white", high = "darkgreen") +
  labs(y = NULL, x = "Wavelength (nm)", fill = "Loading") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.key.height = unit(0.2, "in"))

plot_grid(rotplot, p, nrow = 2, rel_heights = c(1, 3), labels = c("A", "B"))
ggsave("figure_anova.png", width = 7, height = 7, dpi = 400)

