rm(list = ls())

# Here we test for the assumptions of (M)ANOVA across different levels

library(nmgc)
library(dewlap)

data <- read.csv("data/reflectance.csv", header = TRUE)
wl <- paste0("wl", 300:700)
variables <- paste0("PC", 1:4)

#### 1. Pooled data -- across the whole archipelago ####

data <- cbind(data, data.frame(npcomp(data, wl)$x)[, variables])

# 1.1. Multivariate normality
t1 <- test_multinorm(data, variables, grouping = "habitat")[, -1]

# 1.2. Heterogeneity of covariance matrices
t2 <- test_covariance(data, variables, grouping = "habitat")[, -1]

# 1.3. Univariate normality
t3 <- test_multinorm(data, variables, grouping = "habitat", univariate = TRUE)[, -1]

# 1.4. Heterogeneity of variances
t4 <- test_covariance(data, variables, grouping = "habitat", univariate = TRUE)[, -1]

save_table(t1, "table_multinorm_pooled", digits = c(0, 2, 0, 0))
save_table(t2, "table_covariance_pooled", digits = c(1, 0, 4, 0))
save_table(t3, "table_normality_pooled", digits = c(0, 0, 3, 4, 0))
save_table(t4, "table_variance_pooled", digits = c(0, 2, 0, 4))

#### 2. Within islands ####

# 2.1. Multivariate normality
t5 <- test_multinorm(data, variables, grouping = "habitat", nesting = "island")

# 2.2. Homogeneity of covariance matrices
t6 <- test_covariance(data, variables, grouping = "habitat", nesting = "island")

save_table(t5, "table_multinorm", digits = c(0, 0, 2, 4, 0))
save_table(t6, "table_covariance", digits = c(0, 1, 0, 4, 0))

# 2.3. Univariate normality
t7 <- test_multinorm(data, variables, grouping = "habitat", nesting = "island", univariate = TRUE)

# Reduce to the islands with evidence for deviations from multivariate normality
t7 <- t7 %>%
  group_by(island, habitat) %>%
  nest() %>%
  right_join(t5) %>%
  filter(pvalue < 0.05) %>%
  select(island, habitat, data) %>%
  unnest(cols = c(data))

save_table(t7, "table_normality", digits = c(0, 0, 0, 3, 4, 0))

# 2.4. Heterogeneity of variances
t8 <- test_covariance(data, variables, grouping = "habitat", nesting = "island", univariate = TRUE)

# Reduce to the islands with deviations from homogeneity of covariance matrices
t8 <- t8 %>%
  mutate(island = str_replace(island, ".[0-9]$", "")) %>%
  group_by(island) %>%
  nest() %>%
  right_join(t6) %>%
  filter(pvalue < 0.05) %>%
  select(island, data) %>%
  unnest(cols = c(data))

save_table(t8, "table_homoskedasticity", digits = c(0, 0, 3, 0, 4, 0))

# 3. Outliers

test_outliers(data, variables, grouping = "habitat", nesting = "island")

# No outliers detected
