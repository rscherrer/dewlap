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

t1_names <- c("Habitat", "$HZ$", "$P$", "")
t2_names <- c("$\\chi^2$", "df", "$P$", "")
t3_names <- c("Habitat", "Variable", "$W$", "$P$", "")
t4_names <- c("Variable", "$K^2$", "df", "$P$", "")

fname <- "analyses/05-assumptions/table_%s"
t1_fname <- sprintf(fname, "multinorm_pooled")
t2_fname <- sprintf(fname, "covariance_pooled")
t3_fname <- sprintf(fname, "normality_pooled")
t4_fname <- sprintf(fname, "variance_pooled")

save_table(t1, t1_fname, digits = c(0, 2, 0, 0), col.names = t1_names)
save_table(t2, t2_fname, digits = c(1, 0, 4, 0), col.names = t2_names)
save_table(t3, t3_fname, digits = c(0, 0, 3, 4, 0), col.names = t3_names)
save_table(t4, t4_fname, digits = c(0, 2, 0, 4), col.names = t4_names)

#### 2. Within islands ####

# 2.1. Multivariate normality
t5 <- test_multinorm(data, variables, grouping = "habitat", nesting = "island")

# 2.2. Homogeneity of covariance matrices
t6 <- test_covariance(data, variables, grouping = "habitat", nesting = "island")

t5_names <- c("Island", "Habitat", "$HZ$", "$P$", "")
t6_names <- c("Island", "$\\chi^2$", "df", "$P$", "")

t5_fname <- sprintf(fname, "multinorm")
t6_fname <- sprintf(fname, "covariance")

save_table(t5, t5_fname, digits = c(0, 0, 2, 4, 0), col.names = t5_names)
save_table(t6, t6_fname, digits = c(0, 1, 0, 4, 0), col.names = t6_names)

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
