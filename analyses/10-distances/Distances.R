rm(list = ls())

# Here we explore the distances at which differences between sites are detected

library(dewlap)
library(tidyverse)
library(nmgc)
library(RColorBrewer)

data <- read.csv("data/reflectance.csv", header = TRUE)
wl <- paste0("wl", 300:700)
variables <- paste0("PC", 1:4)
data <- cbind(data, npcomp(data, wl, nesting = "island", reduce = 1:4)$x)
data <- data[, -grep("island", colnames(data))[2]]

# Kepp only the most significant islands
keep <- c("Abaco", "Bimini", "Cayman Brac", "Little Cayman", "Long Island")
data <- data %>% filter(island %in% keep)
data <- data %>% select(island, habitat, site, longitude, latitude, PC1, PC2, PC3, PC4)

# Rearrange the data
data <- data %>% gather_("variable", "score", variables)

# Function to compare sites within an island for a given variable
compare_sites <- function(data) {

  data <- data %>%
    group_by(site, habitat, longitude, latitude) %>%
    nest()

  ii <- 1:(nrow(data) - 1)
  names(ii) <- ii
  out <- map_dfr(ii, function(i) {
    jj <- (i + 1):(nrow(data))
    names(jj) <- jj
    map_dfr(jj, function(j) {
      res <- wilcox.test(data$data[[i]]$score, data$data[[j]]$score)
      data.frame(U = res$statistic, pvalue = res$p.value)
    }, .id = "j")
  }, .id = "i")

  out <- out %>% mutate_at(c("i", "j"), as.numeric)
  out <- out %>% mutate(
    habitat_i = map_chr(i, ~ as.character(data$habitat[.x])),
    habitat_j = map_chr(j, ~ as.character(data$habitat[.x]))
  )
  out <- out %>% filter(habitat_i != habitat_j)
  out <- out %>% mutate(
    lon_i = map_dbl(i, ~ data$longitude[.x]),
    lat_i = map_dbl(i, ~ data$latitude[.x]),
    lon_j = map_dbl(j, ~ data$longitude[.x]),
    lat_j = map_dbl(j, ~ data$latitude[.x])
  )
  out <- out %>% mutate(
    distance = pmap_dbl(list(lon_i, lat_i, lon_j, lat_j), function(x1, y1, x2, y2) {
      distm(x = c(x1, y1), y = c(x2, y2))
    })
  )

}

# Apply the function across islands and variables
res_wilcox <- data %>%
  group_by(island, variable) %>%
  nest() %>%
  mutate(test = map(data, compare_sites)) %>%
  select(-data) %>%
  unnest(cols = c(test)) %>%
  select(-i, -j)

# Read in the results of the post-hoc tests
res_posthoc <- read.csv("analyses/07-ANOVA/table_posthoc.csv", header = TRUE)
res_posthoc <- res_posthoc %>% filter(pvalue < 0.05)
res_posthoc <- res_posthoc %>% rename(island = "nesting")

res_postkw <- read.csv("analyses/07-ANOVA/table_postkw.csv", header = TRUE)
res_postkw <- res_postkw %>% filter(pvalue < 0.05)

# Combine parametric and nonparametric significant contrasts
res_posthoc <- res_posthoc %>% group_by(island, variable) %>% nest()
res_postkw <- res_postkw %>% group_by(island, variable) %>% nest()
res_postkw <- res_postkw %>% rename(data2 = "data")

res_posthoc <- res_posthoc %>%
  full_join(res_postkw) %>%
  mutate(data = map2(data, data2, ~ if (is.null(.y)) .x else .y)) %>%
  select(-data2) %>%
  unnest(cols = c(data)) %>%
  rename(habitat_i = "contrast1", habitat_j = "contrast2")

# Keep the site comparisons for the significant contrasts only
res_posthoc <- res_posthoc %>% group_by(island, variable, habitat_i, habitat_j) %>% nest()
res_wilcox <- res_wilcox %>% group_by(island, variable, habitat_i, habitat_j) %>% nest()
res_wilcox <- res_wilcox %>% rename(data2 = "data")

res_wilcox <- res_wilcox %>%
  right_join(res_posthoc) %>%
  select(-data) %>%
  unnest(cols = c(data2))

# Correct P-values for multiple testing and keep the most significant
res_wilcox <- res_wilcox %>% mutate(pvalue = p.adjust(pvalue, method = "BH"))
res_wilcox <- res_wilcox %>% filter(pvalue < 0.05)

# Plot P-values against distances
p <- ggplot(
  data = res_wilcox,
  aes(x = distance / 1000, y = pvalue, shape = island, fill = variable)
) +
  geom_vline(xintercept = 0.5) +
  geom_vline(xintercept = 1, lty = 2) +
  geom_vline(xintercept = 5, lty = 3) +
  geom_vline(xintercept = 10, lty = 4) +
  geom_point(size = 5) +
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  xlab("Distance (km)") +
  ylab("P-value") +
  scale_fill_manual(values = rev(brewer.pal(length(variables), "YlOrRd"))) +
  scale_shape_manual(values = 21:25) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  facet_grid(. ~ paste(habitat_i, habitat_j, sep = " vs ")) +
  labs(fill = "Variable", shape = "Island") +
  theme(legend.position = "bottom", legend.box = "vertical")
p
ggsave("figure_distances.png", p, width = 7, height = 3.5, dpi = 300)

# Plot significant distances across islands
p2 <- p + aes(y = island) + scale_y_discrete() + ylab(NULL)
p2
ggsave("figure_distances2.png", p2, width = 7, height = 3.5, dpi = 300)

# Assess the correlation between distance and P-value
res_spearman <- with(res_wilcox, cor.test(distance, pvalue, method = "spearman"))
res_spearman <- with(res_spearman, list(rho = round(estimate, 3), pvalue = round(p.value, 4)))

label <- 'rho==%s~" P = %s"'
label <- sprintf(label, res_spearman$rho, res_spearman$pvalue)

# Save table
t1 <- res_wilcox %>% select(-U, -pvalue)
t1 <- t1 %>% bind_cols(res_wilcox %>% ungroup %>% select(U, pvalue))
t1 <- t1 %>% add_signif()
t1_fname <- "analyses/10-distances/table_wilcoxon"
t1_names <- c("Island", "Variable", "Contrast", "", "Lon. 1", "Lat. 1", "Lon. 2", "Lat. 2", "Distance (m)", "$U$", "$P$", "")
save_table(t1, t1_fname, digits = c(0, 0, 0, 0, 4, 4, 4, 4, 0, 0, 4), col.names = t1_names)
