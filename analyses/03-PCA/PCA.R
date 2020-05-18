rm(list = ls())

# Exploration of the principal component analysis

library(nmgc)
library(tidyverse)
library(knitr)

# Here we are checking first the correlation between PC1 and brightness
# and second, the variance explained by each PC on each island

data <- read.csv("data/reflectance.csv", header = TRUE)

# Perform PCA on each island
pca <- npcomp(data, paste0("wl", 300:700), nesting = "island", combine = TRUE, reduce = 1:4)

# And on the whole archipelago
pca_pooled <- npcomp(data, paste0("wl", 300:700))

####  1. Correlation between brightness and PC1 ####

# Calculate brightness (mean reflectance)
data$brightness <- data %>% select_at(paste0("wl", 300:700)) %>% apply(1, mean)

# Attach PC1
data <- data %>% select_at(c("island", "habitat", "brightness"))
data$PC1 <- pca$x$PC1
data$PC1pooled <- pca_pooled$x[, 1]
data <- data %>%
  gather_("variable", "PC1", c("PC1", "PC1pooled")) %>%
  mutate(island = ifelse(variable == "PC1pooled", "Archipelago", as.character(island))) %>%
  dplyr::select(-variable)

# Measure correlation
res_cor <- data %>%
  split(f = .[["island"]]) %>%
  map_dfr(
    ~ cor.test(.x$PC1, .x$brightness) %>%
      .[c("estimate", "p.value")] %>%
      data.frame %>%
      rename(r2 = "estimate", pvalue = "p.value") %>%
      mutate(r2 = r2^2),
    .id = "island"
  )

# Reorder
res_cor <- res_cor %>%
  filter(island != "Archipelago") %>%
  rbind(., res_cor %>% filter(island == "Archipelago"))

# Save table
t1 <- res_cor
t1 <- latex %>%
  add_signif() %>%
  mutate(pvalue = ifelse(pvalue < 0.0001, "< 0.0001", pvalue))
colnames(t1)[colnames(t1) == "signif"] <- ""
save_table(t1, "table_brightness", digits = c(0, 3, 0))


# Plot the correlation
p <- data %>%
  filter(island != "Archipelago") %>%
  ggplot(aes(x = brightness, y = PC1)) +
  geom_point() +
  theme_bw() +
  facet_wrap(. ~ island) +
  xlab("Brightness (mean reflectance, %)")

# Add labels
res_cor$r2 <- round(res_cor$r2, 3)
res_cor <- res_cor %>%
  mutate(plabel = round(pvalue, 4)) %>%
  mutate(plabel = ifelse(pvalue < 0.0001, "P < 0.0001", paste("P =", pvalue))) %>%
  mutate(plabel = ifelse(pvalue < 0.05, str_replace(plabel, "$", "*"), plabel))
p + geom_text(
  data = res_cor %>% filter(island != "Archipelago"),
  aes(label = paste0("R\U000B2 = ", r2, "\n", plabel)), x = 40, y = -25, hjust = 1
)

ggsave("figure_brightness.png", width = 6, height = 6, dpi = 300)

# Now make a figure for the whole archipelago
data %>%
  filter(island == "Archipelago") %>%
  ggplot(aes(x = brightness, y = PC1)) +
  geom_point() +
  theme_bw() +
  xlab("Brightness (mean reflectance, %)") +
  geom_text(
    data = res_cor %>% filter(island == "Archipelago"),
    aes(label = paste0("R\U000B2 = ", r2, "\n", plabel)), x = 40, y = -25, hjust = 1
  )

ggsave("figure_brightness_pooled.png", width = 4, height = 4, dpi = 300)

#### 2. Variance explained by the PCs ####

extra_row <- (pca_pooled$sdev / sum(pca_pooled$sdev))[1:4] %>% as.list
extra_row <- c("Archipelago", extra_row)
names(extra_row) <- colnames(pca$sdev)
expvar <- rbind(pca$sdev, extra_row %>% data.frame())

expvar <- expvar %>% mutate(total = expvar %>% dplyr::select(PC1:PC4) %>% apply(1, sum))

# Save table
t2 <- expvar
save_table(t2, "table_expvar", digits = c(0, 3, 3, 3, 3, 3))
