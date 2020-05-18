rm(list = ls())

library(nmgc)
library(tidyverse)
library(cowplot)
library(knitr)

data <- read.csv("data/reflectance.csv", header = TRUE)

#### Analysis within each island ####

# Principal components

variables <- paste0("PC", 1:4)

# SVM-classification on PCs

res <- classify(
  data, variables, grouping = "habitat", nesting = "island",
  to_pcomp = paste0("wl", 300:700), digest = TRUE, test = TRUE, k = 5,
  nrep = 100, seed = 24, method = "SVM", importance = TRUE
)

# Prepare a confusion matrix plot to show on the side

conf <- res$avg$`Long Island`
conf <- conf %>% apply(2, function(x) x / sum(x))
conf <- conf %>%
  data.frame %>%
  rownames_to_column("predicted") %>%
  gather_("true", "freq", colnames(conf))
colnames(conf) <- c("predicted", "true", "freq")

confplot <- ggplot(conf, aes(x = true, y = predicted, fill = freq)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "darkgreen", limits = c(0, 1)) +
  theme_bw() +
  labs(x = "True habitat", y = "Predicted habitat", fill = "Frequency") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

blank <- ggplot() + theme_void()

p <- plot_classif(res, facets = "island")
plot_grid(p, plot_grid(blank, confplot, blank, nrow = 3), ncol = 2, rel_widths = c(2, 1))
ggsave("classif_svm_pca.png", height = 6, width = 9, dpi = 300)

write.csv(res$mean, "table_classif_svm_pca.csv", row.names = FALSE)

# Latex table
latex <- res$mean
head(latex)
latex <- latex %>%
  mutate(signif = ifelse(pvalue < 0.05, "*", "")) %>%
  mutate(signif = ifelse(pvalue < 0.01, "**", signif))
colnames(latex)[colnames(latex) == "signif"] <- ""
latex <- kable(latex, digits = c(0, 3, 0, 1, 0, 4), format = "latex")
texfile <- file("table_classif_svm_pca.tex")
writeLines(latex, texfile)
close(texfile)

res$imp %>%
  gather_("variable", "importance", variables) %>%
  filter(island %in% c("Abaco", "Bimini", "Cayman Brac", "Little Cayman", "Long Island")) %>%
  ggplot(aes(x = variable, y = importance)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = NULL, y = "Relative importance") +
  facet_wrap(. ~ island) +
  geom_hline(yintercept = 0.25, lty = 2)

ggsave("importance_svm_pca.png", height = 4, width = 5, dpi = 300)

# LDA-classification on PCs

res <- classify(
  data, variables, grouping = "habitat", nesting = "island",
  to_pcomp = paste0("wl", 300:700), digest = TRUE, test = TRUE, k = 5,
  nrep = 100, seed = 24, method = "LDA", importance = TRUE
)

p <- plot_classif(res, facets = "island", fill = "coral")
plot_grid(p, plot_grid(blank, confplot, blank, nrow = 3), ncol = 2, rel_widths = c(2, 1))
ggsave("classif_lda_pca.png", height = 6, width = 9, dpi = 300)

write.csv(res$mean, "table_classif_lda_pca.csv", row.names = FALSE)

# Latex table
latex <- res$mean
head(latex)
latex <- latex %>%
  add_signif() %>%
  mutate(pvalue = ifelse(pvalue < 0.0001, "< 0.0001", pvalue))
colnames(latex)[colnames(latex) == "signif"] <- ""
latex <- kable(latex, digits = c(0, 3, 0, 1, 0, 4), format = "latex")
texfile <- file("table_classif_lda_pca.tex")
writeLines(latex, texfile)
close(texfile)

res$imp %>%
  gather_("variable", "importance", variables) %>%
  filter(island %in% c("Abaco", "Bimini", "Cayman Brac", "Little Cayman", "Long Island")) %>%
  ggplot(aes(x = variable, y = importance)) +
  geom_boxplot() +
  theme_bw() +
  labs(x = NULL, y = "Relative importance") +
  facet_wrap(. ~ island) +
  geom_hline(yintercept = 0.25, lty = 2)

ggsave("importance_lda_pca.png", height = 4, width = 5, dpi = 300)

# Reflectance data

# SVM-classification on reflectance

variables <- paste0("wl", seq(300, 700, 50))

res <- classify(
  data, variables, grouping = "habitat", nesting = "island", digest = TRUE, test = TRUE, k = 5,
  nrep = 100, seed = 24, method = "SVM", importance = TRUE
)

p <- plot_classif(res, facets = "island")
plot_grid(p, plot_grid(blank, confplot, blank, nrow = 3), ncol = 2, rel_widths = c(2, 1))
ggsave("classif_svm_refl.png", height = 6, width = 9, dpi = 300)

write.csv(res$mean, "table_classif_svm_refl.csv", row.names = FALSE)

# Latex table
latex <- res$mean
head(latex)
latex <- latex %>%
  add_signif() %>%
  mutate(pvalue = ifelse(pvalue < 0.0001, "< 0.0001", pvalue))
colnames(latex)[colnames(latex) == "signif"] <- ""
latex <- kable(latex, digits = c(0, 3, 0, 1, 0, 4), format = "latex")
texfile <- file("table_classif_svm_refl.tex")
writeLines(latex, texfile)
close(texfile)


res$imp %>%
  gather_("variable", "importance", variables) %>%
  mutate(variable = variable %>% str_replace("wl", "")) %>%
  filter(island %in% c("Abaco", "Bimini", "Cayman Brac", "Little Cayman", "Long Island")) %>%
  ggplot(aes(x = variable, y = importance)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Wavelength (nm)", y = "Relative importance") +
  facet_wrap(. ~ island) +
  geom_hline(yintercept = 1/length(variables), lty = 2)

ggsave("importance_svm_refl.png", height = 4, width = 5, dpi = 300)


# LDA-classification on reflectance

res <- classify(
  data, variables, grouping = "habitat", nesting = "island", digest = TRUE, test = TRUE, k = 5,
  nrep = 100, seed = 24, method = "LDA", importance = TRUE
)

p <- plot_classif(res, facets = "island", fill = "coral")
plot_grid(p, plot_grid(blank, confplot, blank, nrow = 3), ncol = 2, rel_widths = c(2, 1))
ggsave("classif_lda_refl.png", height = 6, width = 9, dpi = 300)

write.csv(res$mean, "table_classif_lda_refl.csv", row.names = FALSE)

# Latex table
latex <- res$mean
head(latex)
latex <- latex %>%
  add_signif() %>%
  mutate(pvalue = ifelse(pvalue < 0.0001, "< 0.0001", pvalue))
colnames(latex)[colnames(latex) == "signif"] <- ""
latex <- kable(latex, digits = c(0, 3, 0, 1, 0, 4), format = "latex")
texfile <- file("table_classif_lda_refl.tex")
writeLines(latex, texfile)
close(texfile)

res$imp %>%
  gather_("variable", "importance", variables) %>%
  mutate(variable = variable %>% str_replace("wl", "")) %>%
  filter(island %in% c("Abaco", "Bimini", "Cayman Brac", "Little Cayman", "Long Island")) %>%
  ggplot(aes(x = variable, y = importance)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Wavelength (nm)", y = "Relative importance") +
  facet_wrap(. ~ island) +
  geom_hline(yintercept = 1/length(variables), lty = 2)

ggsave("importance_lda_refl.png", height = 4, width = 5, dpi = 300)


#### Pooling the islands together ####

# Principal components

variables <- paste0("PC", 1:4)

# SVM-classification of PCs

res <- classify(
  data, variables, grouping = "habitat",
  to_pcomp = paste0("wl", 300:400), digest = TRUE, test = TRUE, k = 5,
  nrep = 100, seed = 24, method = "SVM", importance = TRUE
)

p <- plot_classif(res, bins = 50, ylim = c(0, 100), dfac = 1000, py = 90, add_insets = FALSE)
confplot <- plot_classif(res, type = "confusion") +
  labs(x = "True habitat", y = "Predicted habitat", fill = "Frequency") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
plot_grid(p, plot_grid(blank, confplot, blank, nrow = 3, rel_heights = c(1, 3, 1)), ncol = 2, rel_widths = c(1.1, 1))
ggsave("classif_svm_pca_pooled.png", height = 3, width = 6, dpi = 300)

write.csv(res$mean, "table_classif_svm_pca_pooled.csv", row.names = FALSE)

res$imp %>%
  gather_("variable", "importance", variables) %>%
  ggplot(aes(x = variable, y = importance)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Wavelength (nm)", y = "Relative importance") +
  geom_hline(yintercept = 1/length(variables), lty = 2)

ggsave("importance_svm_pca_pooled.png", height = 3, width = 3, dpi = 300)


# LDA-classification on PCs

res <- classify(
  data, variables, grouping = "habitat",
  to_pcomp = paste0("wl", 300:400), digest = TRUE, test = TRUE, k = 5,
  nrep = 100, seed = 24, method = "LDA", importance = TRUE
)

p <- plot_classif(res, bins = 50, ylim = c(0, 100), dfac = 1000, py = 90, add_insets = FALSE, fill = "coral")
confplot <- plot_classif(res, type = "confusion") +
  labs(x = "True habitat", y = "Predicted habitat", fill = "Frequency") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
plot_grid(p, plot_grid(blank, confplot, blank, nrow = 3, rel_heights = c(1, 3, 1)), ncol = 2, rel_widths = c(1.1, 1))
ggsave("classif_lda_pca_pooled.png", height = 3, width = 6, dpi = 300)

write.csv(res$mean, "table_classif_lda_pca_pooled.csv", row.names = FALSE)

res$imp %>%
  gather_("variable", "importance", variables) %>%
  ggplot(aes(x = variable, y = importance)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Wavelength (nm)", y = "Relative importance") +
  geom_hline(yintercept = 1/length(variables), lty = 2)

ggsave("importance_lda_pca_pooled.png", height = 3, width = 3, dpi = 300)


# Reflectance data

variables <- paste0("wl", seq(300, 700, 50))

# SVM-classification on reflectance

res <- classify(
  data, variables, grouping = "habitat", digest = TRUE, test = TRUE, k = 5,
  nrep = 100, seed = 24, method = "SVM", importance = TRUE
)

p <- plot_classif(res, bins = 50, ylim = c(0, 100), dfac = 1000, py = 90, add_insets = FALSE)
confplot <- plot_classif(res, type = "confusion") +
  labs(x = "True habitat", y = "Predicted habitat", fill = "Frequency") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
plot_grid(p, plot_grid(blank, confplot, blank, nrow = 3, rel_heights = c(1, 3, 1)), ncol = 2, rel_widths = c(1.1, 1))
ggsave("classif_svm_refl_pooled.png", height = 3, width = 6, dpi = 300)

write.csv(res$mean, "table_classif_svm_refl_pooled.csv", row.names = FALSE)

res$imp %>%
  gather_("variable", "importance", variables) %>%
  mutate(variable = variable %>% str_replace("wl", "")) %>%
  ggplot(aes(x = variable, y = importance)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Wavelength (nm)", y = "Relative importance") +
  geom_hline(yintercept = 1/length(variables), lty = 2)

ggsave("importance_svm_refl_pooled.png", height = 4, width = 5, dpi = 300)


# LDA-classification on reflectance

res <- classify(
  data, variables, grouping = "habitat", digest = TRUE, test = TRUE, k = 5,
  nrep = 100, seed = 24, method = "LDA", importance = TRUE
)

p <- plot_classif(res, bins = 50, ylim = c(0, 100), dfac = 1000, py = 90, add_insets = FALSE, fill = "coral")
confplot <- plot_classif(res, type = "confusion") +
  labs(x = "True habitat", y = "Predicted habitat", fill = "Frequency") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
plot_grid(p, plot_grid(blank, confplot, blank, nrow = 3, rel_heights = c(1, 3, 1)), ncol = 2, rel_widths = c(1.1, 1))
ggsave("classif_lda_refl_pooled.png", height = 3, width = 6, dpi = 300)

write.csv(res$mean, "table_classif_lda_refl_pooled.csv", row.names = FALSE)

res$imp %>%
  gather_("variable", "importance", variables) %>%
  mutate(variable = variable %>% str_replace("wl", "")) %>%
  ggplot(aes(x = variable, y = importance)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(x = "Wavelength (nm)", y = "Relative importance") +
  geom_hline(yintercept = 1/length(variables), lty = 2)

ggsave("importance_lda_refl_pooled.png", height = 4, width = 5, dpi = 300)
