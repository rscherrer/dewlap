rm(list = ls())

# Analysis of parallel divergence across islands

library(tidyverse)
library(nmgc)
library(RColorBrewer)

# Are differences in coloration similar on different islands?
# To test this, perform a PCA on the whole archipelago
# Then, measure and plot the pairwise differences between the habitats
# in color space, along each PC
# Re-scale these contrasts if needed, because some islands may be more
# variable than others
# If islands undergo the same divergence process, then they should cluster
# together in this contrast space

data <- read.csv("data/reflectance.csv", header = TRUE)
variables <- paste0("PC", 1:4)
wl <- paste0("wl", 300:700)

# Perform PCA at the scale of the whole archipelago
data <- data <- cbind(data, data.frame(npcomp(data, wl)$x)[, variables])

# Compute contrast scores across islands and variables
newdata <- data %>%
  group_by(island, habitat) %>%
  summarize_at(variables, mean) %>%
  gather_("variable", "score", variables) %>%
  spread_("habitat", "score") %>%
  mutate(
    coa_cop = coastal - coppice,
    coa_man = coastal - mangrove,
    cop_man = coppice - mangrove
  ) %>%
  ungroup() %>%
  group_by(variable) %>%
  mutate_at(
    c("coa_cop", "coa_man", "cop_man"),
    function(x) scale(x, center = FALSE, scale = TRUE)
  ) %>%
  filter(island %in% c("Abaco", "Bimini", "Cayman Brac", "Little Cayman", "Long Island")) %>%
  ungroup()

newdata <- newdata %>% mutate(
  x = coa_cop,
  y = cop_man
)

newdata <- newdata %>% rename(vbl = "variable")

p <- ggplot(newdata, aes(x = x, y = y, fill = vbl, shape = island)) +
  geom_point(size = 5) +
  xlab("Coastal versus coppice") +
  ylab("Coppice versus mangrove") +
  theme_bw() +
  labs(shape = "Island", fill = "Variable") +
  scale_fill_manual(values = rev(brewer.pal(length(variables), "YlOrRd"))) +
  scale_shape_manual(values = 21:25) +
  guides(fill = guide_legend(override.aes = list(shape = 21)))

p

# A MANOVA confirms that variables do not behave in similar ways among islands
X <- as.matrix(newdata[, c("x", "y")])
res_manova <- manova(X ~ vbl, data = newdata)
res_manova <- res_manova %>% summary(test = "Pillai")

# No assumption violated
test_multinorm(newdata, variables = c("x", "y"), grouping = "vbl")
test_covariance(newdata, variables = c("x", "y"), grouping = "vbl")

# Add MANOVA results to the plot
to_get <- c("Pillai", "num Df", "den Df", "approx F", "Pr(>F)")
labs <- map(to_get, ~ res_manova$stats[1, .x])
rounds <- c(3, 2, 0, 0, 4)
labs <- map2(labs, rounds, ~ round(.x, .y))
label <- "MANOVA, Pillai = %s, F(%s, %s) = %s, P = %s"
label <- sprintf(label, labs[[1]], labs[[2]], labs[[3]], labs[[4]], labs[[5]])

p <- p + ggtitle(NULL, label)
p

ggsave("figure_contrasts.png", p, width = 5, height = 4, dpi = 300)
