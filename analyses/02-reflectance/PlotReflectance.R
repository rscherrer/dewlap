# Plot reflectance spectra

rm(list = ls())

library(tidyverse)

colors <- c(coastal = "darkgoldenrod1", coppice = "forestgreen", mangrove = "mediumseagreen")

# Read the data
data <- read.csv("data/reflectance.csv", header = TRUE)

# Reorganize into a long table
smr <- data %>%
  gather_("wl", "reflectance", paste0("wl", 300:700)) %>%
  mutate(wl = as.numeric(gsub("wl", "", wl)))

# Prepare ribbons
smr <- smr %>% group_by(wl, habitat, island)
smr <- smr %>%
  summarize(
    mid = quantile(reflectance, 0.05),
    lower = median(reflectance),
    upper = quantile(reflectance, 0.95)
  )

# Plot the spectra
p <- ggplot(smr, aes(x = wl, color = habitat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = habitat), alpha = 0.2) +
  #geom_line(aes(y = mid)) +
  facet_wrap(. ~ island) +
  xlab("Wavelength (nm)") +
  ylab("Reflectance (%)") +
  theme_bw() +
  labs(fill = "Habitat", color = "Habitat") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors)
p

ggsave("figure_reflectance.png", height = 6, width = 6, dpi = 300)

#####################

# At the scale of the whole archipelago now

# Reorganize into a long table
smr <- data %>%
  gather_("wl", "reflectance", paste0("wl", 300:700)) %>%
  mutate(wl = as.numeric(gsub("wl", "", wl)))

# Prepare ribbons
smr <- smr %>% group_by(wl, habitat)
smr <- smr %>%
  summarize(
    mid = quantile(reflectance, 0.05),
    lower = median(reflectance),
    upper = quantile(reflectance, 0.95)
  )

# Plot the spectra
p <- ggplot(smr, aes(x = wl, color = habitat)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = habitat), alpha = 0.2) +
  #geom_line(aes(y = mid)) +
  xlab("Wavelength (nm)") +
  ylab("Reflectance (%)") +
  theme_bw() +
  labs(fill = "Habitat", color = "Habitat") +
  scale_fill_manual(values = colors) +
  scale_color_manual(values = colors)
p

ggsave("figure_reflectance_pooled.png", height = 2, width = 4, dpi = 300)
