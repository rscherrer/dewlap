## ---------------------------
##
## Script name: 01_process_data.R
##
## Purpose of script: This script extracts, curates, processes and assembles
## the raw data (reflectance profiles and associated metadata) into a dataset
## ready for analysis with statistical tools. Among other things, some outliers
## are removed and the reflectance curves are smoothed using a procedure that
## removes some very localized artifacts we have identified, without denaturing
## the signal. Plots are also made to show how processing affects the reflectance
## curves. Follow the comments for details about each step.
##
## Author: Raphael Scherrer
##
## Date Created: 2021-02-06
##
## This script comes with no guarantee whatsoever.
##
## Find me on GitHub at https://github.com/rscherrer
##
## Email:
## r.scherrer@rug.nl
## raphael.scherrer@evobio.eu
## raph.rjfs@hotmail.fr
##
## ---------------------------

rm(list = ls())

library(tidyverse)
library(pavo)
library(patchwork)

# Read metadata
meta <- read_csv("data/master_database.csv")

# Read spectral data
specs <- map(list.dirs("data/raw")[-1], getspec)

# Read specimen tags for each spectrum
tags <- map(specs, ~ colnames(.x)[-1])

# Remove all double underscores in tags
fieldtags <- map(tags, str_replace, "__", "_")

# Extract their field tags (slightly different for each dataset)
fieldtags <- list(
  str_remove(fieldtags[[1]], "_Refl.*$"),
  str_remove(str_remove(fieldtags[[2]], "_Refl.*$"), "CB_"),
  str_remove(str_remove(fieldtags[[3]], "_Refl.*$"), "LC_"),
  str_remove(fieldtags[[4]], "_Refl.*$"),
  str_remove(fieldtags[[5]], "_Refl.*$"),
  str_remove(fieldtags[[6]], "_Refl.*$"),
  str_remove(str_remove(str_remove(fieldtags[[7]], "_Refl.*$"), "_Relative.*$"), "JASIL_")
)

# Isolate the metadata of interest for each dataset (slightly different for each)
metas <- list(
  meta %>% filter(FieldTag %in% unique(fieldtags[[1]]), Island == "Abaco"),
  meta %>% filter(FieldTag %in% unique(fieldtags[[2]]), Island == "CaymanBrac"),
  meta %>% filter(FieldTag %in% unique(fieldtags[[3]]), Island == "LittleCayman"),
  meta %>% filter(FieldTag %in% unique(fieldtags[[4]]), Island == "LongIsland"),
  meta %>% filter(FieldTag %in% unique(fieldtags[[5]]), Island == "RaggedIsland"),
  meta %>% filter(FieldTag %in% unique(fieldtags[[6]]), Island == "SouthBimini"),
  meta %>% filter(FieldTag %in% unique(fieldtags[[7]]), Island %in% c("Eleuthera", "NorthAndros", "SouthAndros"))
)

# For each dataset...
data <- pmap_dfr(list(specs, metas, fieldtags, tags), function(specs, meta, fieldtags, tags) {

  # Keep only the columns of interest
  meta <- meta %>% select(RGR, FieldTag, Island, Habitat, Locale, lat, long, Sex, notes)

  # Transpose the spectral data
  colnames <- paste0("wl", specs$wl)
  specs <- as_tibble(as.data.frame(t(specs))[-1,])
  colnames(specs) <- colnames

  # Combine data and metadata
  specs %>%
    mutate(tag = tags, FieldTag = fieldtags) %>%
    left_join(meta)

}, .id = "dataset")

# Add spot number (= where on the dewlap reflectance was taken)
data <- data %>% mutate(spot = as.numeric(str_remove(tag, "^.*_")))

# Rename columns
data <- data %>%
  rename(
    code = "RGR", fieldtag = "FieldTag", island = "Island", habitat = "Habitat",
    locale = "Locale", latitude = "lat", longitude = "long", sex = "Sex"
  )

# Add a global individual identifier
data <- data %>% mutate(specimen = paste(island, fieldtag, sep = "_"))

# Rename islands
data <- data %>%
  mutate(
    island = fct_recode(
      island,
      !!!c(
        "Cayman Brac" = "CaymanBrac", "Little Cayman" = "LittleCayman",
        "Long Island" = "LongIsland", "Ragged Island" = "RaggedIsland",
        "Bimini" = "SouthBimini", "South Andros" = "SouthAndros",
        "North Andros" = "NorthAndros"
      )
    )
  )

# Change the order of the islands
data <- data %>%
  mutate(island = factor(island, levels = c("Abaco", "Bimini", "Cayman Brac", "Eleuthera", "Little Cayman", "Long Island", "North Andros", "Ragged Island", "South Andros")))

# Change manually some small males that should be counted as adults (their reflectance could be counted as dewlap center)
data$spot[data$code %in% c("RGR1498", "RGR1499", "RGR1500")] <- 2
data$sex[data$code %in% c("RGR1498", "RGR1499", "RGR1500", "RGR1501")] <- "M"
data$sex[data$code == "RGR1186"] <- "M"
data$spot[data$code == "RGR1186"] <- 2
data$spot[data$tag == "JASIL_3241_Reflection_3"] <- 2
data$spot[data$tag == "2937843_Reflection_3"] <- 2
data$spot[data$tag == "JASIL_3104_Reflection_3"] <- 2
data$spot[data$tag %in% c("JASIL_3106_Reflection_1", "JASIL_3107_Reflection_1", "JASIL_3109_Reflection_1")] <- 2
data$spot[data$tag %in% c("JASIL_3132_Reflection_1", "JASIL_3142_Reflection_1", "JASIL_3144_Reflection_1", "JASIL_3149_Reflection_1")] <- 2
data$spot[data$tag %in% c("JASIL_3157_Reflection_1", "JASIL_3188_Reflection_01", "JASIL_3189_Reflection_01", "JASIL_3193_Reflection_01", "JASIL_3194_Reflection_01")] <- 2
data$spot[data$tag == "2938787__Reflection_03"] <- 2

# Filter
data <- data %>% filter(!is.na(island))
data <- data %>% filter(spot == 2) # 1 = throat, 2 = center, 3 = edge
data <- data %>% filter(habitat %in% c("mangrove", "coppice", "coastal"))
data <- data %>% filter(sex == "M")

theme_set(theme_classic())

# Plot the raw data
plots <- data %>%
  group_by(island) %>%
  nest() %>%
  arrange(island) %>%
  mutate(plot = map2(island, data, function(island, data) {

    # Plot to make for each island
    data %>%
      mutate(
        latitude_lab = paste("lat =", round(as.numeric(latitude), 3)),
        longitude_lab = paste("lon =", round(as.numeric(longitude), 3))
      ) %>%
      pivot_longer(wl300:wl700, names_to = "wl") %>%
      mutate(wl = as.numeric(str_remove(wl, "wl"))) %>%
      ggplot(aes(x = wl, y = value, group = specimen, color = habitat)) +
      geom_line(alpha = 0.6) +
      facet_wrap(~ latitude_lab + longitude_lab, nrow = 1) +
      xlab("Wavelength (nm)") +
      ylab("Reflectance (%)") +
      scale_color_manual(values = c("goldenrod", "forestgreen", "mediumseagreen")) +
      labs(color = "Habitat") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      ggtitle(island)

  })) %>%
  ungroup() %>%
  select(plot) %>%
  c() %>%
  first()

# Combine the plots
P <- wrap_plots(
  wrap_plots(
    plots[[1]],
    wrap_plots(plots[[2]], plot_spacer(), nrow = 1, widths = c(3, 4)),
    wrap_plots(plots[[3]], plot_spacer(), nrow = 1, widths = c(3, 4)),
    wrap_plots(plots[[4]], plot_spacer(), nrow = 1, widths = c(5, 2)),
    wrap_plots(plots[[5]], plot_spacer(), nrow = 1, widths = c(3, 4)),
    ncol = 1
  ),
  wrap_plots(
    plots[[6]],
    wrap_plots(plots[[7]], plot_spacer(), nrow = 1, widths = c(3, 1)),
    plots[[8]],
    wrap_plots(plots[[9]], plot_spacer(), nrow = 1, widths = c(3, 1)),
    plot_spacer(),
    ncol = 1
  ),
  nrow = 1, widths = c(7, 4), guides = "collect"
)

# Save
ggsave("results/processing/raw_data.png", P, height = 10, width = 16, dpi = 300)

# Convert the data to pavo's rspec format
data_rspec <- as.data.frame(t(select(data, wl300:wl700)))
data_rspec$wl <- 300:700
data_rspec <- as.rspec(data_rspec)

# Function to have a look at the data during processing
PLOTFUN <- function(data) {

  data %>%
    pivot_longer(colnames(data)[-1]) %>%
    ggplot(aes(x = wl, y = value, group = name)) +
    geom_line(alpha = 0.2) +
    xlab("Wavelength (nm)") +
    ylab("Reflectance (%)")

}

# Plot the data
p1 <- data_rspec %>% PLOTFUN()

# Function to correct artifacts by linearization
CORRECFUN <- function(data, from, to) {

  data %>%
    as_tibble() %>%
    pivot_longer(colnames(data)[-1]) %>%
    mutate(todo = wl >= from & wl <= to) %>%
    group_by(todo, name) %>%
    nest() %>%
    mutate(data = map2(todo, data, function(todo, data) {

      # Flatten the portion that has the artifact
      if (todo) data$value <- loess.smooth(data$wl, data$value, evaluation = nrow(data), span = 50, degree = 1, family = "gaussian")$y
      return(data)

    })) %>%
    unnest(data) %>%
    ungroup() %>%
    select(-todo) %>%
    pivot_wider(names_from = "name", values_from = "value") %>%
    as.rspec(whichwl = 1)

}

# Correct identified artifacts (spurious spikes in reflectance)
data_rspec <- data_rspec %>% CORRECFUN(from = 540, to = 560)
p2 <- data_rspec %>% PLOTFUN()
data_rspec <- data_rspec %>% CORRECFUN(from = 539, to = 561)
p3 <- data_rspec %>% PLOTFUN()
data_rspec <- data_rspec %>% CORRECFUN(from = 600, to = 630)
p4 <- data_rspec %>% PLOTFUN()
data_rspec <- data_rspec %>% CORRECFUN(from = 599, to = 631)
p5 <- data_rspec %>% PLOTFUN()

# Clamp negative values to zero
data_rspec <- procspec(data_rspec, fixneg = "zero")
p6 <- data_rspec %>% PLOTFUN()

# Smooth
data_rspec <- procspec(data_rspec, opt = "smooth", span = 0.5)
p7 <- data_rspec %>% PLOTFUN()

# Clamp again (smoothing created very slightly negative values)
data_rspec <- procspec(data_rspec, fixneg = "zero")
p8 <- data_rspec %>% PLOTFUN()

# Assemble the plots
P <- wrap_plots(
  p1 + ggtitle("Raw data"),
  p2 + ggtitle("Artifact correction phase 1"),
  p3 + ggtitle("Artifact correction phase 2"),
  p4 + ggtitle("Artifact correction phase 3"),
  p5 + ggtitle("Artifact correction phase 4"),
  p6 + ggtitle("Zero truncation"),
  p7 + ggtitle("Smoothing"),
  p8 + ggtitle("Zero re-truncation"),
  ncol = 2
)

# Save
ggsave("results/processing/processing.png", P, width = 7, height = 10, dpi = 300)

# Replace original with processed data
data[, paste0("wl", 300:700)] <- t(data_rspec[, -1])

# Remove outliers
outliers <- c("RGR1103", "RGR0012", "RGR0446", "RGR1186")
data <- data %>% filter(!code %in% outliers)

# Plot after outliers are removed
p9 <- data %>%
  pivot_longer(wl300:wl700, names_to = "wl") %>%
  mutate(wl = as.numeric(str_remove(wl, "wl"))) %>%
  ggplot(aes(x = wl, y = value, group = specimen)) +
  geom_line(alpha = 0.2) +
  xlab("Wavelength (nm)") +
  ylab("Reflectance (%)")

# Assemble the plots
P <- wrap_plots(
  p8 + ggtitle("With outliers"),
  p9 + ggtitle("Without") + ylim(range(p8$data$value)),
  nrow = 1
)

# Save
ggsave("results/processing/outliers.png", P, width = 7, height = 3, dpi = 300)

# Rearrange the rows
data <- data %>% arrange(habitat, island)

# Plot the processed data
plots <- data %>%
  group_by(island) %>%
  nest() %>%
  arrange(island) %>%
  mutate(plot = map2(island, data, function(island, data) {

    # Plot for each island
    data %>%
      mutate(
        latitude_lab = paste("lat =", round(as.numeric(latitude), 3)),
        longitude_lab = paste("lon =", round(as.numeric(longitude), 3))
      ) %>%
      pivot_longer(wl300:wl700, names_to = "wl") %>%
      mutate(wl = as.numeric(str_remove(wl, "wl"))) %>%
      ggplot(aes(x = wl, y = value, group = specimen, color = habitat)) +
      geom_line(alpha = 0.6) +
      facet_wrap(~ latitude_lab + longitude_lab, nrow = 1) +
      xlab("Wavelength (nm)") +
      ylab("Reflectance (%)") +
      scale_color_manual(values = c("goldenrod", "forestgreen", "mediumseagreen")) +
      labs(color = "Habitat") +
      theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
      ggtitle(island)

  })) %>%
  ungroup() %>%
  select(plot) %>%
  c() %>%
  first()

# Combine the plots
P <- wrap_plots(
  wrap_plots(
    plots[[1]],
    wrap_plots(plots[[2]], plot_spacer(), nrow = 1, widths = c(3, 4)),
    wrap_plots(plots[[3]], plot_spacer(), nrow = 1, widths = c(3, 4)),
    wrap_plots(plots[[4]], plot_spacer(), nrow = 1, widths = c(5, 2)),
    wrap_plots(plots[[5]], plot_spacer(), nrow = 1, widths = c(3, 4)),
    ncol = 1
  ),
  wrap_plots(
    plots[[6]],
    wrap_plots(plots[[7]], plot_spacer(), nrow = 1, widths = c(3, 1)),
    plots[[8]],
    wrap_plots(plots[[9]], plot_spacer(), nrow = 1, widths = c(3, 1)),
    plot_spacer(),
    ncol = 1
  ),
  nrow = 1, widths = c(7, 4), guides = "collect"
)

# Save
ggsave("results/processing/processed_data.png", P, height = 10, width = 16, dpi = 300)

# Save the data
write_csv(data, "data/reflectance.csv")
