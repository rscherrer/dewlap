# Make maps of the archipelago

rm(list = ls())

library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(tmap)
library(tmaptools)

# Load geographical data
world <- ne_countries(scale = "large", returnclass = "sf")

# Filter the Caribbean only
island_map_data <- world %>% filter(subregion == "Caribbean")

# Read the bounding boxes of all islands
island_coordinates <- read_csv("data/island_coordinates.csv")

# Make them a nested list and rremove islands we do not want
island_coordinates <- island_coordinates %>%
  group_by(island) %>%
  filter(!island %in% c("Conception Island", "Rum Cay")) %>%
  nest

islands <- island_coordinates$island

# Intersect the Caribbean map with each ot the islands' bounding boxes to get one map per island
island_map_data <- map(island_coordinates$data, function(curr_island_data) {

  extent <- st_bbox(c(
    xmin = curr_island_data$lonmin,
    xmax = curr_island_data$lonmax,
    ymin = curr_island_data$latmin,
    ymax = curr_island_data$latmax
  ), crs = 4326)

  extent <- st_as_sfc(extent)

  return(st_intersection(island_map_data, extent))

})

island_map_data <- do.call(rbind, island_map_data)

# Load sampling sites
data <- readRDS("data/locations.rds")

# Add island names for facetting
island_map_data$island <- levels(data$island)

# Map the whole archipelago
p <- tm_shape(world, bbox = bb(matrix(c(-85, 18, -70, 30), 2, 2))) +
  tm_style("classic") +
  tm_polygons(col = "palegoldenrod", border.col = "gray50") +
  tm_grid(alpha = 0.5, n.x = 5, labels.size = 1) +
  tm_scale_bar(position = c(0.9, 0.85), just = "right", lwd = 1, text.size = 0.7) +
  tm_compass(position = c(0.05, 0.05), just = "left", size = 2, type = "arrow") +
  tm_shape(island_map_data) +
  tm_polygons(col = "gray10", border.col = "gray10") +
  tm_shape(island_map_data) +
  tm_borders(col = "gray10", lwd = 2)

p

# Save to SVG for further polishing in Inkscape
svg("figures2/map.svg", width = 7, height = 7)
supfig1
dev.off()

#### THE REST DOES NOT WORK AND IS DRIVING ME CRAZY ####

# Map the sampling sites
sampling_sites <- do.call("st_sfc", lapply(seq_len(nrow(data)), function(i) {

  curr_coord <- data[i, c("longitude", "latitude")] %>% as.matrix()
  return(st_point(curr_coord))

}))

sampling_sites <- st_sf(data, geometry = sampling_sites)

# Facetted map of each island
p <- tm_shape(island_map_data) +
  tm_polygons(col = "palegoldenrod", border.col = "gray50") +
  tm_facets(by = "island", nrow = 3) +
  tm_style("classic") +
  tm_scale_bar(position = c(0.1, 0.75), just = "left", lwd = 1, text.size = 2) +
  tm_compass(position = c(0.05, 0.05), just = "left", size = 1, text.size = 0.5, type = "arrow") +
  tm_shape(sampling_sites) +
  tm_dots() +
  tm_symbols(
    title.col = "Habitat",
    col = "habitat",
    shape = "habitat",
    border.col = "black",
    legend.shape.show = FALSE,
    legend.col.is.portrait = FALSE,
    size = 1,
    palette = c("darkgoldenrod1", "forestgreen", "mediumseagreen")
  ) +
  tm_grid(alpha = 0.5, n.x = 4, n.y = 4, labels.size = 0.8) +
  tm_legend(
    legend.outside.position = "bottom",
    legend.outside = TRUE,
    legend.just = "center"
  )
p
