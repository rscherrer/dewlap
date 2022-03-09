## ---------------------------
##
## Script name: 04_median_distance.R
##
## Purpose of script: This script shows the median within-island distance
## between two patches, computed across all islands.
##
## Author: Raphael Scherrer
##
## Date Created: 2022-03-09
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
library(geosphere)

data <- read_csv("metadata/sites.csv")

data %>%
  group_by(island) %>%
  nest() %>%
  mutate(dist = map(data, function(data) {

    y <- with(data, distm(matrix(c(longitude, latitude), ncol = 2), fun = distGeo))
    y[upper.tri(y)]

  })) %>%
  select(-data) %>%
  unnest(dist) %>%
  ungroup() %>%
  select(dist) %>%
  summarize(median = median(dist))


