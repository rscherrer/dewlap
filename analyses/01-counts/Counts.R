# This script was used to generate the table of counts across islands and habitats

rm(list = ls())

library(dewlap)

data <- read.csv("data/reflectance.csv", header = TRUE)
tab <- table(data$island, data$habitat)

save_table(tab, "analyses/01-counts/table_counts")
