rm(list = ls())

library(readr)
library(knitr)

filename <- "results/group_comparisons/machine_learning/PC/lda/summary.csv"

tab <- read_csv(filename)

# Process significance column
colnames(tab)[ncol(tab)] <- ""
tab[[ncol(tab)]][is.na(tab[[ncol(tab)]])] <- ""

# Make a latex table
kable(tab, "latex", booktabs = TRUE, linesep = "", escape = FALSE)
