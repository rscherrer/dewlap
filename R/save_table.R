#' Save tables for this study

save_table <- function(tab, nm, digits = 0) {

  library(knitr)

  write.csv(tab, paste0(nm, ".csv"), row.names = FALSE)
  tab <- knitr::kable(tab, "latex", digits = digits, booktabs = TRUE)
  texfile <- file(paste0(nm, ".tex"))
  writeLines(tab, texfile)
  close(texfile)

}
