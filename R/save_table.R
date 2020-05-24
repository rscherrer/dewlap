#' Save tables for this study
#'
#' @param tab The table
#' @param nm File name (without the extension)
#' @param ... Extra parameters for `knitr::kable`
#'
#' @return Nothing. Just saves.
#'
#' @export

save_table <- function(tab, nm, ...) {

  library(knitr)

  # Save table as CSV
  write.csv(tab, paste0(nm, ".csv"), row.names = FALSE)

  # Save LaTeX table
  if (is.null(col_names)) col_names <- NA
  tab <- knitr::kable(
    tab, "latex", booktabs = TRUE, linesep = "", escape = FALSE, ...
  )
  texfile <- file(paste0(nm, ".tex"))
  writeLines(tab, texfile)
  close(texfile)

}
