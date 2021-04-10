# Function to add a significance column with asterisk symbols to a tibble
# (We are going to need that)
# Use on a tibble that contains P-values from test results in a column
# Arguments:
# D: the data frame
# col: the name of the column containing the P-values
# levels: the significance levels for the three labels *, ** and ***, in that
# order
# dropname: whether to replace the name of the extra column by a blank

mutate_signif <- function(
  D, col = "pvalue", levels = c(0.05, 0.01, 0.001), dropname = TRUE
) {

  D <- D %>%
    mutate(
      .signif = ifelse(get(col) < levels[1], "*", ""),
      .signif = ifelse(get(col) < levels[2], "**", .signif),
      .signif = ifelse(get(col) < levels[3], "***", .signif)
    )

  if (dropname) D <- D %>% rename(" " = ".signif")

  return(D)

}
