#' Classification analysis
#'
#' Perform a replicated classification analysis of a multivariate dataset into categorical labels using machine learning tools and k-fold cross validation
#'
#' @param data A data frame
#' @param variables The variables used to classify
#' @param grouping Name of the grouping variable (the labels)
#' @param nesting Optional nesting variable, if the analysis must be conducted separately on different subsets of the data
#' @param method The data mining model used. Currently supports "SVM" and "LDA".
#' @param k Number of bins for the k-fold cross-validation procedure
#' @param nrep Number of replicate analyses (i.e. number of k-fold cross validations)
#' @param minsize Minimum size required per group for a training data set
#' @param seed Optional random seed to reset at the beginning
#' @param importance Whether to perform sensitivity analysis on the input (takes a while)
#' @param return_machine Whether to return the machines (takes space)
#' @param verbose Whether to display messages
#' @param pb Whether to display progress bars
#' @param to_pcomp Variable to perform PCA on
#' @param center Center the PCA
#' @param scale Scale the PCA
#'
#' @return A list of two tibbles.
#'
#' The first tibble, `full`, contains information about each fitted machine in the analysis, including its nesting level `nesting`, replicate number `repl` and cross-validation bin `kbin` but also classification testing outputs such as its confusion matrix `confmat`, its prediction `accuracy` and the number of observations in the testing set `ntested`. If `return_machine` is TRUE, an extra column is added, that is a list of the actual fitted machines. If `importance` is TRUE, extra columns are added, named after the input variables and containing their respective relative importance scores (these are computed using the `rminer::Importance` function for one-dimensional sensitivity analysis).
#'
#' The second tibble, `digested`, is a summary with, for each `nesting` level: the average confusion matrix `confmat` over all `nrep` replicate total confusion matrices (which are sums of `k` confusion matrices, one per cross-validation bin), the `mean` accuracy and its standard error `stderr` across replicates, Stouffer's Z-statistic `zs` used to combine P-values from binomial tests for each replicate (using the Z-transform method, see Whitlock 2005) and the corresponding combine P-value `pcombined`.
#'
#' @export

classify3 <- function(
  data,
  variables,
  grouping,
  nesting = NULL,
  method = "SVM",
  k = 5,
  nrep = 1,
  minsize = 5,
  seed = NULL,
  importance = FALSE,
  return_machine = FALSE,
  verbose = TRUE,
  pb = TRUE,
  to_pcomp = NULL,
  center = TRUE,
  scale = TRUE
) {

  # Convert the dataset into a tibble if it is not already one
  if (!inherits(data, "tbl")) data <- tibble::as_tibble(data)

  # Compute principal components if needed
  if (!is.null(to_pcomp)) {
    data <- data %>%
      cbind(
        npcomp(
          data, to_pcomp, center, scale, nesting, combine = TRUE,
          reduce = variables
        )$x
      )
  }

  # Random seed
  if (!is.null(seed)) set.seed(seed)

  # Proportion of observations in the testing set
  ptesting <- 1 / k

  assertthat::assert_that(nrow(data) > k)
  assertthat::assert_that(floor(ptesting * nrow(data)) > 0)

  # Define the possible labels
  labels <- unique(data[[grouping]])

  nested <- TRUE

  # Artificial single nesting level if no nesting supplied
  if (is.null(nesting)) {
    nested <- FALSE
    data$nesting <- factor(1)
    nesting <- "nesting"
  }

  # Split the data into nesting levels
  data <- data %>% split(f = .[, nesting])

  # Decide on a looping function depending on whether we want a progress bar
  if (!verbose) pb <- FALSE
  thislapply1 <- thislapply2 <- lapply
  if (pb) {
    if (nested) {
      thislapply1 <- pbapply::pblapply
    } else {
      thislapply2 <- pbapply::pblapply
    }
  }

  # For each nesting level...
  machines <- thislapply1(data, function(data) {

    is_fine <- FALSE

    # Sample the training and testing sets until valid
    while (!is_fine) {

      # Randomly assign data to k testing groups
      groups <- rep(seq_len(k), each = floor(ptesting * nrow(data)))
      if (length(groups) < nrow(data)) {
        groups <- c(groups, seq_len(nrow(data) - length(groups)))
      }
      assertthat::assert_that(length(groups) == nrow(data))
      groups <- sample(groups, replace = FALSE)

      # Check that each label is sufficiently represented within each training
      # group
      is_fine <- all(sapply(seq_len(k), function(j) {

        represented <- table(data[which(groups != j), grouping])
        if (!all(labels %in% names(represented))) return (FALSE)
        return (all(represented > minsize))

      }))

    }

    # For each testing group...
    lapply(seq_len(k), function(j) {

      # Sample indices for a training dataset
      training <- which(groups != j)

      assertthat::assert_that(length(training) < nrow(data))

      # Downsample the training dataset to the size of the least represented
      # label
      targetsize <- min(table(data[training, grouping]))
      training <- do.call(
        "c",
        lapply(
          labels,
          function(lab) {
            sample(
              training[data[training, grouping] == lab], targetsize,
              replace = FALSE
            )
          }
        )
      )

      assertthat::assert_that(
        all(table(data[training, grouping]) == targetsize)
      )

      # Set up model formula
      model <- as.formula(
        paste(grouping, "~", paste(variables, collapse = " + "))
      )

      if (method == "SVM") {

        # Fit a support vector machine to the data
        machine <- rminer::fit(
          model, data = data[training, c(variables, grouping)], model = "svm",
          kernel = "rbfdot", task = "class"
        )

      } else if (method == "LDA") {

        # Or a linear discriminant analysis
        machine <- MASS::lda(
          formula = model, data = data[training, c(variables, grouping)]
        )

      } else stop("unknown method")

      # Sensitivity analysis of the fitted model
      imp <- NULL
      if (importance) {
        if (method == "LDA") {
          imp <- rminer::Importance(
            machine, data = data[training, c(variables, grouping)],
            PRED = function(M, data) predict(M, data)$class
          )
        } else  if (method == "SVM") {
          imp <- rminer::Importance(
            machine, data = data[training, c(variables, grouping)]
          )
        }
        imp <- imp$imp[seq_along(variables)]
        names(imp) <- variables
      }

      # Predict the labels of the remaining data
      predictions <- rminer::predict(
        machine, newdata = data[groups == j, variables]
      )

      if (method == "LDA") predictions <- predictions$class

      # Compare true and predicted labels
      conf <- table(predictions, data[[grouping]][groups == j])

      if (!return_machine) machine <- NULL

      # Return the confusion matrix, the results of the importance analysis
      # and the machine itself, if needed
      return (list(confmat = conf, importance = imp, machine = machine))

    }) # end of cross-validation bin

  }) # end of nesting level

  # Prepare a data frame with results for each machine
  res <- tidyr::expand_grid(nesting = names(data), kbin = seq(k))

  # Fill in that data frame with the output of each machine
  res <- res %>%
    dplyr::group_by(nesting, kbin) %>%
    tidyr::nest() %>%
    dplyr::mutate(

      # Confusion matrix
      confmat = purrr::pmap(
        list(nesting, kbin),
        ~ purrr::pluck(machines, ..1, ..2)$confmat
      ),

      # Vector of importance scores
      importance = purrr::pmap(
        list(nesting, kbin),
        ~ purrr::pluck(machines, ..1, ..2)$importance
      ),

      # Fitted machine
      machine = purrr::pmap(
        list(nesting, kbin),
        ~ purrr::pluck(machines, ..1, ..2)$machine
      )

    ) %>%
    dplyr::select(-data) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      accuracy = purrr::map_dbl(confmat, pdiag),
      ntested = purrr::map_int(confmat, sum)
    )

  # Summarize accuracy across machines for each nesting level
  smr <- res %>%
    dplyr::group_by(nesting) %>%
    tidyr::nest() %>%
    dplyr::mutate(

      # Summed confusion matrix
      confmat = purrr::map(data, ~ Reduce('+', .x$confmat)),

      # Total accuracy
      accuracy = purrr::map_dbl(confmat, pdiag),

      # Number of tested points (= sample size of each nesting level)
      ntested = purrr::map_int(confmat, sum),

      # Is accuracy greater than expected by chance? (one-tailed binomial)
      pvalue = purrr::map2_dbl(
        confmat, ntested,
        ~ binom.test(
          x = sum(diag(.x)),
          p = 1 / length(labels),
          n = .y,
          alternative = "greater"
        )$p.value
      )

    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-data)

  # If the sensitivity analysis was conducted...
  if (importance) {

    # Place importance scores in multiple columns, one per variable
    res <- res %>%
      dplyr::bind_cols(purrr::map_dfr(res$importance, ~ .x))

  }

  # Remove the importance column from the output
  res <- res %>% dplyr::select(-importance)

  # Remove the machine column if the fitted model must not be returned
  if (!return_machine) res <- res %>% dplyr::select(-machine)

  # Return the full and summarized results
  return(list(full = res, digested = smr))

}
