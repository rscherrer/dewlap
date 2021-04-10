#' Classification analysis --- testing versions
#'
#' This function is an old version of `classify`, used only to compare different
#' ways of computing P-values in their tendency to yield false positives.
#' We concluded that comparing a distribution of observed accuracy scores to
#' a null distribution was prone to false discovery because we have total
#' control over the number of replicates, i.e. the number of points in the
#' vector of observed scores.
#'
#' @export

classifyX <- function(
  data,
  variables,
  grouping,
  nesting = NULL,
  method = "SVM",
  k = 5,
  nrep = 1,
  nperm = 0,
  minsize = 5,
  seed = NULL,
  importance = FALSE,
  getmachine = FALSE,
  verbose = TRUE,
  pb = TRUE,
  digest = TRUE,
  topcomp = NULL,
  center = TRUE,
  scale = TRUE
) {

  # Compute principal components if needed
  if (!is.null(topcomp)) {
    data <- data %>%
      cbind(
        npcomp(
          data, topcomp, center, scale, nesting, combine = TRUE,
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
  labels <- unique(data[, grouping])

  nested <- TRUE

  # Artificial single nesting level if no nesting supplied
  if (is.null(nesting)) {
    nested <- FALSE
    data$nesting <- factor(1)
    nesting <- "nesting"
  }

  # Split the data into each nesting level
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
  results <- thislapply1(data, function(data) {

    # For each replicate (and randomized) analysis...
    thislapply2(seq(nrep + nperm), function(i) {

      # Randomize the data if needed
      if (i > nrep) data[[grouping]] <- sample(data[[grouping]])

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

        # Fit the classifier
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

        } else stop("unknown method") # could add other methods here

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
        conf <- table(predictions, data[groups == j, grouping])

        if (!getmachine) machine <- NULL

        return (list(conf = conf, imp = imp, machine = machine))

      }) # for each cross-validation bin
    }) # for each replicate
  }) # for each plot

  # If the results of many machines must be summarized...
  if (digest) {

    # For each plot...
    digested <- purrr::map_dfr(results, function(replicates) {

      # List of confusion matrices per replicate
      confs <- purrr::map(replicates, function(kbins) {

        # For each replicate sum the matrices of each testing set (i.e. k-bin)
        Reduce('+', purrr::map(kbins, ~ .x$conf))

      })

      # Separate the randomized and empirical confusion matrices
      if (nperm > 0) {

        assertthat::assert_that(length(confs) == nrep + nperm)
        randoms <- confs[seq(nrep + 1, length(confs))]
        confs <- confs[seq(nrep)]
        replicates <- replicates[seq(nrep)]

      }

      # Get the average confusion matrix for the plot
      avg <- mavg(confs)

      # Vector of replicate accuracies
      accus <- purrr::map_dbl(confs, pdiag)

      # Mean accuracy (just accuracy if only one replicate)
      mean <- mean(accus)

      prandom <- pttest <- NULL
      prob <- 1 / length(labels)

      # P-value from randomized distribution
      if (nperm > 0) {

        # Null distribution of accuracy scores
        accu0 <- purrr::map_dbl(randoms, pdiag)

        # Compute randomization P-value
        prandom <- 1 - length(which(mean > accu0)) / nperm

        # Parametrize the binomial distribution using the permuted one
        prob <- mean(accu0)

        # Compare entire distributions
        pttest <- t.test(accu0, accus, alternative = "greater")$p.value
        pksmirnov <- ks.test(accu0, accus, alternative = "greater")$p.value

      }

      # Compute binomial P-value
      n <- sum(confs[[1]])
      ntest <- floor(ptesting * n)
      pvalue <- 1 - pbinom(mean * ntest, ntest, 1 / length(labels))
      pvalue2 <- 1 - pbinom(mean * ntest, ntest, prob)

      # Optional data frame with importance scores
      if (importance) {
        imp <- purrr::map_dfr(replicates, function(kbins) {
          purrr::map_dfr(kbins, ~ .x$imp)
        })
      } else imp <- NULL

      # Return a row for the digested data frame over plots
      tibble::tibble(
        accu = mean,
        n = n,
        ptest = ptesting,
        ntest = ntest,
        pvalue = pvalue,
        pvalue2 = pvalue2,
        prob = prob,
        prandom = prandom,
        pttest = pttest,
        pksmirnov = pksmirnov,
        conf = list(avg),
        imp = list(imp)
      )

    })

    digested$nesting <- names(data)
    digested <- digested[, c(ncol(digested), seq(2, ncol(digested) - 1))]

    return(digested)

  }

  return (results)

}
