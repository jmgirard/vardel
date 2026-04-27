# from {agreement} by Jeffrey Girard, PhD
# https://github.com/jmgirard/agreement?tab=readme-ov-file


#' Calculate Chance-Adjusted Agreement
#'
#' Description
#'
#' @param .data *Required.* A matrix or data frame in tall format containing
#'   categorical data where each row corresponds to a single score (i.e.,
#'   assignment of an object to a category) Cells should contain numbers or
#'   characters indicating the discrete category that the corresponding rater
#'   assigned the corresponding object to. Cells should contain \code{NA} if a
#'   particular assignment is missing (e.g., that object was not assigned to a
#'   category by that rater).
#' @param object *Optional.* The name of the variable in \code{.data}
#'   identifying the object of measurement for each observation, in non-standard
#'   evaluation without quotation marks. (default = \code{Object})
#' @param rater *Optional.* The name of the variable in \code{.data} identifying
#'   the rater or source of measurement for each observation, in non-standard
#'   evaluation without quotation marks. (default = \code{Rater})
#' @param score *Optional.* The name of the variable in \code{.data} containing
#'   the categorical score or rating/code for each observation, in non-standard
#'   evaluation without quotation marks. (default = \code{Score})
#' @param approach *Optional.* A string or vector of strings specifying the
#'   chance-adjustment approach(es) to use. Currently, the "alpha", "gamma",
#'   "irsq", "kappa", "pi", and "s" approaches are available. (default =
#'   c("alpha", "gamma", "kappa", "irsq", "pi", "s"))
#' @param categories *Optional.* A vector (numeric, character, or factor)
#'   containing all possible categories that objects could have been assigned
#'   to. When this argument is omitted or set to \code{NULL}, the possible
#'   categories are assumed to be those observed in \code{.data}. However, in
#'   the event that not all possible categories are observed in \code{.data},
#'   this assumption may be misleading and so the possible categories, and their
#'   ordering, can be explicitly specified. (default = NULL)
#' @param weighting *Optional.* A single string specifying the type of weighting
#'   scheme to use. Weighting schemes allow the accommodation of ordered and
#'   unordered categories with the same formulas. Currently, "identity" weights
#'   are available for unordered/nominal categories, both "linear" and
#'   "quadratic" weights are available for ordered categories, and "custom"
#'   weights can be specified via \code{custom_weights}. (default = "identity")
#' @param agreement *Optional.* Either \code{NULL} or a single string specifying
#'   the formula to use in calculating percent observed agreement. Currently,
#'   "objects" is available to calculate agreement averaged across objects,
#'   "pairs" is available to calculate agreement averaged across object-rater
#'   pairs, and "kripp" is available to calculate agreement using Krippendorff's
#'   formula. \code{NULL} sets agreement to the default formula for each
#'   approach (i.e., "kripp" for Krippendorff's alpha, "pairs" for Van Oest's
#'   irsq, and "objects" for all others). (default = NULL)
#' @param bootstrap *Optional.* A single non-negative integer that specifies how
#'   many bootstrap resamplings should be computed (used primarily for
#'   estimating confidence intervals and visualizing uncertainty). To skip
#'   bootstrapping, set this argument to 0. (default = 2000)
#' @param alpha_c *Optional.* Either \code{NULL} or a vector of numbers
#'   corresponding to the alpha_c parameters in Van Oest's formula. If
#'   \code{NULL}, and irsq is estimated, a vector of ones will be used to
#'   implement the uniform prior coefficient. (default = NULL)
#' @param custom_weights *Optional.* Either \code{NULL} or a q-by-q weight
#'   matrix where q is the number of unique categories. Weights must be between
#'   0 (no credit) and 1 (full credit). (default = NULL)
#' @param warnings *Optional.* A single logical value that specifies whether
#'   warnings should be displayed. (default = TRUE).
#' @return An object of type 'cai' containing the results and details.
#'   \describe{\item{approach}{A character vector containing the name of each
#'   approach in order} \item{observed}{A numeric vector containing the raw
#'   observed agreement according to each approach} \item{expected}{A numeric
#'   vector containing the expected chance agreement according to each approach}
#'   \item{adjusted}{A numeric vector containing the chance-adjusted agreement
#'   according to each approach. Note that these values are those typically
#'   named after each approach (e.g., this is the kappa coefficient)}
#'   \item{boot_results}{A list containing the results of the bootstrap
#'   procedure} \item{details}{A list containing the details of the analysis,
#'   such as the formatted \code{codes}, relevant counts, weighting scheme and
#'   weight matrix.} \item{call}{The function call that created these results.}}
#' @references Gwet, K. L. (2014). *Handbook of inter-rater reliability: The
#'   definitive guide to measuring the extent of agreement among raters* (4th
#'   ed.). Gaithersburg, MD: Advanced Analytics.
#' @references van Oest, R. (2019). A new coefficient of interrater agreement:
#'   The challenge of highly unequal category proportions. *Psychological
#'   Methods, 24*(4), 439-451. \url{https://doi.org/10/ggbk3f}
#' @family functions for categorical data
#' @family functions for chance-adjusted agreement
#' @export
cat_vardel_adjusted <- function(.data,
                         object = "ObjectID",
                         rater = "RaterID",
                         score = "Score",
                         approach = c("kappa","s"),
                         categories = NULL,
                         weighting = c("identity"),
                         agreement = NULL,
                         bootstrap = 0,
                         alpha_c = NULL,
                         custom_weights = NULL,
                         warnings = FALSE) {

  # Validate inputs
  assertthat::assert_that(is.data.frame(.data) || is.matrix(.data))
  approach <- match.arg(approach, several.ok = TRUE)
  assertthat::assert_that(is_null(categories) || is_vector(categories))
  #weighting <- match.arg(weighting)
  assertthat::assert_that(
    is.null(agreement) ||
    all(agreement %in% c("objects", "pairs", "kripp"))
  )
  assertthat::assert_that(bootstrap == 0 || assertthat::is.count(bootstrap))
  assertthat::assert_that(assertthat::is.flag(warnings))
  assertthat::assert_that(is_null(alpha_c) || is.numeric(alpha_c))

  assertthat::assert_that(is_null(custom_weights) || is.matrix(custom_weights))

  # Prepare .data for analysis
  d <- prep_data_cat(
    .data = .data,
    object = {{object}},
    rater = {{rater}},
    score = {{score}},
    approach = approach,
    categories = categories,
    weighting = weighting,
    agreement = agreement,
    alpha_c = alpha_c,
    custom_weights = custom_weights,
    bootstrap = bootstrap
  )

  # Prepare empty results in case of errors
  n_approach <- length(approach)
  out <- new_cai(
    approach = approach,
    observed = rep(NA_real_, n_approach),
    expected = rep(NA_real_, n_approach),
    adjusted = rep(NA_real_, n_approach),
    boot_results = list(
      t = matrix(NA, nrow = 1, ncol = 3 * n_approach),
      t0 = rep(NA, 3 * n_approach)
    ),
    details = d,
    call = match.call()
  )

  # Warn about bootstrapping samples with less than 20 objects
  if (d$n_objects < 20 && bootstrap > 0 && warnings == TRUE) {
    warning("With a small number of objects, bootstrap confidence intervals may not be stable.")
  }

  # Warn about bootstrapping with fewer than 1000 resamples
  if (bootstrap > 0 && bootstrap < 1000 && warnings == TRUE) {
    warning("To get stable confidence intervals, consider using more bootstrap resamples.")
  }

  # Warn about there being fewer than 2 categories
  if (d$n_categories < 2) {
    if (warnings == TRUE) {
      warning("Only a single category was observed or requested. Returning NA.\nHint: Try setting the possible categories explicitly with the categories argument")
    }
    return(out)
  }
  if (d$n_raters < 2) {
    if (warnings == TRUE) {
      warning("Only a single rater was observed. Returning NA.")
    }
    return(out)
  }

  # Create function to perform bootstrapping
  boot_function <- function(ratings,
                            index,
                            function_list,
                            categories,
                            weight_matrix,
                            agreement,
                            alpha_c) {

    resample <- ratings[index, , drop = FALSE]
    bsr <- rep(NA_real_, times = length(function_list) * 3)
    # Loop through approaches
    for (i in seq_along(function_list)) {
      bsr[(i * 3 - 2):(i * 3)] <- function_list[[i]](
        codes = resample,
        categories = categories,
        weight_matrix = weight_matrix,
        agreement = agreement,
        alpha_c = alpha_c
      )
    }

    bsr
  }

  # Collect functions into vector to speed up bootstrapping
  expr_list <- parse(text = paste0("calc_", approach))
  function_list <- NULL
  for (i in 1:length(expr_list)) {
    function_list <- c(function_list, eval(expr_list[[i]]))
  }

  # Calculate the bootstrap results
  boot_results <-
    boot::boot(
      data = d$ratings,
      statistic = boot_function,
      R = bootstrap,
      function_list = function_list,
      categories = d$categories,
      weight_matrix = d$weight_matrix,
      agreement = d$agreement,
      alpha_c = d$alpha_c
    )

  # Construct cai class output object
  # out <- new_cai(
  #   approach = approach,
  #   observed = boot_results$t0[seq(from = 1, to = length(approach) * 3, by = 3)],
  #   expected = boot_results$t0[seq(from = 2, to = length(approach) * 3, by = 3)],
  #   adjusted = boot_results$t0[seq(from = 3, to = length(approach) * 3, by = 3)],
  #   boot_results = boot_results,
  #   details = d,
  #   call = match.call()
  # )

  kap <- boot_results$t0[seq(from = 3, to = length(approach) * 3, by = 3)][1]
  sbp <- boot_results$t0[seq(from = 3, to = length(approach) * 3, by = 3)][2]

  res_caa <- tibble::tibble(
    method = "kappa",
    estimate = signif(kap, digits = 3),
    sigma_s = NA,
    sigma_r = NA,
    sigma_vsr = NA,  
    vs = NA,
    vr = NA,
    vsr = NA,   
    message = "",
    warning = "",
    error = ""


  ) |> dplyr::add_row(
    method = "s_bp", 
    estimate = signif(sbp, digits = 3),
    sigma_s = NA,
    sigma_r = NA,
    sigma_vsr = NA,  
    vs = NA,
    vr = NA,
    vsr = NA,   
    message = "",
    warning = "",
    error = ""
  )

  return(res_caa)
}


# Prepare categorical data for analysis -----------------------------------
prep_data_cat <- function(.data,
                          object,
                          rater,
                          score,
                          approach = NULL,
                          categories = NULL,
                          weighting = "identity",
                          agreement = NULL,
                          alpha_c = NULL,
                          custom_weights = NULL,
                          bootstrap = 0) {

  out <- list()

  # Ensure df is a tibble
  df <- tibble::as_tibble(.data)

  # Select the important variables
  df <- dplyr::select(df, {{object}}, {{rater}}, {{score}})

  # Add explicit NA rows to missing object-rater combinations
  #df <- tidyr::complete(df, {{object}}, {{rater}})

  # Reorder df by rater and object so that scores fills out matrix properly
  df <- dplyr::arrange(df, {{rater}}, {{object}})

  # Get and count each variable's unique values
  out$objects <- unique(dplyr::pull(df, {{object}}))
  out$raters <- unique(dplyr::pull(df, {{rater}}))
  out$n_objects <- length(out$objects)
  out$n_raters <- length(out$raters)

  # Pull scores, convert NaN to NA, and count NAs
  scores <- dplyr::pull(df, {{score}})
  scores[is.nan(scores)] <- NA
  out$n_missing_scores <- sum(rlang::are_na(scores))

  # Get and count observed categories
  cat_observed <- sort(unique(scores))
  n_cat_observed <- length(cat_observed)

  # If specified, get and count possible categories
  if (is_null(categories)) {
    cat_possible <- cat_observed
    n_cat_possible <- n_cat_observed
  } else {
    if (is.factor(categories)) {
      cat_possible <- levels(categories)
    } else {
      cat_possible <- unique(categories)
    }
    n_cat_possible <- length(cat_possible)
  }
  out$categories <- cat_possible
  out$n_categories <- n_cat_possible

  # Format ratings into an object-by-rater matrix
  out$ratings <- matrix(
    scores,
    nrow = out$n_objects,
    ncol = out$n_raters,
    dimnames = list(out$objects, out$raters)
  )

  # Drop objects and raters that contain only missing values
  out <- remove_uncoded(out)

  # Validate basic counts
  assertthat::assert_that(out$n_objects >= 1,
              msg = "There must be at least 1 valid object in `.data`.")
  assertthat::assert_that(out$n_raters >= 2,
              msg = "There must be at least 2 raters in `.data`.")

  # Validate categories
  cat_unknown <- setdiff(cat_observed, cat_possible)
  assertthat::assert_that(is_empty(cat_unknown),
              msg = "A category not in `categories` was observed in `.data`.")

  # Get weight matrix
  out$weighting <- weighting
  if (weighting == "custom") {
    assertthat::assert_that(ncol(custom_weights) == nrow(custom_weights))
    assertthat::assert_that(ncol(custom_weights) == n_cat_possible)
    out$weight_matrix <- custom_weights
  } else {
    out$weight_matrix <- calc_weights(weighting, cat_possible)
  }

  # Add other information to d
  out$approach <- approach
  out['agreement'] <- list(agreement)
  out$bootstrap <- bootstrap

  # Set up alpha_c
  assertthat::assert_that(
    is.null(alpha_c) ||
    length(alpha_c) == 1 ||
    length(alpha_c) == n_cat_possible
  )
  if ("irsq" %in% approach && is.null(alpha_c)) {
    alpha_c <- rep(1, times = n_cat_possible)
  } else if ("irsq" %in% approach && length(alpha_c) == 1) {
    alpha_c <- rep(alpha_c, times = n_cat_possible)
  }
  alpha_c[alpha_c == Inf] <- 1e6
  out['alpha_c'] <- list(alpha_c)

  out
}

# Drop objects and raters that contain only missing values ----------------
remove_uncoded <- function(x) {
  mat <- x$ratings
  mat <- mat[rowSums(rlang::are_na(mat)) != ncol(mat), , drop = FALSE]
  mat <- mat[, colSums(rlang::are_na(mat)) != nrow(mat), drop = FALSE]
  out <- x
  out$ratings <- mat
  out$n_objects <- nrow(mat)
  out$n_raters <- ncol(mat)
  out
}

# Calculate chance-adjusted index -----------------------------------------
adjust_chance <- function(poa, pea) {
  (poa - pea) / (1 - pea)
}

# Convert from codes to rater counts in object-by-category matrix ---------
raters_obj_cat <- function(codes, categories) {
  table(
    row(codes),
    factor(unlist(codes), levels = categories),
    useNA = "no"
  )
}

# Convert from codes to object counts in rater-by-category matrix ---------
objects_rat_cat <- function(codes, categories) {
  table(
    col(codes),
    factor(unlist(codes), levels = categories),
    useNA = "no"
  )
}

# Calculate weight matrix -------------------------------------------------
#' @export
calc_weights <- function(type = c("identity", "linear", "quadratic"),
                         categories) {

  type <- match.arg(type, several.ok = FALSE)

  # Count the categories
  n_categories <- length(categories)

  # Convert to numeric if necessary
  if (!is.numeric(categories) && type != "identity") {
    category_values <- 1:n_categories
    warning(
      "Numeric categories are required for ", type, " weights.\n",
      "Converting to integers from 1 to the number of categories."
    )
  } else {
    category_values <- categories
  }

  # Start with diagonal matrix
  weight_matrix <- diag(n_categories)
  rownames(weight_matrix) <- categories
  colnames(weight_matrix) <- categories

  # If categories are ordered, calculate weights
  if (type != "identity") {
    max_distance <- diff(range(category_values))
    for (i in seq_along(categories)) {
      for (j in seq_along(categories)) {
        obs_distance <- category_values[[i]] - category_values[[j]]
        if (type == "linear") {
          weight_matrix[i, j] <- 1 - abs(obs_distance) / max_distance
        } else if (type == "quadratic") {
          weight_matrix[i, j] <- 1 - obs_distance^2 / max_distance^2
        }
      }
    }
  }

  weight_matrix
}

# safe_boot.ci ------------------------------------------------------------

safe_boot.ci <- function(boot.out, level, type, index, ...) {

  # Determine the location of results in bootci object
  if (type == "bca") {
    field <- "bca"
    elems <- 4:5
  } else if (type == "perc") {
    field <- "percent"
    elems <- 4:5
  } else if (type == "basic") {
    field <- "basic"
    elems <- 4:5
  } else if (type == "norm") {
    field <- "normal"
    elems <- 2:3
  }

  # Check for constant or completely missing results and replace with NAs
  stat_var <- stats::var(boot.out$t[, index], na.rm = TRUE)
  if (is.na(stat_var) || dplyr::near(stat_var, 0)) {
    out <- c(NA, NA)
  } else {
    out <-
      boot::boot.ci(
        boot.out = boot.out,
        conf = level,
        type = type,
        index = index,
        ...
      )[[field]][elems]
  }
  out
}

# new_cai -----------------------------------------------------------------

# S3 Constructor for cai class
new_cai <- function(approach = character(),
                    observed = double(),
                    expected = double(),
                    adjusted = double(),
                    boot_results = list(),
                    details = list(),
                    call = rlang::quo(),
                    ...) {

  new_s3_scalar(
    approach = approach,
    observed = observed,
    expected = expected,
    adjusted = adjusted,
    boot_results = boot_results,
    details = details,
    call = call,
    ...,
    class = "agreement_cai"
  )

}

# summary.agreement_cai ---------------------------------------------------

# Summary method for objects of cai class
#' @method summary agreement_cai
#' @export
summary.agreement_cai <- function(object,
                                  digits = 3,
                                  ci = TRUE,
                                  ...) {

  # Validate inputs
  assertthat::assert_that(digits == 0 || assertthat::is.count(digits))
  assertthat::assert_that(assertthat::is.flag(ci))

  # Print function call and header
  cat(
    "\nCall:\n",
    paste(deparse(object$call), sep = "\n", collapse = "\n"),
    "\n\n",
    "Objects = \t", object$details$n_objects,
    "\nRaters = \t", object$details$n_raters,
    "\nCategories = \t{", paste(object$details$categories, collapse = ", "), "}",
    "\nWeighting = \t", object$details$weighting,
    "\n\n",
    "Chance-Adjusted Categorical Agreement with Bootstrapped CIs\n\n",
    sep = ""
  )

  # Create matrix containing basic results
  v <- c(object$observed, object$expected, object$adjusted)
  m <- round(matrix(v, ncol = 3L), digits)
  rownames(m) <- object$approach
  colnames(m) <- c("Observed", "Expected", "Adjusted")

  # Append Confidence Intervals for Adjusted Column if requested
  if (ci == TRUE) {
    ci <- round(
      stats::confint(
        object,
        ...
        )[seq(3, length(object$approach) * 3, by = 3), ],
      digits
    )
    m <- cbind(m, ci)
  }

  # Print matrix
  print(m, print.gap = 3L, na.print = "")
  cat("\n")

}

# confint.agreement_cai ---------------------------------------------------

#' Calculate confidence intervals for chance-adjusted agreement indexes
#'
#' Calculate confidence intervals for chance-adjusted agreement indexes using
#' one of several approaches.
#'
#' @param object A chance-adjusted agreement or "cai" object.
#' @param level A single real number between 0 and 1 that indicates the
#'   confidence level or width of the confidence intervals.
#' @param type A single string indicating which type of bootstrap confidence
#'   interval to calculate. Options currently include "bca" for bias-corrected
#'   and accelerated bootstrap CIs, "perc" for bootstrap percentile CIs, "basic"
#'   for basic bootstrap CIs, and "norm" for normal approximation CIs.
#' @param ... Further arguments to be passed to \code{boot::boot.ci()}.
#' @method confint agreement_cai
#' @export
confint.agreement_cai <- function(object,
                                  level = 0.95,
                                  type = c("bca", "perc", "basic", "norm"),
                                  ...) {

  # Validate inputs
  assertthat::assert_that(is.number(level), level > 0, level < 1)
  type <- match.arg(type, several.ok = FALSE)

  # Calculate CI bounds and store in list
  ci_list <-
    lapply(
      1:ncol(object$boot_results$t),
      function (i) safe_boot.ci(
        object$boot_results,
        level = level,
        type = type,
        index = i,
        ...
      )
    )

  # Convert list to matrix and add dimnames
  out <- matrix(unlist(ci_list), ncol = 2, byrow = TRUE)
  rownames(out) <-
    paste0(
      rep(object$approach, each = 3),
      " ",
      c("Observed", "Expected", "Adjusted")
    )
  colnames(out) <- c("lower", "upper")

  out
}

# Create a new S3 class from a list
new_s3_lst <- function(x, ..., class) {
  stopifnot(is.list(x))
  stopifnot(is.character(class))
  structure(x, ..., class = class)
}

# Create a new S3 class from a scalar
new_s3_scalar <- function(..., class) {
  new_s3_lst(list(...), class = class)
}

# Redirect to the desired formula for percent observed agreement
calc_agreement <- function(codes,
                           categories,
                           weight_matrix,
                           formula = c("objects", "pairs", "kripp")) {

  formula <- match.arg(formula, several.ok = FALSE)

  if (formula == "objects") {
    calc_agreement_objects(codes, categories, weight_matrix)
  } else if (formula == "pairs") {
    calc_agreement_pairs(codes, categories, weight_matrix)
  } else if (formula == "kripp") {
    calc_agreement_kripp(codes, categories, weight_matrix)
  }
}


# Calculate percent observed agreement averaged over objects
# Gwet (2014)
calc_agreement_objects <- function(codes, categories, weight_matrix) {

  # How many raters assigned each object to each category?
  r_oc <- raters_obj_cat(codes, categories)

  # How many raters assigned each object to any category?
  r_o <- rowSums(r_oc)

  # How much agreement was observed for each object-category combination?
  obs_oc <- r_oc * (t(weight_matrix %*% t(r_oc)) - 1)

  # How much agreement was observed for each object across all categories?
  obs_o <- rowSums(obs_oc)

  # How much agreement was maximally possible for each object?
  max_o <- r_o * (r_o - 1)

  # What was the percent observed agreement for each object?
  poa_o <- obs_o[r_o >= 2] / max_o[r_o >= 2]

  # What was the percent observed agreement across all objects?
  poa <- mean(poa_o)

  poa
}



############ CAA approach #1: Bennet's S/brennen prediger 

#' @export
cat_s <- function(.data, ...) {
  cat_vardel_adjusted(.data, approach = "s", agreement = "objects", ...)
}

# Worker function to calculate the S score and its components
calc_s <- function(codes, categories, weight_matrix, agreement, ...) {

  # Default to agreement averaged over objects
  if (is.null(agreement)) agreement <- "objects"

  # Calculate percent observed agreement
  poa <- calc_agreement(codes, categories, weight_matrix, agreement)

  # Calculate percent expected agreement
  pea <- calc_chance_s(codes, categories, weight_matrix)

  # Calculate chance-adjusted index
  cai <- adjust_chance(poa, pea)

  # Create and label output vector
  out <- c(POA = poa, PEA = pea, CAI = cai)

  out
}

# Worker function to calculate expected agreement using the S model of chance
calc_chance_s <- function(codes, categories, weight_matrix) {

  # How many categories were possible?
  n_categories <- length(categories)

  # How much chance agreement is expected for each combination of categories?
  pea_cc <- weight_matrix / n_categories^2

  # How much chance agreement is expected across all combinations of categories?
  pea <- sum(pea_cc)

  pea
}



######## CAA apporach #2

#' @export
cat_kappa <- function(.data, ...) {
  cat_vardel_adjusted(.data, approach = "kappa", agreement = "objects", ...)
}

# Calculate the kappa coefficient and its components
calc_kappa <- function(codes, categories, weight_matrix, agreement, ...) {

  # Default to agreement averaged over objects
  if (is.null(agreement)) agreement <- "objects"

  # Calculate percent observed agreement
  poa <- calc_agreement(codes, categories, weight_matrix, agreement)

  # Calculate percent expected agreement
  pea <- calc_chance_kappa(codes, categories, weight_matrix)

  # Calculate chance-adjusted index
  cai <- adjust_chance(poa, pea)

  # Create and label output vector
  out <- c(POA = poa, PEA = pea, CAI = cai)

  out
}

# Calculate expected agreement using the kappa model of chance
calc_chance_kappa <- function(codes, categories, weight_matrix) {

  n_raters <- ncol(codes)
  n_categories <- length(categories)

  # How many objects did each rater assign to each category?
  o_rc <- objects_rat_cat(codes, categories)

  # How many objects did each rater assign to any category?
  o_r <- rowSums(o_rc)

  # How many objects could each rater have assigned to each category?
  omax_rc <- o_r %*% matrix(1, ncol = n_categories)

  # What was the prevalence of each category for each rater?
  prev_rc <- o_rc / omax_rc

  # What was the prevalance of each category, averaged across all raters?
  prev_c <- colMeans(prev_rc)

  # TODO: Add interpretations and informative variable names
  x <- t(prev_rc) %*% prev_rc #dot-products for each rater-category combination
  y <- prev_c %*% t(prev_c) #dot-products for the average rater
  z <- (x - n_raters * y) / (n_raters - 1) #scaling or correction?

  # What is the probability of two categories being assigned at random?
  exp_cc <- y - z / n_raters

  # How much chance agreement is expected for each combination of categories?
  pea_cc <- weight_matrix * exp_cc

  # How much chance agreement is expected across all combinations of categories?
  pea <- sum(pea_cc, na.rm = TRUE)

  pea
}
