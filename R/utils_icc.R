# One-way ICCs ------------------------------------------------------------

## nested design, single rating
icc_1 <- function(vs, vrins) {
  vs / (vs + vrins)
}

## nested design, average ratings, balanced number of raters
icc_k <- function(vs, vrins, k) {
  vs / (vs + vrins / k)
}

## nested design, average ratings, unbalanced number of raters
icc_khat <- function(vs, vrins, khat) {
  vs / (vs + vrins / khat)
}

# Helper functions --------------------------------------------------------

create_srm <- function(.data,
                       subject = "subject",
                       rater = "rater",
                       score = "score") {
  # TODO: Allow subject, rater, and score to be specified using NSE
  assertthat::assert_that(is.data.frame(.data) || is.matrix(.data))
  cn <- colnames(.data)
  assertthat::assert_that(rlang::is_character(subject, n = 1), subject %in% cn)
  assertthat::assert_that(rlang::is_character(rater, n = 1), rater %in% cn)
  assertthat::assert_that(rlang::is_character(score, n = 1), score %in% cn)

  # Remove missing and non-finite scores
  na_index <- is.na(.data[[score]]) | !is.finite(.data[[score]])
  .data <- .data[!na_index, ]

  # Check if each subject was scored by each rater
  srm <- table(.data[[subject]], .data[[rater]], useNA = "no") > 0

  new_srm(srm)
}

create_parse <- function(.data, subject, rater, scores, v){
  # TODO: Check that is_twoway() works correctly for mv models
  # TODO: Check that everything works when scores contain weird characters
  twoway <- is_twoway(.data, subject, rater)
  if (twoway && v == 1) {
    formula <- brms::bf(paste0(
      scores, " ~ 1 + (1 | ", subject, ") + (1 | ", rater, ")"
    ))
  } else if (!twoway && v == 1) {
    formula <- brms::bf(paste0(
      scores, " ~ 1 + (1 | ", subject, ")"
    ))
  } else if (twoway && v > 1) {
    formula <- brms::bf(paste0(
      "mvbind(", paste(scores, collapse = ", "), ") | mi() ~ ",
      "1 + (1 | ", subject, ") + (1 | ", rater, ")"
    )) + brms::set_rescor(FALSE)
  } else if (!twoway && v > 1) {
    formula <- brms::bf(paste0(
      "mvbind(", paste(scores, collapse = ", "), ") | mi() ~ ",
      "1 + (1 | ", subject, ")"
    )) + brms::set_rescor(FALSE)
  } else {
    stop("Error determining model type")
  }
  return(formula)

}



create_parse_glmmtmb <- function(.data, subject, rater, scores, v){
  # TODO: Check that is_twoway() works correctly for mv models
  # TODO: Check that everything works when scores contain weird characters
  twoway <- is_twoway(.data, subject, rater)
  if (twoway && v == 1) {
    formula <- as.formula(paste0(
      scores, " ~ 1 + (1 | ", subject, ") + (1 | ", rater, ")"
    ))
  } else if (!twoway && v == 1) {
    formula <- as.formula(paste0(
      scores, " ~ 1 + (1 | ", subject, ")"
    ))
  } else if (twoway && v > 1) {
    formula <- as.formula(paste0(
      "mvbind(", paste(scores, collapse = ", "), ") | mi() ~ ",
      "1 + (1 | ", subject, ") + (1 | ", rater, ")"
    )) 
  } else if (!twoway && v > 1) {
    formula <- as.formula(paste0(
      "mvbind(", paste(scores, collapse = ", "), ") | mi() ~ ",
      "1 + (1 | ", subject, ")"
    )) 
  } else {
    stop("Error determining model type")
  }
  return(formula)

}

create_parseaov <- function(.data, subject, rater, scores, v){
  # TODO: Check that is_twoway() works correctly for mv models
  # TODO: Check that everything works when scores contain weird characters
  twoway <- is_twoway(.data, subject, rater)
  if (twoway && v == 1) {
    formula <- paste0("",scores," ~ ",subject, " + ", rater, "")

  } else if (!twoway && v == 1) {
    formula <- paste0("(",scores," ~ ",subject, ")")

  } else if (twoway && v > 1) {
    stop("Multiple scores in AOV? Need to do...")
  } else if (!twoway && v > 1) {
    stop("Multiple scores in AOV? Need to do...")
  } else {
    stop("Error determining model type")
  }
  return(formula)

}


check_convergence <- function(model) {
  warn <- model@optinfo$conv$lme4$messages
  is.null(warn) || !grepl('failed to converge', warn)
}

is_balanced <- function(.data, subject, rater) {
  # How many raters scored each subject?
  ks <- rowSums(table(.data[[subject]], by = .data[[rater]]))

  # Were all subjects scored by the same number of raters?
  length(unique(ks)) == 1
}

is_complete <- function(.data, subject, rater) {
  # How many subjects did each rater score?
  nk <- rowSums(table(.data[[rater]], by = .data[[subject]]))

  # How many unique subjects were there?
  n <- length(unique(.data[[subject]]))

  # Did all raters score all subjects?
  all(nk == n)
}

is_twoway <- function(.data, subject, rater) {
  # How many subjects did each rater score?
  nk <- rowSums(table(.data[[rater]], by = .data[[subject]]))

  # Did all raters NOT score just one subject?
  !all(nk == 1)
}

get_terms <- function(.data, subject, rater, k, error) {
  # Is the design balanced or unbalanced?
  bal <- is_balanced(.data, subject, rater)

  # Is the design complete or incomplete?
  com <- is_complete(.data, subject, rater)

  # Is the design two-way/crossed or one-way/nested?
  two <- is_twoway(.data, subject, rater)

  if (two && error == "Relative" && k == 1 && com) {
    "ICC(C,1)"
  } else if (two && error == "Relative" && k == 1 && !com) {
    "ICC(Q,1)"
  } else if (two && error == "Relative" && k > 1 && com) {
    "ICC(C,k)"
  } else if (two && error == "Relative" && k > 1 && !com) {
    "ICC(Q,khat)"
  } else if (two && error == "Absolute" && k == 1) {
    "ICC(A,1)"
  } else if (two && error == "Absolute" && k > 1 && bal) {
    "ICC(A,k)"
  } else if (two && error == "Absolute" && k > 1 && !bal) {
    "ICC(A,khat)"
  } else if (!two && k == 1) {
    "ICC(1)"
  } else if (!two && k > 1 && bal) {
    "ICC(k)"
  } else if (!two && k > 1 && !bal) {
    "ICC(khat)"
  } else {
    NA_character_
  }
}

get_estimates <- function(m, method = ggdist::mode_qi, ci = 0.95) {
  apply(X = m, MARGIN = 2, FUN = method, .width = ci) |>
    dplyr::bind_rows(.id = "term")
}

bname <- function(x) {
  gsub(pattern = "_", replacement = "", x)
}

get_lmer_ints <- function(m) {
  l <- stats::coef(m)
  l2 <- lapply(l, \(x) datawizard::rownames_as_column(x, var = "id"))
  d <- Reduce(rbind, l2)
  # dd <- as.data.frame(lme4::ranef(m))
  #fixef<-
  # dd_ci <- transform(dd, lower = condval - 1.96*condsd, upper = condval + 1.96*condsd)
  out <- data.frame(
    component = rep(names(l), times = sapply(l2, nrow)),
    term = rep("Intercept", times = nrow(d)),
    id = d$id,
    estimate = d$`(Intercept)`,
    #lower = dd_ci$lower,
    #upper = dd_ci$upper
    lower = NA_real_,
    upper = NA_real_
  )
  out
}



# S3 Generics -------------------------------------------------------------

#' @export
calc_khat <- function(x, ...) {
  UseMethod("calc_khat")
}

#' @export

calc_ks<- function(x, ...) {
  UseMethod("calc_ks")
}

#' @export
calc_q <- function(x, ...) {
  UseMethod("calc_q")
}

#' @export
calc_icc <- function(x, ...) {
  UseMethod("calc_icc")
}

new_ks <- function(x, ...) {
  structure(x, ..., class = "varde_ks")
}

new_srm <- function(x, ...) {
  structure(x, ..., class = "varde_srm")
}

# S3 Methods --------------------------------------------------------------

## calc_ks methods

#' @method calc_ks varde_srm 
#' @export
calc_ks.varde_srm <- function(srm) {

  # Calculate ks from subject-by-rater matrix
  ks <- rowSums(srm)

  # Assign the varde_ks class
  ks <- new_ks(ks)

  # Return ks
  ks

}

#' @method calc_ks data.frame
#' @export
calc_ks.data.frame <- function(.data,
                               subject = "subject",
                               rater = "rater",
                               score = "score") {

  # Create subject-by-rater matrix from .data
  srm <- create_srm(.data, subject = subject, rater = rater, score = score)

  # Calculate ks from the subject-by-rater matrix
  ks <- calc_ks(srm)

  # Return ks
  ks

}

## calc_khat methods

#' @method calc_khat varde_ks
#' @export
calc_khat.varde_ks <- function(ks) {

  # Calculate harmonic mean
  khat <- length(ks) / sum(1 / ks)

  # Return khat
  khat

}

#' @method calc_khat varde_srm
#' @export
calc_khat.varde_srm <- function(srm) {

  # Count raters per subject from subject-by-rater matrix
  ks <- calc_ks(srm)

  # Calculate khat from raters per subject
  khat <- calc_khat(ks)

  # Return khat
  khat

}

#' @method calc_khat data.frame
#' @export
calc_khat.data.frame <- function(.data,
                                 subject = "subject",
                                 rater = "rater",
                                 score = "score",
                                 ...) {

  # Create subject-by-rater matrix from .data
  srm <- create_srm(.data, subject = subject, rater = rater, score = score, ...)

  # Calculate khat from subject-by-rater matrix
  khat <- calc_khat(srm)

  # Return khat
  khat

}

## calc_q methods

#' @method calc_q varde_srm
#' @export
calc_q.varde_srm <- function(srm) {

  # How many raters per subject?
  ks <- calc_ks(srm = srm)

  # How many subjects?
  n <- nrow(srm)

  # What is the harmonic mean of raters per subject?
  khat <- calc_khat(ks = ks)

  # Generate all unique pairs of subject indexes
  spairs <- combn(n, 2)

  # Function to calculate the proportion of overlap for a pair of subjects
  pair_overlap <- function(spair, srm) {
    # What is the index of first subject?
    s1 <- spair[[1]]

    # What is the index of second subject?
    s2 <- spair[[2]]

    # How many raters for first subject?
    k_s1 <- sum(srm[s1, ])

    # How many raters for second subject?
    k_s2 <- sum(srm[s2, ])

    # How many raters shared between subjects?
    k_s1s2 <- sum(colSums(srm[c(s1, s2), ]) > 1)

    # What is the proportion of rater overlap for this pair of subjects?
    (2 * k_s1s2) / (k_s1 * k_s2)

    # NOTE: Because we will iterate over unique pairs rather than all pairs...
    # NOTE: We double the numerator to capture both orderings (A-B and B-A)
    # NOTE: This saves a little time by halving the number of iterations needed
  }

  # Apply function to all unique pairs of subjects and sum across pairs
  total_overlap <- sum(apply(X = spairs, MARGIN = 2, FUN = pair_overlap, srm))

  # Calculate the proportion of non-overlap across subjects and raters
  q <- (1 / khat) - (total_overlap / (n * (n - 1)))

  # Return q
  q

}

#' @method calc_q data.frame
#' @export
calc_q.data.frame <- function(.data,
                              subject = "subject",
                              rater = "rater",
                              score = "score",
                              ...) {

  # Create subject-by-rater matrix from .data
  srm <- create_srm(.data, subject = subject, rater = rater, score = score, ...)

  # Calculate q from subject-by-rater matrix
  q <- calc_q(srm = srm)

  # Return q
  q

}

#' @export print.varde_icc
#' @export
print.varde_icc <- function(x, variances = TRUE, intercepts = TRUE, ...) {
  cat(crayon::blue("# iccs_samples Estimates\n"))
  print(x$iccs_summary, ...)
  if (variances) {
    cat(crayon::blue("\n# Variance Estimates\n"))
    print(x$vars_summary, ...)
  }
  if (intercepts) {
    cat(crayon::blue("\n# Intercept Estimates\n"))
    print(x$ints_summary, ...)
  }
}

#' @export summary.varde_icc
#' @export
summary.varde_icc <- function(x,
                              which = "iccs_samples",
                              ...) {

  match.arg(which, choices = c("iccs_samples", "variances", "intercepts", "model"))
  if (which == "iccs_samples") {
    out <- x$iccs_summary
  } else if (which == "variances") {
    out <- x$vars_summary
  } else if (which == "intercepts") {
    out <- x$ints_summary
  } else if (which == "model") {
    out <- summary(x$model, ...)
  }
  out
}


# calc_icc() --------------------------------------------------------------



#' @param .data dataframe
#' @param subject char
#' @param rater char
#' @param scores Int
#' @param k int (number of raters for ICC(A,1))
#' @return variances/effects
#' @export
calc_vardel_icc <- function(.data,
  subject = "ObjectID",
  rater = "RaterID",
  scores = "Score",
  k = NULL
 # ci = 0.95
  #varde = matrix(),
  ){

  assertthat::assert_that(
    rlang::is_null(k) || rlang::is_integerish(k, n = 1)
  )

  # assertthat::assert_that(
  #   rlang::is_double(ci, n = 1, finite = TRUE),
  #   ci > 0, ci < 1
  # )
  # assertthat::assert_that(
  #   rlang::is_integerish(chains, n = 1, finite = TRUE),
  #   chains >= 1
  # )


  # How many score variables were provided?
  v <- length(scores)

  # Create logical subject-rater matrices
  srm <- lapply(
    X = scores,
    FUN = create_srm,
    .data = .data,
    subject = subject,
    rater = rater
  )
  names(srm) <- scores

  # Count the number of raters who scored each subject
  ks <- lapply(X = srm, FUN = rowSums)

  # Count the number of subjects scored by each rater
  nk <- lapply(X = srm, FUN = colSums)

  # # Remove all subjects that had no raters
  keep <- lapply(ks, function(x) names(x[x > 0])) |>
    unlist() |>
    unique()
  .data <- .data[.data[[subject]] %in% keep, ]

  # # Remove all raters that had no subjects
  keep <- lapply(nk, function(x) names(x[x > 0])) |>
    unlist() |>
    unique()
  .data <- .data[.data[[rater]] %in% keep, ]

  #check design of data
  balanced <- is_balanced(.data, subject,rater)
  complete <- is_complete(.data, subject,rater)
  twoway <- is_twoway(.data, subject, rater)


  # If not specified, set k as the number of unique raters
  if (is.null(k)) {
    k <- length(unique(.data[[rater]]))
  }

  # Calculate the harmonic mean of the number of raters per subject
  khat <- lapply(srm, calc_khat)

  # Calculate the proportion of non-overlap for raters and subjects
  q <- lapply(srm, calc_q)

  if(twoway == FALSE){
      q <- 1/k #since sigmaR cannot be distinguished
    } 

  # Construct mixed-effects formula
  formula <- create_parse(.data, subject, rater, scores, v)
  

  # model_fit <- lme4::lmer(formula = formula,
  #                   data = .data,
  #                 REML=TRUE)
  
  safe_lme <- purrr::quietly(purrr::possibly(lme4::lmer, otherwise = NULL))

  model_fit <- safe_lme(
    formula = formula,
    data = .data 
  )


  if ( is.null(model_fit$result ) ) {
  # we had an error!
    iccs = NULL
    message = length(model_fit$message) > 0
    warning = length(model_fit$warning) > 0
    error = TRUE
  } else {
    #  icc_est <- computeICC_random(fit, subject, k, khat, q, v)
    #iccs_CIs <- ICC_CIs_LME(fit, subject, rater, ci, k, khat, q)
    #obtain confidence intervals of ICCs
    ## TODO: What if one way model fit?

    # In vardel, we will not calculate CIs by simulation:
    # instead we will use ggdist() to estimate them
    #level <- ci

    #index lme4 model
    model_fitted <- model_fit$result

    khat <- khat$Score
    Q <- q$Score
    #two way random effects
    lme_vars <- lme4::VarCorr(model_fitted)
    vs <- lme_vars[[subject]][1] #obtain object name

    #get not specified random effects variances
    ran_eff <- attr(model_fitted@flist,"names")
    ran_eff <- ran_eff[ran_eff != subject]
    vr <- lme_vars[[ran_eff]][1]

    #residual/interaction variances
    #vsr <- (attr(lme4::VarCorr(model_fit), "sc"))^2
    vsr <- sigma(model_fitted)^2 


  
    #only interested in ICC(A,1)
    #iccs <- signif((vs / (vs + vr + vsr)), digits = 3)
    iccs <- c(
    vs / (vs + vr + vsr),
    vs / (vs + (vr + vsr) / khat),
    vs / (vs + vsr),
    vs / (vs + vsr / khat)
    )
    iccs <- sapply(iccs, signif, digits = 3
    )

    
    
    

    message = length(model_fit$message)>0 
    warning = length(model_fit$warning)>0  
    error = FALSE

  }
  
  icc_names <- c(
    "ICC(A,1)", "ICC(A,k)",
    "ICC(C,1)", "ICC(C,k)"
  )
  
  #out <- iccs #ICC(A,1)
  #seedNum <-.data$Seed[1]
  subject_var <- .data$OBJ_VAR[1]
  rater_var <- .data$RATER_VAR[1]

  out <- tibble::tibble(
    method = "icc",
    icc = icc_names,
    estimate = iccs,
    sigma_s = subject_var,
    sigma_r = rater_var,
    sigma_vsr = 1,  
    vs = vs,
    vr = vr,
    vsr = vsr,   
    message = as.character(message),
    warning = as.character(warning), 
    error = as.character(error)
    #SeedNum = seedNum
  )
  #attr(out, "seed")

  return(out)
}



#' @param .data dataframe
#' @param subject char
#' @param rater char
#' @param scores Int 
#' @param k int (number of raters for ICC(A,1))
#' @return variances/effects
#' @export
calc_g_binary_icc <- function(.data,
  subject = "ObjectID",
  rater = "RaterID",
  scores = "Score",
  k = NULL
 # ci = 0.95
  #varde = matrix(),
  ){

  assertthat::assert_that(
    rlang::is_null(k) || rlang::is_integerish(k, n = 1)
  )

  # assertthat::assert_that(
  #   rlang::is_double(ci, n = 1, finite = TRUE),
  #   ci > 0, ci < 1
  # )
  # assertthat::assert_that(
  #   rlang::is_integerish(chains, n = 1, finite = TRUE),
  #   chains >= 1
  # )


  # How many score variables were provided?
  v <- length(scores)

  # Create logical subject-rater matrices
  srm <- lapply(
    X = scores,
    FUN = create_srm,
    .data = .data,
    subject = subject,
    rater = rater
  )
  names(srm) <- scores

  # Count the number of raters who scored each subject
  ks <- lapply(X = srm, FUN = rowSums)

  # Count the number of subjects scored by each rater
  nk <- lapply(X = srm, FUN = colSums)

  # # Remove all subjects that had no raters
  keep <- lapply(ks, function(x) names(x[x > 0])) |>
    unlist() |>
    unique()
  .data <- .data[.data[[subject]] %in% keep, ]

  # # Remove all raters that had no subjects
  keep <- lapply(nk, function(x) names(x[x > 0])) |>
    unlist() |>
    unique()
  .data <- .data[.data[[rater]] %in% keep, ]

  #check design of data
  balanced <- is_balanced(.data, subject,rater)
  complete <- is_complete(.data, subject,rater)
  twoway <- is_twoway(.data, subject, rater)


  # If not specified, set k as the number of unique raters
  if (is.null(k)) {
    k <- length(unique(.data[[rater]]))
  }

  # Calculate the harmonic mean of the number of raters per subject
  khat <- lapply(srm, calc_khat)

  # Calculate the proportion of non-overlap for raters and subjects
  q <- lapply(srm, calc_q)

  if(twoway == FALSE){
      q <- 1/k #since sigmaR cannot be distinguished
    } 

  # Construct mixed-effects formula
  formula <- create_parse_glmmtmb(.data, subject, rater, scores, v)
  

  # model_fit <- lme4::lmer(formula = formula,
  #                   data = .data,
  #                 REML=TRUE)
  
  safe_glmmTMB <- purrr::quietly(purrr::possibly(glmmTMB::glmmTMB, otherwise = NULL))

  model_fit <- safe_glmmTMB(
    formula = formula,
    family = binomial(link = "probit"),
    data = .data 
  )


  if ( is.null(model_fit$result ) ) {
  # we had an error!
    iccs = NULL
    message = length(model_fit$message) > 0
    warning = length(model_fit$warning) > 0
    error = TRUE
  } else {
    #  icc_est <- computeICC_random(fit, subject, k, khat, q, v)
    #iccs_CIs <- ICC_CIs_LME(fit, subject, rater, ci, k, khat, q)
    #obtain confidence intervals of ICCs
    ## TODO: What if one way model fit?

    # In vardel, we will not calculate CIs by simulation:
    # instead we will use ggdist() to estimate them
    #level <- ci

    #index lme4 model
    model_fitted <- model_fit$result

    khat <- khat$Score
    Q <- q$Score
    #two way random effects
    glmmTMB_vars <- glmmTMB::VarCorr(model_fitted)
    #vs <- glmmTMB_vars[[subject]][1] #obtain object name
    vs <- glmmTMB_vars$cond$ObjectID[1]

    #get not specified random effects variances
    #ran_eff <- attr(model_fitted@flist,"names")
    #ran_eff <- ran_eff[ran_eff != subject]
    #vr <- lme_vars[[ran_eff]][1]
    vr <- glmmTMB_vars$cond$RaterID[1]

    #residual/interaction variances
    #vsr <- (attr(lme4::VarCorr(model_fit), "sc"))^2
    vsr <- sigma(model_fitted)^2 


  
    iccs <- c(
    vs / (vs + vr + vsr),
    vs / (vs + (vr + vsr) / khat),
    vs / (vs + vsr),
    vs / (vs + vsr / khat)
    )
    
    iccs <- sapply(iccs, signif, digits = 3
    )
    

    message = length(model_fit$message)>0 
    warning = length(model_fit$warning)>0  
    error = FALSE

  }
  
  icc_names <- c(
    "ICC(A,1)", "ICC(A,k)",
    "ICC(C,1)", "ICC(C,k)"
  )
  
  #out <- iccs #ICC(A,1)
  #seedNum <-.data$Seed[1]

  subject_var <- .data$OBJ_VAR[1]
  rater_var <- .data$RATER_VAR[1]
  resid_vsr <- .data$RES_VAR[1]

  out <- tibble::tibble(
    method = "g_icc",
    icc = icc_names,
    estimate = iccs,
    sigma_s = subject_var,
    sigma_r = rater_var,
    sigma_vsr = resid_vsr,  
    vs = vs,
    vr = vr,
    vsr = vsr,   
    message = as.character(message),
    warning = as.character(warning), 
    error = as.character(error)
    #SeedNum = seedNum
  )
 
  return(out)
}






#' @param .data dataframe
#' @param subject char
#' @param rater char
#' @param scores Int 
#' @param k int (number of raters for ICC(A,1))
#' @return variances/effects
#' @export
calc_g_ordinal_icc <- function(.data,
  subject = "ObjectID",
  rater = "RaterID",
  scores = "Score",
  k = NULL
 # ci = 0.95
  #varde = matrix(),
  ){

  assertthat::assert_that(
    rlang::is_null(k) || rlang::is_integerish(k, n = 1)
  )

  # assertthat::assert_that(
  #   rlang::is_double(ci, n = 1, finite = TRUE),
  #   ci > 0, ci < 1
  # )
  # assertthat::assert_that(
  #   rlang::is_integerish(chains, n = 1, finite = TRUE),
  #   chains >= 1
  # )
  
  #make score as a factor
  .data$Score <- as.factor(.data$Score)

  # How many score variables were provided?
  v <- length(scores)

  # Create logical subject-rater matrices
  srm <- lapply(
    X = scores,
    FUN = create_srm,
    .data = .data,
    subject = subject,
    rater = rater
  )
  names(srm) <- scores

  # Count the number of raters who scored each subject
  ks <- lapply(X = srm, FUN = rowSums)

  # Count the number of subjects scored by each rater
  nk <- lapply(X = srm, FUN = colSums)

  # # Remove all subjects that had no raters
  keep <- lapply(ks, function(x) names(x[x > 0])) |>
    unlist() |>
    unique()
  .data <- .data[.data[[subject]] %in% keep, ]

  # # Remove all raters that had no subjects
  keep <- lapply(nk, function(x) names(x[x > 0])) |>
    unlist() |>
    unique()
  .data <- .data[.data[[rater]] %in% keep, ]

  #check design of data
  balanced <- is_balanced(.data, subject,rater)
  complete <- is_complete(.data, subject,rater)
  twoway <- is_twoway(.data, subject, rater)


  # If not specified, set k as the number of unique raters
  if (is.null(k)) {
    k <- length(unique(.data[[rater]]))
  }

  # Calculate the harmonic mean of the number of raters per subject
  khat <- lapply(srm, calc_khat)

  # Calculate the proportion of non-overlap for raters and subjects
  q <- lapply(srm, calc_q)

  if(twoway == FALSE){
      q <- 1/k #since sigmaR cannot be distinguished
    } 

  # Construct mixed-effects formula
  formula <- create_parse_glmmtmb(.data, subject, rater, scores, v)
  

  # model_fit <- lme4::lmer(formula = formula,
  #                   data = .data,
  #                 REML=TRUE)
  
  safe_ordinal <- purrr::quietly(purrr::possibly(ordinal::clmm, otherwise = NULL))

  model_fit <- safe_ordinal(
    formula = formula,
    link = "probit",
    data = .data 
  )


  if ( is.null(model_fit$result ) ) {
  # we had an error!
    iccs = NULL
    message = length(model_fit$message) > 0
    warning = length(model_fit$warning) > 0
    error = TRUE
  } else {
    #  icc_est <- computeICC_random(fit, subject, k, khat, q, v)
    #iccs_CIs <- ICC_CIs_LME(fit, subject, rater, ci, k, khat, q)
    #obtain confidence intervals of ICCs
    ## TODO: What if one way model fit?

    # In vardel, we will not calculate CIs by simulation:
    # instead we will use ggdist() to estimate them
    #level <- ci

    #index lme4 model
    model_fitted <- model_fit$result

    khat <- khat$Score
    Q <- q$Score
    #two way random effects
    glmmTMB_vars <- glmmTMB::VarCorr(model_fitted)
    #vs <- glmmTMB_vars[[subject]][1] #obtain object name
    vs <- glmmTMB_vars$ObjectID[1]

    #get not specified random effects variances
    #ran_eff <- attr(model_fitted@flist,"names")
    #ran_eff <- ran_eff[ran_eff != subject]
    #vr <- lme_vars[[ran_eff]][1]
    vr <- glmmTMB_vars$RaterID[1]

    #residual/interaction variances
    #vsr <- (attr(lme4::VarCorr(model_fit), "sc"))^2
    #vsr <- sigma(model_fitted)^2 
    vsr <- 1 #constant due to probit link 


  
    #only interested in ICC(A,1) # add all four 
    #iccs <- signif((vs / (vs + vr + vsr)), digits = 3)
    iccs <- c(
    vs / (vs + vr + vsr),
    vs / (vs + (vr + vsr) / khat),
    vs / (vs + vsr),
    vs / (vs + vsr / khat))
    iccs <- sapply(iccs, signif, digits = 3
    )

    message = length(model_fit$message)>0 
    warning = length(model_fit$warning)>0  
    error = FALSE

  }
  
  icc_names <- c(
    "ICC(A,1)", "ICC(A,k)",
    "ICC(C,1)", "ICC(C,k)"
  )
  
  #out <- iccs #ICC(A,1)
  #seedNum <-.data$Seed[1]
  subject_var <- .data$OBJ_VAR[1]
  rater_var <- .data$RATER_VAR[1]

  out <- tibble::tibble(
    method = "g_icc",
    icc = icc_names,
    estimate = iccs,
    sigma_s = subject_var,
    sigma_r = rater_var,
    sigma_vsr = 1,  
    vs = vs,
    vr = vr,
    vsr = vsr,   
    message = as.character(message),
    warning = as.character(warning), 
    error = as.character(error)
    #SeedNum = seedNum
  )
  

  return(out)
}


#' @param .data dataframe
#' @param subject char
#' @param rater char
#' @param scores Int 
#' @param k int (number of raters for ICC(A,1))
#' @return variances/effects
#' @export
calc_aov_icc <- function(.data,
  subject = "ObjectID",
  rater = "RaterID",
  scores = "Score",
  k = NULL
 # ci = 0.95
  #varde = matrix(),
  ){

  assertthat::assert_that(
    rlang::is_null(k) || rlang::is_integerish(k, n = 1)
  )

  # assertthat::assert_that(
  #   rlang::is_double(ci, n = 1, finite = TRUE),
  #   ci > 0, ci < 1
  # )
  # assertthat::assert_that(
  #   rlang::is_integerish(chains, n = 1, finite = TRUE),
  #   chains >= 1
  # )
  
  #make score as a factor
  #.data$Score <- as.factor(.data$Score)

  # How many score variables were provided?
  v <- length(scores)

  # Create logical subject-rater matrices
  srm <- lapply(
    X = scores,
    FUN = create_srm,
    .data = .data,
    subject = subject,
    rater = rater
  )
  names(srm) <- scores

  # Count the number of raters who scored each subject
  ks <- lapply(X = srm, FUN = rowSums)

  # Count the number of subjects scored by each rater
  nk <- lapply(X = srm, FUN = colSums)

  # # Remove all subjects that had no raters
  keep <- lapply(ks, function(x) names(x[x > 0])) |>
    unlist() |>
    unique()
  .data <- .data[.data[[subject]] %in% keep, ]

  # # Remove all raters that had no subjects
  keep <- lapply(nk, function(x) names(x[x > 0])) |>
    unlist() |>
    unique()
  .data <- .data[.data[[rater]] %in% keep, ]

  #check design of data
  balanced <- is_balanced(.data, subject,rater)
  complete <- is_complete(.data, subject,rater)
  twoway <- is_twoway(.data, subject, rater)


  # If not specified, set k as the number of unique raters
  if (is.null(k)) {
    k <- length(unique(.data[[rater]]))
  }

  # Calculate the harmonic mean of the number of raters per subject
  khat <- lapply(srm, calc_khat)

  # Calculate the proportion of non-overlap for raters and subjects
  q <- lapply(srm, calc_q)

  if(twoway == FALSE){
      q <- 1/k #since sigmaR cannot be distinguished
    } 

  # Construct mixed-effects formula
  formula <- create_parseaov(.data, subject, rater, scores, v)
  

  # model_fit <- lme4::lmer(formula = formula,
  #                   data = .data,
  #                 REML=TRUE)

## fit anova model as fixed ("what most people do")
.data$Score <- as.numeric(.data$Score)
.data$ObjectID <- as.factor(.data$ObjectID)
.data$RaterID <- as.factor(.data$RaterID)
  
# attempt to do it by hand  
  # dat <- .data |> # make wide
  #   select(c(ObjectID,RaterID,Score)) |> 
  #   pivot_wider(names_from = RaterID, values_from = Score) |> 
  #   select(!ObjectID)

  # Sum of squared deviations for COLUMNS
#col_rater_ssd <- sum(apply(dat2, 2, function(x) sum((x - mean(x))^2)))

# Sum of squared deviations for ROWS
#row_object_ssd <- sum(apply(dat2, 1, function(x) sum((x - mean(x))^2)))


  
safe_aov <- purrr::quietly(purrr::possibly(lm, otherwise = NULL))
#works for binary and ordinal...

model_fit <- safe_aov(
    formula = formula,
    data = .data
  )


  if ( is.null(model_fit$result ) ) {
  # we had an error!
    iccs = NULL
    message = length(model_fit$message) > 0
    warning = length(model_fit$warning) > 0
    error = TRUE
  } else {
    #  icc_est <- computeICC_random(fit, subject, k, khat, q, v)
    #iccs_CIs <- ICC_CIs_LME(fit, subject, rater, ci, k, khat, q)
    #obtain confidence intervals of ICCs
    ## TODO: What if one way model fit?

    # In vardel, we will not calculate CIs by simulation:
    # instead we will use ggdist() to estimate them
    #level <- ci

    #index lme4 model
    model_fitted <- model_fit$result

    khat <- khat$Score
    Q <- q$Score

    #extract MS rows and columns
    
    aov_mod <- as.data.frame(anova(model_fitted))

    MSr <- aov_mod["RaterID",]$`Mean Sq` #MeanSq Rater
    MSo <- aov_mod["ObjectID",]$`Mean Sq` #MeanSq Object
    MSe <- aov_mod["Residuals",]$`Mean Sq` #MeanSq Error
    n_objects <-  unique(nk[[1]])

    #only interested in ICC(A,1) #McGraw & Wong (1996)
    iccs_a1 <- (MSr-MSe) / (MSr + ((khat-1)*MSe) + 
      ((khat/n_objects)*(MSo - MSe)))

    iccs_ak <- (MSr-MSe) / (MSr + ((MSo - MSe)/khat)) 

    iccs_c1 <- (MSr - MSe) / (MSr + ((khat -1)*MSe))

    iccs_ck <- (MSr - MSe)/(MSr)


    iccs <- c(iccs_a1, iccs_ak, iccs_c1, iccs_ck) 


    message = length(model_fit$message)>0 
    warning = length(model_fit$warning)>0  
    error = FALSE

  }
  
  icc_names <- c(
    "ICC(A,1)", "ICC(A,k)",
    "ICC(C,1)", "ICC(C,k)"
  )
  
  subject_var <- .data$OBJ_VAR[1]
  rater_var <- .data$RATER_VAR[1]
  resid_vsr <- .data$RES_VAR[1]
  
  #out <- iccs #ICC(A,1)
  #seedNum <-.data$Seed[1]

 out <- tibble::tibble(
    method = "aov_icc",
    icc = icc_names,
    estimate = iccs,
    sigma_s = subject_var,
    sigma_r = rater_var,
    sigma_vsr = resid_vsr,  
    vs = MSo,
    vr = MSr,
    vsr = MSe,   
    message = as.character(message),
    warning = as.character(warning), 
    error = as.character(error)
    #SeedNum = seedNum
  )

  return(out)
}