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

#' Calculate Inter-Rater ICC
#'
#' Calculate variance component and inter-rater intraclass correlation estimates
#' using a Bayesian generalizability study.
#'
#' @param .data Either a data frame containing at least the variables identified
#'   in `subject`, `rater`, and `score` or a brmsfit object.
#' @param subject A string indicating the column name in `.data` that contains
#'   an identifier for the subject or thing being scored in each row (e.g.,
#'   person, image, or document). (default = `"subject"`)
#' @param rater A string indicating the column name in `.data` that contains an
#'   identifier for the rater or thing providing the score in each row (e.g.,
#'   rater, judge, or instrument). (default = `"rater"`)
#' @param scores A character vector indicating the column names in `.data` that
#'   contain the numerical scores representing the rating of each row's subject
#'   from that same row's rater (e.g., score, rating, judgment, measurement).
#'   (default = `c("score1", "score2")`)
#' @param k Either `NULL` to set the number of raters you would like to estimate
#'   the reliability of to the total number of unique raters observed in `.data`
#'   or an integer specifying the number of raters you would like to estimate
#'   the reliability of (see details below). (default = `NULL`)
#' @param method A function (ideally from [ggdist::point_interval()]) that
#'   returns a data frame containing a point estimate (`y`) and the lower
#'   (`ymin`) and upper (`ymax`) bounds of an interval estimate. (default =
#'   [ggdist::mode_qi()])
#' @param engine A character vector indicating the choice of estimation
#' framework (`"LME"`, or `"BRMS"`). (default = "`"BRMS"`)
#' @param ci A finite number between 0 and 1 that represents the width of the
#'   credible intervals to estimate (e.g., 0.95 = 95% CI). (default = `0.95`)
#' @param chains An integer representing the number of Markov chains to use in
#'   estimation. Forwarded on to [brms::brm()]. (default = `4`)
#' @param iter An integer representing the total number of interations per chain
#'   (including warmup). Forwarded on to [brms::brm()]. (default = `5000`)
#' @param file Either `NULL` to ignore or a string representing the filename to
#'   save the results to. If a file with that name already exists, the results
#'   will instead be read from that file. (default = `NULL`)
#' @param ... Further arguments passed to [brms::brm()].
#' @return A list object of class "varde_icc" that includes three main elements:
#' * `$iccs_summary`: A [tibble::tibble()] containing summary information about
#'   each ICC estimate.
#' * `$vars_summary`: A [tibble::tibble()] containing summary information about
#'   each variance estimate.
#' * `$ints_summary`: A [tibble::tibble()] containing summary information about
#'   each random intercept estimate.
#' * `$iccs_samples`: A matrix where each row is a single posterior sample and
#'   each column is an ICC estimate.
#' * `$vars_samples`: A matrix where each row is a single posterior sample and
#'   each column is a variance estimate.
#' * `$ints_samples`: A matrix where each row is a single posterior sample and
#'   each column is a random intercept estimate.
#' * `$config`: A list containing the specified `method`, `ci`, `k` values.
#' * `$model`: The brmsfit object created by [brms::brm()] containing the full
#'   results of the Bayesian generalizability study.

calc_icc_old <- function(.data,
                     subject = "ObjectID",
                     rater = "RaterID",
                     scores = "Score",
                     k = NULL,
                     method = ggdist::mode_qi,
                     engine = "BRMS",
                     ci = 0.95,
                     chains = 4,
                     iter = 5000,
                     file = NULL
                     #varde = matrix(),
                     ) {

  assertthat::assert_that(
    rlang::is_null(k) || rlang::is_integerish(k, n = 1)
  )
  assertthat::assert_that(
    rlang::is_double(ci, n = 1, finite = TRUE),
    ci > 0, ci < 1
  )
  assertthat::assert_that(
    rlang::is_integerish(chains, n = 1, finite = TRUE),
    chains >= 1
  )

  # if (!is.null(file)) {
  #   # Add rds extension if needed
  #   if (tools::file_ext(file) != "rds") {
  #     file <- paste0(file, ".rds")
  #   }
  #   # If already exists, read it in
  #   if (file.exists(file)) {
  #     message("Reading results from file")
  #     out <- readRDS(file)
  #     return(out)
  #   }
  # }

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

  # TODO
  ## ADD WARNINGS WHEN DATA IS TRIMMED? ####
  ## fix the output plots for LME objects
  ## Fix the CI of estimates in LME objects?

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

  ## Check type of design
  # Balanced or unbalanced
  ##TODO Recall:
  #two-way complete/incomplete and/or balanced/unbalanced designs
  #one-way is special case of incomplete two-way design (balanced no overlap)
  #tell user which engine


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

  # #check if varde_res is provided
  # if (class(varde)[1] == "varde_res") {
  #
  #   out <- calc_icc(varde, rater, subject, scores, k, khat, q, ci, file)
  #   return(out)
  # }

  #Branch engine cases
  

  switch(engine,


         "LME" = {
           # Construct mixed-effects formula
           formula <- create_parse(.data, subject, rater, scores, v)
           
          
           model_fit <- lme4::lmer(formula = formula,
                             data = .data,
                            REML=TRUE)
           

           #res <- varde(fit, subject, rater, ci, k, khat, q) #obtain model summaries

           #res <- varde(fit, ci=ci) #obtain model summaries

           
           if (check_convergence(model_fit)){
             
            
           #  icc_est <- computeICC_random(fit, subject, k, khat, q, v)
           #iccs_CIs <- ICC_CIs_LME(fit, subject, rater, ci, k, khat, q)
           #obtain confidence intervals of ICCs
             ## TODO: What if one way model fit?

          # In vardel, we will not calculate CIs by simulation:
             # instead we will use ggdist() to estimate them
              level <- ci
              khat <- khat$Score
              Q <- q$Score
              #two way random effects
              lme_vars <- lme4::VarCorr(model_fit)
              vs <- lme_vars[[subject]][1] #obtain object name

              #get not specified random effects variances
              ran_eff <- attr(model_fit@flist,"names")
              ran_eff <- ran_eff[ran_eff != subject]
              vr <- lme_vars[[ran_eff]][1]

              #residual/interaction variances
              #vsr <- (attr(lme4::VarCorr(model_fit), "sc"))^2
              vsr <- sigma(model_fit)^2 

              ## List all (and create SDs)
              #sigmas <- c(S_r = S_r, S_s = S_s, S_sr = S_sr)

            
             
              iccs <- cbind(
                vs / (vs + vr + vsr))
                # vs / (vs + (vr + vsr) / k),
                # vs / (vs + (vr + vsr) / khat),
                # vs / (vs + vsr),
                # vs / (vs + (vr )/ k),
                # vs / (vs + Q*vr + vsr/khat)
           
           icc_names <- c("ICC(A,1)")
          #    , "ICC(A,k)", "ICC(A,khat)",
          #    "ICC(C,1)", "ICC(C,k)", "ICC(Q,khat)"
          #  )
           colnames(iccs) <- paste(
             rep(icc_names, each = v),
             colnames(iccs),
             sep = "__"
           )
            
             # Construct ICC output tibble
           #iccs_estimates <- get_estimates(iccs, method = method, ci = ci)


            #  iccs_estimates <- get_estimates(iccs_CIs$Samples, method = method, ci = ci)
            #  iccs_estimates$term <- icc_names

             # for when u use t(icc_est) instead
             # iccs_estimates <- data.frame(
             #   term = icc_names,
             #   y = t(icc_est),
             #   ymin = iccs_CIs$ci.lower,
             #   ymax = iccs_CIs$ci.upper,
             #   .wdith = rep(ci,length(icc_est)),
             #   .point = rep("mean",length(icc_est)),
             #   .interval = rep("CI",length(icc_est))

             #)
             #rownames(iccs_estimates) <- NULL
           #  nam <- c("ICCa1_ci","ICCakhat_ci","ICCak_ci","ICCc1_ci","ICCqkhat_ci","ICCck_ci")
             #iccs <-data.matrix(iccs_CIs$Samples)



           } else {
             stop("Model did not converge")
           }
         },

         "BRMS" = {
           # Construct mixed-effects formula
           formula <- create_parse(.data, subject, rater, scores, v)
           twoway <- is_twoway(.data, subject, rater)
           # Fit Bayesian mixed-effects model
           fit <- brms::brm(
             formula = formula,
             data = .data,
             chains = chains,
             iter = iter,
             init = "random",
             ...
           )

           # Extract posterior draws from model
           res <- varde(fit, ci = ci)

           # Extract posterior draws as matrices
           if (v > 1) {
             vs <- res$vars_samples[, paste(subject, bname(scores), sep = "__")]
             if (twoway) {
               vr <- res$vars_samples[, paste(rater, bname(scores), sep = "__")]
             } else {
               vr <- rep(NA_real_, length(vs))
             }
             vsr <- res$vars_samples[, paste("Residual", bname(scores), sep = "__")]
           } else {
             vs <- res$vars_samples[, subject]
             if (twoway) {
               vr <- res$vars_samples[, rater]
             } else {
               vr <- rep(NA_real_, length(vs))
             }
             vsr <- res$vars_samples[, "Residual"]
           }

           colnames(vs) <- scores
           colnames(vr) <- scores
           colnames(vsr) <- scores

           # # Calculate the harmonic mean of the number of raters per subject
           # khat <- lapply(srm, calc_khat)
           #
           # # Calculate the proportion of non-overlap for raters and subjects
           # q <- lapply(srm, calc_q)

           # Make matrices for k, khat, and q
           kmat <- matrix(rep(k, times = v * nrow(vs)), ncol = v, byrow = TRUE)
           khatmat <- matrix(
             rep(unlist(khat), times = nrow(vs)),
             ncol = v,
             byrow = TRUE
           )

           qmat <- matrix(
             rep(unlist(q), times = nrow(vs)),
             ncol = v,
             byrow = TRUE
           )

           # Calculate posterior for each intraclass correlation coefficient
           iccs <- cbind(
             vs / (vs + vr + vsr),
             vs / (vs + (vr + vsr) / khatmat),
             vs / (vs + (vr + vsr) / kmat),
             vs / (vs + vsr),
             vs / (vs + qmat * vr + vsr / khatmat),
             vs / (vs + vsr / kmat)
           )
           icc_names <- c(
             "ICC(A,1)", "ICC(A,khat)", "ICC(A,k)",
             "ICC(C,1)", "ICC(Q,khat)", "ICC(C,k)"
           )
           colnames(iccs) <- paste(
             rep(icc_names, each = v),
             colnames(iccs),
             sep = "__"
           )

           # Construct ICC output tibble
           iccs_estimates <- get_estimates(iccs, method = method, ci = ci)


         }

         )




  # iccs_summary <-
  #   data.frame(
  #    # term = icc_names,
  #     term = iccs_estimates$term,
  #     estimate = iccs_estimates$y,
  #     lower = iccs_estimates$ymin,
  #     upper = iccs_estimates$ymax,
  #     raters = rep(c(rep(1, v), unlist(khat), rep(k, v)), times = 2),
  #     error = rep(c("Absolute", "Relative"), each = v * 3)
  #   ) |>
  #   tidyr::separate(col = term, into = c("term", "score"), sep = "__") |>
  #   dplyr::relocate(score, .before = 1) |>
  #   dplyr::arrange(score, error, raters)

  # # if (v == 1) {
  # #   colnames(iccs) <- icc_names
  # # }


  # out <-
  #   varde_icc(
  #     iccs_summary = iccs_summary,
  #     vars_summary = res$vars_summary,
  #     ints_summary = res$ints_summary,
  #     iccs_samples = iccs,
  #     vars_samples = res$vars_samples,
  #     ints_samples = res$ints_samples,
  #     config = list(method = method, ci = ci, k = k),
  #     model = fit
  #   )
  
  out <- iccs #ICC(A,1)



  # if (!is.null(file)) {
  #   try(saveRDS(out, file = file), silent = FALSE)
  # }

  return(out)
}

#' @param .data dataframe
#' @param subject char
#' @param rater char
#' @param scores Int [0-inf]
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

    
    icc_names <- c(
    "ICC(A,1)", "ICC(A,k)",
    "ICC(C,1)", "ICC(C,k)"
  )
    
    # icc_names <- c("ICC(A,1)")
      
    # colnames(iccs) <- paste(
    #   rep(icc_names, each = v),
    #   colnames(iccs),
    #   sep = "__"
    # )
    message = length(model_fit$message)>0 
    warning = length(model_fit$warning)>0  
    error = FALSE

  }
  
  
  
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
#' @param scores Int [0-inf]
#' @param k int (number of raters for ICC(A,1))
#' @return variances/effects
#' @export
calc_g_icc <- function(.data,
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


  
    #only interested in ICC(A,1)
    iccs <- signif((vs / (vs + vr + vsr)), digits = 3)
    
    # icc_names <- c("ICC(A,1)")
      
    # colnames(iccs) <- paste(
    #   rep(icc_names, each = v),
    #   colnames(iccs),
    #   sep = "__"
    # )
    message = length(model_fit$message)>0 
    warning = length(model_fit$warning)>0  
    error = FALSE

  }
  
  
  
  #out <- iccs #ICC(A,1)
  #seedNum <-.data$Seed[1]

  out <- list(
    g_icc = iccs,
    g_message = message,
    g_warning = warning, 
    g_error = error
    #SeedNum = seedNum
  )
  #attr(out, "seed")

  return(out)
}






#' @param .data dataframe
#' @param subject char
#' @param rater char
#' @param scores Int [0-inf]
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

    
    icc_names <- c(
    "ICC(A,1)", "ICC(A,k)",
    "ICC(C,1)", "ICC(C,k)"
  )
      
    # colnames(iccs) <- paste(
    #   rep(icc_names, each = v),
    #   colnames(iccs),
    #   sep = "__"
    # )
    message = length(model_fit$message)>0 
    warning = length(model_fit$warning)>0  
    error = FALSE

  }
  
  
  
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
  #attr(out, "seed")

  return(out)
}


#' @param .data dataframe
#' @param subject char
#' @param rater char
#' @param scores Int [0-inf]
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


  
safe_ordinal <- purrr::quietly(purrr::possibly(lm, otherwise = NULL))
#works for binary and ordinal...

model_fit <- safe_ordinal(
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
    iccs <- signif(( (MSr-MSe) / (MSr + ((khat-1)*MSe) + ((khat/n_objects)*(MSo - MSe)))), digits = 3)

    message = length(model_fit$message)>0 
    warning = length(model_fit$warning)>0  
    error = FALSE

  }
  
  
  
  #out <- iccs #ICC(A,1)
  #seedNum <-.data$Seed[1]

  out <- list(
    aov_icc = iccs,
    aov_message = message,
    aov_warning = warning, 
    aov_error = error
    #SeedNum = seedNum
  )
  #attr(out, "seed")

  return(out)
}