get_decay_cuts <- function(k) {
  i <- 1:k
  # Linear decay formula
  p <- (2 * (k - i + 1)) / (k * (k + 1))
  # Cumulative probabilities (excluding the last one which is 1.0)
  cum_p <- cumsum(p)[-k]
  # Convert to Probit Z-scores
  return(qnorm(cum_p))
}


#' @param n_raters Number of Raters
#' @param n_objects Number of Objets
#' @param target_icc ICC
#' @param category Number of categories
#' @param e_category Threshold equality (TRUE) or linear decay (FALSE)
#' @return data
#' @export
simulate_ordinal <- function(
  n_raters,
  n_objects,
  target_icc,
  k_category,
  e_category,
  icc_type
) {
  if (e_category == TRUE) {
    probs <- seq(
      1 / k_category,
      (k_category - 1) / k_category,
      length.out = k_category - 1
    )
    cuts <- qnorm(probs)
  } else {
    cuts <- get_decay_cuts(k_category)
  }
  ORR <- 5
  valid_sample <- FALSE
  while (!valid_sample) {
    dat <- generate_data_ORR(n_raters, n_objects, target_icc, ORR, icc_type)
    total_var <- dat$OBJ_VAR[1] + dat$RATER_VAR[1] + dat$RES_VAR[1]
    scaled_cuts <- cuts * sqrt(total_var)
    dat$eta <- dat$u_i + dat$v_j + dat$Error
    dat$Score <- as.integer(
      cut(dat$eta, breaks = c(-Inf, scaled_cuts, Inf), labels = FALSE)
    )
    if (length(unique(dat$Score)) > 1) {
      valid_sample <- TRUE
    }
  }
  return(dat)
}

ordinal_sim <- function(
  n_raters,
  n_objects,
  target_icc,
  k_category,
  e_category,
  icc_type,
  seed,
  condition,
  filename,
  reps,
  writeFiles
) {
  set.seed(
    seed,
    kind = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
  res <- simhelpers::repeat_and_stack(
    reps,
    {
      dat <- simulate_ordinal(
        n_raters,
        n_objects,
        target_icc,
        k_category,
        e_category,
        icc_type
      )
      combined_mat <- fit_all_models_ordinal(dat)
      combined_mat
    },
    stack = TRUE
  )
  if (writeFiles == TRUE) {
    saveRDS(res, file = file.path(filename))
  }
  return(res)
}

#' @param P Parameter grid
#' @param Iter Int; # of repetitions per condition
#' @return ICCs
#' @export
run_all_ordinal <- function(P, iter, writeFiles) {
  res <- furrr::future_pmap(
    P,
    ordinal_sim,
    reps = iter,
    writeFiles = writeFiles,
    .progress = TRUE,
    .options = furrr::furrr_options(
      seed = NULL,
      packages = c("vardel", "msm", "glmmTMB", "ordinal")
    )
  )
  P$result <- res
  return(P)
}

fit_all_models_ordinal <- function(dat) {
  if (var(dat$Score) == 0) {
    return(dplyr::tibble(
      method = "error",
      error = "ZeroVarianceData",
      icc = NA,
      estimate = NA,
      ci_lower = NA,
      ci_upper = NA
    ))
  }
  srm_matrix <- create_srm(
    dat,
    subject = "ObjectID",
    rater = "RaterID",
    score = "Score"
  )
  res_list <- list(
    t_icc = calc_vardel_icc(dat, srm = srm_matrix),
    g_icc = calc_g_ordinal_icc(dat, srm = srm_matrix),
    aov_icc = calc_aov_icc(dat, srm = srm_matrix),
    caa = cat_vardel_adjusted(dat, weighting = "quadratic")
  )
  return(dplyr::bind_rows(res_list))
}
