#' @param n_raters Number of Raters
#' @param n_objects Number of Objets
#' @param target_icc ICC
#' @param p probability
#' @return data
#' @export
simulate_binary <- function(n_raters, n_objects, target_icc, p, icc_type) {
  ORR <- 5
  valid_sample <- FALSE
  while (!valid_sample) {
    dat <- generate_data_ORR(n_raters, n_objects, target_icc, ORR, icc_type)
    total_var <- dat$OBJ_VAR[1] + dat$RATER_VAR[1] + dat$RES_VAR[1]
    scaled_intercept <- qnorm(p) * sqrt(total_var)
    dat$eta <- scaled_intercept + dat$u_i + dat$v_j + dat$Error
    dat$Score <- as.integer(dat$eta > 0)
    if (length(unique(dat$Score)) > 1) {
      valid_sample <- TRUE
    }
  }
  return(dat)
}

binary_sim <- function(
  n_raters,
  n_objects,
  target_icc,
  p,
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
      dat <- simulate_binary(n_raters, n_objects, target_icc, p, icc_type)
      combined_mat <- fit_all_models_binary(dat)
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
run_all_binary <- function(P, iter, writeFiles) {
  res <- furrr::future_pmap(
    P,
    binary_sim,
    reps = iter,
    writeFiles = writeFiles,
    .progress = TRUE,
    .options = furrr::furrr_options(
      seed = NULL,
      packages = c("vardel", "msm", "glmmTMB")
    )
  )
  P$result <- res
  return(P)
}

fit_all_models_binary <- function(dat) {
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
    g_icc = calc_g_binary_icc(dat, srm = srm_matrix),
    aov_icc = calc_aov_icc(dat, srm = srm_matrix),
    caa = cat_vardel_adjusted(dat, weighting = "quadratic")
  )
  return(dplyr::bind_rows(res_list))
}
