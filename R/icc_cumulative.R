#' @param .data dataframe
#' @param subject char
#' @param rater char
#' @param score Int
#' @return variances/effects
#' @export
calc_g_ordinal_icc <- function(
  .data,
  subject = "ObjectID",
  rater = "RaterID",
  score = "Score",
  srm = NULL
) {
  if (is.null(srm)) {
    srm <- create_srm(
      .data = .data,
      subject = subject,
      rater = rater,
      score = score
    )
  }
  # Count the number of raters who scored each subject
  ks <- rowSums(srm)
  # Count the number of subjects scored by each rater
  nk <- colSums(srm)
  # Remove all subjects that had no raters
  keep_subj <- names(ks[ks > 0])
  .data <- .data[.data[[subject]] %in% keep_subj, ]
  # Remove all raters that had no subjects
  keep_rater <- names(nk[nk > 0])
  .data <- .data[.data[[rater]] %in% keep_rater, ]
  # Coerce scores to a factor
  .data[[score]] <- as.factor(.data[[score]])
  # Calculate the harmonic mean of the number of raters per subject
  khat <- calc_khat(srm)
  # Construct mixed-effects formula
  formula <- create_parse_glmmtmb(.data, subject, rater, score, 1)
  safe_ordinal <- purrr::quietly(purrr::possibly(
    ordinal::clmm,
    otherwise = NULL
  ))
  model_fit <- safe_ordinal(
    formula = formula,
    link = "probit",
    data = .data
  )
  if (is.null(model_fit$result)) {
    # we had an error!
    iccs <- rep(NA_real_, 4)
    vs <- NA_real_
    vr <- NA_real_
    vsr <- NA_real_
    lower_ci <- rep(NA_real_, 4)
    upper_ci <- rep(NA_real_, 4)
    message <- length(model_fit$message) > 0
    warning <- length(model_fit$warning) > 0
    error <- TRUE
  } else {
    model_fitted <- model_fit$result
    ord_vars <- ordinal::VarCorr(model_fitted)
    vs <- as.numeric(ord_vars[["ObjectID"]])
    vr <- as.numeric(ord_vars[["RaterID"]])
    vsr <- 1
    iccs <- c(
      vs / (vs + vr + vsr),
      vs / (vs + (vr + vsr) / khat),
      vs / (vs + vsr),
      vs / (vs + vsr / khat)
    )
    full_hessian <- model_fitted$Hessian
    full_vcov <- tryCatch(solve(full_hessian), error = function(e) NULL)
    se_a1 <- se_ak <- se_c1 <- se_ck <- NA_real_
    if (!is.null(full_vcov)) {
      st_idx <- grep("^ST", rownames(full_hessian))
      if (length(st_idx) == 2) {
        st_vcov <- full_vcov[st_idx, st_idx, drop = FALSE]
        var_corr <- ordinal::VarCorr(model_fitted)
        st_vals <- c(
          attr(var_corr[["ObjectID"]], "stddev"),
          attr(var_corr[["RaterID"]], "stddev")
        )
        form_ak <- as.formula(paste0("~ x1^2 / (x1^2 + ((x2^2 + 1) / ", khat, "))"))
        form_ck <- as.formula(paste0("~ x1^2 / (x1^2 + (1 / ", khat, "))"))
        se_a1 <- tryCatch(msm::deltamethod(~ x1^2 / (x1^2 + x2^2 + 1), st_vals, st_vcov), error = function(e) NA)
        se_ak <- tryCatch(msm::deltamethod(form_ak, st_vals, st_vcov), error = function(e) NA)
        se_c1 <- tryCatch(msm::deltamethod(~ x1^2 / (x1^2 + 1), st_vals, st_vcov), error = function(e) NA)
        se_ck <- tryCatch(msm::deltamethod(form_ck, st_vals, st_vcov), error = function(e) NA)
      }
    }
    se_vector <- c(se_a1, se_ak, se_c1, se_ck)
    safe_iccs <- pmin(pmax(iccs, 1e-5), 1 - 1e-5)
    logit_iccs <- qlogis(safe_iccs)
    se_logit <- se_vector / (safe_iccs * (1 - safe_iccs))
    lower_ci <- plogis(logit_iccs - (1.96 * se_logit))
    upper_ci <- plogis(logit_iccs + (1.96 * se_logit))
    message <- length(model_fit$message) > 0
    warning <- length(model_fit$warning) > 0
    error <- FALSE
  }
  icc_names <- c("ICC(A,1)", "ICC(A,k)", "ICC(C,1)", "ICC(C,k)")
  subject_var <- .data$OBJ_VAR[1]
  rater_var <- .data$RATER_VAR[1]
  out <- tibble::tibble(
    method = "g_icc",
    icc = icc_names,
    estimate = iccs,
    ci_lower = lower_ci,
    ci_upper = upper_ci,
    sigma_s = subject_var,
    sigma_r = rater_var,
    sigma_vsr = 1,
    vs = vs,
    vr = vr,
    vsr = vsr,
    message = as.character(message),
    warning = as.character(warning),
    error = as.character(error)
  )
  return(out)
}
