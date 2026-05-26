#' @param .data dataframe
#' @param subject char
#' @param rater char
#' @param score Int
#' @return variances/effects
#' @export
calc_g_binary_icc <- function(
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
  # Calculate the harmonic mean of the number of raters per subject
  khat <- calc_khat(srm)
  formula <- create_parse_glmmtmb(.data, subject, rater, score, 1)
  safe_glmmTMB <- purrr::quietly(purrr::possibly(
    glmmTMB::glmmTMB,
    otherwise = NULL
  ))
  model_fit <- safe_glmmTMB(
    formula = formula,
    family = binomial(link = "probit"),
    REML = TRUE,
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
    glmmTMB_vars <- glmmTMB::VarCorr(model_fitted)
    vs <- glmmTMB_vars$cond$ObjectID[1]
    vr <- glmmTMB_vars$cond$RaterID[1]
    vsr <- 1
    iccs <- c(
      vs / (vs + vr + vsr),
      vs / (vs + (vr + vsr) / khat),
      vs / (vs + vsr),
      vs / (vs + vsr / khat)
    )
    se_a1 <- se_ak <- se_c1 <- se_ck <- NA_real_
    full_vcov <- tryCatch(
      glmmTMB::vcov(model_fitted, full = TRUE),
      error = function(e) NULL
    )
    if (!is.null(full_vcov)) {
      theta_idx <- grep("^theta", names(model_fitted$fit$par))
      if (length(theta_idx) == 2 && all(theta_idx <= ncol(full_vcov))) {
        theta_vcov <- full_vcov[theta_idx, theta_idx, drop = FALSE]
        theta_vals <- model_fitted$fit$par[theta_idx]
        if (
          nrow(theta_vcov) == 2 &&
            ncol(theta_vcov) == 2 &&
            !any(is.na(theta_vcov))
        ) {
          form_ak <- as.formula(paste0(
            "~ exp(x1)^2 / (exp(x1)^2 + ((exp(x2)^2 + 1) / ",
            khat,
            "))"
          ))
          form_ck <- as.formula(paste0(
            "~ exp(x1)^2 / (exp(x1)^2 + (1 / ",
            khat,
            "))"
          ))
          se_a1 <- tryCatch(
            msm::deltamethod(
              ~ exp(x1)^2 / (exp(x1)^2 + exp(x2)^2 + 1),
              theta_vals,
              theta_vcov
            ),
            error = function(e) NA
          )
          se_ak <- tryCatch(
            msm::deltamethod(form_ak, theta_vals, theta_vcov),
            error = function(e) NA
          )
          se_c1 <- tryCatch(
            msm::deltamethod(
              ~ exp(x1)^2 / (exp(x1)^2 + 1),
              theta_vals,
              theta_vcov
            ),
            error = function(e) NA
          )
          se_ck <- tryCatch(
            msm::deltamethod(form_ck, theta_vals, theta_vcov),
            error = function(e) NA
          )
        }
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
  resid_vsr <- .data$RES_VAR[1]
  out <- tibble::tibble(
    method = "g_icc",
    icc = icc_names,
    estimate = iccs,
    ci_lower = lower_ci,
    ci_upper = upper_ci,
    sigma_s = subject_var,
    sigma_r = rater_var,
    sigma_vsr = resid_vsr,
    vs = vs,
    vr = vr,
    vsr = vsr,
    message = as.character(message),
    warning = as.character(warning),
    error = as.character(error)
  )
  return(out)
}
