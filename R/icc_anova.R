#' @param .data dataframe
#' @param subject char
#' @param rater char
#' @param score Int
#' @return variances/effects
#' @export
calc_aov_icc <- function(
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
  # Construct mixed-effects formula
  formula <- create_parseaov(.data, subject, rater, score, 1)
  .data$Score <- as.numeric(.data$Score)
  .data$ObjectID <- as.factor(.data$ObjectID)
  .data$RaterID <- as.factor(.data$RaterID)
  safe_aov <- purrr::quietly(purrr::possibly(lm, otherwise = NULL))
  model_fit <- safe_aov(
    formula = formula,
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
    aov_mod <- as.data.frame(
      withCallingHandlers(
        anova(model_fitted),
        warning = function(w) {
          if (grepl("essentially perfect fit", w$message)) {
            invokeRestart("muffleWarning")
          }
        }
      )
    )
    MSr <- aov_mod["RaterID", ]$`Mean Sq`
    MSo <- aov_mod["ObjectID", ]$`Mean Sq`
    MSe <- aov_mod["Residuals", ]$`Mean Sq`
    n_objects <- length(unique(.data[[subject]]))
    # Catch perfect fits to avoid division by zero
    if (isTRUE(all.equal(MSe, 0)) || MSe < 1e-10) {
      iccs <- c(1, 1, 1, 1)
      lower_ci <- rep(NA_real_, 4)
      upper_ci <- rep(NA_real_, 4)
    } else {
      iccs_a1 <- (MSo - MSe) /
        (MSo + ((khat - 1) * MSe) + ((khat / n_objects) * (MSr - MSe)))
      iccs_ak <- (MSo - MSe) / (MSo + ((MSr - MSe) / n_objects))
      iccs_c1 <- (MSo - MSe) / (MSo + ((khat - 1) * MSe))
      iccs_ck <- (MSo - MSe) / MSo
      iccs <- c(iccs_a1, iccs_ak, iccs_c1, iccs_ck)
      df_o <- n_objects - 1
      df_e <- (n_objects - 1) * (khat - 1)
      alpha <- 0.05
      F_c <- MSo / MSe
      FL_c <- F_c / qf(1 - alpha / 2, df_o, df_e)
      FU_c <- F_c / qf(alpha / 2, df_o, df_e)
      icc_c1_lower <- (FL_c - 1) / (FL_c + khat - 1)
      icc_c1_upper <- (FU_c - 1) / (FU_c + khat - 1)
      icc_ck_lower <- 1 - (1 / FL_c)
      icc_ck_upper <- 1 - (1 / FU_c)
      n <- n_objects
      k <- khat
      # ICC(A,1) Approximation
      a_1 <- (k * iccs_a1) / (n * (1 - iccs_a1))
      b_1 <- 1 + (k * iccs_a1 * (n - 1)) / (n * (1 - iccs_a1))
      v_a1 <- (a_1 * MSr + b_1 * MSe)^2 /
        ((a_1 * MSr)^2 / (k - 1) + (b_1 * MSe)^2 / ((n - 1) * (k - 1)))
      F_upper_a1 <- qf(1 - alpha / 2, n - 1, v_a1)
      F_lower_a1 <- qf(alpha / 2, n - 1, v_a1)
      icc_a1_lower <- n *
        (MSo - F_upper_a1 * MSe) /
        (F_upper_a1 * (k * MSr + (k * n - k - n) * MSe) + n * MSo)
      icc_a1_upper <- n *
        (MSo - F_lower_a1 * MSe) /
        (F_lower_a1 * (k * MSr + (k * n - k - n) * MSe) + n * MSo)
      # ICC(A,k) Approximation
      a_k <- iccs_ak / (n * (1 - iccs_ak))
      b_k <- 1 + (iccs_ak * (n - 1)) / (n * (1 - iccs_ak))
      v_ak <- (a_k * MSr + b_k * MSe)^2 /
        ((a_k * MSr)^2 / (k - 1) + (b_k * MSe)^2 / ((n - 1) * (k - 1)))
      F_upper_ak <- qf(1 - alpha / 2, n - 1, v_ak)
      F_lower_ak <- qf(alpha / 2, n - 1, v_ak)
      icc_ak_lower <- n *
        (MSo - F_upper_ak * MSe) /
        (F_upper_ak * (MSr - MSe) + n * MSo)
      icc_ak_upper <- n *
        (MSo - F_lower_ak * MSe) /
        (F_lower_ak * (MSr - MSe) + n * MSo)
      lower_ci <- c(icc_a1_lower, icc_ak_lower, icc_c1_lower, icc_ck_lower)
      upper_ci <- c(icc_a1_upper, icc_ak_upper, icc_c1_upper, icc_ck_upper)
    }
    message <- length(model_fit$message) > 0
    warning <- length(model_fit$warning) > 0
    error <- FALSE
  }
  icc_names <- c("ICC(A,1)", "ICC(A,k)", "ICC(C,1)", "ICC(C,k)")
  subject_var <- .data$OBJ_VAR[1]
  rater_var <- .data$RATER_VAR[1]
  resid_vsr <- .data$RES_VAR[1]
  out <- tibble::tibble(
    method = "aov_icc",
    icc = icc_names,
    estimate = iccs,
    ci_lower = lower_ci,
    ci_upper = upper_ci,
    sigma_s = subject_var,
    sigma_r = rater_var,
    sigma_vsr = resid_vsr,
    vs = (MSo - MSe) / khat,
    vr = (MSr - MSe) / n_objects,
    vsr = MSe,
    message = as.character(message),
    warning = as.character(warning),
    error = as.character(error)
  )
  return(out)
}
