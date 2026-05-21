#' @param n_raters Int
#' @param n_objects Int
#' @param target_icc 0-1 inclusive
#' @param ORR object-to-rater ratio: As a general rule of thumb: Whenever you want
#' your true signal (object) to be $X$ times larger than your systematic bias (rater),
#' you just set $k_{ratio} = X$.(default 5)
#' @return variances/effects
#' @export
generate_data_ORR <- function(n_raters, n_objects, target_icc, ORR, icc_type) {
  # Obtain variances
  sigma_sqr <- get_data_ORR(target_icc, ORR, n_raters, icc_type)
  # Generate effects
  obj_effects <- rnorm(
    n_objects,
    mean = 0,
    sd = sqrt(sigma_sqr$var_object)
  )
  rater_effects <- rnorm(
    n_raters,
    mean = 0,
    sd = sqrt(sigma_sqr$var_rater)
  )
  # Create data structure # fully crossed designs only
  dat <- expand.grid(
    ObjectID = 1:n_objects,
    RaterID = 1:n_raters
  )
  # Generate scores (avoiding tidyverse for speed)
  dat$u_i <- obj_effects[dat$ObjectID]
  dat$v_j <- rater_effects[dat$RaterID]
  dat$Error <- rnorm(nrow(dat), mean = 0, sd = sqrt(sigma_sqr$var_residual))
  dat$OBJ_VAR <- sigma_sqr$var_object
  dat$RATER_VAR <- sigma_sqr$var_rater
  dat$RES_VAR <- sigma_sqr$var_residual
  return(dat)
}


# This attempt at setting variance assumes fixed error variance of 1
# and a ratio of object to ratio variance (ORR)
# Default assumes no rater variance (perfect ICC)
# As a general rule of thumb: Whenever you want your true
# signal (object) to be $X$ times larger than your systematic bias (rater),
# you just set $k_{ratio} = X$.

#' @param target_icc Int
#' @param ORR rater object ratio: default 5 (i.e. Object Variance 5x greater than Rater Var)
#' @param n_rater number of raters
#' @param icc_type ICC type
#' @return variances/effects
#' @export
get_data_ORR <- function(target_icc, ORR, n_raters, icc_type) {
  #error check on ICC value
  if (target_icc >= 1 | target_icc <= 0) {
    stop("ICC must be between 0 and 1")
  }

  #Set variance parameter components
  var_error <- 1
  if (icc_type == 1) {
    #utilize ICC(A,1) parameter formula
    var_object <- (target_icc * ORR * var_error) /
      (ORR - (target_icc * (ORR + 1)))
  } else if (icc_type == 2) {
    #utilize ICC(A,K) parameter formula
    var_object <- (target_icc * ORR * (var_error)) /
      (n_raters * ORR - target_icc * (n_raters * ORR + 1))
  } else if (icc_type == 3) {
    #utilize ICC(C,1) parameter formula
    var_object <- (target_icc * var_error) / (1 - target_icc)
  } else {
    #utilize ICC(C,K) parameter formula
    var_object <- (target_icc * var_error) / (n_raters * (1 - target_icc))
  }

  if (var_object <= 0) {
    stop(
      "Target ICC and Ratio are mathematically impossible with this error variance."
    )
  }
  var_rater <- var_object / ORR

  #return variances
  out <- list(
    var_object = var_object,
    var_rater = var_rater,
    var_residual = var_error
  )

  return(out)
}
