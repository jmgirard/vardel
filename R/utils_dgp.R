#' @param n_raters Int
#' @param n_objects Int
#' @param target_icc 0-1 inclusive
#' @param fixed_obj_var default 1
#' @param rater_resid_ratio default 0.2
#' @return variances/effects
#' @export
generate_data_ROR<- function(n_raters,n_objects,
  target_icc, ROR = 0.006, icc_type = 1){

  # Obtain variances 
  sigma_sqr <- get_data_ROR(target_icc, ROR, n_raters, icc_type)

  # sigma_sqr <- list(
  #   var_object = 30,
  #   var_rater = 0.2,
  #   var_residual = 1.0
  # )
  
  # Generate effects 
  obj_effects   <- rnorm(
    n_objects, mean = 0, sd = sqrt(sigma_sqr$var_object)
  )
  rater_effects <- rnorm(
    n_raters, mean =0, sd = sqrt(sigma_sqr$var_rater)
  )

  # Create data structure # fully crossed designs only
  dat <- expand.grid(
    ObjectID = 1:n_objects,
    RaterID= 1:n_raters
  )
 # Generate scores 
  data <- dat |> tibble::tibble() |>
    dplyr::mutate(
      u_i = obj_effects[ObjectID], #object effect
      v_j = rater_effects[RaterID], #rater effect
      Error = rnorm(dplyr::n(), mean=0, sd = sqrt(sigma_sqr$var_residual)),
      OBJ_VAR = sigma_sqr$var_object,
      RATER_VAR = sigma_sqr$var_rater, 
      RES_VAR = sigma_sqr$var_residual, 
      #implicate grand mean such that
      # y_ij = mu + O_i + R_j + e_ij

     # Score = grand_mean + Effect_O + Effect_R + Error
    )

  return(data)
}

# ' @param target_icc Int
# ' @param fixed_obj_var default to 1
# ' @param rater_resid_ratio default 0.2
# ' @return variances/effects
# ' @export
get_data_ORE <- function(target_icc, ROR) {

  #error check on ICC value 
    if(target_icc >= 1 | target_icc <=0) {
      stop("ICC must be between 0 and 1")
    }
  
  
  #calculate noise(combined rater and residual variance)
  # TODO: CHANGE RATIO to OBJECT/RATER VARIANCE

 # total_noise <- fixed_obj_var * ((1/target_icc) - 1) 
  
  #total_var <- (target_icc/(1-target_icc)) + 1

  #rho = sig2_obj/ (sig2_obj + sig2_rater + 1)
  #sig2_obj = (rho * (sig2_rater+1)) / (1-rho)

  var_object = (target_icc * (rater_resid_ratio + 1))/ (1-target_icc)
  
  #obtain rater variance as a function of ratio  
 # var_residual <- total_noise / (1 + rater_resid_ratio)

 # var_rater <- total_noise - var_residual

  #return variances 
  out <- list(
    var_object = var_object, 
    var_rater = rater_resid_ratio, 
    var_residual = 1)
  
  return(out)
}

# This attempt at setting variance assumes fixed error variance of 1
# and a ratio of object to ratio variance (ROR)
# Default assumes no rater variance (perfect ICC)

#' @param target_icc Int
#' @param ROR rater object ratio: default 0.2 (i.e. Object Variance 5x greater than Rater Var)
#' @return variances/effects
#' @export
get_data_ROR <- function(target_icc, ROR = 0.2, n_raters = 1, icc_type = 1) {

  #error check on ICC value 
    if(target_icc >= 1 | target_icc <=0) {
      stop("ICC must be between 0 and 1")
    }
  
  #Set variance parameter components
  var_error <- 1 # probit residual variance
  #formaula: ICC = O / (O + ROR + E)
  if (icc_type == 1) {
    #utilize ICC(A,1) parameter formula
    #var_object <- (target_icc * var_error) / (1-target_icc * (1 + ROR))
    var_object <- (target_icc * ROR * var_error) / (ROR - (target_icc * (ROR+1)))
  } else if (icc_type == 2){
    #utilize ICC(A,K) parameter formula
    var_object <- (target_icc * ROR * (var_error)) / (n_raters * ROR - target_icc * (n_raters*ROR +1))
  } else if (icc_type == 3){
    #utilize ICC(C,1) parameter formula
    var_object <- (target_icc * var_error) / (1-target_icc)
  } else {
    #utilize ICC(C,K) parameter formula 
    var_object <- (target_icc * var_error) / (n_raters * (1-target_icc))
  }


  if (var_object <= 0) {
  stop("Target ICC and Ratio are mathematically impossible with this error variance.")
}
  var_rater <- var_object * ROR


  #return variances 
  out <- list(
    var_object = var_object, 
    var_rater = var_rater, 
    var_residual = var_error)
  
  return(out)
  }

