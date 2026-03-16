#' @param n_raters Number of Raters
#' @param n_objects Number of Objets
#' @param target_icc ICC [0,1]
#' @param p probability
#' @return data
#' @export
simulate_binary <- function(
  n_raters = 30,
  n_objects= 100,
  target_icc = 0.5,
  p=0.5){
  # set seed first
  
  
  # 1) Set binary data hyper-parameters 
  
  #fixed_obj_var <- 1  # assume using rater-residual ratio for now
  intercept <- qnorm(p) # probit transformation 
  #intercept <- p # 50% on the logit scale 

  # Scenario A: Noise is mostly random error (Rater Variance is low)
  # Ratio 0.2 means: Rater Var is 5x less than Object Variance
   ROR <- 0.20

  # 2) First, obtain (object and rater) random effects 
  # must be fully crossed.
 
  dat <- generate_data_ROR(n_raters,n_objects,
    target_icc, ROR)
  
  # 3)) Generate datasets given the binomial model 
  # calculate Linear Predictor and Probabilities
  df <- dat |> tibble::tibble() %>%
    dplyr::mutate(
    # Map random effects to rows
    # u_i = object_effects[Object_ID],
    # v_j = rater_effects[Rater_ID],
    
    # Linear Predictor (probit scale)
    eta = intercept + u_i + v_j + Error, # latent score 
    
    # Probability scale (probit link)
    #prob = pnorm(eta), 
    
    # Generate Binary Rating
    #Score = rbinom(n(), 1, prob) #bernouli (don't do since not stochastic?)

    
    Score = ifelse(eta > 0 , 1, 0) #threshold on latent scale 

    #Score = eta
    #Seed = SEED
  )

  return(df)

}



#custom binary simulation driver 
binary_sim <- function(n_raters, n_objects, target_icc, p, 
  seed,filename, reps, writeFiles){
    #set seed on each iteration
    set.seed(seed, kind = "L'Ecuyer-CMRG", 
    normal.kind = "Inversion", sample.kind = "Rejection") #parallel

    res <- simhelpers::repeat_and_stack(reps, {


      #DGP 
      dat <- simulate_binary(n_raters, n_objects, target_icc, p)

      #Analyze 

      t_icc <- calc_vardel_icc(dat)
      g_icc <- calc_g_icc(dat)
      caa <- cat_vardel_adjusted(dat)

      combined_mat <- c(t_icc,g_icc,caa)


    }, stack = TRUE) 

  
  if(writeFiles == TRUE){
     #w_res <- as.data.frame(res)
     #write_csv(w_res,file = file.path(filename))
    saveRDS(res, file = file.path(filename))
  }
  return(res)
  

}




#' @param P Parameter grid
#' @param Iter Int; # of repetitions per condition
#' @return ICCs
#' @export
run_all_binary <- function(P, iter, writeFiles){
  res <- furrr::future_pmap(P, binary_sim, reps=iter, writeFiles=writeFiles,
      .progress = TRUE,
    .options = furrr::furrr_options(seed = NULL,
   packages = "vardel"))
  
  params$result <- res
  return(params)
}

##########################################
# ORDINAL SIMULATIONS
##########################################


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
#' @param target_icc ICC [0,1]
#' @param category Number of categories
#' @param e_category Threshold equality (TRUE) or linear decay (FALSE)
#' @return data
#' @export
simulate_ordinal <- function(
  n_raters = 30,
  n_objects= 100,
  target_icc = 0.5,
  k_category = 3,
  e_category = TRUE
){
  
  
  # 1) Set ordinal data hyper-parameters 
  # no grand mean intercept due to thresholds

  if (e_category ==TRUE) {
    probs <- seq(1/k_category, (k_category-1)/k_category, length.out = k_category-1)
    cuts <- qnorm(probs)
  } else {
    cuts <- get_decay_cuts(k_category)
  }


  # Scenario A: Noise is mostly random error (Rater Variance is low)
  # Ratio 0.2 means: Rater Var is 5x less than Object Variance
   ROR <- 0.20

  # 2) First, obtain (object and rater) random effects 
  # must be fully crossed.
 
  dat <- generate_data_ROR(n_raters,n_objects,
    target_icc, ROR)
  
  # 3)) Generate datasets given the binomial model 
  # calculate Linear Predictor and Probabilities
  df <- dat |> tibble::tibble () %>%
    dplyr::mutate(
    # Map random effects to rows
    # u_i = object_effects[Object_ID],
    # v_j = rater_effects[Rater_ID],
    
    # Linear Predictor (probit scale)
    eta = u_i + v_j + Error, # latent score 
    
    # Probability scale (probit link)
    #prob = pnorm(eta), 
    
    # Generate ordinal 
    #Score = rbinom(n(), 1, prob) #bernouli (don't do since not stochastic?)

    
    Score = cut(eta, breaks = c(-Inf, cuts, Inf), labels = FALSE) #thresholds on latent scale 

    #Score = eta
    #Seed = SEED
  )

  return(df)

  
  }



#custom binary simulation driver 
ordinal_sim <- function(n_raters, n_objects, target_icc, k_category, 
  e_category, seed,filename, reps, writeFiles){
    #set seed on each iteration
    set.seed(seed, kind = "L'Ecuyer-CMRG", 
    normal.kind = "Inversion", sample.kind = "Rejection") #parallel

    res <- simhelpers::repeat_and_stack(reps, {


      #DGP 
      dat <- simulate_ordinal(n_raters, n_objects, target_icc,
       k_category, e_category)

      #Analyze 

      t_icc <- calc_vardel_icc(dat)
      g_icc <- calc_g_ordinal_icc(dat)
      caa <- cat_vardel_adjusted(dat, weighting = "quadratic")

      combined_mat <- c(t_icc,g_icc,caa)


    }, stack = TRUE) 

  
  if(writeFiles == TRUE){
     #w_res <- as.data.frame(res)
     #write_csv(w_res,file = file.path(filename))
    saveRDS(res, file = file.path(filename))
  }
  return(res)
  

}




#' @param P Parameter grid
#' @param Iter Int; # of repetitions per condition
#' @return ICCs
#' @export
run_all_ordinal <- function(P, iter, writeFiles){
  res <- furrr::future_pmap(P, ordinal_sim, reps=iter, writeFiles=writeFiles,
      .progress = TRUE,
    .options = furrr::furrr_options(seed = NULL,
   packages = "vardel"))
  
  params$result <- res
  return(params)
}