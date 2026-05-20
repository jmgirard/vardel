#################################
# LINEAR SIMULATIONS 
#################################

#' @param n_raters Number of Raters
#' @param n_objects Number of Objets
#' @param target_icc ICC 
#' @return data
#' @export
simulate_linear <- function(n_raters, n_objects, target_icc, icc_type) {
  ORR <- 5
  dat <- generate_data_ORR(n_raters, n_objects, target_icc, ORR, icc_type)
  df <- dat |> 
    tibble::as_tibble() |> 
    dplyr::mutate(Score = u_i + v_j + Error)
  return(df)
}

linear_sim <- function(n_raters, n_objects, target_icc, icc_type,
  seed, filename, reps, writeFiles){
    set.seed(
      seed, 
      kind = "L'Ecuyer-CMRG", 
      normal.kind = "Inversion", 
      sample.kind = "Rejection"
    )
    res <- simhelpers::repeat_and_stack(reps, {
      dat <- simulate_linear(n_raters, n_objects, target_icc, icc_type)
      aov_icc <- calc_aov_icc(dat)
    }, stack = TRUE) 
  if(writeFiles == TRUE){
    saveRDS(res, file = file.path(filename))
  }
  return(res)
}

#' @param P Parameter grid
#' @param Iter Int; # of repetitions per condition
#' @return ICCs
#' @export
run_all_AOVlinear <- function(P, iter, writeFiles){
  res <- furrr::future_pmap(P, linear_sim, reps=iter, writeFiles=writeFiles,
      .progress = TRUE,
    .options = furrr::furrr_options(seed = NULL,
   packages = "vardel"))
  
  P$result <- res
  return(P)
}



#################################
# Binary Simulations 
#################################


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
    total_var <_ dat$OBJ_VAR[1] + dat$RATER_VAR[1] + dat$RES_VAR[1]
    scaled_intercept <- qnorm(p) * sqrt(total_var)
    
    df <- dat |>
      tibble::as_tibble() |>
      dplyr::mutate(
        eta = intercept + u_i + v_j + Error,
        Score = dplyr::if_else(eta <= 0, 0, 1)
      )
      
    if (length(unique(df$Score)) > 1) {
      valid_sample <- TRUE
    }
  }
  return(df)
}



#custom binary simulation driver 
binary_sim <- function(n_raters, n_objects, target_icc, p, 
  icc_type, seed, condition, filename, reps, writeFiles){
    set.seed(
      seed, 
      kind = "L'Ecuyer-CMRG", 
      normal.kind = "Inversion", 
      sample.kind = "Rejection"
    )
    res <- simhelpers::repeat_and_stack(reps, {

      dat <- simulate_binary(n_raters, n_objects, target_icc, p, icc_type)


      #Analyze 
      t_icc <- calc_vardel_icc(dat)
      t_icc_glmm <- calc_vardel_icc_glmm(dat)
      g_icc <- calc_g_binary_icc(dat,
        subject = "ObjectID",
        rater = "RaterID",
        scores = "Score",
        icc_type = icc_type)
      caa <- cat_vardel_adjusted(dat, weighting = "quadratic")
      aov_icc <- calc_aov_icc(dat)
      
 

      combined_mat <- dplyr::bind_rows(t_icc,t_icc_glmm,g_icc,caa,aov_icc)


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
  
  P$result <- res
  return(P)
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
){
  
  if (e_category ==TRUE) {
    probs <- seq(1/k_category, (k_category-1)/k_category, length.out = k_category-1)
    cuts <- qnorm(probs)
  } else {
    cuts <- get_decay_cuts(k_category)
  }
  ORR <- 5
  valid_sample <- FALSE
  while(!valid_sample) {
    dat <- generate_data_ORR(n_raters, n_objects, target_icc, ORR, icc_type)
    total_var <- dat$OBJ_VAR[1] + dat$RATER_VAR[1] + dat$RES_VAR[1]
    scaled_cuts <- cuts * sqrt(total_var)
    df <- dat |>
      tibble::as_tibble() |>
      dplyr::mutate(
        eta = u_i + v_j + Error,
        Score = cut(eta, breaks = c(-Inf, scaled_cuts, Inf), labels = FALSE)
      )
    if (length(unique(df$Score)) > 1) {
      valid_sample <- TRUE
    }
  }
  return(df)
}



#custom ordinal simulation driver 
ordinal_sim <- function(n_raters, n_objects, target_icc, k_category, 
  e_category, icc_type, seed, condition, filename, reps, writeFiles){
    #set seed on each iteration
    set.seed(seed, kind = "L'Ecuyer-CMRG", 
    normal.kind = "Inversion", sample.kind = "Rejection") #parallel

    res <- simhelpers::repeat_and_stack(reps, {


      #DGP 
      dat <- simulate_ordinal(n_raters, n_objects, target_icc,
       k_category, e_category, icc_type)

      #Analyze 

      t_icc <- calc_vardel_icc(dat)
      g_icc <- calc_g_ordinal_icc(dat, 
        subject = "ObjectID",
        rater = "RaterID",
        scores = "Score",
      icc_type = icc_type)
      caa <- cat_vardel_adjusted(dat, weighting = "quadratic")
      aov_icc <- calc_aov_icc(dat)

      combined_mat <- dplyr::bind_rows(t_icc,g_icc,caa,aov_icc)


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
  
  P$result <- res
  return(P)
}


ordinal_sim_onlyglmmtmb <- function(n_raters, n_objects, target_icc, k_category, 
  e_category, icc_type, seed, condition, filename, reps, writeFiles){
    #set seed on each iteration
    set.seed(seed, kind = "L'Ecuyer-CMRG", 
    normal.kind = "Inversion", sample.kind = "Rejection") #parallel
  
    res <- replicate(n = reps, expr = {
      #some sims conditions (9 to be exact)
      #have some simulations that give out less columns
      #and rbind() messesup, so use
      #dplyr::list_rbind() instead

      dat <- simulate_ordinal(n_raters, n_objects, target_icc,
       k_category, e_category, icc_type)
      
      t_icc <- calc_vardel_icc_glmm(dat)

      if (!is.tibble(t_icc)){
        t_icc <- tibble::tibble(
                method = "icc_glmmtmb",
               icc = NULL,
               estimate = NULL,
               sigma_s = NULL,
               sigma_r = NULL,
               sigma_vsr = NULL,  
               vs = NULL,
               vr = NULL,
               vsr = NULL,   
              message = NULL,
              warning = NULL, 
              error = TRUE,
    #SeedNum = seedNum
  )}

      
      

      combined_mat <- t_icc
      return(combined_mat)

    }, simplify = FALSE)

    # res <- simhelpers::repeat_and_stack(reps, {


    #   #DGP 
    #   dat <- simulate_ordinal(n_raters, n_objects, target_icc,
    #    k_category, e_category, icc_type)

    #   #Analyze 

    #   t_icc <- calc_vardel_icc_glmm(dat)
    #   # g_icc <- calc_g_ordinal_icc(dat, 
    #   #   subject = "ObjectID",
    #   #   rater = "RaterID",
    #   #   scores = "Score",
    #   # icc_type = icc_type)
    #   # caa <- cat_vardel_adjusted(dat, weighting = "quadratic")
    #   # aov_icc <- calc_aov_icc(dat)

    #   #combined_mat <- dplyr::bind_rows(t_icc,g_icc,caa,aov_icc)
    #   combined_mat <- t_icc


    # }, stack = TRUE) 

  
  res <- purrr::list_rbind(res) # fix 

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
run_all_ordinal_onlyglmmtmb <- function(P, iter, writeFiles){
  res <- furrr::future_pmap(P, ordinal_sim_onlyglmmtmb, reps=iter, writeFiles=writeFiles,
      .progress = TRUE,
    .options = furrr::furrr_options(seed = NULL,
   packages = "vardel"))
  
  P$result <- res
  return(P)
}


##########################################
# ANOVA SIMULATIONS
##########################################

#custom ordinal ANOVA simulation driver 
ordinal_AOV_sim <- function(n_raters, n_objects, target_icc, k_category, 
  e_category, icc_type, seed, condition, filename, reps, writeFiles){
    #set seed on each iteration
    set.seed(seed, kind = "L'Ecuyer-CMRG", 
    normal.kind = "Inversion", sample.kind = "Rejection") #parallel

    res <- simhelpers::repeat_and_stack(reps, {


      #DGP 
      dat <- simulate_ordinal(n_raters, n_objects, target_icc,
       k_category, e_category, icc_type)

      #Analyze 

      aov_icc <- calc_aov_icc(dat)
      # g_icc <- calc_g_ordinal_icc(dat)
      # caa <- cat_vardel_adjusted(dat, weighting = "quadratic")


      #combined_mat <- c(t_icc,g_icc,caa)


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
run_ANOVA_ordinal <- function(P, iter, writeFiles){
  res <- furrr::future_pmap(P, ordinal_AOV_sim, reps=iter, writeFiles=writeFiles,
      .progress = TRUE,
    .options = furrr::furrr_options(seed = NULL,
   packages = "vardel"))
  
  P$result <- res
  return(P)
}



#custom binary simulation driver 
binary_AOV_sim <- function(n_raters, n_objects, target_icc, p, icc_type,
  seed,filename, reps, writeFiles){
    #set seed on each iteration
    set.seed(seed, kind = "L'Ecuyer-CMRG", 
    normal.kind = "Inversion", sample.kind = "Rejection") #parallel

    res <- simhelpers::repeat_and_stack(reps, {


      #DGP 
      dat <- simulate_binary(n_raters, n_objects, target_icc, p, icc_type)

      #Analyze 

      # t_icc <- calc_vardel_icc(dat)
      # g_icc <- calc_g_icc(dat)
      # caa <- cat_vardel_adjusted(dat)

      # combined_mat <- c(t_icc,g_icc,caa)

       aov_icc <- calc_aov_icc(dat)


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
run_ANOVA_binary <- function(P, iter, writeFiles){
  res <- furrr::future_pmap(P, binary_AOV_sim, reps=iter, writeFiles=writeFiles,
      .progress = TRUE,
    .options = furrr::furrr_options(seed = NULL,
   packages = "vardel"))
  
  P$result <- res
  return(P)
}