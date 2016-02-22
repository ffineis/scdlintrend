### Script to estimate treatment effect across many simulated (AB)^k studies. ###
library(doParallel)
library(foreach)
# install.packages("~/Desktop/single_case_docs/scdlintrend_1.0.tar.gz", repos = NULL, type = "source")
# library(scdlintrend)

##----------------------------------------------------------------------------------------------------------------
## Simulate multiple studies where tau and sigma are known, obtain treatment effects.
##----------------------------------------------------------------------------------------------------------------

#' @title simulate_multi_cases
#' 
#' @description simulate N SCD studies and estimate N*m standardized mean differences.
#' 
#' @param N number of multi-case simulations
#' @param m number of cases within a simulation
#' @param maximum number of observations in a case
#' @param k order of the study, i.e. (AB)^k
#' @param phi first-order autocorrelation
#' @param tau the standard deviation of individual effects
#' @param sigma the inter-case residual variance
#'
#' @export 
#' 
#' @return (N, m) matrix of estimated effect sizes, "d"

simulate_multi_cases <- function(N = 100, m = 10, n = 30, k = 2, phi = 0.2, tau = 1, sigma = 1){
  
  T_stor <- matrix(NA, nrow = N, ncol = m) #simulations x cases matrix
  rho = (tau^2)/(tau^2 + sigma^2)
  
  # library(doMC)
  # registerDoMC(cores=3)
  
#   temp <- foreach(icount(N),
#           .combine = rbind,
# #           .export = "base",
#           packages = "scdlintrend") %dopar% {
  for(sim in 1:N){
    sim_data <- simulate_ABk(m = m, n = n, k = 2,
                             phi = 0.2,
                             tau = 1,
                             min_Ni = NULL,
                             betas = NULL,
                             sigma = sigma)
    
#     temp <- foreach(i = icount(m), .combine = rbind) %dopar% {
    temp <- vector(mode = "numeric", length = m)
    for(i in 1:m){
      trt_vec_i <- sim_data$df[sim_data$df$case == i, "treatment"]
      y_i <- sim_data$df[sim_data$df$case == i, "outcome"]
      X_i <- design_matrix_ABk(trt_vec_i, k)
      time_vec <- seq_len(length(y_i))
      sigma_i <- AR1_matrix(phi, rho, time_vec)
      midpoints <- gather_midpoints(trt_vec_i, k)
      c_vec <- c_vector(trt_vec, k, midpoints)
      beta_vec <- beta_vector(X_i, y_i, sigma_i)
      
      intercepts <- beta_vec[1:(2*k)]; slopes <- beta_vec[((2*k)+1):(4*k)]
      beta_vec_reorder <- vector(mode = "numeric", length = 4*k)
      for(ii in seq(1, 2*k)){
        beta_vec_reorder[(2*ii)-1] <- intercepts[ii]
      	beta_vec_reorder[2*ii] <- slopes[ii]
      }
      temp[i] <- T_bar(c_vec, beta_vec_reorder, k)/sqrt(tau^2 + sigma^2)
    }
#     return(temp)
    T_stor[sim,] <- temp
  }
  return(T_stor)
}


##----------------------------------------------------------------------------------------------------------------
## Estimate effect size when (sigma^2 + tau^2) is estimated
##----------------------------------------------------------------------------------------------------------------

estimate_effect_size <- function(data, k){

  #Provide initial check of quality of supplied data.
  if(!all(c("case", "time", "treatment", "phase", "outcome") %in% colnames(data))) stop("Data.frame must have case, time, treatment, phase, and outcome columns.")

  m <- max(unique(data$case)) #case vector should be 1:num_cases
  N_is <- sapply(c(1:m), FUN = function(x){sum(as.numeric(data$case == x))}) #number of observations in each case

  A_mat <- matrix(0, nrow = nrow(data), ncol = nrow(data)) #store all m (N_i x N_i) matrices that can be calculated from the var_i function case by case
  B_mat<- matrix(0, nrow = m, ncol = nrow(data))#store all m (4k x 1) beta vectors per case
  Sigma_mat <- matrix(0, nrow = nrow(data), ncol = nrow(data)) 
  E_mat <- matrix(0, nrow = nrow(data), ncol = nrow(data)) #store all m (N_i x N_i) matrices per case

  for(i in 1:m){
    if (i == 1){
      start = 1; stop = N_is[i]
      } else{
        start = N_is[i-1]+1; stop = N_is[i]
      }

    trt_vec <- data[data$case == i, 'treatment']
    y_i <- data[data$case == i, "outcome"]
    times <- data[data$case == i, "time"]
    # sigma <- AR1_matrix(phi, rho, times) #what to do because we don't know rho? Need to make function to estimate phi.
    X<- design_matrix_ABk(trt_vec, k = k)
    c_vec <- c_vector(trt_vec, k)

    A_mat[start:stop, start:stop] <- var_i(X, y_i, sigma, k)
    B_mat[start:stop, start:stop] <- beta_vector(X, y_i, sigma)
    E_mat[start:stop, start:stop] <- t(c_vec)%*%sigma%*%c_vec%*%A_mat_list[[i]]
  }

  one_m <- as.vector(rep(1, m))
  C_mat <- diag(m) - outer(one_m, one_m)/m
  y <- data$outcome

  T_bold <- B_mat%*%y

  numerator <- sum(T_bold)/m
  denominator <- sqrt(t(y)%*%(t(B_mat)%*%C_mat%*%B_mat - E_mat + A_mat)%*%y)

  return(numerator/denominator)
}

