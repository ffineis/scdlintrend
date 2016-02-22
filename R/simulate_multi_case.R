### Script to estimate treatment effect across many simulated (AB)^k studies. ###
library(doParallel)
library(foreach)

##----------------------------------------------------------------------------------------------------------------
## Obtain $\beta_{i}$ vector for a given case
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
  
  cl<-makeCluster(3)
  registerDoParallel(cl)
  
  # foreach(sim = 1:N, .combine = rbind, .export = "base") %dopar% {
  for(sim = 1:N){
    sim_data <- simulate_ABk(m = m, n = n, k = 2,
                             phi = 0.2,
                             tau = 1,
                             min_Ni = NULL,
                             betas = NULL,
                             sigma = sigma)
    
#     temp <- foreach(i = 1:m, .combine = rbind) %dopar% {
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
      T_stor[sim, i] <- T_bar(c_vec, beta_vec_reorder, k)/sqrt(tau^2 + sigma^2)
      # T_bar(c_vec, beta_vec_reorder, k)
    }
    
#     stopCluster(cl)
  }
  # stopCluster(cl)
  return(T_stor)
}

