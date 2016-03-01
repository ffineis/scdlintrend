### Script to estimate treatment effect across many simulated (AB)^k studies. ###
library(doParallel)
library(foreach)
library(ggplot2)
# install.packages("~/Desktop/single_case_docs/scdlintrend_1.0.tar.gz", repos = NULL, type = "source")
# library(scdlintrend)


##----------------------------------------------------------------------------------------------------------------
## Helper function:
## Calculate the analytic expected value of T, (lowercase) delta, and variance, 1_{m}^{t} B^{t} Sigma B 1_{m} (sigma^{2} + tau^{2}).
##----------------------------------------------------------------------------------------------------------------

treatment_effect_distribution <- function(sigma, tau, B, Sigma){
  one_m <- as.vector(rep(1, m))
  # delta <- analytic definition of (lowercase) delta for multiple cases... we only have definition of the single case parameter
  variance <- (tau^2 + sigma^2)*t(one_m)%*%B%*%Sigma%*%t(B)%*%one_m #NOTE: t(B)%*%Sigma%*%B is non-conformable, as B is mxN, and Sigma is NxN.
  return(list(variance = variance)) #also need to return delta
}


##----------------------------------------------------------------------------------------------------------------
## DEPRECATED...
## Simulate multiple studies where tau and sigma are known, obtain treatment effects.
##----------------------------------------------------------------------------------------------------------------

#' @title simulate_treatment_effects
#' 
#' @description simulate N between-subjects studies and estimate N*m estimated treatment effects
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
#' @return (N, m) matrix of estimated treatment effects

simulate_treatment_effects <- function(N = 1000, m = 10, n = 30, k = 2, phi = 0.2, tau = 1, sigma = 1,
                                        plot_treatment_effects = FALSE){
  
  T_stor <- matrix(NA, nrow = N, ncol = m) #(simulations x cases) matrix full of ESTIMATED treatment effects
  delta_stor <- matrix(NA, nrow = N, ncol = m) #(simulations x cases) matrix full ACTUAL treatment effects
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
    
    # temp <- foreach(i = icount(m), .combine = rbind) %dopar% {
    # temp <- vector(mode = "numeric", length = m)

    # B_mat<- matrix(0, nrow = m, ncol = m*4*k)
    # Sigma_mat <- matrix(0, nrow = nrow(sim_data$df), ncol = nrow(sim_data$df))

    for(i in 1:m){

      if (i == 1){
        start = 1; stop = sim_data$N_i[i]
      } else{
        start = cumsum(sim_data$N_i[1:i])[i-1]+1; stop = cumsum(sim_data$N_i[1:i])[i]
      }

      trt_vec_i <- sim_data$df[sim_data$df$case == i, "treatment"]
      y_i <- sim_data$df[sim_data$df$case == i, "outcome"]
      X_i <- design_matrix_ABk(trt_vec_i, k)
      time_vec <- seq_len(length(y_i))
      sigma_i <- AR1_matrix(phi, rho, time_vec)
      midpoints <- gather_midpoints(trt_vec_i, k)
      c_vec <- c_vector(trt_vec, k, midpoints)
      beta_hat_vec <- beta_vector(X_i, y_i, sigma_i)

      # B_mat[i, (4*k*(i-1) + 1):(4*k*i)] <- beta_hat_vec
      # Sigma_mat[start:stop, start:stop] <- sigma_i

      intercepts <- beta_hat_vec[1:(2*k)]; slopes <- beta_hat_vec[((2*k)+1):(4*k)]
      beta_hat_vec_reorder <- vector(mode = "numeric", length = 4*k)
      for(ii in seq(1, 2*k)){
        beta_hat_vec_reorder[(2*ii)-1] <- intercepts[ii]
      	beta_hat_vec_reorder[2*ii] <- slopes[ii]
      }
      T_stor[sim, i] <- T_bar(c_vec, beta_hat_vec_reorder, k) #T_bar is average treatment effect estimate over case i's 2k-1 subphase changes.
      delta_stor[sim, i] <- (1/(2*k - 1))*t(c_vec)%*%sim_data$betas[i, ] #Store the actual treatment effect for each case.
    }
#     return(temp)

    #Get analytic parameters of the distribution of T, the case-averaged treatment effect.
    # analytic_params <- treatment_effect_distribution(sigma, tau, B_betas, Sigma_mat)
  }

  rownames(T_stor) <- paste0("simulation ", seq_len(N))
  colnames(T_stor) <- paste0("case ", seq_len(m))

  if(plot_treatment_effects){
    case_avgs <- data.frame(T = c(apply(T_stor, 1, mean),apply(delta_stor, 1, mean)),
                              statistic = c(rep("estimate", N), rep("actual", N)))
    p <- ggplot(case_avgs, aes(x = T, fill = statistic)) +
            geom_histogram(aes(y = ..density..), #plot density not raw bin count
                            alpha = 0.3, position = "identity") +
            ggtitle("Estimated vs. Actual Treatment Effects")
  }

  return(list(T_stor = T_stor, delta_stor = delta_stor, te_plot = p))
}


##----------------------------------------------------------------------------------------------------------------
## Estimate effect size, d, when (sigma^2 + tau^2) is estimated.
## (Issues with the A matrix (individual variance-related matrix) and E matrix remain)
## ...This is a work in progress...
## NOTE: this function is for when phi and rho are known.
##----------------------------------------------------------------------------------------------------------------

estimate_effect_size <- function(data, k, rho, phi){

  #Provide initial check of quality of supplied data.
  if(!all(c("case", "time", "treatment", "phase", "outcome") %in% colnames(data))) stop("Data.frame must have case, time, treatment, phase, and outcome columns.")

  m <- max(unique(data$case)) #case vector should be 1:num_cases
  N_is <- sapply(c(1:m), FUN = function(x){sum(as.numeric(data$case == x))}) #number of observations in each case

  A_mat <- matrix(0, nrow = nrow(data), ncol = nrow(data)) #store all m (N_i x N_i) matrices that can be calculated from the var_i function case by case
  B_mat<- matrix(0, nrow = 4*k*m, ncol = nrow(data))#store all m (4k x 1) beta vectors per case
  Sigma_mat <- matrix(0, nrow = nrow(data), ncol = nrow(data)) 
  E_mat <- matrix(0, nrow = nrow(data), ncol = nrow(data)) #store all m (N_i x N_i) matrices per case

  for(i in 1:m){
    if (i == 1){
      start = 1; stop = N_is[i]
    } else{
      start = cumsum(N_is[1:i])[i-1]+1; stop = cumsum(N_is[1:i])[i]
    }

    trt_vec <- data[data$case == i, 'treatment']
    y_i <- data[data$case == i, "outcome"]
    times <- data[data$case == i, "time"]
    sigma <- AR1_matrix(phi, rho, times) # This is the case where we do know rho and phi.
    X <- design_matrix_ABk(trt_vec, k = k)
    c_vec <- c_vector(trt_vec, k)

    A_mat[start:stop, start:stop] <- var_i(X, y_i, sigma, k) #clean next three lines up (get rid of functions?) to just compute precision once.
    B_mat[start:stop,] <- solve(t(X)%*%solve(sigma)%*%X)%*%t(X)%*%solve(sigma)
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
  






