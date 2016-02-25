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

simulate_treatment_effects <- function(N = 100, m = 10, n = 30, k = 2, phi = 0.2, tau = 1, sigma = 1,
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
## Given an analytic value of the effect size (delta_{R}), known sigma and tau values, N_i,
## find a Beta vector that yields about this effect size. Assumes that the slope coefficients are
## constant across subphase type.
##----------------------------------------------------------------------------------------------------------------

#' @title beta_given_delta
#' 
#' @description derive a beta vector for an (AB)^k SCD study with desired treatment effect equal to delta. Uses
#'              gradient descent for convergence
#' 
#' @param delta the standardized mean difference effect size
#' @param c_vector vector of multipliers for beta vector
#' @param sigma the true within-case variance
#' @param tau the true between-case variance
#' @param tol absolute difference between true and estimated delta required for convergence
#' @param slope_baseline linear trend during baseline subphase
#' @param slope_treatment linear trend during treatment subphase
#'
#' @export 
#' 
#' @return list containing an estimated beta vector and loss values obtained during convergence

beta_given_delta <- function(delta, c_vector, sigma = 1, tau = 1, tol = 10e-6, slope_baseline = 0.5, slope_treatment = 1){
  
  k = length(c_vector)/4
  beta_tmp <- vector(mode = "numeric", length = 4*k)
  var_scale <- (1/(2*k - 1))*(1/sqrt(tau^2 + sigma^2))
  lr <- 10e-4

  beta_tmp[seq(2, 4*k, by = 4)] <- slope_baseline #if parallel time trends, make slopes during treatment and baseline equal to 1
  beta_tmp[seq(4, 4*k, by = 4)] <- slope_treatment
  beta_tmp[seq(1, 4*k, by = 2)] <- 0
  estimate.delta <- t(c_vector)%*%beta_tmp*var_scale
  loss_stor <- 0.5*(delta-(estimate.delta))^2
  itr_ctr <- 1

  #run GD until beta gives a te such that estimated.delta converged to delta.
  while(abs(delta - estimate.delta) > tol){
    for(i in seq(1, 4*k, by = 2)){
      beta_tmp[i] <- beta_tmp[i] - lr*c_vector[i]*var_scale*(estimate.delta-delta) #Gradient Descent step; add lr*(partial_derivative(loss) wrt beta_i)
    }
    estimate.delta <- t(c_vector)%*%beta_tmp*var_scale
    loss_stor <- c(loss_stor, 0.5*(delta-(estimate.delta))^2)
    itr_ctr <- itr_ctr+1
  }

  return(list(beta = beta_tmp, loss = loss_stor))
}


##----------------------------------------------------------------------------------------------------------------
## Estimate effect size, d, when (sigma^2 + tau^2) is estimated.
## (Issues with the A matrix (individual variance-related matrix) and E matrix remain)
## ...This is a work in progress...
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
      start = cumsum(N_is[1:i])[i-1]+1; stop = cumsum(N_is[1:i])[i]
    }

    trt_vec <- data[data$case == i, 'treatment']
    y_i <- data[data$case == i, "outcome"]
    times <- data[data$case == i, "time"]
    # sigma <- AR1_matrix(phi, rho, times) #what to do because we don't know rho? Need to make function to estimate phi.
    X<- design_matrix_ABk(trt_vec, k = k)
    c_vec <- c_vector(trt_vec, k)

    A_mat[start:stop, start:stop] <- var_i(X, y_i, sigma, k)
    B_mat[i, start:stop] <- beta_vector(X, y_i, sigma)
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
  






