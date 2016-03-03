### Script to estimate treatment effect across many simulated (AB)^k studies. ###

# install.packages("~/Desktop/single_case_docs/scdlintrend_1.0.tar.gz", repos = NULL, type = "source")
# library(scdlintrend)

##----------------------------------------------------------------------------------------------------------------
## Estimate effect size, d, when (sigma^2 + tau^2) is estimated.
## (Issues with the A matrix (individual variance-related matrix), B, and E matrix remain)
## ...This is a work in progress...
## NOTE: this function is for when phi and rho are known.
##----------------------------------------------------------------------------------------------------------------

#' @title estimate_effect_size
#' 
#' @description calculate d, the estimated effect size when rho and phi are KNOWN.
#' 
#' @param data data.frame containing case, time, treatment, phase, and outcome columns
#' @param k number of complete sub-phase cycles, as in (AB)^k
#' @param rho (tau^2)/(tau^2 + sigma^2)
#' @param phi first-order autocorrelation parameter
#'
#' @export 
#' 
#' @return d, the estimate of the effect size, delta_{R}
#' 
#' @examples
#' beta_matrix <- matrix(rep(c(1:8), 3), byrow = T, nrow = 3)
#' sim_data <- simulate_ABk(m = 3, n = 20, k = 2, min_Ni = 10, outcome = T, beta_matrix = beta_matrix)
#' d <- estimate_effect_size(sim_data, 2, 0.5, 0.2)

estimate_effect_size <- function(data, k, rho, phi){

  #Provide initial check of quality of supplied data.
  if(!all(c("case", "time", "treatment", "phase", "outcome") %in% colnames(data))) stop("Data.frame must have case, time, treatment, phase, and outcome columns.")

  m <- max(unique(data$case)) #case vector should be 1:num_cases
  N_is <- sapply(c(1:m), FUN = function(x){sum(as.numeric(data$case == x))}) #number of observations in each case
  N <- sum(N_is)

  A_mat <- matrix(0, nrow = N, ncol = N) #store all m (N_i x N_i) matrices that can be calculated from the var_i function case by case
  B_mat<- matrix(0, nrow = 4*k*m, ncol = N)#store all m (4k x 1) beta vectors per case
  Sigma_mat <- matrix(0, nrow = N, ncol = N) 
  E_mat <- matrix(0, nrow = N, ncol = N) #store all m (N_i x N_i) matrices per case

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

    precision <- solve(sigma);
    xt_prec_x_inv <- solve(t(X)%*%precision%*%X)
    A_mat[start:stop, start:stop] <- diag(N) - precision%*%X%*%xt_prec_x_inv%*%t(X)%*%precision #Note precision%*%X, this is altered from paper.
    B_mat[start:stop,] <- xt_prec_x_inv%*%t(X)%*%precision #Note B_mat's dimensions incompatible with C.
    E_mat[start:stop, start:stop] <- t(c_vec)%*%sigma%*%c_vec%*%A_mat[start:stop, start:stop] #Dimensionality issues.
    Sigma_mat[start:stop, start:stop] <- sigma
  }

  one_m <- as.vector(rep(1, m))
  C_mat <- diag(m) - outer(one_m, one_m)/m
  y <- data$outcome

  T_bold <- B_mat%*%y

  numerator <- sum(T_bold)/m
  denominator <- sqrt(t(y)%*%(t(B_mat)%*%C_mat%*%B_mat - E_mat + A_mat)%*%y) #Dimensionality issues.

  return(list(d = numerator/denominator, B = B_mat, C = C_mat, Sigma = Sigma_mat, A = A_mat, E = E_mat))
}
  

##----------------------------------------------------------------------------------------------------------------
## (Helper Fun) Calculate Var(d) and nu per multi-case study.
##----------------------------------------------------------------------------------------------------------------

nu_and_var_d <- function(delta, B, C, Sigma, A, E){
  m <- dim(C)[1]
  one_m <- as.vector(rep(1, m))
  
  var_d_term_1 <- 1/(t(one_m)%*%t(B)%*%Sigma%*%(B)%*%one_m)
  var_d_term_2_denom <- ((B%*%C%*%t(B)) - E + A)%*%Sigma%*%((B%*%C%*%t(B)) - E + A)%*%Sigma
  var_d <- var_d_term_1 + (delta**2)/(2*sum(diag(var_d_term_2_denom)))

  nu <- sum(diag(var_d_term_2_denom))/(sum(diag(((B%*%C%*%t(B)) - E + A)%*%Sigma))**2)
  
  return(list(var_d = var_d, nu = nu))
}


##----------------------------------------------------------------------------------------------------------------
## (Helper Fun) J, the Hedges' correction factor
##----------------------------------------------------------------------------------------------------------------

J <- function(x) 1-3/(4*x - 1)


##----------------------------------------------------------------------------------------------------------------
## Run simulations given Delta, sigma, tau, and N_i for cases where time profiles are identical
##----------------------------------------------------------------------------------------------------------------

#' @title run_equitime_simulations
#' 
#' @description calculate d and g estimates and variances for N simulations of an m x n multi-SCD studies.
#' 
#' @param N number of simulations
#' @param m number of cases in each multi-SCD study
#' @param n number of time observations for each study
#' @param k number of complete baseline/treatment cycles
#' @param delta true effect size parameter
#' @param sigma within-case variance
#' @param tau inter-case variance
#'
#' @export 
#' 
#' @return d, var_d, g, var_g, d_param: N-vectors related to effect size estimates and variances,
#'                                      and d_param is the true delta parameter
#' 
#' @examples
#' sim_1000 <- run_equitime_simulations()
#' hist(sim_1000$d, col = 'grey80')
#' abline(h = sim_1000$d_param, col = 'red')

run_equitime_simulations <- function(N = 1000, m = 6, n = 20, k = 2, delta = 2, sigma = 1, tau = 1){

  d_stor <- vector(mode = "numeric", length = N)
  var_d_stor <- d_stor
  g_stor <- d_stor
  var_g_stor <- d_stor

  trt_vec_inds <- floor(seq(1, n, length.out = 2*k + 1)); trt_vec <- vector(mode = "character", length = 2*k + 1)
  d_trt_vec <- diff(trt_vec_inds)
  trt_vec <- rep("BASELINE", d_trt_vec[1])
  for(i in 2:length(d_trt_vec)){
    if(i%%2 == 1){
      trt_vec <- c(trt_vec, rep("BASELINE", d_trt_vec[i]))
    } else{
      trt_vec <- c(trt_vec, rep("TREATMENT", d_trt_vec[i]))
    }
  }

  c_vec <- c_vector(trt_vec, k = k, x_is = gather_midpoints(trt_vec, k = k))
  print("Estimating beta given calculated c_vector and given delta parameter...")
  beta <- beta_given_delta(delta, c_vec, sigma = sigma, tau = tau, slope_baseline = 0.5, slope_treatment = 3)$beta
  print("Beta successfully calculated...")
  beta_matrix <- matrix(rep(beta, m), byrow = T, nrow = m)
  rho <- (tau^2)/(tau^2 + sigma^2)

  div_el <- floor(N/10)
  for(i in 1:N){
    if(i%%div_el == 0) print(sprintf("On iteration %d...", i))
    sim_data <- simulate_ABk(m = m, n = n, k = k, outcome = T, beta_matrix = beta_matrix)
    ees <- estimate_effect_size(sim_data, k, rho, phi)
    d_stor[i] <- ees$d

    nu_vd <- nu_and_var_d(delta, ees$B, ees$C, ees$Sigma, ees$A, ees$E)

    var_d_stor[i] <- nu_vd$var_d
    g_stor[i] <- d_stor[i]*J(nu_vd$nu)
    var_g_stor[i] <- nu_vd$var_d*(J(nu_vd)**2)
  }

  return(list(d = d_stor, var_d = var_d_stor, g = g_stor, var_g = var_g_stor, d_param = t(c_vec)%*%beta))
}


