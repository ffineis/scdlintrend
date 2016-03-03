#######################################################################################################################
### Functions supporting linear algebra operations to calculate the treatment effect in the single-individual study ###
#######################################################################################################################


##----------------------------------------------------------------------------------------------------------------
## Create design matrix X for a given case
##----------------------------------------------------------------------------------------------------------------

#' @title design_matrix_ABk
#' 
#' @description construct design matrix for given a given case with k complete baseline-treatment phases.
#' 
#' @param trt_vec character vector indicating sub-phase of a particular case, e.g.
#'					c("BASELINE", "BASELINE", "TREATMENT", "TREATMENT")
#' @param k number of complete sub-phase cycles, as in (AB)^k
#' @param treatment_name string indicating name of the treatment sub-phase as recorded in df.
#'
#' @export 
#' 
#' @return a (N_i, 4k)-dimension design matrix
#' 
#' @examples
#' X_1 <- design_matrix_ABk(sim_df[sim_df$case == 1, "treatment"], k = 2)
design_matrix_ABk <- function(trt_vec, k, treatment_name = "TREATMENT"){
	N_i <- length(trt_vec)
	time <- seq(N_i) #index time beginning at t = 1.
	trt_numeric <- ifelse(trt_vec == treatment_name, 1, 0)

	subphase_ind <- vector(length = N_i, mode = "numeric")
	subphase_ind[1] <- 1
	for(i in 2:N_i) subphase_ind[i] <- ifelse(trt_numeric[i-1] == trt_numeric[i], subphase_ind[i-1], subphase_ind[i-1]+1)

	dummy_matrix <- matrix(NA, nrow = (4*k), ncol = N_i)

	for(ii in 1:(2*k)){
		Dii <- ifelse(subphase_ind==ii, 1, 0)
		dummy_matrix[ii, ] <- Dii
		dummy_matrix[(ii + 2*k),] <- Dii*time  
	}

	return(t(dummy_matrix))
}


##----------------------------------------------------------------------------------------------------------------
## Create autocorrelation matrix for a given case, Sigma_i
##----------------------------------------------------------------------------------------------------------------

#' @title AR1_matrix
#' 
#' @description construct the autocorrelation matrix for a given case, Sigma_i. Original author: Pustejovsky from scdhlm package.
#' 
#' @param phi 1-st order autocorrelation
#' @param rho intraclass variation, i.e. the ratio of intra-person variance to total variance, i.e. (tau^2 / (tau^2 + sigma^2))
#' @param times a sequence 1:N_i for a given case
#'
#' @export 
#' 
#' @return a (N_i, N_i)-dimensional (square) matrix
#' 
#' @examples
#' Sigma_1 <- AR1_matrix(sim_df[sim_df$case == 1, "treatment"], k = 2)

AR1_matrix <- function(phi, rho, times) (1-rho)*(phi^as.matrix(dist(times))) + rho


##----------------------------------------------------------------------------------------------------------------
## Obtain beta vector for a given case
##----------------------------------------------------------------------------------------------------------------

#' @title beta_vector
#' 
#' @description construct the 4k x 1 beta vector corresponding to treatment effect and slope coefficients per subphase
#' 
#' @param X the (4k, N_i)-dimensional design matrix for case i.
#' @param y the vector of N_i outcomes for case i
#' @param Sigma the AR(1) covariance matrix
#'
#' @export 
#' 
#' @return a (4k, 1) vector of beta coefficients
#' 
#' @examples
#' #assuming you've run >sim <- simulate_ABk(...)
#' design_1 <- design_matrix_ABk(sim$df[sim$df$case == 1, "treatment"], k = 2)
#' y_1 <- sim$df[sim$df$case == 1, "outcome"]
#' sigma_1 <- AR1_matrix(phi = 0.4, rho = 0.05, times = 1:14)
#' beta_vec <- beta_vector(design_1, y_1, sigma_1)
#' #NOTE: beta_vec is [intercept_1, intercept_2,...,intercept_k, slope_1, slope_2,...,slope_k], and should match corresp values in sim$betas

beta_vector <- function(X, y, Sigma){
	precision <- solve(Sigma)
	return(solve(t(X)%*%precision%*%X)%*%t(X)%*%precision%*%y) #this beta vector might be ordered improperly, should be [intercept_1, slope_1, ...]
}


##----------------------------------------------------------------------------------------------------------------
## Obtain variance about regression line for case i
##----------------------------------------------------------------------------------------------------------------

#' @title var_i
#' 
#' @description calculate the variance about the regression line for case i
#' 
#' @param X the (4k, N_i)-dimensional design matrix for case i.
#' @param y the vector of N_i outcomes for case i
#' @param Sigma the AR(1) covariance matrix
#' @param k number of complete sub-phase cycles, as in (AB)^k
#'
#' @export 
#' 
#' @return sigma_sq_i estimate, A_i an (N_i, N_i)-dimensional matrix (latter value necessary for combining treatment effects across individuals)
#' 
#' @examples
#' #assuming you've run >sim <- simulate_ABk(...)
#' sigma_calcs <- var_i(design_1, y_1, sigma_1, k = 2)

var_i <- function(X, y, Sigma, k){
	N_i <- length(y)
	precision <- solve(Sigma)
									#probably precision%*%X (N_i x 4k) AND NOT
							 		#precision%*%t(X) (non-conformable)
	A_i <- (1/(N_i-4*k))*(diag(1, N_i) - precision%*%X%*%solve((t(X)%*%precision%*%X))%*%t(X)%*%precision)

	sigma_sq_i <- t(y)%*%A_i%*%y

	return(list(sigma_sq_i = sigma_sq_i, A_i = A_i))
}


##----------------------------------------------------------------------------------------------------------------
## Obtain average treatment effect of case i
##----------------------------------------------------------------------------------------------------------------

#' @title T_bar_individual
#' 
#' @description calculate the mean treatment effect for an individual case
#' 
#' @param beta the (4k, 1) vector of coefficients ordered [intercept1, intercept2,...,interceptk, slope1, slope2,...,slopek]
#' @param c_i the (4k, 1) vector of multipliers and midpoints for counting up beta coefficients
#' @param k number of complete sub-phase cycles, as in (AB)^k
#'
#' @export 
#' 
#' @return the average T estimate for one case

T_bar_individual <- function(beta, c_i, k){
	t(c_i)%*%beta/(2*k - 1) #note, paper could be missing 2*k - 1 denominator factor
}


##################################################################################################################
############################################## TEST IMPLEMENTATION ###############################################
##################################################################################################################

# phi <- 0.4; k <- 2; rho = 0.5 
# sim_df <- simulate_ABk_no_outcome(m = 6, n = 20, k = k) #all cases have identical baseline/treatment time profiles

# trt_vec_1 <- sim_df[sim_df$case == 1, "treatment"] #representative of all time profiles
# c_vec <- c_vector(trt_vec_1, k = k, x_is = gather_midpoints(trt_vec_1, k = k))
# beta_data <- beta_given_delta(3.0, c_vec, sigma = 1, tau = 1, slope_baseline = 0.5, slope_treatment = 3)
# beta <- beta_data$beta

# sim_data <- simulate_ABk_outcomes(sim_df, beta, 2, phi = phi, sigma = 1, tau = 1)
# case <- sim_data$df[sim_data$df$case == 1, ]
# y_1 <- as.matrix(case$outcome)
# trt_vec <- case$treatment
# X_1 <- design_matrix_ABk(trt_vec_1, k)
# time_vec <- case$time
# sigma_1 <- AR1_matrix(phi, rho, time_vec)

# beta_est <- beta_vector(X_1, y_1, sigma_1)

# intercepts <- beta_est[1:(2*k)]; slopes <- beta_est[((2*k)+1):(4*k)]
# beta_est_reorder <- vector(mode = "numeric", length = 4*k)
# for(i in seq(1, 2*k)){
# 	beta_est_reorder[(2*i)-1] <- intercepts[i]
# 	beta_est_reorder[2*i] <- slopes[i]
# }

# print("Estimated linear coefficients:")
# print(t(beta_est)) #coefficients are out of order...
# print("Actual linear coefficients:")
# print(t(beta))
# print("Re-ordered linear coefficients:")
# print(beta_est_reorder)

# T_bar_incorrect <- T_bar_individual(c_vec, beta_vec, k)
# T_bar_correct <- T_bar_individual(c_vec, beta_vec_reorder, k)
# print(sprintf("Average treatment effect: %.2f", T_bar_correct))

# sigma_calcs <- var_i(X_1, y_1, sigma_1, k = 2)
# print(sprintf("Estimated variance about the regression line: %.2f", sigma_calcs$sigma_sq_i))
