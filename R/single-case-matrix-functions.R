### Functions supporting linear algebra to calculate d and g, biased and unbiased (respectively) estimates for delta, the effect size. ###

##-----------------------------------------------------------------------------------------------
## Calculate x_i's, 1 <= x_i <= 2k-1, the mid points of each sub-phase after the first baseline
##-----------------------------------------------------------------------------------------------

#helper function gather_midpoints
#trt_vec is character vector indicating sub-phase of a particular case, e.g.
#		c("BASELINE", "BASELINE", "TREATMENT", "TREATMENT")
#k is number of complete sub-phase cycles equivalent to number of instances treatment sub-phase kicked in

gather_midpoints <- function(trt_vec, k, treatment_name = "TREATMENT"){
	time <- seq(length(trt_vec))
	trt_numeric <- ifelse(trt_vec == treatment_name, 1, 0)

	diff_trt_numeric <- diff(trt_numeric)
	begin_trt_inds <- which(diff_trt_numeric == 1) + 1 #first instances when treatment introduced
	end_trt_inds <- c(which(diff_trt_numeric == -1), length(trt_vec)) #last instantces of treatment subphases

	#check
	if(length(begin_trt_inds) != length(end_trt_inds)) stop("begin_trt_inds needs to be the same length as end_trt_inds!")

	midpoints <- vector(mode = "numeric", length = (2*k - 1))
	treatment_ctr <- 1

	for(ii in 1:(2*k -1)){ #iterate over sub-phases
		if(ii %% 2 == 1){ #in a treatment sub-phase
			midpoints[ii] <- mean(begin_trt_inds[treatment_ctr]:end_trt_inds[treatment_ctr])
		}
		else{ #in a baseline sub-phase
			midpoints[ii] <- mean((end_trt_inds[treatment_ctr]+1):(begin_trt_inds[treatment_ctr+1]-1))
			treatment_ctr <- treatment_ctr + 1
		}
	}
	return(midpoints)
}


##----------------------------------------------------------------------------------------------------------------
## Create $\vec{c}$ vector required for calculating $\bar{T} = B\vec{y}$, the vector of treatment effect estimates
##----------------------------------------------------------------------------------------------------------------

#' @title c_vector
#' 
#' @description construct c vector containing multipliers and midpoints for multiplying with beta matrix.
#'				user should make one of these per case.
#' 
#' @param trt_vec character vector indicating sub-phase of a particular case, e.g.
#'					c("BASELINE", "BASELINE", "TREATMENT", "TREATMENT")
#' @param k number of complete sub-phase cycles, as in (AB)^k
#' @param x_is (2k-1)x1 vector of times at which we evaluate the treatment effect, perhaps subphase introduction times or midpoints of subphases.
#' @param treatment_name string indicating name of the treatment sub-phase as recorded in df.
#'
#' @export 
#' 
#' @return vector of -1's, 1's, and midpoints
#' 
#' @examples
#' c_1 <- c_vector(sim_df[sim_df$case == 1, "treatment"], k = 2)

c_vector <- function(trt_vec, k, x_is = NULL, treatment_name = "TREATMENT"){

	c_vec <- vector(mode = "numeric", length = 4*k) #storage

	if(is.null(x_is)){
		x_is <- gather_midpoints(trt_vec, k, treatment_name) #midpoints beginning with first treatment sub-phase
	}

	x_i_ctr <- 1

	c_vec[1] <- -1; c_vec[2] <- -x_is[x_i_ctr] #not sure about x_is indexing? Possible notation issue in paper.

	for (j in seq(1, (2*k-2))){
		c_vec[(2*j + 1)] <- 2*((-1)**(j+1))
		c_vec[(2*j + 2)] <- (x_is[x_i_ctr] + x_is[x_i_ctr+1])*((-1)**(j+1))
		x_i_ctr <- x_i_ctr + 1
	}

	c_vec[(length(c_vec)-1)] <- 1
	c_vec[length(c_vec)] <- x_is[length(x_is)]

	return(c_vec)
}


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
## Obtain $\beta_{i}$ vector for a given case
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
	return(solve(t(X)%*%precision%*%X)%*%t(X)%*%precision%*%y) #this beta vector is ordered incorrectly, should be [intercept_1, slope_1, ...]
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

T_bar <- function(beta, c_i, k){
	t(c_i)%*%beta/(2*k - 1) #note, paper is missing this 2*k - 1 denominator factor...
}





##################################################################################################################
############################################## TEST IMPLEMENTATION ###############################################
##################################################################################################################

# phi <- 0.4; rho <- 0.00; k <- 2
# sim <- simulate_ABk(m = 6, n = 20, k = k, phi = phi)
# sim_df <- sim$df
# betas <- sim$betas
# y_1 <- as.matrix(sim_df[sim_df$case == 1, "outcome"])
# trt_vec <- sim_df[sim_df$case == 1, "treatment"]
# X_1 <- design_matrix_ABk(trt_vec, k)
# time_vec <- seq_len(length(y_1))
# sigma_1 <- AR1_matrix(phi, rho, time_vec)
# midpoints <- gather_midpoints(trt_vec, k)
# c_vec <- c_vector(trt_vec, k, midpoints)

# beta_vec <- beta_vector(X_1, y_1, sigma_1)

# print("Estimated linear coefficients:")
# print(t(beta_vec)) #coefficients are out of order...
# print("Actual linear coefficients:")
# print(betas[1, ])

# intercepts <- beta_vec[1:(2*k)]; slopes <- beta_vec[((2*k)+1):(4*k)]
# beta_vec_reorder <- vector(mode = "numeric", length = 4*k)
# for(i in seq(1, 2*k)){
# 	beta_vec_reorder[(2*i)-1] <- intercepts[i]
# 	beta_vec_reorder[2*i] <- slopes[i]
# }

# T_bar_incorrect <- T_bar(c_vec, beta_vec, k)
# T_bar_correct <- T_bar(c_vec, beta_vec_reorder, k)
# print(sprintf("Average treatment effect: %.2f", T_bar_correct))

# sigma_calcs <- var_i(X_1, y_1, sigma_1, k = 2)
# print(sprintf("Estimated variance about the regression line: %.2f", sigma_calcs$sigma_sq_i))
