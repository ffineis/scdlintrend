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
#' @param treatment_name string indicating name of the treatment sub-phase as recorded in df.
#'
#' @export 
#' 
#' @return vector of -1's, 1's, and midpoints
#' 
#' @examples
#' c_1 <- c_vector(sim_df[sim_df$case == 1, "treatment"], k = 2)

c_vector <- function(trt_vec, k, treatment_name = "TREATMENT"){

	c_vec <- vector(mode = "numeric", length = 4*k) #storage
	x_is <- gather_midpoints(trt_vec, k, treatment_name) #midpoints beginning with first treatment sub-phase
	x_i_ctr <- 1

	c_vec[1] <- -1; c_vec[2] <- x_is[x_i_ctr]
	x_i_ctr <- 1 #not sure about this? Possible notation issue in paper.

	for (j in seq(1, (2*k-2), by = 2)){
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

#' @title design_matrix
#' 
#' @description construct design matrix for given a given case
#' 
#' @param trt_vec character vector indicating sub-phase of a particular case, e.g.
#'					c("BASELINE", "BASELINE", "TREATMENT", "TREATMENT")
#' @param k number of complete sub-phase cycles, as in (AB)^k
#' @param treatment_name string indicating name of the treatment sub-phase as recorded in df.
#'
#' @export 
#' 
#' @return a (2k, N_i)-dimensional design matrix
#' 
#' @examples
#' X_1 <- design_matrix(sim_df[sim_df$case == 1, "treatment"], k = 2)




##----------------------------------------------------------------------------------------------------------------
## Create autocorrelation matrix for a given case, Sigma_i
##----------------------------------------------------------------------------------------------------------------

#' @title AR1_matrix
#' 
#' @description construct the autocorrelation matrix for a given case, Sigma_i
#' 
#' @param phi 1-st order autocorrelation
#' @param times a sequence 1:N_i for a given case
#'
#' @export 
#' 
#' @return a (N_i, N_i)-dimensional (square) matrix
#' 
#' @examples
#' Sigma_1 <- AR_1_matrix(sim_df[sim_df$case == 1, "treatment"], k = 2)

AR1_matrix <- function(phi, times) phi^as.matrix(dist(times))


