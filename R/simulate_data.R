######################################################################################################
################## Simulate multi-individual study given \delta, sigma, tau, #########################
########################### and \vec{beta} consistent with \delta  ###################################
######################################################################################################


##--------------------------------------------------------------------------------------------
## (Helper fun) Generate number of time observations (N_i, or N_a dep. on notation) per case
##--------------------------------------------------------------------------------------------

#' @title gather_Nis
#'
#' @description generate a vector of number of observations for each of the m cases.
#'
#' @param m number of cases
#' @param max_Ni maximum number of observations across all cases
#' @param min_Ni minimum number of observations across all cases.
#'        If min_Ni == max_Ni, all cases have same number of obs.
#'
#' @export 
#'
#' @return vector of length m of lengths of observational study between min_Ni and max_Ni
gather_Nis <- function(m, max_Ni = 12, min_Ni = 12){
  
  if(min_Ni > max_Ni) stop("Cannot have min_Ni greater than max_Ni.")

  if(max_Ni == min_Ni) Nis <- rep(max_Ni, m)
  
  if(max_Ni > min_Ni){
    possible_lengths <- c(min_Ni:max_Ni)
    
    Nis <- vector(length = m, mode = "numeric")
    Nis <- sapply(Nis, FUN = function(x){tmp <- x
                                         while(tmp == 0){
                                           tmp <- rpois(1, max_Ni)
                                           tmp <- ifelse(!(tmp %in% possible_lengths), 0, tmp)}
                                         tmp
    })
  }
  return(Nis) 
}


##-------------------------------------------------------------------------------------------
## Generate AB^k data with potentially variable number of observations per case
## User has option to generate outcome data or not.
##-------------------------------------------------------------------------------------------

#' @title simulate_ABk
#' 
#' @description Simulates data from an $AB^{k}$ study; user has option to simulate obvervational data
#'              when a matrix of beta coefficients has been supplied.
#' 
#' @param m number of cases
#' @param n maximum number of observations in a case
#' @param k number of phases, i.e. complete AB cycles. There are 2k sub-phases.
#' @param min_Ni minimum allowable number of time observations in a case.
#' @param outcome Boolean, does user want simulated outcome data?
#' @param beta_matrix m x 4k matrix with linear beta coefficients
#' @param phi the autocorrelation coefficient
#' @param sigma within-case variance
#' @param tau between-case variance
#'   
#' @export 
#' 
#' @return data.frame with case, time, treatment, phase, and outcome columns
#' 
#' @examples
#' sim_df <- simulate_ABk(m = 6, n = 12, k = 2)

simulate_ABk <- function(m = 6, n = 12, k = 1, min_Ni = NULL,
                         outcome = FALSE, beta_matrix = NULL,
                         phi = 0.2, sigma = 1, tau = 1){
  
  if(is.null(min_Ni)) min_Ni <- n
  
  #investigate this...
  if((2*k - 1) >= min_Ni) stop("Too many phases given smallest number of case observations, min_Ni.")
  
  N_i <- gather_Nis(m, n, min_Ni) #gather number of obs in each case
  # tau_i <- rnorm(n = m, mean = 0, sd = tau) #individual case effects
    
  #get m x (2k+1) matrix indicating baseline/treatment introduction times.
  introduction_time_mat <- t(sapply(N_i, FUN = function(x){floor(seq(1, x, length.out = (2*k) + 1))})) #indicate (start,end] of treatment/baseline sub-phases.
  zero_start_time <- introduction_time_mat; zero_start_time[, 1] <- rep(0, nrow(zero_start_time)) 

  #get m x k matrix indicating ends of phase, i.e. k
  phase_matrix <- introduction_time_mat[,2:ncol(introduction_time_mat)]
  phase_matrix <- as.matrix(phase_matrix[,seq(2,ncol(phase_matrix), by = 2)], ncol = k)
  
  #check dimensions...
  if(!((nrow(introduction_time_mat)==m) & (ncol(introduction_time_mat) == 2*k + 1))) stop("Dimensions of baseline/treatment intro time matrix are incorrect.")

  #storage under !outcome
  simulation_df <- data.frame(matrix(0, nrow = sum(N_i), ncol = 4))
  colnames(simulation_df) <- c("case", "time", "treatment", "phase")

  #storage under outcome == T
  if((outcome) & !(is.null(beta_matrix))){
    if(!all(dim(beta_matrix) == c(m, 4*k))) stop("Dimensions of beta_matrix should be mx4k.")
    sigmas <- rep(sigma, m)
    taus <- rnorm(m, mean = 0, sd = tau)
    simulation_df$outcome <- rep(0, nrow(simulation_df))
    proceed_with_outcomes <- TRUE
  }

  #fill in simulation_df.
  for(case in 1:m){
    
    #row indexing...
    if(case == 1) start <- 1
    if(case > 1) start <- cumsum(N_i[1:case])[case-1] + 1

    end <- cumsum(N_i[1:case])[case]
    
    simulation_df$case[start:end] <- rep(case, N_i[case])
    simulation_df$time[start:end] <- 1:N_i[case]
    
    state <- "BASELINE" #begin study in baseline
    phase <- 1
    intro_break <- 2 #column index of introduction_time_mat corresponding to time when sub-phase (baseline or treatment) ends
    for(ii in 0:(N_i[case]-1)){
      
      simulation_df$treatment[start+ii] <- state
      simulation_df$phase[start+ii] <- phase

      if(((ii+1) == introduction_time_mat[case, intro_break]) & (state == "BASELINE")){
        state <- "TREATMENT"
        intro_break <- intro_break + 1
      }
      if(((ii+1) == introduction_time_mat[case, intro_break]) & (state == "TREATMENT")){
        state <- "BASELINE"
        intro_break <- intro_break + 1
      }
      if((ii+1) == phase_matrix[case, phase]){ #increase the phase, k (i.e. AB^k) count.
        phase <- phase + 1
      }
    }

    #If user wants the outcome and beta_matrix is correct, simulate outcome data.
    if(proceed_with_outcomes){
      tmp <- simulation_df[simulation_df$case == case, ]
      agg_time_in_subphase <- aggregate(1:nrow(tmp), by = list(tmp$treatment, tmp$phase), function(x) length(x))$x
      
      intercept_vec <- as.numeric(rep(beta_matrix[case, seq(1, 4*k, by = 2)], agg_time_in_subphase))
      slope_vec <- as.numeric(rep(beta_matrix[case, seq(2, 4*k, by = 2)], agg_time_in_subphase))
      time_trend <- slope_vec*tmp$time

      err <- as.numeric(arima.sim(list(ar=phi), n=nrow(tmp), sd = sigmas[case]))
      simulation_df$outcome[start:end] <- as.numeric(err + time_trend + intercept_vec + taus[case])
    }
  }

  simulation_df$phase <- as.factor(simulation_df$phase)
  simulation_df$case <- as.factor(simulation_df$case)

  return(simulation_df)
}


##-----------------------------------------------------------------------------------------------
## Calculate x_i's, 1 <= i <= 2k-1, the mid points of each sub-phase after the first baseline
##-----------------------------------------------------------------------------------------------

#' @title gather_midpoints
#'
#' @description from vector of time observation points, determine midpoint of each baseline or treatment subphase
#' @param trt_vec character vector indicating sub-phase of a particular case, e.g.
#           c("BASELINE", "BASELINE", "TREATMENT", "TREATMENT")
#' @param k number of complete baseline, treatment subphases
#' @param treatment_name string contained in simulated data's 'treatment' column
#'
#' @export
#' 
#' @return vector of time midpoints of subphases

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
## Create c vector required for calculating T = t(c)Beta, the average treatment effect for an individual
##----------------------------------------------------------------------------------------------------------------

#' @title c_vector
#' 
#' @description construct c vector containing multipliers and midpoints for multiplying with beta matrix.
#'        user should make one of these per case.
#' 
#' @param trt_vec character vector indicating sub-phase of a particular case, e.g.
#'          c("BASELINE", "BASELINE", "TREATMENT", "TREATMENT")
#' @param k number of complete sub-phase cycles, as in (AB)^k
#' @param x_is (2k-1)x1 vector of times at which we evaluate the treatment effect, perhaps subphase introduction times or midpoints of subphases.
#' @param treatment_name string indicating name of the treatment sub-phase as recorded in df.
#'
#' @export 
#' 
#' @return vector of -1's, 1's, and f(midpoints)
#' 
#' @examples
#' c_1 <- c_vector(sim_df[sim_df$case == 1, "treatment"], k = 2)

c_vector <- function(trt_vec, k, x_is = NULL, treatment_name = "TREATMENT"){

  c_vec <- vector(mode = "numeric", length = 4*k) #storage

  if(is.null(x_is)){
    x_is <- gather_midpoints(trt_vec, k, treatment_name) #midpoints beginning with first treatment sub-phase
  } else{
    if (length(x_is) != (2*k -1)) stop("x_is needs to have 4k values.")
  }

  x_i_ctr <- 1

  c_vec[1] <- -1; c_vec[2] <- -x_is[x_i_ctr] #not sure about x_is indexing? Can i's get indexed along with j in the paper?

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

beta_given_delta <- function(delta, c_vec, sigma = 1, tau = 1, tol = 10e-7, slope_baseline = 0.5, slope_treatment = 1){
  
  k = length(c_vec)/4
  beta_tmp <- vector(mode = "numeric", length = 4*k)
  var_scale <- (1/(2*k - 1))*(1/sqrt(tau^2 + sigma^2))
  lr <- 10e-4

  beta_tmp[seq(2, 4*k, by = 4)] <- slope_baseline #If parallel time trends, make slopes during treatment and baseline equal to 1
  beta_tmp[seq(4, 4*k, by = 4)] <- slope_treatment
  beta_tmp[seq(1, 4*k, by = 2)] <- 0
  estimate.delta <- t(c_vec)%*%beta_tmp*var_scale
  loss_stor <- 0.5*(delta-(estimate.delta))^2
  itr_ctr <- 1

  #run GD until beta gives a te such that estimated.delta converged to delta.
  while(abs(delta - estimate.delta) > tol){
    for(i in seq(1, 4*k, by = 2)){
      beta_tmp[i] <- beta_tmp[i] - lr*c_vec[i]*var_scale*(estimate.delta-delta) #Gradient Descent step; add lr*(partial_derivative(loss) wrt beta_i)
    }
    estimate.delta <- t(c_vec)%*%beta_tmp*var_scale
    loss_stor <- c(loss_stor, 0.5*(delta-(estimate.delta))^2)
    itr_ctr <- itr_ctr+1
  }

  return(list(beta = beta_tmp, loss = loss_stor))
}


##----------------------------------------------------------------------------------------------------------------
## Given a data.frame outlining treatment/baseline subphases, a beta vector to be used across m individual cases, 
## a value for the first-order autocorrelation (phi), within-case and between-case variance (sigma and tau),
## generate outcome vectors for each individual case
##----------------------------------------------------------------------------------------------------------------

#' @title simulate_ABk_outcomes
#' 
#' @description generate outcome vectors for each individual case
#' 
#' @param df data.frame from simulate_ABk_no_outcome
#' @param beta_vec vector of linear coefficients, shared over individuals
#' @param k the number of complete phases
#' @param sigma within-case variance
#' @param tau between-case variance
#' @param treatment_name 
#' @param slope_treatment string indicating name of the treatment sub-phase as recorded in df
#'
#' @export 
#' 
#' @return list containing an estimated beta vector and loss values obtained during convergence

simulate_ABk_outcomes <- function(df, beta_vec, k, phi = 0.2,
                                  sigma = 1, tau = 1, treatment_name = "TREATMENT"){

  if(!all(colnames(df) %in% c("case", "time", "treatment", "phase"))) stop("Column names of df are incorrect.")
  eij <- list()
  df$outcome <- rep(0, nrow(df))

  m <- length(unique(df$case))
  sigmas <- rep(sigma, m)
  taus <- rnorm(m, mean = 0, sd = tau)

  for(i in 1:m){
    case <- df[df$case == i, ]
    agg_time_in_subphase <- aggregate(1:nrow(case), by = list(case$treatment, case$phase), function(x) length(x))$x
    
    intercept_vec <- as.numeric(rep(beta_vec[seq(1, 4*k, by = 2)], agg_time_in_subphase))
    slope_vec <- as.numeric(rep(beta_vec[seq(2, 4*k, by = 2)], agg_time_in_subphase))
    time_trend <- slope_vec*case$time

    eij[[i]] <- as.numeric(arima.sim(list(ar=phi), n=nrow(case), sd = sigmas[i]))
    df[df$case == i, "outcome"] <- as.numeric(eij[[i]] + time_trend + intercept_vec + taus[i])
  }

  return(df)
}


######################################################################################################
######################################## AUXILIARY FUNCTIONS #########################################
######################################################################################################

##----------------------------------------------------------------
## Visualize AB^k data with linear trends
##----------------------------------------------------------------

#' @title visualize_sim
#' 
#' @description plot simulated AB^k data with linear trends.
#' 
#' @param df similar to output of simulate_AB
#' @param treatment_name string indicating name of the treatment sub-phase as recorded in df.
#'
#' @export 
#' 
#' @return NULL
#' 
#' @examples
#' sim <- simulate_AB(m = 6, n = 12, k = 2, phi = 0.15, tau = 1)
#' visualize_sim(sim$df)

visualize_sim <- function(df, treatment_name = "TREATMENT"){
  m <- length(unique(df$case))
  par(mfrow = c(m,1))
  par(mar = rep(1.8, 4))
  for(i in c(1:m)){
    tmp <- df[df$case == i,]
    colors <- ifelse(tmp$treatment == treatment_name, 'red', 'blue')
    plot(1:nrow(tmp), tmp$outcome, pch = 16, col = colors,
         ylab = sprintf("case %d", i),
         xlab = "time",
         cex = 0.9)
  }
}


##----------------------------------------------------------------
## Visualize the ACF of a particular AB^k data case with linear trends
##----------------------------------------------------------------

#' @title visualize_case_acf
#' 
#' @description plot ACF of a particular AB^k test case
#' 
#' @param df similar to output of simulate_AB
#' @param case string of identifier of case in df you wish to plot the ACF for.
#'
#' @export 
#' 
#' @return NULL
#' 
#' @examples
#' sim <- simulate_AB(m = 6, n = 12, k = 2, phi = 0.15, tau = 1)
#' visualize_case_acf(sim$df)

visualize_case_acf <- function(df, case = 1){
  outcome <- df[df$case == case,"outcome"]
  acf(outcome)
}


######################################################################################################
####################################### TEST SIMULATION FUNCS ########################################
######################################################################################################

## Instance where all cases have identical baseline/treatment time profiles
# sim_data <- simulate_ABk_no_outcome(m = 6, n = 20, k = 2) #no outcome data.

# trt_vec_1 <- sim_data[sim_data$case == 1, "treatment"] #representative of all time profiles
# c_vec <- c_vector(trt_vec_1, k = 2, x_is = gather_midpoints(trt_vec_1, k = 2))
# beta_data <- beta_given_delta(3.0, c_vec, sigma = 1, tau = 1, slope_baseline = 0.5, slope_treatment = 3)
# beta <- beta_data$beta

# sim_data <- simulate_ABk_outcomes(sim_data, beta, 2, phi = 0.2, sigma = 1, tau = 1) #append outcome data.
# visualize_sim(sim_data)

## Other instance: cases do not have identical baseline/treatment time profiles
# beta_matrix <- matrix(rep(c(1:8), 3), byrow = T, nrow = 3)
# sim_data <- simulate_ABk(m = 3, n = 20, k = 2, min_Ni = 10, outcome = T, beta_matrix = beta_matrix)
# visualize_sim(sim_data[sim_data$case == 1, ])




