### First stab at simulating SCD study data ###

##----------------------------------------------------------------
## (Helper fun) Generate number of time observations (N_i, or N_a dep. on notation) per case
##----------------------------------------------------------------

# Get N_i counts: number of observations for each of the m cases.
# m is number of cases
# max_Ni is maximum number of observations across all cases
# min_Ni is minimum number of observations across all cases
#
# returns vector of length m of lengths of observational study between min_Ni and max_Ni
gather_Nis <- function(m, max_Ni = 12, min_Ni = 6){
  
  if(min_Ni > max_Ni) stop("Cannot have min_Ni greater than max_Ni")
  
  possible_lengths <- c(min_Ni:max_Ni)
  
  Nis <- vector(length = m, mode = "numeric")
  Nis <- sapply(Nis, FUN = function(x){tmp <- x
                                       while(tmp == 0){
                                         tmp <- rpois(1, max_Ni)
                                         tmp <- ifelse(!(tmp %in% possible_lengths), 0, tmp)}
                                       tmp
  })
  return(Nis) 
}


##----------------------------------------------------------------
## Generate AB^k data with variable number of observations per case
##----------------------------------------------------------------

#' @title simulate_AB
#' 
#' @description Simulates data from an $AB^{k}$ test with linear trends in baseline and treatment sub-phases.
#' 
#' @param m number of cases
#' @param n maximum number of observations in a case
#' @param k number of phases, i.e. complete AB cycles. There are 2k sub-phases.
#' @param phi autocorrelation parameter 
#' @param min_Ni minimum allowable number of time observations in a case.
#' @param betas (m, 4k)-dimensional matrix containing ($\beta_{0}$, $\beta_{1}$) parameters per every sub-phase. If NULL, this has defaults.
#' @param intercept_lambda poisson distribution parameter for picking baseline $\beta$ intercept.
#' @param sigmas vector of inter-case residual variances
#'   
#' @export 
#' 
#' @return data.frame with case, time, treatment, phase, and outcome columns
#' 
#' @examples
#' simulate_AB(m = 6, n = 12, k = 2, phi = 0.15, intercept_lambda = 15)

simulate_AB <- function(m = 6, n = 12, k = 1,
                        phi = 0.2,
                        min_Ni = NULL,
                        betas = NULL,
                        intercept_lambda = 10,
                        sigmas = NULL){
  
  if(is.null(min_Ni)) min_Ni <- ceiling(n/2)
  
  #investigate this...
  if((2*k - 1) >= min_Ni) stop("Too many phases given smallest number of case observations, min_Ni.")
  
  N_i <- gather_Nis(m, n, min_Ni) #gather number of obs in each case
    
  #get m x (2k+1) matrix indicating baseline/treatment introduction times.
  introduction_time_mat <- t(sapply(N_i, FUN = function(x){floor(seq(1, x, length.out = (2*k) + 1))})) #indicate (start,end] of treatment/baseline sub-phases.
  zero_start_time <- introduction_time_mat; zero_start_time[, 1] <- rep(0, nrow(zero_start_time)) 

  #get m x k matrix indicating ends of phase, i.e. k
  phase_matrix <- introduction_time_mat[,2:ncol(introduction_time_mat)]
  phase_matrix <- as.matrix(phase_matrix[,seq(2,ncol(phase_matrix), by = 2)], ncol = k)
  
  #check dimensions...
  if(!((nrow(introduction_time_mat)==m) & (ncol(introduction_time_mat) == 2*k + 1))) stop("Dimensions of baseline/treatment intro time matrix are incorrect.")
  
  # if(is.null(betas)){   #we want m*(2k) mean beta parameters and m*2k slope beta parameters
  #   betas <- matrix(NA, nrow = m, ncol = 4*k) #there are 4mk parameters to fit...
  
  #   for(ii in 1:m){
  #     betas[ii, seq(1, ncol(betas), by = 2*k)] <- rpois(1, intercept_lambda) #baseline intercept ~10
  #     betas[ii, seq(2, ncol(betas), by = 2*k)] <- rnorm(1, 5, 0.001) #baseline slope ~0.5
  #     betas[ii, seq(3, ncol(betas), by = 2*k)] <- rpois(1, intercept_lambda) #treatment intercept ~10
  #     betas[ii, seq(4, ncol(betas), by = 2*k)] <- rnorm(1, -10, 0.001) #treatment slope ~-1.5
  #   }
  # }
  if(is.null(betas)){
    betas <- matrix(NA, nrow = m, ncol = 4*k) #For simulation, assume that slope params are shared across baseline sub-phases, and same for treatment sub-phases
                                              #although intercept parameters will need to be adjusted for time so that baseline/treatment sub-phases look consistent.
    colnames(betas) <- rep(c("intercept", "slope"), times = 2*k)
    baseline_intercept <- 40; treatment_intercept <- 10
    baseline_slope <- 5; treatment_slope <- -10                                        
    for(ii in 1:m){
      betas[ii, seq(1, ncol(betas), by = 2*k)] <- (baseline_intercept+baseline_slope)-(zero_start_time[ii,seq(1,2*k, by = 2)]+1)*baseline_slope #time offset to make baseline sub-phases look consistent.
      betas[ii, seq(2, ncol(betas), by = 2*k)] <- baseline_slope 
      betas[ii, seq(3, ncol(betas), by = 2*k)] <- (treatment_intercept+treatment_slope)-(zero_start_time[ii,seq(1,2*k, by = 2)]+1)*treatment_slope
      betas[ii, seq(4, ncol(betas), by = 2*k)] <- treatment_slope
    }
  }

  if(is.null(sigmas) | length(sigmas) != m) sigmas <- rep(1, m) #gather inter-case residual variances
  
  #check to see beta matrix is completely filled out...
  if(any(is.na(betas))) stop("Beta matrix was not filled appropriately, NA values remain.")
  
  #fill out timeseries data with AR(1) process data + linear trends according to beta matrix.
  simulation_df <- data.frame(matrix(0, nrow = sum(N_i), ncol = 5))
  eij <- list()
  colnames(simulation_df) <- c("case", "time", "treatment", "phase", "outcome")
  
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
    
    #accrue time spent relative to each sub-phase to avoid having to offset intercept coefficients... Should I do this?
    # normalized_time <- c(unlist(sapply(diff(zero_start_time[case,]), FUN = function(x){1:x})))
    
    #add in intercept based on sub-phase
    # intercept_vec <- ifelse(simulation_df[simulation_df$case==case,"treatment"] == "BASELINE", betas[case,1],betas[case,3])
    intercept_vec <- NULL
    agg_time_in_subphase <- aggregate(simulation_df[simulation_df$case==case,"outcome"], by = list(simulation_df[simulation_df$case==case,"treatment"], simulation_df[simulation_df$case==case,"phase"]), function(x) length(x))$x
    beta_col_selector <- seq(1, 4*k, by = 2)
    for(jj in 1:(2*k)){
        intercept_vec <- as.numeric(c(intercept_vec, rep(betas[case, beta_col_selector[jj]],agg_time_in_subphase[jj])))
    }

    #add in time trend based on sub-phase
    slope_vec <- as.numeric(ifelse(simulation_df[simulation_df$case==case,"treatment"] == "BASELINE", betas[case,seq(2, 4*k, by = 4)], betas[case,seq(4, 4*k, by = 4)]))
    # time_trend <- slope_vec * normalized_time
    time_trend <- slope_vec*simulation_df[simulation_df$case==case,"time"]
    
    #First-order auto-regressive process with time trend: 
    # x1 = arima.sim(list(ar=.4), n=100) + 1:100 #-1<ar<1 controls how fast autocorrelation dies
    eij[[case]] <- arima.sim(list(ar=phi), n=N_i[case], sd = sigmas[case])
    simulation_df$outcome[start:end] <- eij[[case]] + time_trend + intercept_vec
  }

  simulation_df$phase <- as.factor(simulation_df$phase)
  simulation_df$case <- as.factor(simulation_df$case)

  return(list(df = simulation_df, eij = eij, betas = betas, N_i = N_i))
}


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
#' sim <- simulate_AB(m = 6, n = 12, k = 2, phi = 0.15, intercept_lambda = 15)
#' visualize_sim(sim$df)

visualize_sim <- function(df, treatment_name = "TREATMENT"){
  m <- length(unique(df$case))
  dt <- as.data.table(df)
  par(mfrow = c(m,1))
  par(mar = rep(1.8, 4))
  for(i in c(1:m)){
    tmp <- dt[case == i]
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
#' sim <- simulate_AB(m = 6, n = 12, k = 2, phi = 0.15, intercept_lambda = 15)
#' visualize_case_acf(sim$df)

visualize_case_acf <- function(df, case = 1){
  outcome <- df[df$case == case,"outcome"]
  acf(outcome)
}


##----------------------------------------------------------------
## TEST SIMULATION FUNCS
##----------------------------------------------------------------
# sim <- simulate_AB(m = 6, n = 20, k = 2,
#                    phi = 0.2,
#                    min_Ni = NULL,
#                    betas = NULL,
#                    intercept_lambda = 10)
# par.original <- par()
# visualize_sim(sim$df)
# par(par.original)
# visualize_case_acf(sim$df)

