library(purrr)

myalgorithm_simulation <- sienaAlgorithmCreate(projname = 'one.step.condF', cond = FALSE,
                                               useStdInits = FALSE, nsub = 0, n3 = 2, simOnly = TRUE)

forward_simulation <- function(n_samples = 1, thetatilde, effects, dataset) {
  ans <- siena07(myalgorithm_simulation, data = dataset, effects = effects, thetaValues = thetatilde,
                 batch = TRUE, silent = FALSE, verbose = TRUE)
  return(list(sf = ans$sf[1:n_samples,], sf2 = ans$sf2[1:n_samples,,]))
}

RM2 <- function(input_forward_simulation = NULL, which_parameters, n_iterations, 
                theta_start, Sigma_start, fixed_Sigma = FALSE, fixed_theta = FALSE,
                learning_rates, preconditioning, min_Sigma = 0.0001) {
  
  n_samples <- input_forward_simulation$n_samples
  effects   <- input_forward_simulation$effects
  dataset   <- input_forward_simulation$dataset
  
  # positions parameters in thetaValues
  which_lambda <- which_parameters$lambda
  which_theta  <- which_parameters$theta
  which_random <- which_parameters$random
  
  if(is.null(which_random)) qqq <- 0 
  else {
    if(is.numeric(which_random)) which_random <- list(which_random) # convert to list
    qqq <- length(which_random) # number variance parameters
  }
  
  epsilon_theta <- learning_rates$theta
  if(qqq>0) epsilon_Sigma <- learning_rates$Sigma
  iterations_decrease_learning_rate <- learning_rates$iterations_decrease_learning_rate
  number_iterations_average <- learning_rates$number_iterations_average
  stopifnot(length(iterations_decrease_learning_rate) == length(number_iterations_average))
  reduce_gradient_factor <- learning_rates$reduce_gradient_factor
  reduced <- 1
  
  thetatilde <- matrix(NA_real_, 1, length(which_theta) + length(unlist(which_random)))
  thetatilde[,which_theta] <- theta <- theta_start
  if(qqq>0) Sigma <- Sigma_start
  
  chains_parameters <- chains_statistics <- list(theta = NULL, Sigma = NULL)
  chains_parameters$theta <- chains_statistics$theta <- matrix(NA_real_, length(which_theta), n_iterations)
  chains_parameters$Sigma <- chains_statistics$Sigma <- matrix(NA_real_, qqq, n_iterations)
  
  for(i in 1:n_iterations) {
    if(i %% 50 == 0) print(paste0("iteration ", i,"/",n_iterations))
    
    # draw random parameters:
    #for(j in seq_along(which_random)) thetatilde[,which_random[[j]]] <- sqrt(Sigma[j])*rnorm(length(which_random[[j]]))
    for(j in seq_along(which_random)) thetatilde[,which_random[[j]]] <- Sigma[j]*rnorm(length(which_random[[j]]))
    
    # simulate network:
    invisible(capture.output(sim <- forward_simulation(n_samples, thetatilde, effects, dataset)))
    
    # store simulated statistics theta
    chains_statistics$theta[,i] <- sim$sf2[which_theta]
    if(!fixed_theta) { # update theta
      gradient_theta = sim$sf[which_theta]
      theta = theta - epsilon_theta * apply(preconditioning, 1, function(x) sum(x*gradient_theta))
      if(!is.null(which_lambda)) theta[which_lambda] <- map_dbl(theta[which_lambda], ~ max(.x,0.01))
    }
    # store parameter and change thetatilde for next iterations
    thetatilde[,which_theta] <- chains_parameters$theta[,i] <- theta
    
    if(i == 1) { # compute observed statistic for variance
      observed_variance  <- rep(0, qqq)
      simulated_variance <- rep(0, qqq)
      for(j in 1:qqq) observed_variance[j] = var(sim$sf2[which_random[[j]]] - sim$sf[which_random[[j]]]) 
    }
    # store simulated statistics variance
    for(j in 1:qqq) chains_statistics$Sigma[j,i] <- simulated_variance[j] <- var(sim$sf2[which_random[[j]]]) 
    if(!fixed_Sigma) { # update Sigma
      gradient_Sigma = simulated_variance - observed_variance
      Sigma = map_dbl(Sigma - epsilon_Sigma*gradient_Sigma, ~ max(.x,min_Sigma))
    }
    # store parameters
    chains_parameters$Sigma[,i] <- Sigma
    
    if(i %in% iterations_decrease_learning_rate) {
      epsilon_theta <- epsilon_theta * reduce_gradient_factor
      thetatilde[,which_theta] <- theta <- rowMeans(chains_parameters$theta[,i-(0:(number_iterations_average[reduced]-1))])
      if(qqq > 0) {
        epsilon_Sigma <- epsilon_Sigma * reduce_gradient_factor
        Sigma <- rowMeans(chains_parameters$Sigma[,i-(0:(number_iterations_average[reduced]-1)), drop=F])
      } 
      reduced <- reduced + 1
    }
    
    if(max(Sigma)>1000) stop(paste("Sigma too big, iteration",i))
    if(max(abs(theta))>1000) stop(paste("max(abs(theta)) too big, iteration",i))
    
  }
  
  return(list(parameters = chains_parameters, statistics = chains_statistics))
}
