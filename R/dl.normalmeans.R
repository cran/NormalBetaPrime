###############################################
# Function to implement Gibbs sampler for the #
# Dirichlet-Laplace prior for normal means    #
###############################################

######################
# FUNCTION ARGUMENTS #
######################
# x = n*1 noisy vector
#
# tau.est = method for selecting the hyperparmeter tau in Dir(a) scale parameter
#    Option #1: 'fixed'. User-specified (not recommended)
#    Option #1: 'est.sparsity'
#    Option #2: 'reml'
#    NOTE: Default is MMLE if user does not specify either
#
# tau = parameter in Dir(tau,..., tau) scale parameter
#     Default is 1/n. Ignored if 'est.sparsity', 'reml', 'uniform', 
#     or 'halfCauchy' is used
#
# sigma2 = variance term
#      User specifies this. Default is 1
#
# var.select = method for selecting variables
#      Can select "threshold" for thresholding rule. Defaults to this 
#      Can also use "intervals" to select based on marginal credible intervals
#      NOTE: Default is "threshold" if user does not specify either
#
# max.steps = # of iterations to run MCMC. Default is 10,000
#
# burnin = # of samples in burn-in. Default is 5,000
#
#
##################
# RETURNS A LIST #
##################
# theta.hat = posterior mean estimate for theta
# theta.med = posterior median estimate for theta
# theta.var = posterior varince estimate for theta
# theta.intervals = 95% posterior credible intervals for each component
# theta.classifications = binary vector of classifications, according to 
#                         classification method chosen by user
# tau.est = empirical Bayes estimate of the hyperparameter 'tau'

dl.normalmeans = function(x, tau.est=c("fixed", "est.sparsity", "reml", "uniform", "truncatedCauchy"), 
                          tau=1/length(x), sigma2=1, var.select = c("threshold", "intervals"), 
                          max.steps=10000, burnin=5000){
  
  #####################################
  # Check that burnin < max.steps and #
  # that vector is non-empty          #
  #####################################
  if (burnin > max.steps)
    stop("Burnin cannot be greater than # of iterations.")
  if (length(x) == 0)
    stop("Please enter a vector length greater than 0.")
  
  #####################################
  # Number of samples in noisy vector #
  #####################################
  n <- length(x)
  
  ################################################
  # For the hyperparameter a, renamed as a.param #
  ################################################
  a.param <- tau
  
  if( tau.est=="fixed" )            # If 'fixed' method is used.
    a.param <- tau 
  
  if( tau.est=="est.sparsity" )     # if 'est.sparsity' method is used.
    a.param <- est.sparsity(x, sigma2)
  
  if( tau.est=="reml" )             # If 'reml' method is used.
    a.param <- dl.MMLE(x, sigma2)
  
  # If user specified a different 'a,' it must be greater than 0.
  if ( tau <= 0 )
    stop("ERROR: 'tau' should be positive. \n")
  
  # Check that sigma2 is greater than 0.
  if (sigma2 <=0 )
    stop("ERROR: sigma2 should be greater than 0. \n")
  
  if ( tau.est=="uniform" ){          # if 'uniform' method is used.
    a.param <- 1/n                  # starting value of the random walk
    a.samples <- rep(NA, max.steps) # Create a vector for samples of 'a'
    a.samples[1] <- a.param
  }
  
  if ( tau.est=="truncatedCauchy" ){  # if 'truncatedCauchy' method is used.
    a.param <- 1/n                  # starting value of the random walk
    a.samples <- rep(NA, max.steps) # Create a vector for samples of 'a'
    a.samples[1] <- a.param
  }
  
  ############################################
  # Initial guesses for theta, psi, phi, tau #
  ############################################
  theta <- rep(mean(x), n)
  psi <- rep(1, n)
  phi <- rep(1, n)
  tau <- rep(1, n)
  aux.t <- rep(NA, n) # For the auxiliary T_1, ..., T_n
  
  ###########################
  # For storing the samples #
  ###########################
  theta.samples <- rep(list(rep(NA,n)), max.steps)
  kappa.samples <- rep(list(rep(NA,n)), max.steps)
  
  ###########################
  # Start the Gibbs sampler #
  ###########################
  j <- 0
  while (j < max.steps) {
    j <- j + 1
    
    ###########
    # Counter #
    ###########
    if (j %% 1000 == 0) {
      cat("Iteration:", j, "\n")
    }
    
    ####################
    # Sample theta_i's #
    ####################
    kappa <- 1/(1+psi*phi^2*tau^2)  
    theta.means <- (1-kappa)*x     # conditional means for thetas
    theta.sds <- sqrt(sigma2*(1-kappa))   # conditional variances for thetas
    theta <- rnorm(n, mean=theta.means, sd=theta.sds)  # sample thetas as a block
    
    # Store theta samples
    theta.samples[[j]] <- theta
    
    # Store kappa samples
    kappa.samples[[j]] <- kappa
    
    
    # Function to block-sample from the GIG density
    gig.sample <- function(x){
      GIGrvg::rgig(1, lambda=u, chi=x, psi=w)
    }
    gig.sample.vec <- Vectorize(gig.sample) # handles a vector as input
    
    ##################
    # Sample psi_i's #
    ##################
    u <- 1/2
    v <- (1/(sigma2*tau^2))*(theta^2/phi^2) # nx1 vector
    v <- pmax(v, .Machine$double.eps) # For numerical stability
    w <- 1
    
    # Update psi's as a block
    psi <- gig.sample.vec(v)
    
    ########################
    # Sample tau and phi's #
    ########################
    
    # Sample auxiliary T_i's
    u <- a.param-1
    v <- (2/sqrt(sigma2))*abs(theta) #nx1 vector
    v <- pmax(v, .Machine$double.eps) # For numerical stability
    w <- 1
    # Update auxiliary T_i's as a block
    aux.t <- gig.sample.vec(v)
    
    # Set tau = sum of T_i's and phi_i = T_i/tau
    tau <- sum(aux.t)
    phi <- aux.t/tau
    
    # If 'uniform' was the method of 'tau.est'
    if ( (tau.est == "uniform") & (j >= 2) ){
      
      # Draw from proposal density
      prop.sd <- 1e-4
      a.star <- rtruncnorm(1, a=1/n, b=1, mean=a.samples[j-1], sd=prop.sd)
      
      # Ratio q(a.star | a)/q(a | a.star)
      q.ratio <- (pnorm((1-a.samples[j-1])/prop.sd)-pnorm((1/n-a.samples[j-1])/prop.sd))/(pnorm((1-a.star)/prop.sd)-pnorm((1/n-a.star)/prop.sd))
      
      # Ratio pi(a.star | rest)/pi(a | rest)
      pi.ratio <- 2^(n*(a.samples[j-1]-a.star))*(gamma(a.samples[j-1])/gamma(a.star))^n * prod(psi^(a.star-a.samples[j-1]))
      
      # Acceptance probability
      alpha <- min(1, (q.ratio*pi.ratio) )
      
      # Accept/reject algorithm
      if ( runif(1) < alpha){
        a.samples[j] <- a.star
      } else {
        a.samples[j] <- a.samples[j-1]
      }
      # Update a.param
      a.param <- a.samples[j]
    }
    
    # If 'truncatedCauchy' was the method of 'tau.est'
    if ( (tau.est == "truncatedCauchy") & (j >= 2) ){
      
      # Draw from proposal density
      prop.sd <- 1e-4
      a.star <- rtruncnorm(1, a=1/n, b=1, mean=a.samples[j-1], sd=prop.sd)
      
      # Ratio q(a.star | a)/q(a | a.star)
      q.ratio <- (pnorm((1-a.samples[j-1])/prop.sd)-pnorm((1/n-a.samples[j-1])/prop.sd))/(pnorm((1-a.star)/prop.sd)-pnorm((1/n-a.star)/prop.sd))
      
      # Ratio pi(a.star | rest)/pi(a | rest)
      pi.ratio <- 2^(n*(a.samples[j-1]-a.star))*(gamma(a.samples[j-1])/gamma(a.star))^n * prod(psi^(a.star-a.samples[j-1]))
      
      # Acceptance probability
      alpha <- min(1, (q.ratio*pi.ratio) )
      
      # Accept/reject algorithm
      if ( runif(1) < alpha){
        a.samples[j] <- a.star
      } else {
        a.samples[j] <- a.samples[j-1]
      }
      # Update a.param
      a.param <- a.samples[j]
    }
    
  }
  
  ###################
  # Discard burn-in #
  ################### 
  theta.samples <- tail(theta.samples,max.steps-burnin)
  kappa.samples <- tail(kappa.samples,max.steps-burnin)  
  
  # Return estimate of tau if fully Bayes method of estimating it.
  if( (tau.est == "uniform") | (tau.est == "truncatedCauchy") ){
    a.param <- mean(tail(a.samples,max.steps-burnin))
  }
  
  #######################################
  # Extract the posterior mean, median, #
  # 2.5th, and 97.5th percentiles       #
  #######################################
  theta.sample.mat <- simplify2array(theta.samples)
  rm(theta.samples)
  # Posterior mean
  theta.hat <- rowMeans(theta.sample.mat)
  # Posterior median
  theta.med <- apply(theta.sample.mat, 1, median)
  # Posterior var
  theta.sd <- apply(theta.sample.mat, 1, sd)
  theta.var <- theta.sd^2
  # endpoints of 95% posterior credible intervals
  theta.intervals <- apply(theta.sample.mat, 1, function(x) quantile(x, prob=c(.025,.975)))
  
  ##############################
  # Perform variable selection #
  ##############################
  # Initialize vector of binary entries: 0 for inactive variable b_i, 1 for active b_j
  dl.classifications <- rep(0,n)
  
  if(var.select=="threshold"){
    # Estimate the shrinkage factor kappa_i's from the MCMC samples
    kappa.sample.mat <- simplify2array(kappa.samples)
    rm(kappa.samples)
    kappa.estimates <- rowMeans(kappa.sample.mat)
    
    # Return indices of the signals according to our classification rule
    signal.indices <- which((1-kappa.estimates)>=0.5)
    # Reset classified signals as 1
    dl.classifications[signal.indices] <- 1
  }
  
  if(var.select=="intervals"){
    # Find the active covariates
    for(k in 1:n){
      if(theta.intervals[1,k] < 0 && theta.intervals[2,k] < 0){
        dl.classifications[k] <- 1
      }
      else if(theta.intervals[1,k] > 0 && theta.intervals[2,k] > 0){
        dl.classifications[k] <- 1
      }
    }
  }
  
  ###############################
  # Return list of theta.hat,   #
  # theta.med, theta.intervals, #
  # and dl.classifications      #
  ###############################
  # theta.hat = posterior mean point estimator
  # theta.med = posterior median point estimator
  # theta.var = posterior variance estimate
  # theta.intervals = endpoints of 95% posterior credible intervals
  # dl.classifications = selected variables
  # tau.estimate = estimate of 'tau'
  
  dl.output <- list(theta.hat = theta.hat,
                    theta.med = theta.med,
                    theta.var = theta.var,
                    theta.intervals = theta.intervals,
                    dl.classifications = dl.classifications,
                    tau.estimate = a.param)
  # Return list
  return(dl.output)
}

