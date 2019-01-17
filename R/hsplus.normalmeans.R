###########################################
# Function to implement Gibbs sampler for #
# the horseshoe+ prior for normal means   #
###########################################

######################
# FUNCTION ARGUMENTS #
######################
# x = n*1 noisy vector
#
# tau.est = method for selecting the hyperparmeter tau in scale parameter 
#    Option #1: 'fixed.' User-specified (not recommended). Default is 1
#    Option #1: 'est.sparsity'
#    Option #2: 'reml'
#    NOTE: Default is 'reml' if user does not specify either
#
# tau = global shrinkage parameter
#       Default is 1. Ignored if 'est.sparsity' or 'reml' is used.
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
# theta.var = posterior variance estimate for theta
# theta.intervals = 95% posterior credible intervals for each component
# theta.classifications = binary vector of classifications, according to 
#                         classification method chosen by user
# tau.est = empirical Bayes estimate of the global parameter tau

hsplus.normalmeans = function(x, tau.est=c("fixed", "est.sparsity", "reml", "uniform", "truncatedCauchy"), 
                              tau=1/length(x), sigma2=1, var.select=c("threshold", "intervals"), 
                              max.steps = 10000, burnin=5000){
  
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
  
  ##############################
  # For the hyperparameter tau #
  ##############################
  
  if( tau.est=="fixed" )             # If 'fixed' method is used.
    tau <- tau
  
  if( tau.est=="est.sparsity" )      # if 'est.sparsity' method is used.
    tau <- est.sparsity(x, sigma2)
  
  if( tau.est=="reml" )              # If 'reml' method is used.
    tau <- hsplus.MMLE(x, sigma2)
  
  # If user specified a different 'tau,' it must be greater than 0.
  if ( tau <= 0 )
    stop("ERROR: tau should be positive. \n")
  
  # Check that sigma2 is greater than 0.
  if (sigma2 <=0 )
    stop("ERROR: sigma2 should be greater than 0. \n")
  
  if ( tau.est=="uniform" ){              # if 'uniform' method is used.
    tau <- 1/n                      # starting value of the random walk
    tau.samples <- rep(NA, max.steps)     # Create a vector for samples of 'a'
    tau.samples[1] <- tau
  }
  
  
  if ( tau.est=="truncatedCauchy" ){      # if 'truncatedCauchy' method is used.
    tau <- 1/n                      # starting value of the random walk
    tau.samples <- rep(NA, max.steps)     # Create a vector for samples of 'a'
    tau.samples[1] <- tau
  }
  
  #######################################
  # Initial guesses for hyperparameters #
  #######################################
  lambda2 <- rep(1, n)
  nu <- rep(1, n)
  xi <- 1
  eta2 <- rep(1, n)
  phi <- rep(1, n)
  
  ###########################
  # For storing the samples #
  ###########################
  theta.samples <- rep(list(rep(NA,n)), max.steps)
  kappa.samples <- rep(list(rep(NA,n)), max.steps)
  
  # Start the Gibbs sampler
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
    kappa <- 1/(1+lambda2*tau)
    theta.means <- (1-kappa)*x     # conditional means for thetas
    theta.sds <- sqrt(sigma2*(1-kappa))   # conditional variances for thetas
    theta <- rnorm(n, mean=theta.means, sd=theta.sds)  # sample thetas
    
    # Store theta samples
    theta.samples[[j]] <- theta
    
    # Store kappa samples
    kappa.samples[[j]] <- kappa
    
    alpha <- 1 # Shape parameter is 1 in IG densities
    # Function for block-updating draws from inverse-gamma density
    ig.sample <- function(x){
      rigamma(1, alpha, x)
    } 
    ig.sample.vec <- Vectorize(ig.sample) # handles a vector as input
    
    #####################
    # Sample lambda_i's #
    #####################
    beta <- 1/nu + (theta^2)/(2*tau*sigma2) # nx1 vector
    beta <- pmax(beta, 1e-10) # For numerical stability
    # Update the lambda_i's as a block
    lambda2 <- ig.sample.vec(beta)
    
    #################
    # Sample nu_i's #
    #################
    beta <- 1/eta2 + 1/lambda2 # nx1 vector
    beta <- pmax(beta, 1e-10) # For numerical stability
    # Update the nu_i's as a block
    nu <- ig.sample.vec(beta)
    
    #############
    # Sample xi #
    #############
    beta <- 1 + 1/tau # scalar
    xi <- rigamma(1, alpha, beta)
    
    ###################
    # Sample eta2_i's #
    ###################
    beta <- 1/nu + 1/phi # nx1 vector
    beta <- pmax(beta, 1e-10) # For numerical stability
    # Update the eta2_i's as a block
    eta2 <- ig.sample.vec(beta)
    
    ##################
    # Sample phi_i's #
    ##################
    beta <- 1+1/eta2 #nx1 vector
    beta <- pmax(beta, 1e-10) # For numerical stability
    # Update the phi_i's as a block
    phi <- ig.sample.vec(beta)
    
    # If 'uniform' was the method of 'tau.est'
    if ( (tau.est == "uniform") & (j >= 2) ){
      
      # Draw from proposal density
      prop.sd <- 0.1
      tau.star <- rtruncnorm(1, a=1/n, b=1, mean=tau.samples[j-1], sd=prop.sd)
      
      # Ratio q(tau.star | tau)/q(tau | tau.star)
      q.ratio <- ( pnorm((1-tau.samples[j-1])/prop.sd)-pnorm((1/n-tau.samples[j-1])/prop.sd) )/( pnorm((1-tau.star)/prop.sd)-pnorm((1/n-tau.star)/prop.sd) )
      
      # Ratio pi(a.star | rest)/pi(a | rest)
      pi.ratio <- (tau.samples[j-1]/tau.star)^n * prod( (1+(lambda2/(eta2*tau^2)))/(1+(lambda2/(eta2*tau^2))) )
      
      # Acceptance probability
      alpha <- min(1, q.ratio*pi.ratio)
      
      # Accept/reject algorithm
      if ( runif(1,0,1) < alpha){
        tau.samples[j] <- tau.star
      } else {
        tau.samples[j] <- tau.samples[j-1]
      }
      # Update tau
      tau <- tau.samples[j]
    }
    
    # If 'truncatedCauchy' was the method of 'tau.est'
    if ( (tau.est == "truncatedCauchy") & (j >= 2) ){
      
      # Draw from proposal density
      prop.sd <- 1e-2
      tau.star <- rtruncnorm(1, a=1/n, b=1, mean=tau.samples[j-1], sd=prop.sd)
      
      # Ratio q(tau.star | tau)/q(tau | tau.star)
      q.ratio <- ( pnorm((1-tau.samples[j-1])/prop.sd)-pnorm((1/n-tau.samples[j-1])/prop.sd) )/( pnorm((1-tau.star)/prop.sd)-pnorm((1/n-tau.star)/prop.sd) )
      
      # Ratio pi(a.star | rest)/pi(a | rest)
      pi.ratio <- ( ((1+tau.samples[j-1])*tau.samples[j-1])/((1+tau.star)*tau.star) )^n * prod( (1+(lambda2/(eta2*tau^2)))/(1+(lambda2/(eta2*tau^2))) )
      
      # Acceptance probability
      alpha <- min(1, q.ratio*pi.ratio)
      
      # Accept/reject algorithm
      if ( runif(1,0,1) < alpha){
        tau.samples[j] <- tau.star
      } else {
        tau.samples[j] <- tau.samples[j-1]
      }
      # Update tau
      tau <- tau.samples[j]
    }
    
  }
  
  ###################
  # Discard burn-in #
  ################### 
  theta.samples <- tail(theta.samples,max.steps-burnin)
  kappa.samples <- tail(kappa.samples,max.steps-burnin)  
  
  # Return estimate of tau if fully Bayes method of estimating it.
  if( (tau.est == "uniform") | (tau.est == "truncatedCauchy") ){
    tau.samples <- tail(tau.samples,max.steps-burnin)
    tau.param <- mean(tau.samples)
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
  hsplus.classifications <- rep(0,n)
  
  if(var.select=="threshold"){
    # Estimate the shrinkage factor kappa_i's from the MCMC samples
    kappa.sample.mat <- simplify2array(kappa.samples)
    rm(kappa.samples)
    kappa.estimates <- rowMeans(kappa.sample.mat)
    
    # Return indices of the signals according to our classification rule
    signal.indices <- which((1-kappa.estimates)>=0.5)
    # Reset classified signals as 1
    hsplus.classifications[signal.indices] <- 1
  }
  
  if(var.select=="intervals"){
    # Find the active covariates
    for(k in 1:n){
      if(theta.intervals[1,k] < 0 && theta.intervals[2,k] < 0){
        hsplus.classifications[k] <- 1
      }
      else if(theta.intervals[1,k] > 0 && theta.intervals[2,k] > 0){
        hsplus.classifications[k] <- 1
      }
    }
  }
  
  
  ########################################
  # Return list of theta.est, theta.med, #
  # lower.endpoints, upper.endpoints,    #
  # and hsplus.classifications           #
  ########################################
  # theta.hat = posterior mean point estimator
  # theta.med = posterior median point estimator
  # theta.var = posterior variance estimate
  # theta.intervals = endpoints of 95% posterior credible intervals
  # hsplus.classifications = selected variables
  # tau.estimate = empirical Bayes estimate of global parameter tau
  
  hsplus.output <- list(theta.hat = theta.hat,
                        theta.med = theta.med,
                        theta.var = theta.var,
                        theta.intervals = theta.intervals,
                        hsplus.classifications = hsplus.classifications,
                        tau.estimate = tau)
  # Return list
  return(hsplus.output)
}
