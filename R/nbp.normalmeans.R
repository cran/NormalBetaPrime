###################################
# Function to implement NBP Gibbs # 
# sampler for normal means        #
###################################

######################
# FUNCTION ARGUMENTS #
######################
# x = n*1 noisy vector
# a = first hyperparameter in the Beta Prime(a,b) scale parameter. 
#     Defaults to 1/n. Ignored if 'est.sparsity', 'reml', 'uniform',
#     or 'halfCauchy' is used.
#
# a.est = method for selecting the hyperparmeter a in Beta Prime(a,b)
#    Option #1: 'fixed.' User-specified (not recommended)
#    Option #2: 'est.sparsity'
#    Option #3: 'reml'
#    Option #4: 'uniform'
#    Option #5: 'halfCauchy'
#    NOTE: Default is "halfCauchy" if user does not specify either
# 
# b = second hyperparameter in Beta Prime(a,b) scale parameter.
#     Default is 1/2+1/n. 
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
##################
# RETURNS A LIST #
##################
# theta.hat = posterior mean estimate for theta
# theta.med = posterior median estimate for theta
# theta.var = posterior variance estimate for theta
# theta.intervals = 95% posterior credible intervals for each component
# theta.classifications = binary vector of classifications, according to 
#                         classification method chosen by user
# a.est = empirical Bayes estimate of the hyperparameter b

nbp.normalmeans = function(x, a.est=c("fixed", "est.sparsity", "reml", "uniform", "truncatedCauchy"), 
                           a=1/length(x), b=1/2+1/length(x), sigma2=1, var.select = c("threshold", "intervals"), 
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
  
  ############################
  # For the hyperparameter b #
  ############################
  # If user specified a different 'b,' it must be greater than 0.
  if ( b <= 0 )
    stop("ERROR: 'b' should be positive. \n")
  
  ############################
  # For the hyperparameter a #
  ############################
  a.param <- a
  
  if( a.est=="fixed" )              # If 'fixed' method is used.
    a.param <- a
  
  if( a.est=="est.sparsity" )       # if 'est.sparsity' method is used.
    a.param <- est.sparsity(x, sigma2)
  
  if( a.est=="reml" )               # If 'reml' method is used.
    a.param <- nbp.MMLE(x, a, sigma2)
  
  # If user specified a different 'a,' it must be greater than 0.
  if ( a <= 0 )
    stop("ERROR: 'a' should be positive. \n")
  
  # Check that sigma2 is greater than 0.
  if (sigma2 <=0 )
    stop("ERROR: sigma2 should be greater than 0. \n")
  
  if ( a.est=="uniform" ){            # if 'uniform' method is used.
    a.param <- 1/n                    # starting value of the random walk
    a.samples <- rep(NA, max.steps)   # Create a vector for samples of 'a'
    a.samples[1] <- a.param
  }
  
  if ( a.est=="truncatedCauchy" ){    # if 'truncatedCauchy' method is used.
    a.param <- 1/n                    # starting value of the random walk
    a.samples <- rep(NA, max.steps)   # Create a vector for samples of 'a'
    a.samples[1] <- a.param
  }
  
  ## REPARAMETRIZE Beta'(a,b) as lambda*xi, 
  #  WHERE lambda~G(a,1), xi~IG(b,1).
  ## 
  
  ##########################################
  # Initial guesses for theta, lambda, tau #
  ##########################################
  theta <- rep(mean(x), n)
  lambda <- rep(1, n)
  xi <- rep(1, n)
  kappa <- 1/(1+lambda*xi)
  
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
    theta.means <- (1-kappa)*x   # conditional means for thetas
    theta.sds <- sqrt(sigma2*(1-kappa))   # conditional variances for thetas
    theta <- rnorm(n, mean=theta.means, sd=theta.sds)  # sample thetas as a block
    
    # Store theta samples
    theta.samples[[j]] <- theta
    
    ################################
    # Sample lambda_i's and xi_i's #
    ################################
    # Parameters for lambda_i's
    u <- a.param - 1/2
    v <- theta^2/(sigma2*xi) # nx1 vector
    v <- pmax(v, 1e-10) # For numerical stability
    w <- 2
    
    # Update the lambda_i's as a block
    gig.sample <- function(x){
      GIGrvg::rgig(1, lambda=u, chi=x, psi=w)
    }
    gig.sample.vec <- Vectorize(gig.sample) # handles a vector as input
    lambda <- gig.sample.vec(v)
    
    # Rate and scale parameters for xi_i's
    alpha <- b + 1/2
    beta <- theta^2/(2*sigma2*lambda)+1 # nx1 vector
    beta <- pmax(beta, 1e-10) # For numerical stability
    
    # Update the xi_i's as a block
    ig.sample <- function(x){
      rigamma(1, alpha, x)
    } 
    ig.sample.vec <- Vectorize(ig.sample) # handles a vector as input
    xi <- ig.sample.vec(beta)
    
    # If 'uniform' was the method of 'a.est'
    if ( (a.est == "uniform") & (j >= 2) ){
      
      # Draw from proposal density
      prop.sd <- 1e-3
      a.star <- rtruncnorm(1, a=1/n, b=1, mean=a.samples[j-1], sd=prop.sd)
      
      # Ratio pi(a.star | rest)/pi(a | rest)
      omega <- lambda*xi
      pi.ratio <- ((gamma(a.samples[j-1])*(gamma(a.star+b)))/(gamma(a.star)*gamma(a.samples[j-1]+b)))^n * prod((omega/(1+omega))^(a.star-a.samples[j-1])) 
      
      # Ratio q(a.star | a)/q(a | a.star)
      q.ratio <- (pnorm((1-a.samples[j-1])/prop.sd)-pnorm((1/n-a.samples[j-1])/prop.sd))/(pnorm((1-a.star)/prop.sd)-pnorm((1/n-a.star)/prop.sd))
      
      # Acceptance probability
      alpha <- min(1, pi.ratio * q.ratio)
      
      # Accept/reject algorithm
      if ( runif(1) < alpha ){
        a.samples[j] <- a.star
      } else {
        a.samples[j] <- a.samples[j-1]
      }
      # Update a.param
      a.param <- a.samples[j]
    }
    
    # If 'truncatedCauchy' was the method of 'a.est'
    if ( (a.est == "truncatedCauchy") & (j >= 2) ){
      
      # Draw from proposal density
      prop.sd <- 1e-3
      a.star <- rtruncnorm(1, a=1/n, b=1, mean=a.samples[j-1], sd=prop.sd)
      
      # Ratio q(a.star | a)/q(a | a.star)
      q.ratio <- ( pnorm((1-a.samples[j-1])/prop.sd)-pnorm((1/n-a.samples[j-1])/prop.sd) )/( pnorm((1-a.star)/prop.sd)-pnorm((1/n-a.star)/prop.sd) )
      
      # Ratio pi(a.star | rest)/pi(a | rest)
      omega <- lambda*xi
      pi.ratio <- (((1+a.samples[j-1])*gamma(a.samples[j-1])*(gamma(a.star+b)))/((1+a.star)*gamma(a.star)*gamma(a.samples[j-1]+b)))^n * prod((omega/(1+omega))^(a.star-a.samples[j-1])) 
      
      # Acceptance probability
      alpha <- min(1, q.ratio*pi.ratio)
      
      # Accept/reject algorithm
      if ( runif(1) < alpha ){
        a.samples[j] <- a.star
      } else {
        a.samples[j] <- a.samples[j-1]
      }
      # Update a.param
      a.param <- a.samples[j]
    }
    
    ###################################
    # Compute the posterior shrinkage #
    # factors and store them          #
    ###################################
    kappa <- 1/(1+lambda*xi)    # shrinkage factors
    kappa.samples[[j]] <- kappa
    
  }
  
  ###################
  # Discard burn-in #
  ################### 
  theta.samples <- tail(theta.samples,max.steps-burnin)
  kappa.samples <- tail(kappa.samples,max.steps-burnin)  
  
  # Return estimate of 'a' if fully Bayes method of estimating it.
  if( (a.est == "uniform") | (a.est == "truncatedCauchy") ){
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
  nbp.classifications <- rep(0,n)
  
  if(var.select=="threshold"){
    # Estimate the shrinkage factor kappa_i's from the MCMC samples
    kappa.sample.mat <- simplify2array(kappa.samples)
    rm(kappa.samples)
    kappa.estimates <- rowMeans(kappa.sample.mat)
    
    # Return indices of the signals according to our classification rule
    signal.indices <- which((1-kappa.estimates)>=0.5)
    # Reset classified signals as 1
    nbp.classifications[signal.indices] <- 1
  }
  
  if(var.select=="intervals"){
    # Find the active covariates
    for(k in 1:n){
      if(theta.intervals[1,k] < 0 && theta.intervals[2,k] < 0){
        nbp.classifications[k] <- 1
      }
      else if(theta.intervals[1,k] > 0 && theta.intervals[2,k] > 0){
        nbp.classifications[k] <- 1
      }
    }
  }
  
  
  ########################################
  # Return list of theta.hat, theta.med, #
  # theta.var, theta.intervals, and      #
  # nbp.classifications                  #
  ########################################
  # theta.hat = posterior mean point estimator
  # theta.med = posterior median point estimator
  # theta.var = posterior variance estimate 
  # theta.intervals = endpoints of 95% posterior credible intervals
  # nbp.classifications = selected variables
  # a.estimate = estimated hyperparameter a
  
  nbp.output <- list(theta.hat = theta.hat,
                     theta.med = theta.med,
                     theta.var = theta.var,
                     theta.intervals = theta.intervals,
                     nbp.classifications = nbp.classifications,
                     a.estimate = a.param)
  # Return list
  return(nbp.output)
}