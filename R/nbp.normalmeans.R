###################################
# Function to implement NBP Gibbs # 
# sampler for normal means        #
###################################

######################
# FUNCTION ARGUMENTS #
######################
# x = n*1 noisy vector
# a = first hyperparameter in the Beta Prime(a,b) scale parameter. 
#     Defaults to 1/2+1/n
#
# b.est = method for selecting the hyperparmeter b in Beta Prime(a,b)
#    Option #1: 'fixed.' User-specified (not recommended)
#    Option #2: 'est.sparsity'
#    Option #3: 'reml'
#    NOTE: Default is "reml" if user does not specify either
# 
# b = second hyperparameter in Beta Prime(a,b) scale parameter.
#     Defaults is 1/n. Ignored if 'est.sparsity' or 'reml' is used. 
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
# theta.intervals = 95% posterior credible intervals for each component
# theta.classifications = binary vector of classifications, according to 
#                         classification method chosen by user
# b.est = empirical Bayes estimate of the hyperparameter b

nbp.normalmeans = function(x, a=1/2+1/length(x), b.est=c("fixed", "est.sparsity", "reml"), 
                   b=1/length(x), sigma2=1, var.select = c("threshold", "intervals"), 
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

  #######################################################
  # Make sure that the default selection method is the  #
  # thresholding method.                                #
  #######################################################
  if((var.select != "threshold") & (var.select != "intervals")){
      stop("Please enter either 'threshold' or 'intervals' for selection.")
   }
  
	############################
	# For the hyperparameter a #
	############################
  # If user specified a different 'a,' it must be greater than 0.
  if ( a<= 0 )
    stop("ERROR: 'a' should be positive. \n")
  
	############################
	# For the hyperparameter b #
	############################
	
  if( b.est=="fixed" )              # If 'fixed' method is used.
    b <- b
    
  if( b.est=="est.sparsity" )       # if 'est.sparsity' method is used.
    b <- est.sparsity(x, sigma2)
    
  if( b.est=="reml" )               # If 'reml' method is used.
    b <- nbp.MMLE(x, a, sigma2)

  # If user specified a different 'b,' it must be greater than 0.
  if ( b<= 0 )
    stop("ERROR: 'b' should be positive. \n")
  
  # Check that sigma2 is greater than 0.
  if (sigma2 <=0 )
    stop("ERROR: sigma2 should be greater than 0. \n")
  
  
  ## REPARAMETRIZE Beta'(a,b) AS lambda*xi, 
  #  WHERE lambda~IG(a,1), xi~G(b,1).
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
		if (j %% 1000 == 0) 
		   cat("Iteration:", j, "\n")
		
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
	  
	  # Rate and scale parameters for lambda_i's
	  alpha <- a + 1/2
	  beta <- theta^2/(2*sigma2*xi)+1 #nx1 vector
	  beta[beta < .Machine$double.eps] <- .Machine$double.eps # For numerical stability
	  
	  # Update the lambda_i's as a block
	  ig.sample <- function(x){
	    rigamma(1, alpha, x)
	  } 
	  lambda <- sapply(beta, ig.sample)
	  
	  # Parameters for xi_i's
    u <- b - 1/2
    v <- theta^2/(sigma2*lambda) #nx1 vector
    v[v < .Machine$double.eps] <- .Machine$double.eps # For numerical stability
	  w <- 2
	  
	  # Update the xi_i's as a block
	  gig.sample <- function(x){
	    rgig(1, lambda=u, chi=x, psi=w)
	  }
	  xi <- sapply(v, gig.sample)
	  
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

	#######################################
	# Extract the posterior mean, median, #
	# 2.5th, and 97.5th percentiles       #
	#######################################
	theta.sample.mat <- simplify2array(theta.samples)
	# Posterior mean
	theta.hat <- rowMeans(theta.sample.mat)
	# Posterior median
	theta.med <- apply(theta.sample.mat, 1, function(x) median(x))
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
	# Return list of theta.est, theta.med, #
	# lower.endpoints, upper.endpoints,    #
	# and nbp.classifications              #
	########################################
	# theta.hat = posterior mean point estimator
	# theta.med = posterior median point estimator
	# theta.intervals = endpoints of 95% posterior credible intervals
	# nbp.classifications = selected variables
	# b.estimate = estimated hyperparameter b
	
	nbp.output <- list(theta.hat = theta.hat,
	                   theta.med = theta.med,
	                   theta.intervals = theta.intervals,
	                   nbp.classifications = nbp.classifications,
	                   b.estimate = b)
	# Return list
	return(nbp.output)
}