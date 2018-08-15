#########################################################
# Function to implement EM/Gibbs sampling alogirhtm for #
# linear regression Using the normal beta prime prior   #
#########################################################

######################
# FUNCTION ARGUMENTS #
######################
# X = design matrix (n*p). It should already be centered
# y = response vector (n*1)
# method.hyperparameters = method for estimating 'a' and 'b' in the Beta Prime (a, b) prior.
#     'fixed' = user-specified hyperparameters 'a' and 'b'. Not recommended
#     'mml' = estimate 'a' and 'b' from maximum marginal likelihood.
#     DEFAULT is 'mml.'
# a = hyperparameter in BP(a,b) prior. Ignored if 'mml' method is used. 
# b = hyperparameter in BP(a,b) prior. Ignored if 'mml' method is used. 
# c = hyperparameter in the IG(c,d) prior for the variance sigma2. Defaults to 1/2.
# d = hyperparameter in the IG(c,d) prior for the variance sigma2. Defaults to 1/2.
# max.steps = # of times to run MCMC
# burnin = # of samples in burn-in
# selection = method for variable selection. Default is "dss."
#        DSS = 'decoupled shrinkage and selection' method advocated by
#               Hanh and Carvalho, which solves an optimization problem
#        intervals = select using the 95% credible intervals
# DEFAULT is 'DSS.'
#

##################
# RETURNS A LIST #
##################
# beta.hat = posterior mean estimate for beta (p*1 vector)
# beta.med = posterior median estimate for beta (p*1 vector)
# beta.intervals = 95% posterior credible intervals for each coefficient
# beta.classifications = binary vector of classifications
# sigma2.hat = estimated variance sigma2
# a.estimate = MML estimate of shape parameter 'a'
# b.estimate = MML estimate of shape parameter 'b'

nbp = function(X, y, method.hyperparameters=c("mml", "fixed"), a=0.5, b=0.5, 
               c=1e-5, d=1e-5, max.steps = 20000, burnin=10000, 
               selection=c("dss", "intervals")){
  
  # Burnin time should not exceed the total number of iterations.
  if (burnin > max.steps)
    stop("ERROR: Burn-in cannot be greater than # of iterations. \n")
  
  # If method for choosing 'a' is not either "fixed" or "mml," print error.
  if ( (method.hyperparameters != "fixed") & (method.hyperparameters != "mml") )
    stop("ERROR: Specify either 'fixed' or 'mml' for setting shape parameters (a,b) in beta prime(a,b).")
  
  # Check that 'a' and 'b' are both strictly greater than zero.
  if ( a<= 0 )
    stop("ERROR: Shape parameter 'a' should be positive. \n")
  if (b <= 0)
    stop("ERROR: Shape parameter 'b' should be positive. \n")
  
  # Extract dimensions n and p
  n <- nrow(X)
	p <- ncol(X)

	# Time saving
	XtX <- crossprod(X, X)	
	Xty <- crossprod(X, y)
	# List to hold the draws of beta
	beta.samples <- rep(list(rep(NA,p)), max.steps)
	
	# To store the lambda and xi samples
	lambda.samples <- matrix(0, max.steps, p)
	xi.samples <- matrix(0, max.steps, p)
	
	#######################################
	# Initial guesses for beta and sigma2 #
	#######################################
	beta <- rep(0, p)
	sigma2 <- 1
	
	# To store the sigma2 samples
	sigma2.samples <- rep(0, max.steps)
	
	###############################
	# Initial guess for lambda_i, #
	# xi_i, a and b               # 
	###############################
  
  lambda <- rep(1, p) # Initial guesses for lambda_i's
  xi <- rep(1, p)     # Initial guesses for xi_i's
  lambdaxi <- lambda*xi
  
  ########################################
  # If method.hyperparameters is "fixed" #
  ########################################
  if ( method.hyperparameters == "fixed"){
    a <- a
    b <- b
  }
  #####################################
  # If method.hyperpameters is "mml" #
  #####################################
  # Initialize for EM Monte Carlo if MML approach is selected
  if ( method.hyperparameters == "mml" ){
    # Initialize a and b
    a <- 0.01
    b <- 0.01
    
    # Initialize
    ab.current <- c(a,b)
    ab.update <- c(a,b)
    # b.current <- b
    # b.update <- b
    
    # Maximum number of iterations for EM algorithm
    EM.maxiter <- burnin
    EM.tol <- 1e-06 # Tolerance for the EM algorithm
    EM.dif <- 1     # Initialize

    # To store the samples of hyperparameters 'a' and 'b'
    a.samples <- vector(mode="numeric", length=0)
    b.samples <- vector(mode="numeric", length=0)
  }


  ###########################
	# Start the Gibbs sampler #
  ###########################
	j <- 0
	while (j < max.steps) {
		j <- j + 1

		if (j %% 1000 == 0) {
			cat("Iteration:", j, "\n")
		}
		
		###############
		# Sample beta #
		###############
		if (p <= n){
		  lambdaxi <- lambda*xi
		  ridge <- XtX + diag(1/(lambdaxi))
		  
		  # Conditional posterior mean
		  inv.ridge <- solve(ridge)
		  cond.M <- inv.ridge %*% Xty
		  
		  # Draw from MVN density
		  beta <- mvrnorm(1, cond.M, sigma2*inv.ridge)
		  
		} else if (p>n){
		  # Use the more efficient sampling algorithm if p>n
		  # from Bhattacharya et al. (2016)
		  lambdaxi <- lambda*xi
		  
		  L <- chol((1/sigma2)*(XtX+diag(1/as.numeric(lambdaxi),p,p)))
		  v <- solve(t(L), Xty/sigma2)
		  mu <- solve(L,v)
		  u <- solve(L,stats::rnorm(p))
		  beta <- mu + u
		}
		
		# Save the most recent estimate of beta to the list
		beta.samples[[j]] <- beta
		
		################################
		# Sample lambda_i's and xi_i's #
		################################
		
		# Shape and scale for lambda_i's
		ig.shape <- a + 1/2
		ig.scale <- (beta^2/(2*sigma2*xi))+1 # Scale parameters for i=1,..., p.
		
		# Parameters for xi_i's
		u <- b - 1/2
		w <- 2
		
		# Update the lambda_i's as a block
		ig.sample <- function(x){
		  rigamma(1, ig.shape, x)
		} 
		lambda <- sapply(ig.scale, ig.sample)
		
		# Store the lambda samples
		lambda.samples[j, ] <- lambda
		
		# Update the xi_i's as a block
		v <- beta^2/(sigma2*lambda)   # chi parameter for i=1, ..., p.
		# For numerical stability
		v <- pmax(v, 1e-10) # For numerical stability
		
		gig.sample <- function(x){
		    rgig(1, lambda=u, chi=x, psi=w)
    }
		xi <- sapply(v, gig.sample)
		
		# Store the xi samples
		xi.samples[j, ] <- xi
		
		#################		
		# Sample sigma2 #
		#################
		ig.shape <- (n+p+2*c)/2
		resid <- t(y-X%*%beta)%*%(y-X%*%beta)
		scaleterm.2 <- t(beta) %*% diag(1/lambdaxi) %*% beta
		ig.scale <- min((resid+scaleterm.2+2*d)/2, 1e10) # For numerical stability
		sigma2 <- rigamma(1, ig.shape, ig.scale) 
		
		# Store sigma2 samples
		sigma2.samples[j] <- sigma2
		
		#################################################
		# Update hyperparameters 'a' and 'b' in the EM  #
		# algorithm if method.hyperparameters == "mml" #
		#################################################
	
		if (method.hyperparameters == "mml"){
		  
		  if ( (j %% 100 == 0) & (j <= EM.maxiter) & ( EM.dif >= EM.tol) ) {
		  
		    ab.current <- ab.update  # Previous ab update gets saved here
		    # b.current <- b.update
		    low <- j - 99
		    high <- j
		
		    # Estimate E(ln(lambda_i)), i=1,...,p, from the Gibbs sampler
		    ln.lambda.terms <- colMeans(log(lambda.samples[low:high, ]))
		  
		    # Update a
		    fa <- function(a){
		       -p*digamma(a) - sum(ln.lambda.terms)  
		    }
		    a <- uniroot(f=fa, lower=.Machine$double.eps, upper=.Machine$double.xmax, maxiter=2000)$root

		    # Store the a value in samples
		    a.samples[as.integer(j/100)] <- a
		  
		    # Estimate E(ln(xi_i)), i=1,...,p, from the Gibbs sampler
		    ln.xi.terms <- colMeans(log(xi.samples[low:high, ]))
		  
		    # Update b
		    fb <- function(b){
		       -p*digamma(b) + sum(ln.xi.terms) 
		    }
		    b <- uniroot(f=fb, lower=.Machine$double.eps, upper=.Machine$double.xmax, maxiter=2000)$root

		    # Store the b value in samples
		    b.samples[as.integer(j/100)] <- b
		  
		    # Store the new a and b in ab.update
		    # Once it loops back up, ab.update gets stored in ab.current
		    ab.update <- c(a,b)
		    # b.update <- b
		  
		    # Update the difference
		    EM.dif <- sum((ab.update-ab.current)^2)
		    # EM.dif <- abs(b.update-b.current)
		  }
		}
	}
	
  ###################
  # Discard burn-in #
  ###################
  beta.samples <- tail(beta.samples,max.steps-burnin)
  sigma2.samples <- tail(sigma2.samples, max.steps-burnin)
  
  ##################################################################
  # Extract the posterior mean, median, 2.5th, 97.5th percentiles, #
  # as well as the estimates of sigma2 and MML estimates of (a, b) #
  ##################################################################
  beta.sample.mat <- simplify2array(beta.samples)
  # Posterior mean
	beta.hat <- rowMeans(beta.sample.mat)
	# Posterior median
	beta.med <- apply(beta.sample.mat, 1, function(x) median(x))
	# endpoints of 95% posterior credible intervals
	beta.intervals <- apply(beta.sample.mat, 1, function(x) quantile(x, prob=c(.025,.975)))
	# MML estimates of shape parameters (a, b)
	a.estimate <- a
	b.estimate <- b
	sigma2.hat <- mean(sigma2.samples)

  ##############################
  # Perform variable selection #
  ##############################
	# Initialize vector of binary entries: 0 for inactive variable b_i, 1 for active b_j
	nbp.classifications <- rep(0,p)
	
	################################
	# If selection is based on DSS #
	################################
	# Use local linear approximation to solve an optimization problem
	#
	if (selection == "dss"){
	  y1 <- X %*% beta.hat # Replace with this
	  weights <- abs(beta.hat) # Weights are based on the posterior mean
	                           # Posterior mean is never exactly zero, so no risk of equaling zero.
	 
	  # Fit LARS object using 10-fold CV, with design matrix reweighted by abs(weights)
	  model_cv <- cv.glmnet(x = X, y = y, type.measure = "mse", nfolds = 10, alpha = 1,
	                         penalty.factor = 1 / abs(weights), keep = TRUE)
	  
	  # Find indices of nonzero coefficients
	  dss.estimates <- coef(model_cv, s = model_cv$lambda.min)[-1]
	  active.set <- which(dss.estimates != 0)
	  
	  # Reset the entries in the classification vector to be one
	  nbp.classifications[active.set] <- 1
	}
	
	################################
	# If selection is based on the # 
	# 95% credible intervals       #
	################################
	if(selection == "intervals"){
	  # Find the active covariates
  	for(k in 1:p){
	    if(beta.intervals[1,k] < 0 && beta.intervals[2,k] < 0)
	      nbp.classifications[k] <- 1
	    else if(beta.intervals[1,k] > 0 && beta.intervals[2,k] > 0)
	      nbp.classifications[k] <- 1
	  }
	}
	
	
  ######################################
  # Return list of beta.est, beta.med, #
	# lower.endpoints, upper.endpoints,  #
	# and igg.classifications            #
  ######################################
	# beta.hat = posterior mean point estimator
	# beta.med = posterior median point estimator
	# beta.intervals = endpoints of 95% posterior credible intervals
	# nbp.classifications = results of variable selection
	# sigma2.hat = estimated variance
	# a.est = MML estimate of shape parameter 'a'
	# b.est = MML estimate of shape parameter 'b'
	
  nbp.output <- list(beta.hat = beta.hat,
                     beta.med = beta.med,
                     beta.intervals = beta.intervals,
                     nbp.classifications = nbp.classifications,
                     sigma2.hat = sigma2.hat,
                     a.estimate = a.estimate,
                     b.estimate = b.estimate)
  # Return list
	return(nbp.output)
}
