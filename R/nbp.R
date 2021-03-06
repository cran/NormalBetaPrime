#######################################################
# Function to implement Monte Carlo EM alogirhtm for  #
# linear regression Using the normal beta prime prior #
#######################################################

######################
# FUNCTION ARGUMENTS #
######################
# X = design matrix (n*p). It should already be centered
# y = response vector (n*1)
# method.hyperparameters = method for estimating 'a' and 'b' in the Beta Prime (a, b) prior.
#     'fixed' = user-specified hyperparameters 'a' and 'b'. Not recommended
#     'mml' = estimate 'a' and 'b' from maximum marginal likelihood.
# a = hyperparameter in BP(a,b) prior. Ignored if 'mml' method is used. 
# b = hyperparameter in BP(a,b) prior. Ignored if 'mml' method is used. 
# c = hyperparameter in the IG(c,d) prior for the variance sigma2. Defaults to 1e-5.
# d = hyperparameter in the IG(c,d) prior for the variance sigma2. Defaults to 1e-5.
# selection = method for variable selection. Default is "dss."
#        DSS = 'decoupled shrinkage and selection' method advocated by
#               Hanh and Carvalho, which solves an optimization problem
#        intervals = select using the 95% credible intervals
# max.steps = # of times to run MCMC
# burnin = # of samples in burn-in

##################
# RETURNS A LIST #
##################
# beta.hat = posterior mean estimate for beta (p*1 vector)
# beta.med = posterior median estimate for beta (p*1 vector)
# beta.var = posterior variance estimate for beta (p*1 vector)
# beta.intervals = 95% posterior credible intervals for each coefficient
# beta.classifications = binary vector of classifications
# sigma2.estimate = estimate of unknown error variance sigma2
# a.estimate = MML estimate of 'a'
# b.estimate = MML estimate of 'b'

nbp = function(X, y, method.hyperparameters=c("fixed","mml"), a=0.5, b=0.5, c=1e-5, d=1e-5, 
                     selection=c("dss","intervals"), max.steps = 15000, burnin=10000) {
  
  # Burnin time should not exceed the total number of iterations.
  if (burnin > max.steps)
    stop("ERROR: Burn-in cannot be greater than # of iterations. \n")
  
  # If method for choosing 'a' is not either "fixed" or "mml," print error.
  if ( (method.hyperparameters != "fixed") & (method.hyperparameters != "mml") )
    stop("ERROR: Specify method for setting 'a' and 'b' in beta prime(a,b).")
  
  # If selection method is set as other than "dss" or "intervals,"
  # print error.
  if ( (selection != "dss") & (selection !="intervals") )
    stop("ERROR: Specify method for variable selection.")
  
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
  # To store the draws of sigma2
  sigma2.samples <- rep(0, max.steps)
  
  # To store the lambda and xi samples
  lambda.samples <- matrix(0, max.steps, p)
  xi.samples <- matrix(0, max.steps, p)
  
  #######################################
  # Initial guesses for beta and sigma2 #
  #######################################
  # Initial guess for beta and sigma2
  beta <- rep(0, p)
  sigma2 <- 1
  
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
  ####################################
  # If method.hyperpameters is "mml" #
  ####################################
  # Initialize for EM Monte Carlo if MML approach is selected
  if ( method.hyperparameters == "mml" ){
    # Initialize a and b
    a <- 0.01
    b <- 0.01
    
    # Initialize
    ab.current <- c(a,b)
    ab.update <- c(a,b)
    
    # Maximum number of iterations for EM algorithm
    EM.maxiter <- burnin
    EM.tol <- 1e-03 # Tolerance for the EM algorithm
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
    
    ## REPARAMETRIZE Beta'(a,b) as lambda*xi, 
    #  WHERE lambda~G(a,1), xi~IG(b,1).
    ################################
    # Sample lambda_i's and xi_i's #
    ################################
    
    # Parameters for lambda_i's
    u <- a-1/2
    v <- beta^2/(sigma2*xi)       # chi parameter for i=1, ..., p.
    v <- pmax(v, .Machine$double.eps)   # For numerical stability
    w <- 2
    
    # To sample from the GIG density
    gig.sample <- function(x){
      GIGrvg::rgig(1, lambda=u, chi=x, psi=w)
    }
    gig.sample.vec <- Vectorize(gig.sample) # handles a vector as input
    
    # Update the lambda_i's as a block
    lambda <- gig.sample.vec(v)
    # Store the lambda samples
    lambda.samples[j, ] <- pmax(lambda, 1e-8) # for numerical stability
    
    # Scale parameter for xi_i's
    ig.shape <- b+1/2
    ig.scale <- (beta^2/(2*sigma2*lambda))+1 # Scale parameters for i=1,..., p.
    
    # To sample from the IG density
    ig.sample <- function(x){
      rigamma(1, ig.shape, x)
    } 
    ig.sample.vec <- Vectorize(ig.sample) # handles a vector as input
    
    # Update the xi_i's as a block
    xi <- ig.sample.vec(ig.scale)
    # Store the xi samples
    xi.samples[j, ] <- pmax(xi, 1e-8) # For numerical stability
    
    #################		
    # Sample sigma2 #
    #################
    ig.shape <- (n+p+2*c)/2
    
    # Efficiently compute sum of squared residuals
    resid <- sum((y-X%*%beta)^2)
    # Efficiently compute t(beta)%*%diag(1/lambdaxi)%*%beta
    scaleterm.2 <- (1/lambdaxi)*beta
    scaleterm.2 <- sum(beta*scaleterm.2)   
    ig.scale <- min((resid+scaleterm.2+2*d)/2, 1e4) # For numerical stability
    
    # Sample sigma2
    sigma2 <- rigamma(1, ig.shape, ig.scale) 
    
    # Store sigma2 samples
    sigma2.samples[j] <- sigma2
    
    #################################################
    # Update hyperparameters 'a' and 'b' in the EM  #
    # algorithm if method.hyperparameters == "mml"  #
    #################################################
    
    if (method.hyperparameters == "mml"){
      
      if ( (j %% 100 == 0) & (j <= EM.maxiter) & ( EM.dif >= EM.tol) ) {
        
        ab.current <- ab.update  # Previous ab update gets saved here
        
        low <- j - 99
        high <- j
        
        # Estimate E(ln(lambda_i)), i=1,...,p, from the Gibbs sampler
        ln.lambda.terms <- colMeans(log(lambda.samples[low:high, ]))
        
        # Update a
        fa <- function(a){
          -p*digamma(a) + sum(ln.lambda.terms)  
        }
        a <- uniroot(f=fa, lower=.Machine$double.eps, upper=.Machine$double.xmax, maxiter=2000)$root
        
        # Store the a value in samples
        a.samples[as.integer(j/100)] <- a
        
        # Estimate E(ln(xi_i)), i=1,...,p, from the Gibbs sampler
        ln.xi.terms <- colMeans(log(xi.samples[low:high, ]))
        
        # Update b
        fb <- function(b){
          p*digamma(b) + sum(ln.xi.terms) 
        }
        b <- uniroot(f=fb, lower=.Machine$double.eps, upper=.Machine$double.xmax, maxiter=2000)$root
        
        # Store the b value in samples
        b.samples[as.integer(j/100)] <- b
        
        # Store the new a and b in ab.update
        # Once it loops back up, ab.update gets stored in ab.current
        ab.update <- c(a,b)
        
        # Update the difference
        EM.dif <- sqrt(sum((ab.update-ab.current)^2))
      }
    }
  }
  
  ###################
  # Discard burn-in #
  ###################
  beta.samples <- tail(beta.samples,max.steps-burnin)
  sigma2.samples <- tail(sigma2.samples, max.steps-burnin)
  # Remove samples of latent variables
  rm(lambda.samples, xi.samples)
  
  ##############################################
  # Extract the posterior mean, median,        #
  # variance, &  2.5th, and 97.5th percentiles #
  ##############################################
  beta.sample.mat <- simplify2array(beta.samples)
  rm(beta.samples)
  
  # Posterior mean
  beta.hat <- rowMeans(beta.sample.mat)
  # Posterior median
  beta.med <- apply(beta.sample.mat, 1, median)
  # Posterior var
  beta.sd <- apply(beta.sample.mat, 1, sd)
  beta.var <- beta.sd^2
  # endpoints of 95% posterior credible intervals
  beta.intervals <- apply(beta.sample.mat, 1, function(x) quantile(x, prob=c(.025,.975)))
  
  
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
  
  ##################################
  # Store other summary statistics #
  # of the interest                #
  ##################################
  # Estimate of unknown variance, sigma2
  sigma2.estimate <- mean(sigma2.samples)
  # MML estimate of 'a'
  a.estimate <- tail(a.samples, 1)
  # MML estimate of 'b'
  b.estimate <- tail(b.samples, 1)
  
  ######################################
  # Return list of beta.hat, beta.med, #
  # beta.var, beta.intervals, and      #
  # nbp.classifications,               #
  ######################################
  # beta.hat = posterior mean point estimator
  # beta.med = posterior median point estimator
  # beta.var = posterior variance estimates
  # beta.intervals = endpoints of 95% posterior credible intervals
  # nbp.classifications = variable selection
  # sigma2.estimate = posterior estimate of unknown variance sigma2
  # a.estimate = MML estimate of 'a'
  # b.estimate = MML estimate of 'b'
  
  nbp.output <- list(beta.hat = beta.hat,
                     beta.med = beta.med,
                     beta.var = beta.var,
                     beta.intervals = beta.intervals,
                     nbp.classifications = nbp.classifications,
                     sigma2.estimate = sigma2.estimate,
                     a.estimate = a.estimate,
                     b.estimate = b.estimate
  )
  # Return list
  return(nbp.output)
}
