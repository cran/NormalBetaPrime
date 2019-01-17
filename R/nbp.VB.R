#######################################################
# Function to implement variational EM alogirhtm for  #
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
# n.iter = # of iterations to run of VB 
# tol = threshold 

##################
# RETURNS A LIST #
##################
# beta.hat = posterior mean estimate for beta (p*1 vector)
# beta.post.var = posterior variance estimate for beta (p*1 vector)
# beta.intervals = 95% posterior credible intervals for each coefficient
# beta.classifications = binary vector of classifications
# sigma2.est = estimate of unknown error variance sigma2
# a.estimate = final EB estimate of 'a'
# b.estimate = final EB estimate of 'b'

nbp.VB = function(X, y, method.hyperparameters=c("fixed", "mml"), a=0.5, b=0.5, 
                  c=1e-5, d=1e-5, selection=c("dss","intervals"), tol=0.001, n.iter=1000){
  
  # If method for choosing 'a' is not either "fixed" or "mmle," print error.
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

  # Saves time
  XtX <-  t(X) %*% X
  Xty <- t(X) %*% y
  
  # Initialize dif and delta
  dif <- 1
  delta <- tol
  elbo.current <- 0
  elbo.update <- 0
  # For counter
  t <- 1      
  J <- n.iter
  
  # Initialize a.star, b.star
  if(method.hyperparameters=="fixed"){
    a.star <- a
    b.star <- b
    a.samples <- rep(a, n.iter+1)
    b.samples <- rep(b, n.iter+1)
  } 
  if(method.hyperparameters=="mml"){
    # If MML is used, initialize a.star=b.star=0.001
    a.star <- .001
    b.star <- .001
    a.samples <- numeric(0)
    b.samples <- numeric(0)
    a.samples[1] <- a.star
    b.samples[1] <- b.star
  }
  ###################################  
  # Initialize values of parameters #
  # in variational density          #
  ###################################
  l.star <- 2                   # Stays fixed for all updates
  c.star <- (n+p+2*c)/2         # Stays fixed for all updates
  d.star <- 1                   # Initial d.star
  k.star <- rep(1, p)           # Initial k_i, i=1,...,p
  v.star <- rep(1, p)           # Initial v_i, i=1,...,p
  beta.star <- rep(0, p)            
  Phi.star <- matrix(0, p, p)
  Sigma.star <- matrix(0, p, p)
  
  
  # Loop for Variational EM algorithm
  while ( (t <= J) & (dif > delta) ){
      
      elbo.current <- elbo.update # Previous ELBO update gets saved here
    
      ##########################################
      # E-STEP: Update variational parameters. #
      ##########################################
      #Update m.star
      m.star <- a.star - 1/2  
      # Update u.star
      u.star <- b.star + 1/2  
      
      # Time saving (compute expectations)
      EV.lambda <- EX.gig.vec(k.star, l.star, m.star) # Expected values of lambda_i^2, i=1,...,p
      EV.lambda.inv <- EXinv.gig.vec(k.star, l.star, m.star) # Expected values of lambda_i^(-2), i=1,...,p
      EV.xi.inv <- EXinv.invGamma.vec(u.star, v.star) # Expected values of xi_i^(-2), i=1,...,p
      EV.sigma2 <- EX.invGamma(c.star, d.star) # Expected value of sigma^2
      EV.sigma2.inv <- EXinv.invGamma(c.star, d.star) # Expected value of sigma^(-2)
      EV.loglambda <- ElogX.gig.vec(k.star, l.star, m.star) # Expected values of log(lambda_i^2), i=1,...,p
      EV.logxi <- ElogX.invGamma.vec(u.star, v.star) # Expected values of log(xi_i^2), i = 1,.,,,p
      
      # Update D.star
      D.star <- diag(pmax(EV.lambda.inv*EV.xi.inv,1e-10)) # For numerical stability
      D.star.inv <- solve(D.star) # Inverse of D.star
    
      # Time saving
      Phi.star <- D.star.inv - D.star.inv %*% t(X) %*% solve(diag(n)+X %*% D.star.inv %*% t(X)) %*% X %*% D.star.inv
      
      # Update variational Sigma.star
      Sigma.star <- EV.sigma2 * Phi.star
    
      # Update beta.star (px1 vector)
      beta.star <- Phi.star %*% Xty
      
      # Find E(beta_i^2), i=1,...,p
      EV.beta.sq <- rep(0,p) # Vector to store E(beta_i^2)'s
      EV.beta.sq <- beta.star^2 + diag(Sigma.star)
      
      # Update k.star (px1 vector)
      k.star <- EV.beta.sq*EV.sigma2.inv*EV.xi.inv 
      # Update v.star (px1 vector)
      v.star <- 0.5*EV.beta.sq*EV.sigma2.inv*EV.lambda.inv+1
    
      # Update d.star
      EV.resid <- sum(diag(XtX%*%Sigma.star))+sum((y-X%*%beta.star)^2)
      EV.bDb <- sum(EV.lambda.inv*EV.xi.inv*EV.beta.sq) 
      d.star <- 0.5*(EV.resid+EV.bDb+2*d)
      
      ##################################
      # M-STEP: Update hyperparameters #
      ##################################
      if(method.hyperparameters=="mml"){
        # Update a
        a.star <- a.update(a.star, EV.loglambda) 
        a.samples[t+1] <- a.star
        
        # Update b
        b.star <- b.update(b.star, EV.logxi)
        b.samples[t+1] <- b.star
      }
      
    #######################################
    # Compute ELBO and calculate new diff #
    #######################################
    # log-determinant term
    logdet.term <- determinant(Sigma.star, logarithm=TRUE)$modulus[1]  
    
    # Split up terms in sum for readability
    term1 <- -(n/2)*log(2*pi)+p/2+p*log(2)+p*lgamma(u.star)-p*lgamma(a.star)-p*lgamma(b.star)
    term2 <- c*log(d)-c.star*log(d.star)+lgamma(c.star)-lgamma(c)+0.5*logdet.term
    term3 <- -(m.star/2)*sum(log(k.star/l.star)) + sum(log(besselK(sqrt(k.star*l.star),m.star))) - u.star*sum(log(v.star))
    term4 <- sum((k.star/2-1)*EV.lambda)+(l.star/2)*sum(EV.lambda.inv)+sum((v.star-1)*EV.xi.inv)
    # Compute new ELBO
    elbo.update <- term1+term2+term3+term4
    
    # Calculate dif and iterate
    dif <- abs(elbo.update-elbo.current)
    t <- t+1
  }
  # End while loop 
  
  #####################################
  # Extract posterior mean, variance, #
  # and 95% credible intervals        #
  #####################################
  beta.hat <- beta.star
  rm(beta.star)
  beta.var <- diag(Sigma.star)
  rm(Sigma.star)
  # Extract 95% credible intervals
  lower.lim <- qnorm(0.025, mean=beta.hat, sd=sqrt(beta.var))
  upper.lim <- qnorm(0.975, mean=beta.hat, sd=sqrt(beta.var))
  beta.intervals <- cbind(lower.lim, upper.lim)
  # Estimate of unknown variance sigma2
  sigma2.estimate <- EV.sigma2
  
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
      if(beta.intervals[k,1] < 0 && beta.intervals[k,2] < 0)
        nbp.classifications[k] <- 1
      else if(beta.intervals[k,1] > 0 && beta.intervals[k,2] > 0)
        nbp.classifications[k] <- 1
    }
  }

  
  #############################
  # Return list of beta.hat,  #
  # beta.var, beta.intervals, #
  # and nbp.classifications   #
  #############################
  # beta.hat = posterior mean/mode estimate for beta
  # beta.var = posterior med
  # beta.intervals = endpoints of 95% posterior credible intervals
  # sigma2.estimate = posterior estimate of unknown variance sigma2
  # a.estimate = final estimate of a
  # b.estimate = final estimate of b
  
  nbp.output <- list(beta.hat = beta.hat,
                     beta.var = beta.var,
                     beta.intervals = beta.intervals,
                     nbp.classifications = nbp.classifications,
                     sigma2.estimate = sigma2.estimate,
                     a.estimate = a.star,
                     b.estimate = b.star)
  # Return list
  return(nbp.output)
}