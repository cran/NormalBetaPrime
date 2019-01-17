#########################################
# Helper functions for NBP normal means #
#########################################

# This file contains the functions for obtaining Empirical Bayes estimates
# of the hyperparameter 'a' in the sparse normal means problem.

################################################
# Function to find empirical Bayes estimate of #
# parameter based on estimated sparsity level  #            
################################################

######################
# FUNCTION ARGUMENTS #
######################
# x = n*1 noisy vector
# sigma2 = user-specified variance

###################
# RETURNS A VALUE #
###################
# est.signalprop = empirical Bayes estimate based on estimated sparsity level
#' @keywords internal
est.sparsity = function(x, sigma2){
  n <- length(x)
  
  # universal threshold
  univ.threshold <- sqrt(2*log(n)) 
  
  # normalize the vector x based on sigma2
  x <- x/sqrt(sigma2)
  
  # Count the # of entries in x that are greater than the universal threshold
  q <- length(x[(x>univ.threshold)==TRUE])
  
  # Reset est.signalprop as the EB estimator max(1/n, q/n)
  est.signalprop <- max(1/n, q/n)
  
  # Return est.signalprop
  return(est.signalprop) 
}


####################################################
# Function to find empirical Bayes estimate of     #
# hyperparameter 'a' in the NBP prior based on     #
# restricted marginal maximum likelihood (MMLE)    #                 
####################################################

######################
# FUNCTION ARGUMENTS #
######################
# y = n*1 noisy vector
# a = user-specified hyperparameter. Defaults to a=1/2+1/n
# Sigma2 = user-specified variance. Defaults to 1

#############################
# RETURNS ESTIMATE OF       #
# HYPERPARAMETER 'a' IN NBP #
#############################
# a = empirical Bayes estimate of a based on MMLES
#' @keywords internal
nbp.MMLE <- function(y, b=1/2+1/length(y), Sigma2=1){
  
  a <- stats::optimize(f = MMLE.nbp, b=b, data = y, data.var = Sigma2, 
                       lower = (1 / length(y)), upper = 1, maximum = TRUE)$maximum
  return(a)
}

#################
# Term.nbpMMLE: #
#################
# This function first takes the joint density f(y, lambda) for a single y,
# which is given by integrating out theta from the joint f(y, theta, lambda). 
# Then we integrate out lambda to approximate f(y).
#
# u=lambda, t=b.
#
#' @keywords internal
Term.nbpMMLE <- function(y, b, t, data.var){ 
  
  # Joint of f(y, u), where u ~ BetaPrime(a,b)
  f <- function(u) exp( -y^2 / (2*(1+u)*data.var) ) * u^(t-1) * (1+u)^(-(t+b+1/2))
  
  # Integrate out u
  return(log( pracma::quadinf(f, 0, Inf, tol=1e-8)$Q ) ) 
}

# Helper functions for estimating the NBP's restricted MMLE
#' @keywords internal
Term.nbpMMLE.vec <- Vectorize(Term.nbpMMLE) #handles a vector as input

# Helper functions for estimating the NBP's restricted MMLE
#' @keywords internal
MMLE.nbp <- function(data, b, t, data.var){ #the quantity to be maximized
  sum(Term.nbpMMLE.vec(data, b, t, data.var)) - length(data)*lbeta(t,b)
}


####################################################
# Function to find empirical Bayes estimate of tau #
# in the HS+ prior based on estimated on           #
# restricted marginal maximum likelihood (MMLE)    #                 
####################################################

######################
# FUNCTION ARGUMENTS #
######################
# y = n*1 noisy vector
# Sigma2 = user-specified variance

#############################
# RETURNS ESTIMATE OF       #
# HYPERPARAMETER tau IN HS+ #
#############################
# tau = empirical Bayes estimate of a based on MMLE
#' @keywords internal
hsplus.MMLE <- function(y, Sigma2){
  tau <- stats::optimize(f = MMLE.hsplus, data = y, data.var = Sigma2, 
                         lower = (1 / length(y)), upper = 1, maximum = TRUE)$maximum
  
  # Return the estimate of tau
  return(tau)
}


####################
# Term.hsplusMMLE: #
####################
# This funciton first takes the joint density of f(y, lambda) which is obtained
# by integrating out theta from f(y, theta, lambda). We then integrate out lambda
# to approximate the marginal f(y). 
# 
# u = lambda, t = tau.
#
#' @keywords internal
Term.hsplusMMLE <- function(y, t, data.var){ 
  
  # Joint of f(y, u), where u = lambda.
  f <- function(u){   
    exp( -y^2/(2*(1+u)*data.var) ) * (1+u)^(-1/2) * (1/t) * ( log(u/t)/( (u/t)^2-1) )  
  }
  
  # Integrate out u
  return(log( pracma::quadinf(f, 0, Inf, tol=1e-8)$Q ) ) 
}

## Helper functions for estimating the MMLE
#' @keywords internal
Term.hsplusMMLE.vec <- Vectorize(Term.hsplusMMLE) #handles a vector as input

## Helper functions for estimating the MMLE
#' @keywords internal
MMLE.hsplus <- function(t, data, data.var){ 
  # the quantity to be maximized
  sum(Term.hsplusMMLE.vec(data, t, data.var))
}


#######################################################
# Function to find empirical Bayes estimate of 'a' in #
# the Dirichlet-Laplace prior based on restricted     # 
#marginal maximum likelihood (MMLE)                   #                 
#######################################################

######################
# FUNCTION ARGUMENTS #
######################
# y = n*1 noisy vector
# Sigma2 = user-specified variance

###########################
# RETURNS ESTIMATE OF     #
# HYPERPARAMETER a  in DL #
###########################
# a = empirical Bayes estimate of a based on MMLE
#' @keywords internal
dl.MMLE <- function(y, Sigma2){
  a <- stats::optimize(f = MMLE.dl, data = y, data.var = Sigma2, 
                       lower = (1 / length(y)), upper = 1, maximum = TRUE)$maximum
  
  # Return the estimate of a
  return(a)
}

##############
# Term.MMLE: #
###############
# This functon first takes the joint density of f(y, theta), where f(y|theta)
# follows N(theta, sigma^2) and f(theta) is given by Eq. (11) in Bhattacharya et al.
# (2015). Then we integrate out theta to approximate f(y).
#
# u = theta, t=a. 
#
#' @keywords internal
Term.dlMMLE <- function(y, t, data.var){ 
  
  # Joint of f(y, u), where u = theta.
  f <- function(u){   
    2*exp(-(y-u)^2 / (2*data.var)) * (1/gamma(t)) * (1/(2^((1+t)/2)) ) * (u^((t-1)/2)) * besselK(sqrt(2*u), 1-t)
  }
  
  # Integrate out u
  return(log( pracma::quadinf(f, 0, Inf, tol=1e-8)$Q ) ) 
}

# Helper functions for estimating the DL's restricted MMLE
#' @keywords internal
Term.dlMMLE.vec <- Vectorize(Term.dlMMLE) #handles a vector as input

## Helper functions for estimating the DL's restricted MMLE
#' @keywords internal
MMLE.dl <- function(t, data, data.var){ 
  # the quantity to be maximized
  sum(Term.dlMMLE.vec(data, t, data.var))
}


##############################################
##############################################
## EXPECTED VALUES IN VARIATIONAL DENSITIES ## 
##############################################
##############################################
# Compute the expected values of X, X^-1, and log(X) for X~IG(alpha, beta),
# i.e. the inverse gamma distribution. 
# Takes the shape parameters as arguments.

# E(X) for X~IG(alpha, beta)
#' @keywords internal
EX.invGamma <- function(alpha, beta){
  if( (alpha <= 0) | (beta <= 0) ){
    stop("Invalid shape parameters.")
  }
  
  if ((alpha>0) & (alpha <= 1)){
    stop("The mean is undefined.")
  } else{
    return(beta/(alpha-1))
  }
}
#' @keywords internal
EX.invGamma.vec <- Vectorize(EX.invGamma, "beta") # Vectorize

#E(X^(-1)) for X~IG(alpha, beta)
#' @keywords internal
EXinv.invGamma <- function(alpha, beta){
  if( (alpha <= 0) | (beta <= 0) ){
    stop("Invalid shape parameters.")
  } else {
    return(alpha/beta)
  }
}
#' @keywords internal
EXinv.invGamma.vec <- Vectorize(EXinv.invGamma, "beta") # Vectorize

#E[log(X)] for X~IG(alpha, beta)
#' @keywords internal
ElogX.invGamma <- function(alpha, beta){
  if( (alpha <= 0) | (beta <= 0) ){
    stop("Invalid shape parameters.")
  } else {
    return(log(beta)-digamma(alpha))
  }
}
#' @keywords internal
ElogX.invGamma.vec <- Vectorize(ElogX.invGamma, "beta") #Vectorize

# Compute the expected values of X, X^-1, and log(X) for X~GIG(u,v,p), 
# i.e. the inverse gamma distribution. 
# Takes the shape parameters as arguments.

#E(X) for X~GIG(u, v, p)
#' @keywords internal
EX.gig <- function(u, v, p){
  if( (u<=0) | (v<=0) ){
    stop("Invalid shape parameters.")
  } else {
    ratio <- max(HyperbolicDist::besselRatio(sqrt(u*v), p, orderDiff=1),1e-8) # for numerical stability
    return(ratio*sqrt(v)/sqrt(u))
  }
}
#' @keywords internal
EX.gig.vec <- Vectorize(EX.gig, "u") # Vectorize

#E(X^(-1)) for X~GIG(u, v, p)
#' @keywords internal
EXinv.gig <- function(u, v, p){
  if( (u<=0) | (v<=0) ){
    stop("Invalid shape parameters.")
  } else {
    ratio <- max(HyperbolicDist::besselRatio(sqrt(u*v), p, orderDiff=1),1e-8) # for numerical stability
    return(ratio*sqrt(u)/sqrt(v)-(2*p)/v)
    
  }
}
#' @keywords internal
EXinv.gig.vec <- Vectorize(EXinv.gig, "u") # Vectorize

#E[log(X)] for X~GIG(u, v, p)
#' @keywords internal
ElogX.gig <- function(u, v, p){
  if( (u<=0) | (v<=0) ){
    stop("Invalid shape parameters.")
  } else {
    # Numerical differential term approximation
    den <- max(besselK(sqrt(u*v),p), 1e-8) # For numerical stability
    num <- pracma::numderiv(function(x){besselK(sqrt(u*v), nu=x)}, x0=p, tol=1e-8)$df
    
    return(0.5*log(v)-0.5*log(u)+(num/den))
  }
}
#' @keywords internal
ElogX.gig.vec <- Vectorize(ElogX.gig, "u") # Vectorize


#########################################
#########################################
## VARIATIONAL EM ALGORITHM UPDATES OF ##
## THE SHAPE PARAMETERS (a,b) IN NBP   ## 
#########################################  
#########################################

# For empirical Bayes updates of shape parameter 'a'.
# Arguments: 'a'=current iteration of 'a', 'EV.loglambda'=px1-dimensional vector
#' @keywords internal
a.update <- function(a, EV.loglambda){
  p <- length(EV.loglambda)
  fa <- function(a){ p*digamma(a) - sum(EV.loglambda) }   
  
  # Find the root of f(a)
  a.new <- uniroot(f=fa, lower=.Machine$double.eps, upper=.Machine$double.xmax, tol=1e-3, maxiter=2000)$root
  return(min(a.new,1)) # Bessel function computations are unstable for large 'm.star' (large 'a.star')
}

# For empirical Bayes updates of shape parameter 'b'.
# Arguments: 'b'=current iteration of 'b', 'EV.logxi'=px1 dimensional vector
#' @keywords internal
b.update <- function(b, EV.logxi){
  p <- length(EV.logxi)
  fb <- function(b){ -p*digamma(b) - sum(EV.logxi) }   
  
  # Find the root of f(b)
  b.new <- uniroot(f=fb, lower=.Machine$double.eps, upper=.Machine$double.xmax, tol=1e-3, maxiter=2000)$root
  
  return(b.new)
}