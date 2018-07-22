#########################################
# Helper functions for NBP normal means #
#########################################

# This file contains the functions for obtaining Empirical Bayes estimates
# of the hyperparameter 'b' in the sparse normal means problem.

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
# hyperparameter 'b' in the NBP prior based on     #
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
# HYPERPARAMETER 'b' IN NBP #
#############################
# b = empirical Bayes estimate of a based on MMLES
#' @keywords internal
nbp.MMLE <- function(y, a=1/2+1/length(y), Sigma2=1){
  
  b <- stats::optimize(f = MMLE.nbp, data = y, a=a, data.var = Sigma2, 
                       lower = (1 / length(y)), upper = 1, maximum = TRUE)$maximum
  return(b)
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
Term.nbpMMLE <- function(y, a, t, data.var){ 
  
  # Joint of f(y, u), where u = lambda 
  f <- function(u) exp( -y^2 / (2*(1+u)*data.var) ) * u^(a-1) * (1+u)^(-a-t-1)
  
  # Integrate out u
  return(log( pracma::quadinf(f, 0, Inf, tol=1e-8)$Q ) ) 
}

# Helper functions for estimating the NBP's restricted MMLE
#' @keywords internal
Term.nbpMMLE.vec <- Vectorize(Term.nbpMMLE) #handles a vector as input

# Helper functions for estimating the NBP's restricted MMLE
#' @keywords internal
MMLE.nbp <- function(t, a, data, data.var){ #the quantity to be maximized
  sum(Term.nbpMMLE.vec(data, t, a, data.var)) - length(Term.nbpMMLE.vec)*lbeta(a,t)
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