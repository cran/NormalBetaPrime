\name{nbp}
\alias{nbp}

\title{Normal-Beta Prime Regression}

\description{This function implements the Monte Carlo EM (Gibbs sampling) approach for the normal-beta prime (NBP) model in the standard linear regression model,

\deqn{y = X \beta + \epsilon,} 

where \eqn{\epsilon \sim N_n (0, \sigma^2 I_n)}. This is achieved by placing the normal-beta prime (NBP) prior of Bai and Ghosh (2019) on the coefficients of \eqn{\beta}. In the case where \eqn{p > n}, we utilize an efficient sampler from Bhattacharya et al. (2016) to reduce the computational cost of sampling from the full conditional density of \eqn{\beta} to \eqn{O(n^2 p)}. The hyperparameters can be set deterministically by the user or they may be automatically selected by marginal maximum likelihood (MML).

Note that the Gibbs sampling implementation is slower than the variational Bayes (VB) function, \code{nbp.VB}, but \code{nbp} tends to give much more accurate point esimates and posterior approximations. It is recommended that the user use \code{nbp}.
}

\usage{
nbp(X, y, method.hyperparameters=c("fixed","mml"), a=0.5, b=0.5, c=1e-5, d=1e-5,  
    selection=c("dss", "intervals"), max.steps = 15000, burnin=10000)}

\arguments{
  \item{X}{       \eqn{n \times p} design matrix. Should be centered.}

  \item{y}{       \eqn{n \times 1} response vector. Should be centered.}
  
  \item{method.hyperparameters}{ The method for estimating the shape parameters \eqn{(a,b)}. If \code{"mml"} is chosen, the function estimates the marginal maximum likelihood (MML) estimates of \eqn{(a,b)} by the EM algorithm described in Bai and Ghosh (2019). If \code{"fixed"} is chosen, then \eqn{(a,b)} are not estimated from the data. }

  \item{a}{ Shape parameter for \eqn{\beta'(a,b)}. The default is 0.5. The user may specify a different value for \code{a} (\eqn{a>0}). This is ignored if the method for estimating hyperparameters is \code{"mml"}.}

  \item{b}{ Shape parameter for \eqn{\beta'(a,b)}. The default is 0.5. The user may specify a different value for \code{b} (\eqn{b>0}). This is ignored if the method for estimating hyperparameters is \code{"mml"}.}
  
  \item{c}{ The shape parameter for the \eqn{IG(c,d)} prior on unknown variance parameter, \eqn{\sigma^2}. The default is \eqn{10^{-5}}. }

\item{d}{ The rate parameter for the \eqn{IG(c,d)} prior on the unknown variance parameter, \eqn{\sigma^2}. The default is \eqn{10^{-5}}. }

  \item{selection}{ The method of variable selection. \code{"dss"} implements the decoupled selection and shrinkage (DSS) method of Hahn and Carvalho (2015) to select variables. \code{"intervals"} performs variable selection by examining the 95 percent posterior credible intervals for each coefficient, \eqn{\beta_i, i = 1, \ldots, p}.}

\item{max.steps}{ The total number of iterations to run in the Gibbs sampler. Defaults to 15,000.}

  \item{burnin}{ The number of burn-in iterations for the Gibbs sampler. Defaults to 10,000.}

}

\details{
The function implements the normal-beta prime (NBP) model of Bai and Ghosh (2019) using Gibbs sampling. The posterior variances and 95 percent credible intervals for each of the \eqn{p} covariates are also returned so that the user may assess uncertainty quantification. The full model is:

\deqn{Y | (X, \beta ) \sim N_n(X \beta, \sigma^2 I_n),}
\deqn{\beta_i | \omega_i^2 \sim N(0, \sigma^2 \omega_i^2), i = 1, ..., p,}
\deqn{\omega_i^2 \sim \beta'(a,b), i = 1, ..., p,}
\deqn{\sigma^2 \sim IG(c,d),}

where \eqn{\beta'(a,b)} denotes the beta prime density,

\deqn{\pi(\omega_i^2) = \frac{\Gamma(a+b)}{\Gamma(a) \Gamma(b)} (\omega_i^2)^{a-1} (1+\omega_i^2)^{-a-b}. }

}

\value{The function returns a list containing the following components:
\item{beta.hat}{ The posterior mean estimate of \eqn{\beta}. }
\item{beta.med}{ The posterior median estimate of \eqn{\beta}. }
\item{beta.var}{ The posterior variance estimates for each \eqn{\beta_i, i=1, \ldots, p.} }
\item{beta.intervals}{ The 95 percent credible intervals for all \eqn{p} estimates in \eqn{\beta}. }
\item{nbp.classifications}{ A \eqn{p}-dimensional binary vector with "1" if the covariate is selected and "0" if it is deemed irrelevant.}
\item{sigma2.esimate}{ Estimate of unknown variance component \eqn{\sigma^2}.}
\item{a.estimate}{ MML estimate of shape parameter \eqn{a}. If \eqn{a} was fixed \emph{a priori}, returns fixed \eqn{a}.}
\item{b.estimate}{ MML estimate of shape parameter \eqn{b}. If \eqn{b} was fixed \emph{a prior}, returns fixed \eqn{b}.}
}

\references{
Bai, R. and Ghosh, M. (2019). "On the beta prime prior for scale parameters in high-dimensional Bayesian regression models." Pre-print, arXiv:1807.06539.

Bhattacharya, A., Chakraborty, A., and Mallick, B.K. (2016). "Fast sampling with Gaussian scale mixture priors in high-dimensional regression." \emph{Biometrika}, \bold{69}(2): 447-457.

Hahn, P. R. and Carvalho, C. M. (2015). "Decoupling shrinkage and selection in Bayesian linear models: A posterior summary perspective." \emph{Journal of the American Statistical Association}, \bold{110}(509):435-448.
}

\author{
Ray Bai and Malay Ghosh
}

\examples{
###############################
## Example on diabetes data. ## 
###############################
data(diabetes)
attach(diabetes)
X <- scale(diabetes$x) # Center and scale X
y <- scale(diabetes$y) # Center and scale y

################################
# Fit the NBP regression model #
################################
# Should use default of 15,000 for max.steps and 10,000 for burnin
nbp.model <- nbp(X=X, y=y, method.hyperparameters="mml", 
                max.steps=5000, burnin=2500, selection="dss")

nbp.model$beta.med     # posterior median estimates
nbp.model$a.estimate   # MML estimate of shape parameter 'a'
nbp.model$b.estimate   # MML estimate of shape parameter 'b'

nbp.model$beta.intervals         # 95 percent posterior credible intervals
nbp.model$nbp.classifications    # Variables selected

\donttest{
#
#
####################################################
# TRIM32 gene expression data analysis.            #
# Running this code will allow you to reproduce    #
# the results in Section 7 of Bai and Ghosh (2019) #
####################################################

# Load the data
data(eyedata)

# Set seed
set.seed(1)

# Center design matrix X and response vector y
X <- scale(genes, center=TRUE, scale=TRUE) # gene expression data (covariates)
y <- scale(trim32, center=TRUE, scale=TRUE) # levels of TRIM32 (response)
  
###########################
# Implement the NBP model #
###########################
nbp.model = nbp(X,y, method.hyperparameters="mml", selection="dss")

# Variables selected
active.indices <- which(nbp.model$nbp.classifications != 0) # active genes
active.estimates <- nbp.model$beta.med[active.indices]
active.CIs <- nbp.model$beta.intervals[, active.indices]


###################################
# Evaluate predictive performance #
###################################
k <- 5  # Number of folds

# Randomly select indices for the 5 different folds 
folds <- split(sample(nrow(X), nrow(X), replace=FALSE), as.factor(1:k)) 

# To store the mean square prediction error
mspe.nbp <- rep(NA, k)

for(i in 1:k){
  
  # Split data into training set and test set based on folds
  test_ind = folds[[i]]
  
  # 80 percent training set
  y.train = y[-test_ind]
  X.train = X[-test_ind, ]
  
  # Rest is test set 
  y.test = y[test_ind] 
  X.test = X[test_ind, ]
  
  # Run NBP model on training set
   nbp.train = nbp(X=X.train, y=y.train, method.hyperparameters="mml", selection="dss")
   beta.train <- nbp.train$beta.med
   
  # Obtain predictions on test set
  nbp.pred = crossprod(t(X.test), beta.train)
  
  # Compute the MSPE on test set
  mspe.nbp[i] = mean((nbp.pred - y.test)^2)
}

mean.nbp.mspe <- mean(mspe.nbp)
mean.nbp.mspe

#
#
}

}
