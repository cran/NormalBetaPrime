\name{hsplus.normalmeans}
\alias{hsplus.normalmeans}

\title{Normal Means Estimation and Hypothesis Testing with the Horseshoe+ Prior}

\description{This function implements the horseshoe+ model of Bhadra et al. (2017) for obtaining a sparse estimate of \eqn{\theta = (\theta_1, ..., \theta_n)} in the normal means problem,

\deqn{X_i = \theta_i + \epsilon_i,} 

where \eqn{\epsilon_i \sim N(0, \sigma^2)}. This is achieved by placing the horseshoe+ (HS+) prior on the individual \eqn{\theta_i}'s. The sparsity parameter \eqn{\tau} can be specified \emph{a priori}, or it can be estimated from the data, either by: 1) using the estimate of sparsity level by van der Pas et al. (2014), or 2) by taking a restricted marginal maximum likelihood (REML) estimate. Multiple testing can also be performed by either thresholding the shrinkage factor in the posterior mean, or by examining the marginal 95 percent credible intervals.
}

\usage{
hsplus.normalmeans(x, tau.est=c("fixed", "est.sparsity", "reml"), tau=1, 
                  sigma2=1, var.select=c("threshold", "intervals"), 
                  max.steps = 10000, burnin=5000)}

\arguments{
  \item{x}{an \eqn{n \times 1} multivariate normal vector.}

  \item{tau.est}{The method for estimating the sparsity parameter \eqn{\tau}. If \code{"fixed"} is chosen, then \eqn{\tau} is not estimated from the data. If \code{"est.sparsity"} is chosen, the empirical Bayes estimate of sparsity level from van der Pas et al. (2014) is used. If \code{"reml"} is chosen, \eqn{\tau} is estimated from restricted marginal maximum likelihood.}

  \item{tau}{The global parameter \eqn{\tau} in the HS+ prior. Controls the sparsity of the model. Defaults to \eqn{1}. User may specify a different value for \code{tau} (\eqn{\tau > 0}). This is ignored if the method \code{"est.sparsity"} or \code{"reml"} is used to estimate \eqn{\tau}.}
  
  \item{sigma2}{The variance parameter. Defaults to 1. User needs to specify the noise parameter if it is different from 1 (\eqn{\sigma^2 > 0}). }
  
  \item{var.select}{The method of variable selection. \code{"threshold"} selects variables by thresholding the shrinkage factor in the posterior mean. \code{"intervals"} will classify \eqn{\theta_i}'s as either signals (\eqn{\theta_i \neq 0}) or as noise (\eqn{\theta_i = 0}) by examining the 95 percent posterior credible intervals.  }

\item{max.steps}{The total number of iterations to run in the Gibbs sampler. Defaults to 10,000.}

  \item{burnin}{The number of burn-in iterations for the Gibbs sampler. Defaults to 5,000.}

}

\details{
The function implements sparse estimation and multiple hypothesis testing on a multivariate normal mean vector, \eqn{\theta = (\theta_1, \ldots, \theta_n)} with the horseshoe+ prior of Bhadra et al. (2015). The full model is:

\deqn{X | \theta \sim N_n( \theta, \sigma^2 I_n),}
\deqn{\theta_i | \lambda_i \sim N(0, \sigma^2 \lambda_i^2), i = 1, \ldots, n,}
\deqn{\lambda_i \sim C^{+}(0, \tau \eta_i), i = 1, \ldots, n,}
\deqn{\eta_i \sim C^{+}(0, 1), i = 1, \ldots, n.}

 \eqn{\tau} is the global parameter that controls the sparsity of the solution. It can be estimated by the empirical Bayes estimate of the estimated sparsity (\code{"est.sparsity"}) given in van der Pas et al. (2014), or it can be estimated by restricted marginal maximum likelihood (\code{"reml"}) in the interval \eqn{[1/n, 1]}.

The posterior mean is of the form \eqn{[E(1-\kappa_i | X_i)] X_i, i = 1, \ldots, n}. The "threshold" method for variable selection is to classify \eqn{\theta_i} as signal (\eqn{\theta_i \neq 0}) if \eqn{E(1 - \kappa_i | X_i) > 1/2}.

}

\value{The function returns a list containing the following components:
\item{theta.hat}{        The posterior mean estimate of \eqn{\theta}. }
\item{theta.med}{        The posterior median estimate of \eqn{\theta}. }
\item{theta.intervals}{        The 95 percent credible intervals for all \eqn{n} components of \eqn{\theta}. }
\item{hsplus.classifications}{        An \eqn{n}-dimensional binary vector with "1" if the covariate is selected and "0" if it is deemed irrelevant.}
\item{tau.estimate}{        The estimate of the sparsity level. If user specified \code{"fixed"} for \code{tau.est}, then it simply returns the fixed \eqn{\tau}.}
}

\references{
Bhadra, A., Datta, J., Polson, N. G., and Willard, B. (2017). "The Horseshoe+ Estimator for Ultra-Sparse Signals." \emph{Bayesian Analysis}, \bold{12}(4):1105-1131.

Makalic, E., Schmidt, D. F., and Hopper, J. L. (2016). "Bayesian Robust Regression with the Horseshoe+ Estimator." \emph{AI 2016: Advances in Artificial Intelligence}, 429-440.

van der Pas, S. L., Kleijn, B. J. K., and van der Vaart, A. W. (2014). "The Horseshoe Estimator: Posterior Concentration Around Nearly Black Vectors." \emph{Electronic Journal of Statistics}, \bold{8}(2):2585-2618.

van der Pas, S. L., Szabo, B. T., and van der Vaart, A. (2017). "Adaptive posterior contraction rates for the horseshoe." \emph{Electronic Journal of Statistics}, \bold{11}(2):3196-3225.

}

\author{
Ray Bai and Malay Ghosh
}

\examples{
###################################################
## Example on synthetic data.                    ## 
## 5 percent of entries in a sparse vector theta ##
## are set equal to signal value A =7.           ##
###################################################
n <- 100
sparsity.level <- 5
A <- 7

# Initialize theta vector of all zeros
theta.true <- rep(0,n)
# Set (sparsity.level) percent of them to be A
q <- floor(n*(sparsity.level/100))
# Pick random indices of theta.true to equal A
signal.indices <- sample(1:n, size=q, replace=FALSE)

###################
# Generate data X #
###################
theta.true[signal.indices] <- A
X <- theta.true + rnorm(n,0,1)

#################################
# Run the horseshoe+ model on X #
#################################
# For optimal performance, should set max.steps=10,000 and burnin=5000.

# Fix tau at 1/n
hsplus.model <- hsplus.normalmeans(X, tau.est="fixed", tau=1/length(X), 
                                  sigma2=1, var.select="threshold", 
                                  max.steps=1000, burnin=500)

hsplus.model$theta.intervals           # 95 percent credible intervals
hsplus.model$hsplus.classifications    # Classifications
hsplus.model$tau.estimate              # Estimate of sparsity level

\donttest{
######################################################
######################################################
## Prostate cancer data analysis.                   ##
## Running this code will allow you to reproduce    ##
## the results in Section 6 of Bai and Ghosh (2018) ##
######################################################
######################################################

# Load the data
data(singh2002)
attach(singh2002)

# First 50 rows are the cancer patients,
# and the last 52 rows are the control subjects.d

prostate.data <- singh2002$x

# Perform 2-sample t-test and obtain z-scores
n <- ncol(prostate.data)
test.stat <- rep(NA,n)
z.scores <- rep(NA, n)

# Fill in the vectors
for(i in 1:n){
  test.stat[i] <- t.test(prostate.data[51:102,i], 
                        prostate.data[1:50,i])$statistic
  z.scores[i] <- qnorm(pt(test.stat[i],100))
}

###########################################
# Apply horseshoe+ model on the z-scores. #
# Here sigma2 is known with sigma2= 1     #
###########################################

hsplus.model <- hsplus.normalmeans(z.scores, tau.est="reml", var.select="threshold")

##########################################
# How many genes flagged as significant? #
##########################################
num.sig <- sum(hsplus.model$hsplus.classifications != 0)
num.sig

####################################################### 
# Estimated effect size for 10 most significant genes #
#######################################################
most.sig <- c(610,1720,332,364,914,3940,4546,1068,579,4331)
hsplus.model$theta.hat[most.sig]
}

}
