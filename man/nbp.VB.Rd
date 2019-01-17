\name{nbp.VB}
\alias{nbp.VB}

\title{Variational EM Algorithm for Normal-Beta Prime Regression}

\description{This function implements the variational EM approach for the normal-beta prime (NBP) model in the standard linear regression model,

\deqn{y = X \beta + \epsilon,} 

where \eqn{\epsilon \sim N_n (0, \sigma^2 I_n)}. This is achieved by placing the normal-beta prime (NBP) prior of Bai and Ghosh (2019) on the coefficients of \eqn{\beta}. Mean field variational Bayes (MFVB) is used to approximate the posterior \eqn{\pi(\beta | y)} with an appropriate variational density \eqn{q (\beta)}. The hyperparameters can be set deterministically by the user or they may be automatically selected by marginal maximum likelihood (MML).

It is recommended that the user use the Monte Carlo implementation, \code{nbp}, rather than the variational Bayes (VB) function, \code{nbp.VB}, for more accurate point esimates and posterior approximations.
}

\usage{
nbp.VB(X, y, method.hyperparameters=c("fixed","mml"), a=0.5, b=0.5, c=1e-5, d=1e-5,  
    selection=c("dss", "intervals"), tol=0.001, n.iter=1000)}

\arguments{
  \item{X}{       \eqn{n \times p} design matrix. Should be centered.}

  \item{y}{       \eqn{n \times 1} response vector. Should be centered.}
  
  \item{method.hyperparameters}{ The method for estimating the shape parameters \eqn{(a,b)}. If \code{"mml"} is chosen, the function estimates the marginal maximum likelihood (MML) estimates of \eqn{(a,b)} by the EM algorithm described in Bai and Ghosh (2019). If \code{"fixed"} is chosen, then \eqn{(a,b)} are not estimated from the data. }

  \item{a}{ Shape parameter for \eqn{\beta'(a,b)}. The default is 0.5. The user may specify a different value for \code{a} (\eqn{a>0}). This is ignored if the method for estimating hyperparameters is \code{"mml"}.}

  \item{b}{ Shape parameter for \eqn{\beta'(a,b)}. The default is 0.5. The user may specify a different value for \code{b} (\eqn{b>0}). This is ignored if the method for estimating hyperparameters is \code{"mml"}.}
  
  \item{c}{ The shape parameter for the \eqn{IG(c,d)} prior on unknown variance parameter, \eqn{\sigma^2}. The default is \eqn{10^{-5}}. }

\item{d}{ The rate parameter for the \eqn{IG(c,d)} prior on the unknown variance parameter, \eqn{\sigma^2}. The default is \eqn{10^{-5}}. }

  \item{selection}{ The method of variable selection. \code{"dss"} implements the decoupled selection and shrinkage (DSS) method of Hahn and Carvalho (2015) to select variables. \code{"intervals"} performs variable selection by examining the 95 percent posterior credible intervals for each coefficient, \eqn{\beta_i, i = 1, \ldots, p}.}

\item{tol}{ The convergence criterion. If the absolute value of the difference between the current ELBO and the previous ELBO falls below \code{tol}, then the variational EM algorithm terminates.}

  \item{n.iter}{ The maximum number of coordinate ascent iterations to run. Defaults to 1000.}
  

}

\details{
The function implements the normal-beta prime (NBP) model of Bai and Ghosh (2019) using mean field variational Bayes (MFVB). The posterior variances and 95 percent credible intervals for each of the \eqn{p} covariates are also returned so that the user may assess uncertainty quantification. The full model is:

\deqn{Y | (X, \beta ) \sim N_n(X \beta, \sigma^2 I_n),}
\deqn{\beta_i | \omega_i^2 \sim N(0, \sigma^2 \omega_i^2), i = 1, ..., p,}
\deqn{\omega_i^2 \sim \beta'(a,b), i = 1, ..., p,}
\deqn{\sigma^2 \sim IG(c,d),}

where \eqn{\beta'(a,b)} denotes the beta prime density,

\deqn{\pi(\omega_i^2) = \frac{\Gamma(a+b)}{\Gamma(a) \Gamma(b)} (\omega_i^2)^{a-1} (1+\omega_i^2)^{-a-b}. }

}

\value{The function returns a list containing the following components:
\item{beta.hat}{ The posterior mean estimate of \eqn{\beta}. }
\item{beta.var}{ The posterior variance estimates for each \eqn{\beta_i, i=1, \ldots, p.} }
\item{beta.intervals}{ The 95 percent credible intervals for all \eqn{p} estimates in \eqn{\beta}. }
\item{nbp.classifications}{ A \eqn{p}-dimensional binary vector with "1" if the covariate is selected and "0" if it is deemed irrelevant.}
\item{sigma2.esimate}{ Estimate of unknown variance component \eqn{\sigma^2}.}
\item{a.estimate}{ MML estimate of shape parameter \eqn{a}. If \eqn{a} was fixed \emph{a priori}, returns fixed \eqn{a}.}
\item{b.estimate}{ MML estimate of shape parameter \eqn{b}. If \eqn{b} was fixed \emph{a prior}, returns fixed \eqn{b}.}
}

\references{
Bai, R. and Ghosh, M. (2019). "On the beta prime prior for scale parameters in high-dimensional Bayesian Regression models." Pre-print, arXiv:1807.06539.

Blei, D. M., Kucukelbir, A., and McAuliffe, J. D. (2017). "Variational inference: A review for statisticians." \emph{Journal of the American Statistical Association}, \bold{112}(518):859-877.

Hahn, P. R. and Carvalho, C. M. (2015). "Decoupling shrinkage and selection in Bayesian linear models: A posterior summary perspective." \emph{Journal of the American Statistical Association}, \bold{110}(509):435-448.
}

\author{
Ray Bai and Malay Ghosh
}

\examples{

#################################
## Example on synethetic data. ## 
#################################
n <- 50
p <- 80
X <- matrix(rnorm(n*p,-3,3), nrow=n, ncol=p)
beta.true <- c(rep(2,5), rep(0,p-5)) # True beta has five entries of '2' and rest '0'

X <- scale(X) # Center and scale X
y <- crossprod(t(X), beta.true) + rnorm(n)

################################
# Fit the NBP regression model #
# using variational Bayes      #
################################
nbp.model <- nbp.VB(X=X, y=y, method.hyperparameters="mml", selection="dss")

nbp.model$beta.hat     # posterior median estimates
nbp.model$beta.var     # posterior variance estimates
nbp.model$nbp.classifications    # Variables selected

}
