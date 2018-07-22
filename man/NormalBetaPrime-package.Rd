\name{NormalBetaPrime-package}
\alias{NormalBetaPrime-package}
\alias{NormalBetaPrime}
\docType{package}
\title{Normal-Beta Prime}
\description{
This package contains the \code{nbp} and \code{nbp.normalmeans} functions for implementing Bayesian linear regression, variable selection, sparse normal means estimation, and multiple hypothesis testing with the normal-beta prior (NBP) prior, as introduced by Bai and Ghosh (2018) <arXiv:1807.02421> and Bai and Ghosh (2018) <arXiv:1807.06539>. Functions for normal means estimation and hypothesis testing are also available for the Dirichlet-Laplace (DL) prior of Bhattacharya et al. (2015) and the horseshoe+ (HS+) prior of Bhadra et al. (2017): \code{dl.normalmeans} and \code{hsplus.normalmeans} respectively.
}

\details{

The DESCRIPTION file:
\packageDESCRIPTION{Normal-Beta Prime}
\packageIndices{Normal-Beta Prime}

This package implements the NBP model for Bayesian linear regression, variable selection, sparse normal means estimation, and multiple hypothesis testing. Our package also contains functions for sparse normal means estimation and multiple hypothesis testing for the DL and HS+ priors.
}
\author{
Ray Bai and Malay Ghosh

Maintainer: Ray Bai <raybai07@ufl.edu>
}

\references{
Bai, R. and Ghosh, M. (2018). "Large-Scale Multiple Hypothesis Testing With the Normal Beta-Prime Prior." Pre-print, arXiv:1807.02421.

Bai, R. and Ghosh, M. (2018). "On the Beta Prime Prior for Scale Parameters in High-Dimensional Regression Models." Pre-print, arXiv:1807.06539.

Bhadra, A., Datta, J., Polson, N. G., and Willard, B. (2017). "The Horseshoe+ Estimator for Ultra-Sparse Signals." \emph{Bayesian Analysis}, \bold{12}(4):1105-1131.

Bhattacharya, A., Pati, D., Pillai, N. S., and Dunson, D. B. (2015). "Dirichlet-Laplace Priors for Optimal Shrinkage." \emph{Journal of the American Statistical Association}, \bold{110}(512):1479-1490.

van der Pas, S., Scott, J., Chakraborty, A., and Bhattacharya, A. horseshoe: Implementation of the Horseshoe Prior. https://CRAN.R-project.org/package=horseshoe.

van der Pas, S. L., Kleijn, B. J. K., and van der Vaart, A. W. (2014). "The Horseshoe Estimator: Posterior Concentration Around Nearly Black Vectors." \emph{Electronic Journal of Statistics}, \bold{8}(2):2585-2618.

van der Pas, S. L., Szabo, B. T., and van der Vaart, A. (2017). "Adaptive posterior contraction rates for the horseshoe." \emph{Electronic Journal of Statistics}, \bold{11}(2):3196-3225.

}