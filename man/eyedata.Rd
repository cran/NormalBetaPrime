\name{eyedata}
\alias{eyedata}

\title{The Bardet-Biedl syndrome gene expression data from Scheetz et al.
(2006)}
\description{ Gene expression data from the microarray study by Scheetz et al. (2006). This data is also available in the \code{flare} package.
 }
\usage{data(eyedata)}
\format{
This data is a list that consists of the following:
\describe{
    \item{\code{x}:}{is a design matrix with 120 rows and 200 columns.}
    \item{\code{y}:}{is a numeric vector with 120 samples of gene expression levels of the TRIM32 gene.}
  }
}

\details{The data set contains gene expression data (200 genes for 120 samples) and gene expression levels of TRIM32 from the microarray experiments of mammalianeye tissue samples of Scheetz et al. (2006).
}
\source{ 
https://cran.r-project.org/web/packages/flare/
}
\references{
Scheetz, T. E., Kim K.-Y. A., Swiderski R. E., Philp A. R., Braun, T. A., Knudtson, K. L., Dorrance, A. M., DiBona, G. F., Huang, J., Casavant, T. L., Sheffield, V. C., and Stone, E. M. (2006). "Regulation of gene expression in the mammalian eye and its relevance to eye disease." \emph{Proceedings of the National Academy of Sciences}, \bold{103}(39):14429-14434.
}