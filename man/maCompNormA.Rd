\name{maCompNormA}
\alias{maCompNormA}
\alias{maCompNormEq}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Weights for composite normalization}
\description{
This function is used for composite normalization with intensity dependent weights. The function should be used as an argument to the main normalization function \code{\link{maNormMain}}. It only applies when two normalization procedures are combined. 
}
\usage{
maCompNormA()
maCompNormEq()
}
%- maybe also `usage' for other objects documented here.

\value{
  A function which takes as arguments \code{x} and \code{n}, the spot average log-intensities A and the number of normalization procedures. This latter function returns a matrix of weights for combining two normalization procedures, rows correspond to spots and columns to normalization procedures. The weights for the first procedure are given by the empirical cumulative distribution function of the spot average log-intensities A. Note that when performing composite normalization as described in Yang et al. (2002), the first normalization procedure is the global fit and the second procedure is the within-print-tip-group fit.\cr \cr
For \code{maCompEq}, equal weights are given for each procedure.

}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, D. M. Lin, V. Peng, J. Ngai, and T. P. Speed (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research}, Vol. 30, No. 4.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}, Yee Hwa (Jean) Yang.}

\seealso{\code{\link{maNormMain}}, \code{\link{maNormLoess}}, \code{\link{ecdf}}.}

\examples{
# See examples for maNormMain
}
\keyword{smooth}% at least one, from doc/KEYWORDS

