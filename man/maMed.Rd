\name{maMed}
\alias{maMed}
\title{Stratified median calculation}
\description{
This function computes the median of values in \code{y} separately within values of \code{x}. It is used by \code{\link{maNormMed}} for median location normalization.
}
\usage{
maMed(x, y, subset=TRUE)
}
\arguments{
  \item{x}{Variables used to stratify the data.}
  \item{y}{A numeric vector.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of points used to compute the median. }
}

\value{
A numeric vector of median values. 
}


\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, and T. P. Speed (2001). Normalization for cDNA microarray data. In M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), \emph{Microarrays: Optical Technologies and Informatics}, Vol. 4266 of \emph{Proceedings of SPIE}.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, D. M. Lin, V. Peng, J. Ngai, and T. P. Speed (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research}, Vol. 30, No. 4.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}


\seealso{\code{\link{maNormMain}}, \code{\link{maNormMed}}, \code{\link{median}}.}


\examples{
# See examples for maNormMain.

}

\keyword{univar}
\keyword{robust}
