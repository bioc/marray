\name{maMAD}
\alias{maMAD}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Stratified MAD calculation}
\description{
This function computes the median absolute deviation (MAD) of values in \code{y} separately within values of \code{x}. It is used by \code{\link{maNormMAD}} for MAD scale normalization.
}
\usage{
maMAD(x, y, geo=TRUE, subset=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{Variables used to stratify the data.}
  \item{y}{A numeric vector.}
  \item{geo}{If \code{TRUE}, the MAD of each group is divided by the
  geometric mean of the MADs across groups (cf. Yang et al. (2002)). This allows observations to
  retain their original units.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of points used to compute the MAD. }
}

\value{
A numeric vector of MAD values. 
}


\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, and T. P. Speed (2001). Normalization for cDNA microarray data. In M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), \emph{Microarrays: Optical Technologies and Informatics}, Vol. 4266 of \emph{Proceedings of SPIE}.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, D. M. Lin, V. Peng, J. Ngai, and T. P. Speed (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research}, Vol. 30, No. 4.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}


\seealso{\code{\link{maNormMain}}, \code{\link{maNormMAD}}, \code{\link{mad}}.}

\examples{
# See examples for maNormMain.

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{robust}% __ONLY ONE__ keyword per line
