\name{maNormMAD}
\alias{maNormMAD}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{MAD scale normalization function}
\description{
This function is used for scale normalization using the median absolute deviation (MAD) of intensity log-ratios for a group of spots. It can be used for within or between array normalization. The function should be used as an argument to the main normalization function \code{\link{maNormMain}}.
}

\usage{
maNormMAD(x=NULL, y="maM", geo=TRUE, subset=TRUE)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{Name of accessor function for spot statistic used to stratify the data, usually a layout parameter, e.g. \code{maPrintTip} or \code{maPlate}. If \code{x} is not a character, e.g. NULL, the data are not stratified.}
  \item{y}{Name of accessor function for spot statistics, usually
  \code{maM}.}
  \item{geo}{If \code{TRUE}, the MAD of each group is divided by the
  geometric mean of the MADs across groups (cf. Yang et al. (2002)). This allows observations to
  retain their original units.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of points used to compute the scale normalization values. }
}

\value{
A function with bindings for the above arguments. This latter function takes as argument an object of class \code{"\link{marrayRaw}"}  (or possibly \code{"\link{marrayNorm}"}), and returns a vector of values used to scale the location normalized log-ratios. It calls the function \code{\link{maMAD}}, which is not specific to microarray objects.
}


\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, and T. P. Speed (2001). Normalization for cDNA microarray data. In M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), \emph{Microarrays: Optical Technologies and Informatics}, Vol. 4266 of \emph{Proceedings of SPIE}.
\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, D. M. Lin, V. Peng, J. Ngai, and T. P. Speed (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research}, Vol. 30, No. 4.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}


\seealso{\code{\link{maNormMain}}, \code{\link{maMAD}}, \code{\link{mad}}.}

\examples{
# See examples for maNormMain.

}

\keyword{univar}% at least one, from doc/KEYWORDS
\keyword{robust}% __ONLY ONE__ keyword per line
