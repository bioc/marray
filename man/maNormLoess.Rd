\name{maNormLoess}
\alias{maNormLoess}
\title{Intensity dependent location normalization function}
\description{
This function is used for intensity dependent location normalization, using the robust local regression function \code{\link{loess}}.  It should be used as an argument to the main normalization function \code{\link{maNormMain}}.
}
\usage{
maNormLoess(x="maA", y="maM", z="maPrintTip", w=NULL, subset=TRUE,
span=0.4, ...)
}
\arguments{
  \item{x}{Name of accessor method for spot statistics, usually \code{maA}.}
  \item{y}{Name of accessor method for spot statistics, usually \code{maM}.}
  \item{z}{Name of accessor method for spot statistic used to stratify the data, usually a layout parameter, e.g. \code{maPrintTip} or \code{maPlate}. If \code{z} is not a character, e.g. NULL,  the data are not stratified.}
  \item{w}{An optional numeric vector of weights.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of points used to compute the fits.}
  \item{span}{The argument \code{span} which controls the degree of
    smoothing in the  \code{\link{loess}} function.}
  \item{...}{Misc arguments}
}

\value{
 A function with bindings for the above arguments. This latter function takes as argument an object of class \code{"\link{marrayRaw}"}  (or possibly \code{"\link{marrayNorm}"}), and returns a vector of fitted values to be substracted from the raw log-ratios. It calls the function \code{\link{maLoess}}, which is not specific to microarray objects.
}

\references{
S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, and T. P. Speed (2001). Normalization for cDNA microarray data. In M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), \emph{Microarrays: Optical Technologies and Informatics}, Vol. 4266 of \emph{Proceedings of SPIE}.
\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, D. M. Lin, V. Peng, J. Ngai, and T. P. Speed (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research}, Vol. 30, No. 4.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}


\seealso{\code{\link{maNormMain}}, \code{\link{maLoess}}, \code{\link{loess}}.}

\examples{
# See examples for maNormMain.

}

\keyword{smooth}

