\name{maNorm2D}
\alias{maNorm2D}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{2D spatial location normalization function}
\description{
This function is used for 2D spatial location normalization, using the robust local regression function \code{\link{loess}}.  It should be used as an argument to the main normalization function \code{\link{maNormMain}}.
}

\usage{
maNorm2D(x="maSpotRow", y="maSpotCol", z="maM", g="maPrintTip", w=NULL,
subset=TRUE, span=0.4, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{Name of accessor method for spot row coordinates, usually \code{maSpotRow}.}
  \item{y}{Name of accessor method for spot column coordinates, usually \code{maSpotCol}.}
  \item{z}{Name of accessor method for spot statistics, usually the log-ratio \code{maM}.}
  \item{g}{Name of accessor method for print-tip-group indices, usually \code{maPrintTip}.}
  \item{w}{An optional numeric vector of weights.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of points used to
  compute the fits. }
\item{span}{The argument \code{span} which controls the degree of
  smoothing in the  \code{\link{loess}} function.}
  \item{...}{Misc arguments}
}

\details{
The spot statistic named in \code{z} is regressed on spot row and column coordinates, separately within print-tip-group, using the \code{\link{loess}} function.
}

\value{
 A function with bindings for the above arguments. This latter function takes as argument an object of class \code{"\link{marrayRaw}"} (or possibly \code{"\link{marrayNorm}"}), and returns a vector of fitted values to be substracted from the raw log-ratios. It calls the function \code{\link{ma2D}}, which is not specific to microarray objects.
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}


\seealso{\code{\link{maNormMain}}, \code{\link{ma2D}}, \code{\link{loess}}.}

\examples{
# See examples for maNormMain.

}
\keyword{smooth}% at least one, from doc/KEYWORDS

