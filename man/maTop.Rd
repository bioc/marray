\name{maTop}
\alias{maTop}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Identify extreme values}
\description{
This function determines which values in a numeric vector are above or below user supplied cut-offs.
}
\usage{
maTop(x, h=1, l=1)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A "numeric" vector.}
  \item{h}{A "numeric", upper cut-off.}
  \item{l}{A "numeric", lower cut-off.}
}

\value{
A "logical" vector indicating which entries are above or below the cut-offs.}
}

\references{
S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
 }


\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{maPlot}}, \code{\link{maImage}}, \code{\link{quantile}}.}

\examples{
# See examples for maPlot.

}
\keyword{misc}% at least one, from doc/KEYWORDS
