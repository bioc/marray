\name{maLoessLines}
\alias{maLoessLines}
\alias{maLowessLines}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Add smoothed fits to a plot}
\description{
This function may be used to compute and plot loess or lowess fits for an existing plot. The plot can be produced by \code{\link{plot}}, \code{\link{maPlot}}, or \code{\link{maPlot.func}}.
}

\usage{
maLoessLines(subset=TRUE, weights=NULL, loess.args=list(span = 0.4,
degree=1, family="symmetric", control=loess.control(trace.hat =
"approximate", iterations=5, surface="direct")), col=2, lty=1, lwd=2.5, ...)

maLowessLines(subset = TRUE, f = 0.3, col = 2, lty = 1, lwd = 2.5, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{subset}{A "logical" or "numeric" vector indicating the subset of points used to compute the fits. }
  \item{weights}{Optional "numeric" vector of weights -- for \code{maLoessLines} only.} 
  \item{loess.args}{List of optional arguments for the \code{\link{loess}} functions -- for \code{maLoessLines} only.}
  \item{f}{The smoother span for the \code{\link{lowess}} function -- for \code{maLowessLines} only.}
  \item{col}{The fitted line colors.}
  \item{lty}{The fitted line types.}
  \item{lwd}{The fitted line widths.}
  \item{\dots}{Optional graphical parameters, see \code{\link{par}}.}
}

\value{
A function with bindings for \code{subset}, \code{weights}, \code{loess.args},  \code{col}, \code{lty}, \code{lwd}, and \code{\ldots}. This latter function takes as arguments \code{x} and \code{y}, the abscissa and ordinates of points on the plot, and \code{z} a vector of discrete values used to stratify the points. Loess (or lowess) fits are performed separately within values of \code{z}.
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}


\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{loess}}, \code{\link{lowess}}, \code{\link{maPlot}}, \code{\link{maPlot.func}}.}

\examples{
# See examples for maPlot.

}

\keyword{aplot}% at least one, from doc/KEYWORDS
\keyword{smooth}% __ONLY ONE__ keyword per line
