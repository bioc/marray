\name{maLegendLines}
\alias{maLegendLines}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Add a legend to a plot}
\description{
This function may be used to add a legend for lines in plots such as those produced by \code{\link{plot}}, \code{\link{maPlot}}, or \code{\link{maPlot.func}}.
}

\usage{
maLegendLines(legend="", col=2, lty=1, lwd=2.5, ncol=1, ...)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{legend}{A vector of "character" strings to appear in the legend.}
  \item{col}{Line colors for the legend.}
  \item{lty}{Line types for the legend.}
  \item{lwd}{Line widths for the legend.}
  \item{ncol}{The number of columns in which to set the legend items (default is 1, a vertical legend).}
  \item{\dots}{Optional graphical parameters, see \code{\link{par}}.}
}


\value{
A function with bindings for \code{legend}, \code{col}, \code{lty}, \code{lwd}, \code{ncol}, and \code{\ldots}. This latter function takes as arguments \code{x} and \code{y}, the coordinates for the location of the legend on the plot, and it adds the legend to the current plot.
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}


\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{legend}}, \code{\link{maPlot}}, \code{\link{maPlot.func}}.}

\examples{
# See examples for maPlot.

}
\keyword{aplot}% at least one, from doc/KEYWORDS

