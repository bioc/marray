\name{maText}
\alias{maText}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Highlight points on a plot}
\description{
This function may be used to highlight a subset of points on an existing plot, such as a plot produced by \code{\link{plot}}, \code{\link{maPlot}}, or \code{\link{maPlot.func}}.
}

\usage{
maText(subset=NULL, labels=as.character(1:length(subset)), ...)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{subset}{A "logical" or "numeric" vector indicating the subset of points to
  highlight.}
\item{labels}{One or more character strings or expressions specifying the text to be written.}
  \item{\dots}{Optional graphical parameters, see \code{\link{par}}.}
}

\value{
A function with bindings for \code{subset}, \code{labels}, and \code{\ldots}. This latter function takes as arguments \code{x} and \code{y}, the absissa and ordinates of points on the plot.
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{text}}, \code{\link{maPlot}}, \code{\link{maPlot.func}}.}

\examples{
# See examples for maPlot.


}
\keyword{aplot}% at least one, from doc/KEYWORDS

