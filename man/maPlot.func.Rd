\name{maPlot.func}
\alias{maPlot.func}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Scatter-plots with fitted curves and text}
\description{
 This function produces scatter-plots of \code{x} vs. \code{y}. It also allows the user to highlight and annotate subsets of points on the plot, and display fitted curves from robust local regression or other smoothing procedures.
}

\usage{
maPlot.func(x, y, z, 
	lines.func = maLowessLines(subset = TRUE, f = 0.3, col = 1:length(unique(z)), lty = 1, lwd = 2.5),
    	text.func = maText(), 
	legend.func = maLegendLines(legend = as.character(unique(z)), col = 1:length(unique(z)), lty = 1, lwd = 2.5, ncol = 1),
	 ...)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{A "numeric" vector for the abscissa.}
  \item{y}{A "numeric" vector for the ordinates.}
  \item{z}{A vector of statistic used to stratify the data, smoothed curves are fitted separately within values of \code{z}}
  \item{lines.func}{A function for computing and plotting smoothed fits of \code{y} as a function of \code{x}, separately within values of \code{z}, e.g. \code{\link{maLoessLines}}.}
  \item{text.func}{A function for highlighting a subset of points, e.g., \code{maText}.}
  \item{legend.func}{A function for adding a legend to the plot, e.g. \code{\link{maLegendLines}}.}
  \item{\dots}{Optional graphical parameters, see \code{\link{par}}.}
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}


\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{maPlot}}, \code{\link{maLoessLines}}, \code{\link{maLegendLines}}, \code{\link{maText}}, \code{\link{plot}}, \code{\link{lowess}}, \code{\link{loess}}, \code{\link{legend}}.}

\examples{
# See examples for maPlot.


}
\keyword{hplot}% at least one, from doc/KEYWORDS

