\name{maPlot}
\alias{maPlot}


%- Also NEED an `\alias' for EACH other topic documented here.
\title{Scatter-plots for cDNA microarray spot statistics}
\description{The function \code{maPlot} produces scatter-plots of
  microarray spot statistics for the classes \code{"\link{marrayRaw}"} and \code{"\link{marrayNorm}"}. It also allows the user to highlight and annotate subsets of points on the plot, and display fitted curves from robust local regression or other smoothing procedures.
}

\usage{
maPlot(m, x="maA", y="maM", z="maPrintTip", lines.func, text.func, legend.func, ...)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{m}{Microarray object of class \code{"\link{marrayRaw}"} and \code{"\link{marrayNorm}"}.}
  \item{x}{Name of accessor function for the abscissa spot statistic, typically a slot name for the microarray object \code{m}, such as \code{maA}.}
  \item{y}{Name of accessor function for the ordinate spot statistic, typically a slot name for the microarray object \code{m}, such as \code{maM}.}
  \item{z}{Name of accessor method for the spot statistic used to stratify the data, typically a slot name for the microarray layout object (see \code{"\link{marrayLayout}"}) such as \code{maPlate} or a method such as \code{maPrintTip}. If \code{z} is NULL, the data are not stratified.}
  \item{lines.func}{Function for computing and plotting smoothed fits of \code{y} as a function of \code{x}, separately within values of \code{z}, e.g. \code{\link{maLoessLines}}. If \code{lines.func} is NULL, no fitting is performed.}
  \item{text.func}{Function for highlighting a subset of points, e.g., \code{\link{maText}}. If \code{text.func} is NULL, no points are highlighted.}
  \item{legend.func}{Function for adding a legend to the plot, e.g. \code{\link{maLegendLines}}. If \code{legend.func} is NULL, there is no legend.}
  \item{\dots}{Optional graphical parameters, see \code{\link{par}}.}
}

\details{
This function calls the general function \code{\link{maPlot.func}}, which is not specific to microarray data. If there are more than one array in the batch, the plot is done for the first array, by default. Default graphical parameters are chosen for convenience using the function \code{\link{maDefaultPar}} (e.g. color palette,  axis labels,  plot title) but the user has the option to overwrite these parameters at any point.}
}

\references{
S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{maPlot.func}}, \code{\link{maDefaultPar}}, \code{\link{maLoessLines}}, \code{\link{maLegendLines}}, \code{\link{maText}}, \code{\link{plot}}, \code{\link{lowess}}, \code{\link{loess}}, \code{\link{legend}}.}

\examples{

# To see the demo type demo(marrayPlots)

# Examples use swirl dataset, for description type ? swirl
data(swirl)

# - Default arguments
maPlot(swirl)

# Lowess fit using all spots
maPlot(swirl, z=NULL, legend.func=NULL)

# Loess fit using all spots
maPlot(swirl, z=NULL, legend.func=maLegendLines(legend="All spots",col="green"), lines.func=maLoessLines(loess.args=list(span=0.3),col="green"))


# Pre-normalization MA-plot for the Swirl 81 array, with the lowess fits for 
# individual grid columns and 1\% tails of M highlighted
defs <- maDefaultPar(swirl[, 1], x = "maA", y = "maM", z = "maGridCol")
legend.func <- do.call("maLegendLines", defs$def.legend)
lines.func <- do.call("maLowessLines", c(list(TRUE, f = 0.3), defs$def.lines))
text.func<-maText(subset=maTop(maM(swirl)[,1],h=0.01,l=0.01), labels="o", col="violet")
maPlot(swirl[, 1], x = "maA", y = "maM", z = "maGridCol", lines.func=lines.func, text.func = text.func, legend.func=legend.func, main = "Swirl array 81: pre-normalization MA-plot")

}

\keyword{hplot}% at least one, from doc/KEYWORDS

