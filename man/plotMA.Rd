\name{plot}
\docType{methods}
\alias{plot.marrayRaw}
\alias{plot.marrayNorm}
\alias{addText}
\alias{addLines}
\alias{addPoints}
\alias{text,marrayRaw-method}
\alias{text,marrayNorm-method}
\alias{lines,marrayRaw-method}
\alias{lines,marrayNorm-method}
\alias{points,marrayRaw-method}
\alias{points,marrayNorm-method}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Scatter-plots for cDNA microarray spot statistics}
\description{The function \code{maPlot} or \code{plot} produces scatter-plots of
  microarray spot statistics for the classes \code{"\link{marrayRaw}"},
  \code{"\link{marrayNorm}"}. It also allows the user to highlight and
  annotate subsets of points on the plot, and display fitted curves from
  robust local regression or other smoothing procedures.}
\usage{
\method{plot}{marrayRaw}(x, xvar = "maA", yvar = "maM", zvar="maPrintTip", lines.func,text.func,legend.func, \dots)
\method{plot}{marrayNorm}(x, xvar = "maA", yvar = "maM", zvar="maPrintTip", lines.func,text.func,legend.func, \dots)
addText(object, xvar="maA", yvar="maM", subset=NULL, labels=as.character(1:length(subset)), \dots)
addPoints(object, xvar="maA", yvar="maM", subset=TRUE, \dots)
addLines(object, xvar="maA", yvar="maM", zvar="maPrintTip", subset=TRUE, \dots)
\S4method{text}{marrayRaw}(x, xvar = "maA", yvar = "maM", \dots)
\S4method{text}{marrayNorm}(x, xvar = "maA", yvar = "maM", \dots)
\S4method{lines}{marrayRaw}(x, xvar = "maA", yvar = "maM", zvar = "maPrintTip", \dots)
\S4method{lines}{marrayNorm}(x, xvar = "maA", yvar = "maM", zvar = "maPrintTip",\dots)
\S4method{points}{marrayRaw}(x, xvar = "maA", yvar = "maM", \dots)
\S4method{points}{marrayNorm}(x, xvar = "maA", yvar = "maM", \dots)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{Microarray object of class \code{"\link{marrayRaw}"}, \code{"\link{marrayNorm}"}.}
  \item{object}{Microarray object of class \code{"\link{marrayRaw}"}, \code{"\link{marrayNorm}"}.}
  \item{xvar}{Name of accessor function for the abscissa spot statistic, typically a slot name for the microarray object \code{x}, such as \code{maA}.}
  \item{yvar}{Name of accessor function for the ordinate spot statistic, typically a slot name for the microarray object \code{x}, such as \code{maM}.}
  \item{zvar}{Name of accessor method for the spot statistic used to stratify the data, typically a slot name for the microarray layout object (see \code{"\link{marrayLayout}"}) such as \code{maPlate} or a method such as \code{maPrintTip}. If \code{zvar} is \code{NULL}, the data are not stratified.}
  \item{lines.func}{Function for computing and plotting smoothed fits of \code{y} as a function of \code{x}, separately within values of \code{zvar}, e.g. \code{\link{maLoessLines}}. If \code{lines.func} is \code{NULL}, no fitting is performed.}
  \item{text.func}{Function for highlighting a subset of points, e.g., \code{\link{maText}}. If \code{text.func} is \code{NULL}, no points are highlighted.}
  \item{legend.func}{Function for adding a legend to the plot,  e.g. \code{\link{maLegendLines}}. If \code{legend.func} is \code{NULL}, there is no legend.}
  \item{subset}{logical vector or numeric values indicating the subset of points to be plotted.}
  \item{labels}{One or more character strings or expressions specifying the text to be written.}
  \item{\dots}{Optional graphical parameters, see \code{\link{par}}.}
}
\details{
This function calls the general function \code{\link{maPlot.func}},
which is not specific to microarray data. If there are more than one
array in the batch, the plot is done for the first array, by
default. Default graphical parameters are chosen for convenience using
the function \code{\link{maDefaultPar}} (e.g. color palette,  axis
labels,  plot title) but the user has the option to overwrite these
parameters at any point.
}
\references{
S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}
\author{Jean Yee Hwa Yang}
\seealso{\code{\link{maPlot.func}}, \code{\link{maDefaultPar}}, \code{\link{maLoessLines}}, \code{\link{maLegendLines}}, \code{\link{maText}}, \code{\link{plot}}, \code{\link{lowess}}, \code{\link{loess}}, \code{\link{legend}}.}
\examples{
# To see the demo type demo(marrayPlots)

# Examples use swirl dataset, for description type ? swirl
data(swirl)

# Pre-normalization MA-plot for the Swirl 93 array, with the lowess fits for 
# individual print-tip-groups. 
# - Default arguments
plot(swirl[,3])

# Lowess fit using all spots
plot(swirl[,3], zvar=NULL, legend.func=NULL)

# Loess fit using all spots
plot(swirl[,3], zvar=NULL, legend.func=maLegendLines(legend="All spots",col="green"), lines.func=maLoessLines(loess.args=list(span=0.3),col="green"))
}
\keyword{hplot}% at least one, from doc/KEYWORDS
