\name{maDefaultPar}
\alias{maDefaultPar}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Default graphical parameters for microarray objects}
\description{
This function returns default graphical parameters for microarray objects. The parameters may be passed as arguments to the functions \code{\link{maBoxplot}} and \code{\link{maPlot}}.
}

\usage{
maDefaultPar(m, x, y, z)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{m}{Microarray object of class \code{"\link{marrayRaw}"}, \code{"\link{marrayNorm}"}, or \code{"\link{marrayTwo}"}.}
  \item{x}{Name of accessor method for the abscissa spot statistic, typically a slot name for the microarray object \code{m}, such as \code{maA}.}
  \item{y}{Name of accessor method for the ordinate spot statistic, typically a slot name for the microarray object \code{m}, such as \code{maM}.}
  \item{z}{Name of accessor method for the spot statistic used to stratify the data, typically a slot name for the microarray layout object (see \code{"\link{marrayLayout}"}) such as \code{maPlate} or a method such as \code{maPrintTip}.}

}

\value{
A list with elements
  \item{def.box}{default graphical parameters for \code{\link{maBoxplot}}.}
  \item{def.plot}{default graphical parameters for \code{\link{maPlot}}.}
  \item{def.lines}{default graphical parameters for functions such as \code{\link{maLoessLines}} used in \code{\link{maPlot}}.}  
  \item{def.legend}{default graphical parameters for functions such as \code{\link{maLegendLines}} used in \code{\link{maPlot}}.}   
  \item{def.text}{default graphical parameters for functions such as \code{\link{maText}} used in \code{\link{maPlot}}.} 

}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}


\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{maBoxplot}}, \code{\link{maPlot}}, \code{\link{maLegendLines}}, \code{\link{maLoessLines}}, \code{\link{maText}}, \code{\link{maDotsDefaults}}.}

\examples{
# See examples for maPlot.

}
\keyword{iplot}% at least one, from doc/KEYWORDS
\keyword{dplot}% __ONLY ONE__ keyword per line
\keyword{color}% __ONLY ONE__ keyword per line
