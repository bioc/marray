\name{maBoxplot}
\alias{maBoxplot}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Boxplots for cDNA microarray spot statistics}

\description{ The function \code{maBoxplot} produces boxplots of microarray spot
statistics for the classes  \code{\link{marrayRaw}} and
\code{\link{marrayNorm}}.We encourage users to use "boxplot"  rather than "maBoxplot".  The name of the arguments have changed.}

\usage{
maBoxplot(m, x="maPrintTip", y="maM", ...)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{m}{Microarray object of class \code{"\link{marrayRaw}"} and \code{"\link{marrayNorm}"}}
  \item{x}{Name of accessor method for the spot statistic used to stratify the data, typically a slot name for the microarray layout object (see \code{"\link{marrayLayout}"}) such as \code{maPlate} or a method such as \code{maPrintTip}. If \code{x} is NULL, the data are not stratified.}
  \item{y}{Name of accessor method for the spot statistic of interest, typically a slot name for the microarray object \code{m}, such as \code{maM}.}
  \item{\dots}{Optional graphical parameters, see \code{\link{par}}.}
}

\details{
If there are more than one array in the batch, the function produces a boxplot for each array in the batch. Such plots are useful when assessing the need for between array normalization, for example, to deal with scale differences among different arrays. Default graphical parameters are chosen for convenience using the function \code{\link{maDefaultPar}} (e.g. color palette,  axis labels,  plot title) but the user has the option to overwrite these parameters at any point.
}

\references{
S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{boxplot}}, \code{\link{maDefaultPar}}.}

\examples{
## see example in boxplot
}
\keyword{hplot}% at least one, from doc/KEYWORDS

