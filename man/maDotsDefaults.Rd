\name{maDotsDefaults}
\alias{maDotsDefaults}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Replace graphical default parameters by user supplied parameters}

\description{
This function may be used to compare default graphical parameters for microarray diagnostic plots to user supplied parameters given in \code{...}. User supplied parameters overwrite the defaults. It is used in \code{\link{maBoxplot}}, \code{\link{maPlot}}, and \code{\link{maImage}}.}


\usage{
maDotsDefaults(dots, defaults)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{dots}{List of user supplied parameters, e.g. from \code{list(...)}.}
  \item{defaults}{List of default parameters, e.g. from the function \code{\link{maDefaultPar}}.}
}

\value{
  \item{args}{List of graphical parameters.}
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{maDefaultPar}}, \code{\link{maBoxplot}}, \code{\link{maPlot}}, \code{\link{maImage}}.}

\examples{
dots<-list(xlab="X1", ylab="Y1")
defaults<-list(xlab="X1", ylab="Y2", col=2)
pars<-maDotsDefaults(dots, defaults)

do.call("plot",c(list(x=1:10), pars))

}
\keyword{misc}% at least one, from doc/KEYWORDS

