\name{ma2D}
\alias{ma2D}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Stratified bivariate robust local regression}
\description{
This function performs robust local regression of a variable \code{z} on predictor variables \code{x} and \code{y}, separately within values of a fourth variable \code{g}. It is used by \code{\link{maNorm2D}} for 2D spatial location normalization.
}
\usage{
ma2D(x, y, z, g, w=NULL, subset=TRUE, span=0.4, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector of predictor variables.}
  \item{y}{A numeric vector of predictor variables.}
  \item{z}{A numeric vector of responses.}
  \item{g}{Variables used to stratify the data.}
  \item{w}{An optional numeric vector of weights.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of points used to compute the fits. }
  \item{span}{The argument \code{span} which controls the degree of
    smoothing in the  \code{\link{loess}} function.}
  \item{...}{Misc arguments}
}
\details{
\code{z} is regressed on \code{x} and \code{y}, separately within values of \code{g} using the \code{\link{loess}} function.
}
\value{
A numeric vector of fitted values.
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}


\seealso{\code{\link{maNormMain}}, \code{\link{maNorm2D}}, \code{\link{loess}}.}


\examples{
# See examples for maNormMain.
}
\keyword{smooth}% at least one, from doc/KEYWORDS

