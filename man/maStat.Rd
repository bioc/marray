\name{maStat}
\alias{maStat}
\alias{widget.Stat}

\title{A function to calculates various statistics}
\description{
   \code{maStat} calculates various statistcs for genes in the array
   `expr' using the functions in `funNames'.It returns an matrix of
   values.
 }
\usage{
maStat(expr, funNames, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{expr}{An object of class "matrix", "dataframe",
    "marrayRaw", "marrayNorm" or "exprSet" that contains gene expression
    data.}    
  \item{funNames}{Either a character strings containning the names of
    functions or a list of statistics functions.}
  \item{\dots}{further arguments to the functions is 'funNames'}.
}

\value{
  A matrix of values.
}

\author{Yee Hwa (Jean) Yang}

\seealso{\code{\link{meanFun}}, \code{\link{ttestFun}}}

\examples{
 X <- matrix(rnorm(1000), 100,10)
 maStat(X, funNames=c("meanFun", "ttestFun"))
 maStat(X, funNames=list(mean=meanFun(), t=ttestFun()))
}
\keyword{htest}% at least one, from doc/KEYWORDS

