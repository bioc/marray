\name{maSelectGnames}
\alias{maSelectGnames}

\title{Select genes according to the values of a few different statistics}
\description{
  Select genes by considering the \code{\link{union}} or
  \code{\link{intersect}} of multiple statistics.
}

\usage{
maSelectGnames(statdata, crit1 = 50, crit2 = crit1, sub = TRUE, selectstat, operate = c("intersect", "union"))
}

\arguments{
  \item{statdata}{A numerical matrix where the rows corresponds to genes
    and the columns corresponds to various statistics corresponding to a
    particular gene.}
  \item{crit1}{The number of points to be selected.
   If crit1 < 1, the crit1*100\% spots with the smallest M values
   will be selected. If crit1 >= 1, the crit spots  with the
   smallest M values are selected.}  
 \item{crit2}{Similar to "crit1".   If crit2 < 1, the crit2*100\%
   spots with the largest M values will be selected. If crit2 >= 1,
   the crit2 spots with the largest M values are selected.}  
 \item{sub}{A "logical" or "numeric" vector indicating the subset of
   genes to be consider.}
 \item{selectstat}{A integer value indicating the statistics where the
   final ranking is based on.}
 \item{operate}{The operation used to combined different rankings}
}
\details{
  This functions calls \code{\link{stat.gnames}} to select say the 100
  most extreme genes from various statistics and combined the different
  gene lists by either union or intersection.
}
\value{
  A vector of numeric values.
}

\author{Jean Yee Hwa Yang}

\seealso{\code{\link{stat.gnames}}, \code{\link{order}}}

\examples{
X <- matrix(rnorm(1000), 100,10)
Xstat <- maStat(X, funNames=c("meanFun", "ttestFun"))
maSelectGnames(Xstat[,1:2], crit1=50)
}
\keyword{manip}% at least one, from doc/KEYWORDS

