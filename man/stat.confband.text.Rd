\name{stat.confband.text}

\alias{stat.confband.text}

\title{Rank genes according to the value of a statistic.}
  
\description{
  Select values based on intensities binning.
}

\usage{
stat.confband.text(M, A, crit1=0.025, crit2=crit1, nclass=5)
}

\arguments{
 \item{A}{a vector giving the x-coordinates of the points in the scatter
          plot. In the microarray context, this could be a vector of
          average log intensities. ie A} 
 \item{M}{a vector giving the y-coordinates of the points in the scatter
   plot. In the microarray context, this could be a vector of log
   intensity ratios.} 
 \item{crit1}{The number of points to be selected.
   If crit1 < 1, the crit1*100\% spots with the smallest M values
   will be selected. If crit1 >= 1, the crit spots  with the
   smallest M values are selected.}  
 \item{crit2}{Similar to "crit1".   If crit2 < 1, the crit2*100\%
   spots with the largest M values will be selected. If crit2 >= 1,
   the crit2 spots with the largest M values are selected.}  
 \item{nclass}{A single number giving the approximate number of
   intensity depedent groups to consider.} 
% \item{labs}{ one or more character strings or expressions specifying the
%   text to be written.  If this string is not specified, by
%   default the index of the vector `M' will be used.}
% \item{output}{logical, defaulting to `FALSE'. If `TRUE' a vector
%   containning the index to the vector `M' that are  
%   highlighted.} 
}


\value{
  A vector of selected spot index.
}

\examples{
library(marray)
data(swirl)
aveA <- apply(maA(swirl), 1, mean.na)
aveM <- apply(maM(swirl), 1, mean.na)
stat.confband.text(aveM, aveA, crit1=20, crit2=50, nclass=5)
}

\seealso{ \code{\link{stat.gnames}}}

\keyword{manip}
