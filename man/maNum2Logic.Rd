\name{maNum2Logic}
\alias{maNum2Logic}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Convert a numeric vector of indices to a logical vector}
\description{
This function converts a numeric vector of indices to a logical
vector. It is used for subsetting purposes.
}
\usage{
maNum2Logic(n=length(subset), subset=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{the length of the logical vector to be produced.}
  \item{subset}{a numeric vector of indices. A logical vector may also
  be supplied, in which case it is also the value of the function.}
}

\value{
 a logical vector.
}


\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}


\examples{
maNum2Logic(10, 1:3)
}
\keyword{manip}% at least one, from doc/KEYWORDS
