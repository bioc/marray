\name{rm.na}
\alias{rm.na}

\title{Remove missing values}
\description{
  Remove NA's, NAN's and INF's from a vector.
}

\usage{
rm.na(x)
}

\arguments{
  \item{x}{A numeric vector}
}

\value{
  A vector with all NA's remove.
}

\author{Jean Yang}

\examples{
x <- round(rnorm(10), 2)
x[c(2,4,5)] <- NA
x
rm.na(x)
}
\keyword{manip}% at least one, from doc/KEYWORDS

