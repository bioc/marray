\name{numNAFun}
\alias{numNAFun}

\title{A function that calculates the number of NA of each gene in a
  multi-slide microarray experiment. }

\description{
  Returns a function that calculates the number of NA of each gene across
  multi-slide microarray experiment.
}

\usage{
numNAFun()
}

\arguments{
}


\value{
  Returns a function.
}


\author{Yee Hwa (Jean) Yang}

\seealso{\code{\link{is.na}}}

\examples{
x <- rnorm(10)
x[8] <- NA
y <- numNAFun()
y(x)
}
\keyword{univar}% at least one, from doc/KEYWORDS


