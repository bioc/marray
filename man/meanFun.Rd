\name{meanFun}
\alias{meanFun}

\title{A function that calculates the average of each gene in a
  multi-slide microarray experiment}

\description{
  Returns a function that calculate \code{mean}.
}

\usage{
meanFun(y = NULL, na.rm = TRUE)
}

\arguments{
  \item{y}{vector of class labels. Must consist of integers 1 and 2.}
  \item{na.rm}{If TRUE, the NA's are omit from the analysis.}
}

\value{
   Returns a function with bindings for 'y' (class labels) that will
   calculate average across slides.
}


\author{Yee Hwa (Jean) Yang }
\seealso{\code{\link{mean}}}

\examples{
x <- rnorm(10)
y <- meanFun()
y(x)
}

\keyword{univar}% at least one, from doc/KEYWORDS

