\name{ttestFun}
\alias{ttestFun}

\title{A function that calculates the t-statistics
  of each gene in a multi-slide microarray experiment.}

\description{
  Returns a function that calculate \code{t.test}.
}
\usage{
ttestFun(y = NULL, var.equal = FALSE, alternative = c("two.sided", "less", "greater"), na.rm = TRUE)
}

\arguments{
  \item{y}{vector of class labels. Must consist of integers 1 and 2.}
  \item{var.equal}{ a logical variable indicating whether to treat the two
          variances as being equal. If `TRUE' then the pooled variance
          is used to estimate the variance otherwise the Welch
          approximation to the degrees of freedom is used.}
  \item{alternative}{ a character string specifying the alternative hypothesis,
          must be one of `"two.sided"' (default), `"greater"' or
          `"less"'.  You can specify just the initial letter.}
  \item{na.rm}{If TRUE, the NA's are omit from the analysis.}

}

\value{
  Returns a function that performs \code{t.test}.
}


\author{Yee Hwa (Jean) Yang}

\seealso{\code{\link{t.test}}}

\examples{
x <- rnorm(10)
y <- ttestFun()
y(x)
}
\keyword{univar}% at least one, from doc/KEYWORDS

