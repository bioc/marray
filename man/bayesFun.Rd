\name{bayesFun}
\alias{bayesFun}

\title{A function that calculates an odds ratio of each gene in a
  multi-slide microarray experiment}

\description{
  This is a wrapper function around stat.bay.est which takes
  independent, sufficient estimates of the  effect and its variance for
  each gene in a multi-slide microarray experiment and returns an odds
  ratio for each gene: log( Pr(the gene is differentially expressed) /
  Pr(the gene is not differentially expressed) ). 
}

\usage{
bayesFun(...)
}

\arguments{
  \item{...}{arguments for the fuction \code{lm.series} and \code{ebayes}}
}

\details{
}

\value{
  Returns a function.
}

\references{I. Lönnstedt and T. P. Speed.  Replicated Microarray Data.
  Statistical Sinica, Accepted, see
  http://www.stat.berkeley.edu/users/terry/zarray/Html/papersindex.html} 

\author{Yee Hwa (Jean) Yang}

\seealso{\code{lm.series},\code{ebayes}}

\examples{
  bayesFun()
}

\keyword{models}
