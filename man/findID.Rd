\name{findID}
\alias{findID}

\title{Find ID when given an accession number}
\description{
  Search gene ID with a vector of accession number from gene names or ID
  values.
}

\usage{
findID(text, Gnames = gnames, ID = "Name")
}

\arguments{
  \item{text}{A character strings of gene names or id names}
  \item{Gnames}{An objects of \code{marrayRaw}, \code{marrayNorm},
    \code{exprSet} or data.frame of gene names information.}
  \item{ID}{The column of ID corresponding to 'text'}
}

\value{
  A numeric vector the gene ID.
}
\author{Yee Hwa (Jean) Yang}

\seealso{\code{\link{grep}}}

\examples{
library(marrayInput)
data(swirl)
findID("fb24a09", swirl, ID="ID")
findID("geno1", swirl)
}

\keyword{manip}% at least one, from doc/KEYWORDS
