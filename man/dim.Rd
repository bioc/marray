\name{dim}
\alias{dim.marrayRaw}
\alias{dim.marrayNorm}
\alias{dim.marrayInfo}

\title{Retrieve the Dimensions of an marrayRaw, marrayNorm or marrayInfo Object}
\description{
  Retrieve the number of rows (genes) and columns (arrays) for an
  marrayRaw, marrayNorm or marrayInfo object.
}
\usage{
\method{dim}{marrayRaw}(x)
}
\arguments{
  \item{x}{an object of class \code{marrayRaw}, \code{marrayNorm} or \code{marrayInfo}}
}
\details{
Microarray data objects share many analogies with ordinary matrices in
which the rows correspond to spots or genes and the columns to arrays. 
These methods allow one to extract the size of microarray data objects
in the same way that one would do for ordinary matrices. 

A consequence is that row and column commands \code{nrow(x)}, \code{ncol(x)} and so on also work.
}
\value{
Numeric vector of length 2.
The first element is the number of rows (genes) and the second is the number of columns (arrays).
}
\author{modified from Gordon Smyth's function}

\seealso{
  \code{\link[base]{dim}} in the base package.
}
\examples{
M <- A <- matrix(11:14,4,2)
rownames(M) <- rownames(A) <- c("a","b","c","d")
colnames(M) <- colnames(A) <- c("A1","A2")
MA <- new("marrayNorm", maM=M,maA=A)
dim(MA)
dim(M)
}
\keyword{array}
