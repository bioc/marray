\name{write.list}
\alias{write.list}

\title{Data Output}

\description{
  Writes information from a list into a text file.
}

\usage{
write.list(x, filename = "data", append = FALSE, closefile = TRUE, outfile)
}
\arguments{
  \item{x}{the list object to be written.}
  \item{filename}{a character string representing the file name.}
  \item{append}{logical; if true, the data \code{x} is appended to file
    \code{filename}.}
  \item{closefile}{logical indicating if the file connection should be closed.}
  \item{outfile}{file name or connections.}
}
\details{
  This function may be called recursively if there exists list structure
  within a list.
}

\author{Jean Yee Hwa Yang}

\seealso{\code{\link{write.table}}, \code{\link{write}}}

\examples{
data(swirl)
test <- list(A = 1:10, B= maM(swirl)[1:10,], C=list(x=1:10, y=1:4),
             D = summary(maA(swirl[,1])))
write.list(test, filename="test.txt")
}
\keyword{file}
