\name{write.marray}
\alias{write.marray}

\title{Data Output}
\description{
  Calls the function write.table with predefine argument.   The entries
  in each line (row) are separated by tab. 
}

\usage{
write.marray(mraw, file="maRawResults.xls", val="maM", ...)
}

\arguments{
  \item{mraw}{the object to be written, either a \code{marrayRaw} or
    \code{marrayNorm} object.}
  \item{file}{a character string representing the file name.}
  \item{val}{a character string representing the slotNames to be written.}
  \item{\dots}{further arguments to \code{\link{write.table}}.}
}

\details{
  see \code{\link{write.table}}
}

\author{Jean Yee Hwa Yang}

\seealso{\code{\link{write.table}}, \code{\link{write.list}}}

\examples{
data(swirl)
write.marray(swirl[1:10,])
}

\keyword{file}% at least one, from doc/KEYWORDS

