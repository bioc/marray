\name{write.xls}
\alias{write.xls}

\title{Data Output}
\description{
  Calls the function write.table with predefine argument.   The entries
  in each line (row) are separated by tab. 
}

\usage{
write.xls(res, file = "test.xls", ...)
}

\arguments{
  \item{res}{the object to be written, typically a data frame.  If not, it
          is attempted to coerce \code{x} to a data frame.}
  \item{file}{a character string representing the file name.}
  \item{\dots}{further arguments to \code{\link{write.table}}.}
}
\details{
  see \code{\link{write.table}}
}

\author{Jean Yee Hwa Yang}

\seealso{\code{\link{write.table}}, \code{\link{write.list}}}

\examples{
data(swirl)
write.xls(maM(swirl)[1:10,], "normM.xls")
}

\keyword{file}% at least one, from doc/KEYWORDS

