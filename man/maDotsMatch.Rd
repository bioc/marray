\name{maDotsMatch}
\alias{maDotsMatch}

\title{Replace default arguments of a function by user supplied values}

\description{
  This function may be used to replace default arguements for any
  functions to user supplied parameters.}


\usage{
maDotsMatch(dots, defaults)
}

\arguments{
  \item{dots}{List of user supplied argements, e.g. from \code{list(...)}.}
  \item{defaults}{List of  formal arguments of a function, e.g. from the
    function \code{\link{formals}}.} 
}

\value{
  \item{args}{List of argument of a function.}
}

\author{Jean Yee Hwa Yang}

\seealso{\code{\link{maDefaultPar}}, \code{\link{maDotsDefaults}}}

\examples{
dots<-list(x=1:10, y=11:20)
argsfun <- maDotsMatch(dots, formals(args(plot)))
do.call("plot", argsfun)
}
\keyword{misc}% at least one, from doc/KEYWORDS

