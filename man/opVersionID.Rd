\name{opVersionID}
\alias{opVersionID}

\title{Determine the operon oligo set ID}
\description{
  This functions looks the operon ID and determine whether it belongs to 
  "Human Genome Oligo Set V1", "Human Genome Oligo Set V2", "Mouse
  Genome Oligo Set V1" or "Mouse Genome Oligo Set V2". 
}

\usage{
opVersionID(opID)
}

\arguments{
  \item{opID}{A character strings representing operon ID}
}

\value{
  A value "operonh1", "operonh2", "operonm1" or "operonm2" to represents
  "Human Genome Oligo Set V1", "Human Genome Oligo Set V2", "Mouse
   Genome Oligo Set V1" or "Mouse Genome Oligo Set V2". 
}

\references{\url{http://oparray.operon.com/}}
\author{Jean Yee Hwa Yang}

\seealso{\code{\link{URLstring}}, \code{\link{htmlPage}}}

\examples{
 opVersionID("M000205_01")
 URLstring[opVersionID("M000205_01")]
}

\keyword{manip}% at least one, from doc/KEYWORDS
