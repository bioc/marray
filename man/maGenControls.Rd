\name{maGenControls}
\alias{maGenControls}
\alias{controlCode}

\title{Generating a vector recording the control status of the spotted probe sequences. }
\description{
  ControlCode is a matrix representing certain regular expression
  pattern and the control status of the spotted probe sequences.
  This function uses `grep' searches for matches to `pattern' (its first argument)
  within the character vector `x' (second argument). 
}

\usage{
maGenControls(Gnames, controlcode, id = "ID")
}

\arguments{
  \item{Gnames}{An object of class \code{matrix}, \code{data.frame} or \code{marrayInfo}
    which contains  description of spotted probe sequences.}
  \item{controlcode}{A character matrix of n by 2 columns. The first
    column contains a few regular expression of spotted probe sequences
    and the second column contains the corresponding control status.}
  \item{id}{the column number of column name in \code{Gnames} that
    contains description of each spot on the array.}
}

\value{
  A vector of characters recording the control status of the spotted probe sequences.
}

\author{Jean Yee Hwa Yang}

\seealso{\code{\link{grep}}}

\examples{
data(swirl)
maControls(swirl) <- maGenControls(maGnames(swirl), id="Name")
table(maControls(swirl))
}
\keyword{manip}% at least one, from doc/KEYWORDS

