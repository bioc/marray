\name{checkTargetInfo}
\alias{checkTargetInfo}
\title{Verifying the order between intensities matrix and target file information}
\description{
  Check that the foreground and backgruond intensities are stored in the
  same order as provided in the first column of target file. 
}
\usage{
checkTargetInfo(mraw)
}
\arguments{
  \item{mraw}{Object of class \code{marrayRaw} or \code{marryNorm}.}
}
\value{
  A logical value.  This function returns "TRUE" if the first column from
  the Target information is the same order as the foreground and
  backgruond intensities.  
}

\author{Yee Hwa (Jean) Yang}

\examples{
data(swirl)
tmp <- cbind(colnames(maGf(swirl)), swirl@maTargets@maInfo)
colnames(tmp)[1] <- "FileNames"
swirl@maTargets@maInfo <- tmp
checkTargetInfo(swirl)

## reorder
swirl@maTargets@maInfo <- tmp[, c(2:4, 1)]
checkTargetInfo(swirl)

## Missing Target information
checkTargetInfo(new("marrayRaw"))
}

\keyword{file}% at least one, from doc/KEYWORDS
