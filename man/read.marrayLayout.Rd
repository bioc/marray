\name{read.marrayLayout}
\alias{read.marrayLayout}
\alias{widget.marrayLayout}

\title{Create objects of class marrayLayout}
\description{
  This function creates objects of class \code{\link{marrayLayout}} to store layout parameters for  two-color cDNA microarrays.
}
\usage{
read.marrayLayoutfname = NULL, ngr, ngc, nsr, nsc, pl.col = NULL,
                 ctl.col = NULL, sub.col = NULL, notes = fname, skip,
                 sep = "\t", quote = "\"", ...)
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{fname}{the name of the file that stores plate and control
    information. This is usually a file obtained from a database.}
  \item{ngr}{the number of rows of grids per image.}
  \item{ngc}{the number of columns of grids per image.}
  \item{nsr}{the number of rows of spots per grid.}
  \item{nsc}{the number of columns of spots per grid.}
  \item{pl.col}{the column number in \code{fname} that contains plate information.}
  \item{ctl.col}{the column number in \code{fname} that contains control information.}
    \item{sub.col}{the column number in \code{fname} that contains full ID information.}
  \item{notes}{object of class character, vector of explanatory text.}
  \item{skip}{the number of lines of the data file to skip before
    beginning to read data.} 
  \item{sep}{the field separator character.  Values on each line of the
    file are separated by this character. The
    default is to read a tab delimited file.}
  \item{quote}{the set of quoting characters. By default, this is disable
  by setting `quote="\""'.}
  \item{\dots}{further arguments to \code{\link{scan}}.}
}

\details{}
\value{
  An object of class \code{\link{marrayLayout}}.
}

\references{http://www.bioconductor.org/ }
\author{Jean Yang \email{yeehwa@stat.berkeley.edu}}
\note{}

\examples{
datadir <- system.file("swirldata", package="marray")

### Reading in control information from file
skip <-  grep("Row", readLines(file.path(datadir,"fish.gal"), n=100)) - 1
swirl.layout <- read.marrayLayout(fname=file.path(datadir,"fish.gal"), ngr=4, ngc=4,
nsr=22, nsc=24, ctl.col=4, skip=skip)

### Setting control information.
swirl.gnames <- read.marrayInfo(file.path(datadir,"fish.gal"), info.id=4:5, labels=5, skip=21)
x <-  maInfo(swirl.gnames)[,1]
y <- rep(0, maNspots(swirl.layout))
y[x == "control"] <- 1
slot(swirl.layout, "maControls") <- as.factor(y)
}

\keyword{file}% at least one, from doc/KEYWORDS
\keyword{connection}
