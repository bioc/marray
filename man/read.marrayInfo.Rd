\name{read.marrayInfo}
\alias{read.marrayInfo}
\alias{widget.marrayInfo}

\title{Create objects of class marrayInfo}
\description{
  This function creates objects of class \code{\link{marrayInfo}}. The \code{\link{marrayInfo}} class is used to store
  information regarding the target mRNA samples co-hybridized on the arrays or
  the spotted probe sequences (e.g. data frame of gene names, annotations, and other identifiers).
}
\usage{
read.marrayInfo(fname, info.id=NULL, labels=NULL, notes=fname, sep="\t",skip=0, quote="\"", ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{fname}{the name of the file that stores information on target samples or probe sequences. This is usually a file obtained from a database.}
  \item{info.id}{the column numbers in \code{fname} that contain the
    required information.}
  \item{labels}{the column number in \code{fname} which contains the names that the user would like to use to label spots or arrays (e.g. for default titles in \code{\link{maPlot}}.}
  \item{notes}{object of class character, vector of explanatory text}
    \item{sep}{the field separator character.  Values on each line of the
    file are separated by this character. The
    default is to read a tab delimited file.}
\item{skip}{the number of lines of the data file to skip before
  beginning to read data.}
\item{quote}{the set of quoting characters. By default, this is disable
  by setting `quote="\""'.}
  \item{\dots}{further arguments to \code{\link{scan}}.}
}
\value{
  An object of class \code{\link{marrayInfo}}.
}
\references{http://www.bioconductor.org/}
\author{Jean Yang, \email{yeehwa@stat.berkeley.edu}}

\examples{
datadir <- system.file("swirldata", package="marray")

## Reading target sample information
swirl.samples <- read.marrayInfo(file.path(datadir, "SwirlSample.txt"))

## Reading probe sequence information
swirl.gnames <- read.marrayInfo(file.path(datadir, "fish.gal"),
                                info.id=4:5, labels=5, skip=21)
}

\keyword{file}% at least one, from doc/KEYWORDS
\keyword{connection}
