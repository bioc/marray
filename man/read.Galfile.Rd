\name{read.Galfile}
\alias{read.Galfile}

\title{Reading GenePix Gal file}
\description{
  Reading a standard Gal file containing gene information.
}
\usage{
read.Galfile(galfile, path = ".", info.id = c("Name", "ID"),
labels = "ID", notes = "", sep = "\t", skip = NULL,   ncolumns=4, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{galfile}{a character string representing the Gal file.}
  \item{path}{a character string representing the data directory. By
    default this is set to the current working directory ("."). }
  \item{info.id}{the column numbers in `fname' that contain the required
    information.}
  \item{labels}{the column number in \code{fname} which contains the names that the user would like to use to label spots or arrays (e.g. for default titles in \code{\link{maPlot}}.}
  \item{notes}{object of class character, vector of explanatory text}
  \item{sep}{the field separator character.  Values on each line of the
    file are separated by this character. The
    default is to read a tab delimited file.}
  \item{skip}{the number of lines of the data file to skip before
    beginning to read data.}
  \item{ncolumns}{an integer representing the number of columns of
    sub-array (print-tips) on a slides.}
  \item{\dots}{further arguments to \code{\link{scan}}.}
}

\value{
  \item{gnames}{An object of class \code{\link{marrayInfo}}.}
  \item{layout}{An object of class \code{\link{marrayLayout}}.}
}

\author{Yee Hwa (Jean) Yang}

\seealso{\code{\link{read.marrayInfo}}, \code{\link{read.marrayLayout}}}

\examples{
library(marrayInput)
datadir <- system.file("data", package="marrayInput")
try <- read.Galfile(galfile="fish.gal", path=datadir)
names(try)
try$layout
try$gnames
}

\keyword{file}% at least one, from doc/KEYWORDS
