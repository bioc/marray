\name{read.marrayRaw}
\alias{read.marrayRaw}
\alias{read.Spot}
\alias{read.GenePix}
\alias{read.SMD}
\alias{read.Agilent}
\alias{widget.marrayRaw}

\title{Create objects of class "marrayRaw"}
\description{
  This function reads in cDNA microarray data from a directory and creates objects of
  class \code{"\link{marrayRaw}"} from spot quantification data files obtained from image analysis software or databases.
}

\usage{
read.marrayRaw(fnames, path=".", name.Gf=NULL, name.Gb=NULL, name.Rf=NULL,
name.Rb=NULL,name.W=NULL, layout=NULL, gnames=NULL, targets=NULL,
notes=NULL, skip=NULL, sep="	", quote="\"", DEBUG=FALSE, ...)

read.GenePix(fnames = NULL, path = NULL, name.Gf = "F532 Median",
name.Gb ="B532 Median", name.Rf = "F635 Median", name.Rb = "B635 Median",
name.W ="Flags", layout = NULL, gnames = NULL, targets = NULL,
notes = NULL, skip=NULL, sep = "	", quote = "\"", DEBUG=FALSE, ...)

read.SMD(fnames = NULL, path = NULL, name.Gf = "Ch1 Intensity (Median)",
name.Gb = "Ch1 Background (Median)", name.Rf = "Ch2 Intensity (Median)",
name.Rb = "Ch2 Background (Median)", name.W = NULL, info.id = c("Name",
"Clone ID"), layout = NULL, gnames = NULL, targets = NULL, notes = NULL, skip = NULL, sep = "\t", quote = "\"", DEBUG=FALSE, ...)

read.Spot(fnames = NULL, path = ".", name.Gf = "Gmean", name.Gb =
"morphG", name.Rf = "Rmean", name.Rb = "morphR",name.W = NULL, layout =
NULL, gnames = NULL, targets = NULL, notes = NULL, skip = NULL, sep = "\t", quote = "\"", DEBUG=FALSE, ...)

read.Agilent(fnames = NULL, path=NULL, name.Gf = "gMedianSignal", name.Gb = "gBGMedianSignal", name.Rf = "rMedianSignal", name.Rb = "rBGMedianSignal", name.W= NULL, layout = NULL, gnames = NULL, targets = NULL, notes=NULL, skip=NULL, sep="\t", quote="\"", DEBUG=FALSE, info.id=NULL, ...)

widget.marrayRaw(ext = c("spot", "xls", "gpr"), skip = 0, sep = "\t",  quote = "\"", ...)
}

\arguments{
  \item{fnames}{a vector of character strings containing the file names of
    each spot quantification data file.  These typically end in \code{.spot}
    for the software \code{Spot} or \code{.gpr} for the software \code{GenePix}.}
  \item{path}{a character string representing the data directory. By
    default this is set to the current working directory (".").  In the
    case where \code{fnames} contains the full path name, path should be set
    to NULL.}
  \item{name.Gf}{character string for the column header for green foreground intensities.}
  \item{name.Gb}{character string for the column header for green background intensities.}
  \item{name.Rf}{character string for the column header for red foreground intensities.}
  \item{name.Rb}{character string for the column header for red background intensities.}
  \item{name.W}{character string for the column header for spot quality weights.}
  \item{layout}{object of class \code{"\link{marrayLayout}"}, containing microarray layout parameters.}
  \item{gnames}{object of class \code{"\link{marrayInfo}"} containing probe sequence information.}
  \item{targets}{object of class \code{"\link{marrayInfo}"} containing target sample information.}
  \item{notes}{object of class \code{"character"}, vector of explanatory
    text.}
  \item{info.id}{object of class \code{"character"}, vector containing
    the name of the colums of the SMD file containing oligo information
    you want to retrieve. By default, this is set to read Homo sapiens
    data. You may need to modify this argument if your are working on
    another genome.}
  \item{skip}{the number of lines of the data file to skip before beginning to read in data.}
  \item{sep}{the field separator character.  Values on each line of the
    file are separated by this character. The
    default is to read a tab delimited file.}
  \item{quote}{the set of quoting characters. By default, this is
    disabled by setting \code{quote="\""}.}
  \item{ext}{a characters string representing suffix of different image
    analysis output files.}
  \item{DEBUG}{a logical value, if TRUE, a series of echo statements will be printed.}
  \item{\dots}{further arguments to \code{\link{scan}}.}
}
\value{
  An object of class \code{"\link{marrayRaw}"}.
}
\references{\url{http://www.bioconductor.org/}.}
\author{Jean Yang, \email{yeehwa@stat.berkeley.edu}}

\seealso{\code{\link{scan}}, \code{\link{read.marrayLayout}},
  \code{\link{read.marrayInfo}} }

\examples{
datadir <- system.file("swirldata", package="marray")

## Quick guide
swirl.targets <- read.marrayInfo(file.path(datadir, "SwirlSample.txt"))
data <- read.Spot(path=datadir, targets=swirl.targets)

## Alternate commands
skip <-  grep("Row", readLines(file.path(datadir,"fish.gal"), n=100)) - 1

swirl.layout <- read.marrayLayout(ngr=4, ngc=4, nsr=22, nsc=24)

swirl.targets <- read.marrayInfo(file.path(datadir, "SwirlSample.txt"))

swirl.gnames <- read.marrayInfo(file.path(datadir, "fish.gal"),
                                info.id=4:5, labels=5, skip=skip)

x <-  maInfo(swirl.gnames)[,1]
y <- rep(0, maNspots(swirl.layout))
y[x == "control"] <- 1
slot(swirl.layout, "maControls") <- as.factor(y)

fnames <- dir(path=datadir,pattern="spot")
swirl<- read.Spot(fnames, path=datadir,
                       layout = swirl.layout,
                       gnames = swirl.gnames,
                       targets = swirl.targets)

}

\keyword{file}% at least one, from doc/KEYWORDS
