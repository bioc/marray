\name{read.marrayRaw}
\alias{read.marrayRaw}
\alias{read.Spot}
\alias{read.GenePix}
\alias{read.SMD}
\alias{widget.marrayRaw}

\title{Create objects of class "marrayRaw"}
\description{
  This function reads in cDNA microarray data from a directory and creates objects of
  class \code{"\link{marrayRaw}"} from spot quantification data files obtained from image analysis software or databases.
}

\usage{
read.marrayRaw(fnames, path=".", name.Gf, name.Gb=NULL, name.Rf,
name.Rb=NULL,name.W=NULL, layout=NULL, gnames=NULL, targets=NULL,
notes=NULL, skip=0, sep="	", quote="\"",...)

read.GenePix(fnames = NULL, path = ".", name.Gf = "F532 Mean",
name.Gb ="B532 Median", name.Rf = "F635 Mean", name.Rb = "B635 Median",
name.W =NULL, layout = NULL, gnames = NULL, targets = NULL,
notes = NULL, skip=0, sep = "	", quote = "\"", ...)

read.SMDfnames = NULL, path = ".", name.Gf = "CH1I_MEAN",
                 name.Gb = "CH1B_MEDIAN", name.Rf = "CH2I_MEAN",
                 name.Rb = "CH2B_MEDIAN", name.W = NULL, layout =
                 NULL, gnames = NULL, targets = NULL, notes = NULL,
                 skip = 0, sep = "\t", quote = "", ...)

read.Spot(fnames = NULL, path = ".", name.Gf = "Gmean",
name.Gb = "morphG", name.Rf = "Rmean", name.Rb = "morphR",
name.W = NULL, layout = NULL, gnames = NULL, targets = NULL,
notes = NULL, skip = 0, sep = "	", quote = "\"", ...)

widget.marrayRawext = c("spot", "xls", "gpr"), skip = 0, sep = "\t",
                 quote = "", ...)
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
  \item{notes}{object of class \code{"character"}, vector of explanatory text.}
  \item{skip}{the number of lines of the data file to skip before beginning to read in data.}
  \item{sep}{the field separator character.  Values on each line of the
    file are separated by this character. The
    default is to read a tab delimited file.}
  \item{quote}{the set of quoting characters. By default, this is
    disabled by setting \code{quote="\""}.}
  \item{\dots}{further arguments to \code{\link{scan}}.}
}
\details{
  
}
\value{
  An object of class \code{"\link{marrayRaw}"}.
}
\references{\url{http://www.bioconductor.org/}.}
\author{Jean Yang, \email{yeehwa@stat.berkeley.edu}}

\seealso{\code{\link{scan}}, \code{\link{read.marrayLayout}},
  \code{\link{read.marrayInfo}} }

\examples{
datadir <- system.file("data", package="marrayInput")

skip <-  grep("Row", readLines(file.path(datadir,"fish.gal"), n=100)) - 1

swirl.layout <- read.marrayLayout(ngr=4, ngc=4, nsr=22, nsc=24)

swirl.targets <- read.marrayInfo(file.path(datadir, "SwirlSample.txt"))

swirl.gnames <- read.marrayInfo(file.path(datadir, "fish.gal"),
                                info.id=4:5, labels=5, skip=skip)

x <-  maInfo(swirl.gnames)[,1]
y <- rep(0, maNspots(swirl.layout))
y[x == "control"] <- 1
slot(swirl.layout, "maControls") <- as.factor(y)

fnames <- dir(path=datadir,pattern=paste("*", "spot", sep="\."))
swirl<- read.Spot(fnames, path=datadir,
                       layout = swirl.layout,
                       gnames = swirl.gnames,
                       targets = swirl.targets)

}

\keyword{file}% at least one, from doc/KEYWORDS
