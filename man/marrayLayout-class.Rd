\name{marrayLayout-class}
\docType{class}

\alias{marrayLayout-class}
\alias{marrayLayout}

\alias{maNgr}
\alias{maNgr<-}
\alias{maNgc}
\alias{maNgc<-}
\alias{maNsr}
\alias{maNsr<-}
\alias{maNsc}
\alias{maNsc<-}
\alias{maNspots}
\alias{maNspots<-}
\alias{maSub}
\alias{maSub<-}
\alias{maPlate}
\alias{maPlate<-}
\alias{maControls}
\alias{maControls<-}

\alias{[,marrayLayout-method}
\alias{maControls<-,marrayLayout-method}
\alias{maControls,marrayLayout-method}
\alias{maGridCol,marrayLayout-method}
\alias{maGridRow,marrayLayout-method}
\alias{maLayout<-,marrayRaw,marrayLayout-method}
\alias{maLayout<-,marrayNorm,marrayLayout-method}
\alias{maNgc<-,marrayLayout,numeric-method}
\alias{maNgc,marrayLayout-method}
\alias{maNgr<-,marrayLayout,numeric-method}
\alias{maNgr,marrayLayout-method}
\alias{maNotes<-,marrayLayout,character-method}
\alias{maNotes,marrayLayout-method}
\alias{maNsc<-,marrayLayout,numeric-method}
\alias{maNsc,marrayLayout-method}
\alias{maNspots<-,marrayLayout,numeric-method}
\alias{maNspots,marrayLayout-method}
\alias{maNsr<-,marrayLayout,numeric-method}
\alias{maNsr,marrayLayout-method}
\alias{maPlate<-,marrayLayout-method}
\alias{maPlate,marrayLayout-method}
\alias{maPrintTip,marrayLayout-method}
\alias{maSpotCol,marrayLayout-method}
\alias{maSpotRow,marrayLayout-method}
\alias{maSub<-,marrayLayout,logical-method}
\alias{maSub<-,marrayLayout,numeric-method}
\alias{maSub,marrayLayout-method}
\alias{print,marrayLayout-method}
\alias{show,marrayLayout-method}

%\alias{maNotes}

\alias{maPrintTip}
\alias{maGridRow}
\alias{maGridCol}
\alias{maSpotRow}
\alias{maSpotCol}


\title{Class "marrayLayout", classes and methods for layout parameters of cDNA microarrays}
\description{This class is used to keep track of important layout
  parameters for two-color cDNA microarrays. It contains slots for: the
  total number of spotted probe sequences on the array, the dimensions
  of the spot and grid matrices, the plate origin of the probes,
  information on spotted control sequences (e.g. probe sequences which
  should have equal abundance in the two target samples, such as
  housekeeping genes). The terms \emph{print-tip-group}, \emph{grid},
  \emph{spot matrix}, and \emph{sector} are used interchangeably and
  refer to a set of spots printed using the same print-tip.} 


\section{Objects from the Class}{
Objects can be created by calls of the form \code{  new('marrayLayout',}\cr
\code{    maNgr      = ...., # Object of class numeric}\cr
\code{    maNgc      = ...., # Object of class numeric}\cr
\code{    maNsr      = ...., # Object of class numeric}\cr
\code{    maNsc      = ...., # Object of class numeric}\cr
\code{    maNspots   = ...., # Object of class numeric}\cr
\code{    maSub      = ...., # Object of class logical}\cr
\code{    maPlate    = ...., # Object of class factor}\cr
\code{    maControls = ...., # Object of class factor}\cr
\code{    maNotes    = ...., # Object of class character}\cr
\code{  )}
}

\section{Slots}{
  \describe{
    \item{\code{maNgr}:}{Object of class \code{"numeric"}, number of rows for the grid matrix.}
    \item{\code{maNgc}:}{Object of class \code{"numeric"}, number of columns for the grid matrix.}
    \item{\code{maNsr}:}{Object of class \code{"numeric"}, number of rows for the spot matrices.}
    \item{\code{maNsc}:}{Object of class \code{"numeric"}, number of columns for the spot matrices.}
    \item{\code{maNspots}:}{Object of class \code{"numeric"}, total number of spots on the array, equal to \eqn{maNgr x maNgc x maNsr x maNsc}.}
    \item{\code{maSub}:}{Object of class \code{"logical"},  indicating which spots are currently being considered.}
    \item{\code{maPlate}:}{Object of class \code{"factor"}, recording the plate origin of the spotted probe sequences.}
    \item{\code{maControls}:}{Object of class \code{"factor"}, recording the control status of the spotted probe sequences.}
    \item{\code{maNotes}:}{Object of class \code{"character"}, any notes concerning the microarray layout, e.g., printing conditions.}
  }
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "marrayLayout")}: subsetting operator for spots on the array, ensures that all slots are subset properly.}
    \item{maControls<-}{\code{signature(object = "marrayLayout")}: slot assignment method.} 
    \item{maControls}{\code{signature(object = "marrayLayout")}: slot accessor method.}
    \item{maGridCol}{\code{signature(object = "marrayLayout")}: method which computes a vector of grid column coordinates for each spot.}
    \item{maGridRow}{\code{signature(object = "marrayLayout")}: method which computes a vector of grid row coordinates for each spot.}
    \item{maLayout<-}{\code{signature(object = "marrayRaw", value = "marrayLayout")}: slot assignment method.} 
    \item{maLayout<-}{\code{signature(object = "marrayNorm", value = "marrayLayout")}: slot assignment method.} 
    \item{maNgc}{\code{signature(object = "marrayLayout")}: slot accessor method.}
 \item{maNgc<-}{\code{signature(object = "marrayLayout", value = "numeric")}: slot assignment method.}
    \item{maNgr}{\code{signature(object = "marrayLayout")}: slot accessor method.}
\item{maNgr<-}{\code{signature(object = "marrayLayout", value = "numeric")}: slot assignment method.}
    \item{maNotes}{\code{signature(object = "marrayLayout")}: slot accessor method.}
    \item{maNotes<-}{\code{signature(object = "marrayLayout", value = "character")}: slot assignment method.} 
    \item{maNsc}{\code{signature(object = "marrayLayout")}: slot accessor method.}
    \item{maNsc<-}{\code{signature(object = "marrayLayout", value = "numeric")}: slot assignment method.}
    \item{maNspots}{\code{signature(object = "marrayLayout")}: slot accessor method.}
 \item{maNspots<-}{\code{signature(object = "marrayLayout", value = "numeric")}: slot assignment method.}
    \item{maNsr}{\code{signature(object = "marrayLayout")}: slot accessor method.}
\item{maNsr<-}{\code{signature(object = "marrayLayout", value = "numeric")}: slot assignment method.}
    \item{maPlate}{\code{signature(object = "marrayLayout")}: slot accessor method.}
    \item{maPlate<-}{\code{signature(object = "marrayLayout")}: slot assignment method.} 
    \item{maPrintTip}{\code{signature(object = "marrayLayout")}: method which computes a vector of print-tip-group indices for each spot.}
    \item{maSpotCol}{\code{signature(object = "marrayLayout")}: method which computes a vector of spot column coordinates for each spot.}
    \item{maSpotRow}{\code{signature(object = "marrayLayout")}: method which computes a vector of spot row coordinates for each spot.}
    \item{maSub}{\code{signature(object = "marrayLayout")}: slot accessor method.}
    \item{maSub<-}{\code{signature(object = "marrayLayout", value = "logical")}: slot assignment method.} 
    \item{maSub<-}{\code{signature(object = "marrayLayout", value = "numeric")}: slot assignment method.} 
    \item{print}{\code{signature(x = "marrayLayout")}: print method for \code{"\link{marrayLayout}"} class.}
  }
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{marrayRaw}}, \code{\link{marrayNorm}},
  \code{\link{marrayInfo}} and \code{\link{[-methods}}}.} 

\examples{
## See marrayRaw

}

\keyword{classes}
