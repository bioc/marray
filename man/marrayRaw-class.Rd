\name{marrayRaw-class}
\docType{class}

\alias{marrayRaw-class}
\alias{marrayRaw}

\alias{maRf}
\alias{maRf<-}
\alias{maGf}
\alias{maGf<-}
\alias{maRb}
\alias{maRb<-}
\alias{maGb}
\alias{maGb<-}
\alias{maW}
\alias{maW<-}
\alias{maLayout}
\alias{maLayout<-}

\alias{[,marrayRaw-method}
\alias{cbind,marrayRaw-method}
%\alias{coerce,marrayRaw,marrayNorm-method}
\alias{coerce,marrayRaw,exprSet-method}
\alias{coerce,marrayRaw,marraySpots-method}
\alias{maA,marrayRaw-method}
\alias{maControls<-,marrayRaw-method}
\alias{maControls,marrayRaw-method}
\alias{maGb<-,marrayRaw,matrix-method}
\alias{maGb<-,marrayRaw,NULL-method}
\alias{maGb,marrayRaw-method}
\alias{maGf<-,marrayRaw,matrix-method}
\alias{maGf,marrayRaw-method}
%\alias{maGnames<-,marrayRaw,marrayInfo-method}
\alias{maGnames,marrayRaw-method}
\alias{maGridCol,marrayRaw-method}
\alias{maGridRow,marrayRaw-method}
%\alias{maLayout<-,marrayRaw,marrayLayout-method}
\alias{maLayout,marrayRaw-method}
\alias{maLG,marrayRaw-method}
\alias{maLR,marrayRaw-method}
\alias{maM,marrayRaw-method}
\alias{maNgc<-,marrayRaw,numeric-method}
\alias{maNgc,marrayRaw-method}
\alias{maNgr<-,marrayRaw,numeric-method}
\alias{maNgr,marrayRaw-method}
\alias{maNotes<-,marrayRaw,character-method}
\alias{maNotes,marrayRaw-method}
\alias{maNsamples,marrayRaw-method}
\alias{maNsc<-,marrayRaw,numeric-method}
\alias{maNsc,marrayRaw-method}
\alias{maNspots<-,marrayRaw,numeric-method}
\alias{maNspots,marrayRaw-method}
\alias{maNsr<-,marrayRaw,numeric-method}
\alias{maNsr,marrayRaw-method}
\alias{maPlate<-,marrayRaw-method}
\alias{maPlate,marrayRaw-method}
\alias{maPrintTip,marrayRaw-method}
\alias{maRb<-,marrayRaw,matrix-method}
\alias{maRb<-,marrayRaw,NULL-method}
\alias{maRb,marrayRaw-method}
\alias{maRf<-,marrayRaw,matrix-method}
\alias{maRf,marrayRaw-method}
\alias{maSpotCol,marrayRaw-method}
\alias{maSpotRow,marrayRaw-method}
\alias{maSub<-,marrayRaw-method}
\alias{maSub,marrayRaw-method}
%\alias{maTargets<-,marrayRaw,marrayInfo-method}
\alias{maTargets,marrayRaw-method}
\alias{maW<-,marrayRaw,matrix-method}
\alias{maW,marrayRaw-method}
\alias{print,marrayRaw-method}
\alias{show,marrayRaw-method}
\alias{summary,marrayRaw-method}

%\alias{maNgr}
%\alias{maNgc}
%\alias{maNsr}
%\alias{maNsc}
%\alias{maNspots}
%\alias{maSub}
%\alias{maPlate}
%\alias{maControls}

%\alias{maPrintTip}
%\alias{maGridRow}
%\alias{maGridCol
%\alias{maSpotRow}
%\alias{maSpotCol}

\alias{maGnames}
\alias{maGnames<-}
\alias{maTargets}
\alias{maTargets<-}
%\alias{maNotes}

\alias{maLR}
\alias{maLG}

%\alias{maA}
%\alias{maM}

\alias{maNsamples}


\title{Class "marrayRaw", classes and methods for pre-normalization cDNA microarray intensity data}

\description{This class represents pre-normalization intensity data for a batch of cDNA microarrays. A \emph{batch of arrays} consists of a collection of arrays with the same layout  (\code{"\link{marrayLayout}"}). The class contains slots for the green (Cy3) and red (Cy5) foreground and background intensities, the layout of the arrays, and descriptions of the target samples hybridized to the arrays and probe sequences spotted onto the arrays.}

\section{Objects from the Class}{
Objects can be created by calls of the form 
\code{  new('marrayRaw',}\cr
\code{    maRf      = ...., # Object of class matrix}\cr
\code{    maGf      = ...., # Object of class matrix}\cr
\code{    maRb      = ...., # Object of class matrix}\cr
\code{    maGb      = ...., # Object of class matrix}\cr
\code{    maW       = ...., # Object of class matrix}\cr
\code{    maLayout  = ...., # Object of class marrayLayout}\cr
\code{    maGnames  = ...., # Object of class marrayInfo}\cr
\code{    maTargets = ...., # Object of class marrayInfo}\cr
\code{    maNotes   = ...., # Object of class character}\cr
\code{  )}
}

\section{Slots}{
  \describe{
    \item{\code{maRf}:}{Object of class \code{"matrix"}, red foreground intensities, rows correspond to spotted probe sequences, columns to arrays in the batch.}
    \item{\code{maGf}:}{Object of class \code{"matrix"}, green foreground intensities, rows correspond to spotted probe sequences, columns to arrays in the batch. }
    \item{\code{maRb}:}{Object of class \code{"matrix"}, red background intensities, rows correspond to spotted probe sequences, columns to arrays in the batch. }
    \item{\code{maGb}:}{Object of class \code{"matrix"}, green background intensities, rows correspond to spotted probe sequences, columns to arrays in the batch. }
    \item{\code{maW}:}{Object of class \code{"matrix"}, spot quality weights, rows correspond to spotted probe sequences, columns to arrays in the batch.}
    \item{\code{maLayout}:}{Object of class \code{"\link{marrayLayout}"}, layout parameters for the cDNA microarrays.}
    \item{\code{maGnames}:}{Object of class \code{"\link{marrayInfo}"}, description of spotted probe sequences.}
    \item{\code{maTargets}:}{Object of class \code{"\link{marrayInfo}"},  description of target samples hybridized to the arrays.}
    \item{\code{maNotes}:}{Object of class \code{"character"},  any notes concerning the microarray experiments, e.g. hybridization or scanning conditions.}
  }
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "marrayRaw")}: subsetting operator for spots on the array and arrays in the batch, ensures that all slots are subset properly.}
    \item{coerce}{\code{signature(from = "marrayRaw", to = "marrayNorm")}: coerce an object of class \code{"\link{marrayRaw}"} into an object of class \code{"\link{marrayNorm}"}. }
    \item{coerce}{\code{signature(from = "marrayRaw", to = "exprSet")}:  coerce an object of class \code{"\link{marrayRaw}"} into an object of class \code{"\link{exprSet-class}"}. }
\item{coerce}{\code{signature(from = "marrayRaw", to = "marraySpots")}:  coerce an object of class \code{"\link{marrayRaw}"} into an object of class \code{"\link{marraySpots}"}. }
    \item{maA}{\code{signature(object = "marrayRaw")}: function which computes average log-intensities (base 2) A for an object of class \code{"\link{marrayRaw}"}.}
    \item{maControls<-}{\code{signature(object = "marrayRaw")}: slot assignment method.} 
    \item{maControls}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maGb}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maGb<-}{\code{signature(object = "marrayRaw", value = "matrix")}: slot assignment method.} 
    \item{maGb<-}{\code{signature(object = "marrayRaw", value = "NULL")}: slot assignment method.} 
    \item{maGf}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maGf<-}{\code{signature(object = "marrayRaw", value = "matrix")}: slot assignment method.}
    \item{maGnames}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maGnames<-}{\code{signature(object = "marrayRaw", value = "marrayInfo")}: slot assignment method.}
    \item{maGridCol}{\code{signature(object = "marrayRaw")}: method which computes a vector of grid column coordinates for each spot.}
    \item{maGridRow}{\code{signature(object = "marrayRaw")}: method which computes a vector of grid row coordinates for each spot.}
    \item{maLayout}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maLayout<-}{\code{signature(object = "marrayRaw", value = "marrayLayout")}: slot assignment method. }
    \item{maLG}{\code{signature(object = "marrayRaw")}: method which computes green log-intensities (base 2) for an object of class \code{"\link{marrayRaw}"}.}
    \item{maLR}{\code{signature(object = "marrayRaw")}: method  which computes red log-intensities (base 2) for an object of class \code{"\link{marrayRaw}"}.}
    \item{maM}{\code{signature(object = "marrayRaw")}: method which computes intensity log-ratios (base 2) M for an object of class \code{"\link{marrayRaw}"}.}
    \item{maNgc}{\code{signature(object = "marrayRaw")}: slot accessor method.}
 \item{maNgc<-}{\code{signature(object = "marrayRaw", value = "numeric")}:  slot assignment method.}
    \item{maNgr}{\code{signature(object = "marrayRaw")}: slot accessor method.}
 \item{maNgr<-}{\code{signature(object = "marrayRaw", value = "numeric")}:  slot assignment method.}
    \item{maNotes}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maNotes<-}{\code{signature(object = "marrayRaw", value = "character")}: slot assignment method.}
    \item{maNsamples}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maNsc}{\code{signature(object = "marrayRaw")}: slot accessor method.}
\item{maNsc<-}{\code{signature(object = "marrayRaw", value = "numeric")}:  slot assignment method.}
    \item{maNspots}{\code{signature(object = "marrayRaw")}: slot accessor method.}
\item{maNspots<-}{\code{signature(object = "marrayRaw", value = "numeric")}:  slot assignment method.}
    \item{maNsr}{\code{signature(object = "marrayRaw")}: slot accessor method.}
\item{maNsr<-}{\code{signature(object = "marrayRaw", value = "numeric")}:  slot assignment method.}
    \item{maPlate}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maPlate<-}{\code{signature(object = "marrayRaw")}: slot assignment method.}
    \item{maPrintTip}{\code{signature(object = "marrayRaw")}: method which computes a vector of print-tip-group indices for each spot.}
    \item{maRb}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maRb<-}{\code{signature(object = "marrayRaw", value = "matrix")}: slot assignment method.}
    \item{maRb<-}{\code{signature(object = "marrayRaw", value = "NULL")}: slot assignment method. }
    \item{maRf}{\code{signature(object = "marrayRaw")}: slot accessor method. }
    \item{maRf<-}{\code{signature(object = "marrayRaw", value = "matrix")}: slot assignment method. }
    \item{maSpotCol}{\code{signature(object = "marrayRaw")}: method which computes a vector of spot column coordinates for each spot.}
    \item{maSpotRow}{\code{signature(object = "marrayRaw")}: method which computes a vector of spot row coordinates for each spot.}
    \item{maSub}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maSub<-}{\code{signature(object = "marrayRaw")}: slot assignment method.}
    \item{maTargets}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maTargets<-}{\code{signature(object = "marrayRaw", value = "marrayInfo")}: slot assignment method. }
    \item{maW}{\code{signature(object = "marrayRaw")}: slot accessor method.}
    \item{maW<-}{\code{signature(object = "marrayRaw", value = "matrix")}: slot assignment method.}
    \item{print}{\code{signature(x = "marrayRaw")}:  print method for \code{"\link{marrayRaw}"} class.}
  }
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{marrayLayout}}, \code{\link{marrayNorm}}, \code{\link{marrayInfo}}, \code{\link{marraySpots}}, \code{\link{marrayTwo}}, \code{\link{exprSet-class}}, \code{\link{print-methods}}, \code{\link{[-methods}}, \code{\link{coerce-methods}}.}


\examples{
# Examples use swirl dataset, for description type ? swirl

data(swirl)

# Object of class marrayRaw for the 4 swirl arrays
swirl

# Object of class marrayLayout
maLayout(swirl)

# Access only the first 100 spots of the third array
swirl[1:100,3]

# Accessor methods -- How many spots on the array
maNspots(swirl)

# Density plot of log-ratios M for third array
plot(density(maM(swirl[,3])))

# Assignment methods -- Replace maNotes slot
maNotes(swirl)
maNotes(swirl)<-"This is a zebrafish microarray"
maNotes(swirl)

# Coerce methods -- Convert marrayRaw object into exprSet
e<-as(swirl, "exprSet")
e

}
\keyword{classes}
