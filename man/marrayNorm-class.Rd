\name{marrayNorm-class}
\docType{class}

\alias{marrayNorm-class}
\alias{marrayNorm}

\alias{maA}
\alias{maA<-}
\alias{maM}
\alias{maM<-}
\alias{maMloc}
\alias{maMloc<-}
\alias{maMscale}
\alias{maMscale<-}

\alias{[,marrayNorm-method}
\alias{cbind,marrayNorm-method}
\alias{coerce,marrayRaw,marrayNorm-method}
%\alias{coerce,marrayNorm,exprSet-method}
%\alias{coerce,marrayNorm,marraySpots-method}
\alias{maA<-,marrayNorm,matrix-method}
\alias{maA,marrayNorm-method}
\alias{maControls<-,marrayNorm-method}
\alias{maControls,marrayNorm-method}
%\alias{maGnames<-,marrayNorm,marrayInfo-method}
\alias{maGnames,marrayNorm-method}
\alias{maGridCol,marrayNorm-method}
\alias{maGridRow,marrayNorm-method}
%\alias{maLayout<-,marrayNorm,marrayLayout-method}
\alias{maLayout,marrayNorm-method}
\alias{maMloc<-,marrayNorm,matrix-method}
\alias{maMloc,marrayNorm-method}
\alias{maM<-,marrayNorm,matrix-method}
\alias{maM,marrayNorm-method}
\alias{maMscale<-,marrayNorm,matrix-method}
\alias{maMscale,marrayNorm-method}
\alias{maNgc<-,marrayNorm,numeric-method}
\alias{maNgc,marrayNorm-method}
\alias{maNgr<-,marrayNorm,numeric-method}
\alias{maNgr,marrayNorm-method}
\alias{maNormCall,marrayNorm-method}
\alias{maNotes<-,marrayNorm,character-method}
\alias{maNotes,marrayNorm-method}
\alias{maNsamples,marrayNorm-method}
\alias{maNsc<-,marrayNorm,numeric-method}
\alias{maNsc,marrayNorm-method}
\alias{maNspots<-,marrayNorm,numeric-method}
\alias{maNspots,marrayNorm-method}
\alias{maNsr<-,marrayNorm,numeric-method}
\alias{maNsr,marrayNorm-method}
\alias{maPlate<-,marrayNorm-method}
\alias{maPlate,marrayNorm-method}
\alias{maPrintTip,marrayNorm-method}
\alias{maSpotCol,marrayNorm-method}
\alias{maSpotRow,marrayNorm-method}
\alias{maSub<-,marrayNorm-method}
\alias{maSub,marrayNorm-method}
%\alias{maTargets<-,marrayNorm,marrayInfo-method}
\alias{maTargets,marrayNorm-method}
\alias{maW<-,marrayNorm,matrix-method}
\alias{maW,marrayNorm-method}
\alias{print,marrayNorm-method}
\alias{show,marrayNorm-method}
\alias{summary,marrayNorm-method}
\alias{maNormCall}

\title{Class "marrayNorm", classes and methods for post-normalization cDNA microarray intensity data}


\description{
This class represents post-normalization intensity data for a batch of cDNA microarrays. A \emph{batch of arrays} consists of a collection of arrays with the same layout  (\code{"\link{marrayLayout}"}). The class contains slots for the average log-intensities A, the normalized log-ratios M, the location and scale normalization values, the layout of the arrays, and descriptions of the target samples hybridized to the arrays and probe sequences spotted onto the arrays.}


\section{Objects from the Class}{
Objects can be created by calls of the form 
\code{  new('marrayNorm',}\cr
\code{    maA        = ...., # Object of class matrix}\cr
\code{    maM        = ...., # Object of class matrix}\cr
\code{    maMloc     = ...., # Object of class matrix}\cr
\code{    maMscale   = ...., # Object of class matrix}\cr
\code{    maW        = ...., # Object of class matrix}\cr
\code{    maLayout   = ...., # Object of class marrayLayout}\cr
\code{    maGnames   = ...., # Object of class marrayInfo}\cr
\code{    maTargets  = ...., # Object of class marrayInfo}\cr
\code{    maNotes    = ...., # Object of class character}\cr
\code{    maNormCall = ...., # Object of class call}\cr
\code{  )}
}


\section{Slots}{
  \describe{
    \item{\code{maA}:}{Object of class \code{"matrix"}, average log-intensities (base 2) A, rows correspond to spotted probe sequences, columns to arrays in the batch.}
    \item{\code{maM}:}{Object of class \code{"matrix"}, intensity log-ratios (base 2) M, rows correspond to spotted probe sequences, columns to arrays in the batch.}
    \item{\code{maMloc}:}{Object of class \code{"matrix"}, location normalization values, rows correspond to spotted probe sequences, columns to arrays in the batch.}
    \item{\code{maMscale}:}{Object of class \code{"matrix"},  scale normalization values, rows correspond to spotted probe sequences, columns to arrays in the batch.}
    \item{\code{maW}:}{Object of class \code{"matrix"}, spot quality weights, rows correspond to spotted probe sequences, columns to arrays in the batch.}
    \item{\code{maLayout}:}{Object of class \code{"\link{marrayLayout}"}, layout parameters for cDNA microarrays.}
    \item{\code{maGnames}:}{Object of class \code{"\link{marrayInfo}"}, description of spotted probe sequences.}
    \item{\code{maTargets}:}{Object of class \code{"\link{marrayInfo}"}, description of target samples hybridized to the arrays.}
    \item{\code{maNotes}:}{Object of class \code{"character"}, any notes concerning the microarray experiments, e.g. hybridization or scanning conditions.}
    \item{\code{maNormCall}:}{Object of class \code{"call"}, function call for normalizing the batch of arrays.}
  }
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "marrayNorm")}:  subsetting operator for spots on the array and arrays in the batch, ensures that all slots are subset properly.}
    \item{coerce}{\code{signature(from = "marrayRaw", to = "marrayNorm")}: coerce an object of class \code{"\link{marrayRaw}"} into an object of class \code{\link{marrayNorm}}. }
%    \item{coerce}{\code{signature(from = "marrayNorm", to = "exprSet")}: coerce an object of class \code{"\link{marrayNorm}"} into an object of class \code{\link{exprSet-class}}. }
% \item{coerce}{\code{signature(from = "marrayNorm", to = "marraySpots")}: coerce an object of class \code{"\link{marrayNorm}"} into an object of class \code{"\link{marraySpots}"}. }
    \item{maA}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maA<-}{\code{signature(object = "marrayNorm", value = "matrix")}: slot assignment method.} 
    \item{maControls<-}{\code{signature(object = "marrayNorm")}: slot assignment method.} 
    \item{maControls}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maGnames}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maGnames<-}{\code{signature(object = "marrayNorm", value = "marrayInfo")}: slot assignment method.} 
    \item{maGridCol}{\code{signature(object = "marrayNorm")}: method which computes a vector of grid column coordinates for each spot.}
    \item{maGridRow}{\code{signature(object = "marrayNorm")}: method which computes a vector of grid row coordinates for each spot.}
    \item{maLayout}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maLayout<-}{\code{signature(object = "marrayNorm", value = "marrayLayout")}: slot assignment method.} 
    \item{maM}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maM<-}{\code{signature(object = "marrayNorm", value = "matrix")}: slot assignment method.} 
    \item{maMloc}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maMloc<-}{\code{signature(object = "marrayNorm", value = "matrix")}: slot assignment method.} 
    \item{maMscale}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maMscale<-}{\code{signature(object = "marrayNorm", value = "matrix")}: slot assignment method.} 
    \item{maNgc}{\code{signature(object = "marrayNorm")}: slot accessor method.}
 \item{maNgc<-}{\code{signature(object = "marrayNorm", value = "numeric")}: slot assignment method.} 
    \item{maNgr}{\code{signature(object = "marrayNorm")}: slot accessor method.}
\item{maNgr<-}{\code{signature(object = "marrayNorm", value = "numeric")}: slot assignment method.} 
    \item{maNormCall}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maNotes}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maNotes<-}{\code{signature(object = "marrayNorm", value = "character")}: slot assignment method.} 
    \item{maNsamples}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maNsc}{\code{signature(object = "marrayNorm")}: slot accessor method.}
\item{maNsc<-}{\code{signature(object = "marrayNorm", value = "numeric")}: slot assignment method.} 
    \item{maNspots}{\code{signature(object = "marrayNorm")}: slot accessor method.}
\item{maNspots<-}{\code{signature(object = "marrayNorm", value = "numeric")}: slot assignment method.} 
    \item{maNsr}{\code{signature(object = "marrayNorm")}: slot accessor method.}
\item{maNsr<-}{\code{signature(object = "marrayNorm", value = "numeric")}: slot assignment method.} 
    \item{maPlate}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maPlate<-}{\code{signature(object = "marrayNorm")}: slot assignment method.} 
    \item{maPrintTip}{\code{signature(object = "marrayNorm")}: method which computes a vector of print-tip-group indices for each spot.}
    \item{maSpotCol}{\code{signature(object = "marrayNorm")}: method which computes a vector of spot column coordinates for each spot.}
    \item{maSpotRow}{\code{signature(object = "marrayNorm")}: method which computes a vector of spot row coordinates for each spot.}
    \item{maSub}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maSub<-}{\code{signature(object = "marrayNorm")}: slot assignment method.} 
    \item{maTargets}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maTargets<-}{\code{signature(object = "marrayNorm", value = "marrayInfo")}: slot assignment method.}  
    \item{maW}{\code{signature(object = "marrayNorm")}: slot accessor method.}
    \item{maW<-}{\code{signature(object = "marrayNorm", value = "matrix")}: slot assignment method.} 
    \item{print}{\code{signature(x = "marrayNorm")}: print method for \code{"\link{marrayNorm}"} class.}
  }
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{marrayLayout}}, \code{\link{marrayRaw}},
  \code{\link{marrayInfo}}}

\examples{
# Examples use swirl dataset, for description type ? swirl

data(swirl)

# Median normalization
mnorm<-maNorm(swirl[,2:3],norm="m")

# Object of class marrayNorm for the second and third swirl arrays
mnorm

# Function call
maNormCall(mnorm)

#  Object of class marrayInfo -- Probe sequences
maGnames(mnorm)

#  Object of class marrayInfo -- Target samples
 maTargets(mnorm)

# Density plot of log-ratios M for third array
plot(density(maM(mnorm[,2])), lwd=2, col=2, main="Density plots of log-ratios M")
lines(density(maM(swirl[,3])), lwd=2)
abline(v=0)
legend(2,1,c("Pre-normalization","Post-normalization"))
}

\keyword{classes}
