\name{marrayInfo-class}
\docType{class}

\alias{marrayInfo-class}
\alias{marrayInfo}

\alias{maLabels}
\alias{maLabels<-}
\alias{maInfo}
\alias{maInfo<-}
\alias{maNotes}
\alias{maNotes<-}

\alias{[,marrayInfo-method}
\alias{maGnames<-,marrayRaw,marrayInfo-method}
\alias{maGnames<-,marrayNorm,marrayInfo-method}
\alias{maGnames<-,marraySpots,marrayInfo-method}
\alias{maInfo<-,marrayInfo,data.frame-method}
\alias{maInfo,marrayInfo-method}
\alias{maLabels<-,marrayInfo,character-method}
\alias{maLabels<-,marrayInfo,numeric-method}
\alias{maLabels,marrayInfo-method}
\alias{maNotes<-,marrayInfo,character-method}
\alias{maNotes,marrayInfo-method}
\alias{maTargets<-,marrayRaw,marrayInfo-method}
\alias{maTargets<-,marrayNorm,marrayInfo-method}
\alias{print,marrayInfo-method}
\alias{rbind,marrayInfo-method}
\alias{summary,marrayInfo-method}

\title{Class "marrayInfo", description of target samples or spotted probe sequences}

\description{This class is used to store information on target samples hybridized to a batch of arrays or probe sequences spotted onto these arrays. It is not specific to the microarray context.}


x\section{Objects from the Class}{
Objects can be created by calls of the form \code{  new('marrayInfo',}\cr
\code{    maLabels = ...., # Object of class character}\cr
\code{    maInfo   = ...., # Object of class data.frame}\cr
\code{    maNotes  = ...., # Object of class character}\cr
\code{  )}
}

\section{Slots}{
  \describe{
    \item{\code{maLabels}:}{Object of class \code{"character"}, vector of spot or array labels.}
    \item{\code{maInfo}:}{Object of class \code{"data.frame"}. If the object of class \code{"marrayInfo"} is used to describe probe sequences, rows of \code{maInfo} correspond to spots and columns to various gene identifiers and annotations. 
If the object of class \code{"marrayInfo"} is used to describe target samples hybridized to the arrays, rows of \code{maInfo} correspond to arrays and columns to various descriptions of the hybridizations, e.g., names of Cy3 and Cy5 samples, labels for the arrays etc.}
    \item{\code{maNotes}:}{Object of class \code{"character"}, any notes on the target samples or spotted probe sequences.}
  }
}

\section{Methods}{
  \describe{
    \item{[}{\code{signature(x = "marrayInfo")}: subsetting operator for spots on the array or arrays in the batch, ensures that all slots are subset properly.}
    \item{maGnames<-}{\code{signature(object = "marrayRaw", value = "marrayInfo")}: slot assignment method.}
    \item{maGnames<-}{\code{signature(object = "marrayNorm", value = "marrayInfo")}: slot assignment method.}
    \item{maGnames<-}{\code{signature(object = "marraySpots", value = "marrayInfo")}: slot assignment method.}
    \item{maInfo}{\code{signature(object = "marrayInfo")}: slot accessor method.}
    \item{maInfo<-}{\code{signature(object = "marrayInfo", value = "data.frame")}: slot assignment method.}
    \item{maLabels}{\code{signature(object = "marrayInfo")}: slot accessor method.}
    \item{maLabels<-}{\code{signature(object = "marrayInfo", value = "character")}: slot assignment method.}
    \item{maLabels<-}{\code{signature(object = "marrayInfo", value = "numeric")}: slot assignment method.}
    \item{maNotes}{\code{signature(object = "marrayInfo")}: slot accessor method.}
    \item{maNotes<-}{\code{signature(object = "marrayInfo", value = "character")}: slot assignment method.}
    \item{maTargets<-}{\code{signature(object = "marrayRaw", value = "marrayInfo")}: slot assignment method.}
    \item{maTargets<-}{\code{signature(object = "marrayNorm", value = "marrayInfo")}: slot assignment method.}
    \item{print}{\code{signature(x = "marrayInfo")}: print method for \code{"\link{marrayInfo}"} class.}
  }
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}


\author{Jean Yang and Sandrine Dudoit}

\seealso{\code{\link{marrayLayout}}, \code{\link{marrayRaw}}, \code{\link{marrayNorm}}.}

\examples{
## See marrayRaw
}
\keyword{classes}
