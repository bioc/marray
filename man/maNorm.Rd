\name{maNorm}
\alias{maNorm}
\title{Simple location and scale normalization function}
\description{
This function is a simple wrapper function around the main normalization function \code{\link{maNormMain}}. It allows the user to choose from a set of six basic location and scale normalization procedures. The function operates on an object of class \code{"\link{marrayRaw}"} (or possibly \code{"\link{marrayNorm}"}, if normalization is performed in several steps) and returns an object of class \code{"\link{marrayNorm}"}.
}
\usage{
maNorm(mbatch, norm=c("printTipLoess", "none", "median", "loess",
"twoD", "scalePrintTipMAD"), subset=TRUE, span=0.4, Mloc=TRUE,
Mscale=TRUE, echo=FALSE, ...)
}
\arguments{
  \item{mbatch}{Object of class \code{\link{marrayRaw}}, containing intensity 
data for the batch of arrays to be normalized. 
An object of class \code{"\link{marrayNorm}"} may also be passed if 
normalization is performed in several steps.}
\item{norm}{Character string specifying the normalization procedures: 
\describe{
\item{none}{no normalization}
\item{median}{for global median location normalization}
\item{loess}{for global intensity or A-dependent location normalization using 
the \code{\link{loess}} function}
\item{twoD}{for 2D spatial location normalization using the 
\code{\link{loess}} function}
\item{printTipLoess}{for within-print-tip-group intensity dependent location 
normalization using the \code{\link{loess}} function}
\item{scalePrintTipMAD}{for within-print-tip-group intensity dependent 
location normalization followed by within-print-tip-group scale normalization 
using the median absolute deviation (MAD). \cr
}
This argument can be specified using the first letter of each method.}}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of points used to
  compute the  normalization values.}
 \item{span}{The argument \code{span} which controls the degree of smoothing in the \code{\link{loess}} function.}
\item{Mloc}{If \code{TRUE}, the location normalization values are stored in the slot \code{maMloc} of the object of class \code{"\link{marrayNorm}"} returned by the function, if \code{FALSE}, these values are not retained.}
  \item{Mscale}{If \code{TRUE}, the scale normalization values are stored in the slot \code{maMscale} of the object of class \code{"\link{marrayNorm}"} returned by the function, if \code{FALSE}, these values are not retained.}
\item{echo}{If \code{TRUE}, the index of the array currently being
  normalized is printed.}
\item{...}{Misc arguments}
}
\details{
See \code{\link{maNormMain}} for details and also more general procedures.
}
\value{
 \item{mnorm}{An object of class \code{"\link{marrayNorm}"}, containing the normalized intensity data.}
}


\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, and T. P. Speed (2001). Normalization for cDNA microarray data. In M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), \emph{Microarrays: Optical Technologies and Informatics}, Vol. 4266 of \emph{Proceedings of SPIE}.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, D. M. Lin, V. Peng, J. Ngai, and T. P. Speed (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research}, Vol. 30, No. 4.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{maNormMain}}, \code{\link{maNormScale}}.}

\examples{
# Examples use swirl dataset, for description type ? swirl
data(swirl)

# Global median normalization for swirl arrays 2 and 3
mnorm<-maNorm(swirl[,2:3], norm="median", echo=TRUE)

# Within-print-tip-group loess location normalization for swirl array 1
mnorm<-maNorm(swirl[,1], norm="p", span=0.45)

}




\keyword{smooth}% at least one, from doc/KEYWORDS
\keyword{methods}% __ONLY ONE__ keyword per line
