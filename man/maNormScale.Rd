\name{maNormScale}
\alias{maNormScale}
\title{Simple scale normalization function}
\description{
This function is a simple wrapper function around the main normalization
function \code{\link{maNormMain}}. It allows the user to choose from a
set of two basic scale normalization procedures. The function operates
on an object of class \code{"\link{marrayRaw}"} (or possibly
\code{"\link{marrayNorm}"}, if normalization is performed in several
steps) and returns an object of class \code{"\link{marrayNorm}"}. This
function can be used to conormalize a batch of arrays 
(\code{norm="globalMAD"} option).
}
\usage{
maNormScale(mbatch, norm=c("globalMAD", "printTipMAD"), subset=TRUE, geo=TRUE,  Mscale=TRUE, echo=FALSE)
}

\arguments{
  \item{mbatch}{An object of class \code{"\link{marrayRaw}"}, containing 
intensity data for the batch of arrays to be normalized. An object of class 
\code{\link{marrayNorm}} may also be passed if normalization is performed 
in several steps.}
\item{norm}{A character string specifying the normalization
  procedures: 
\describe{
\item{globalMAD}{for global scale
  normalization using the median absolute deviation (MAD), this allows between 
slide scale normalization}
\item{printTipMAD}{for within-print-tip-group scale normalization
  using the median absolute deviation (MAD).}
This argument can be
  specified using the first letter of each method.}}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of points used to
  compute the  normalization values.}
 \item{geo}{If \code{TRUE}, the MAD of each group is divided by the
  geometric mean of the MADs across groups (cf. Yang et al. (2002)). This allows observations to
  retain their original units.}
  \item{Mscale}{If \code{TRUE}, the scale normalization values are stored in the slot \code{maMscale} of the object of class \code{"\link{marrayNorm}"} returned by the function, if \code{FALSE}, these values are not retained.}
\item{echo}{If \code{TRUE}, the index of the array currently being
  normalized is printed.}
}
\details{
See \code{\link{maNormMain}} for details and more general procedures.
}
\value{
\item{mnorm}{An object of class \code{"\link{marrayNorm}"}, containing the normalized intensity data.}
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, and T. P. Speed (2001). Normalization for cDNA microarray data. In M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), \emph{Microarrays: Optical Technologies and Informatics}, Vol. 4266 of \emph{Proceedings of SPIE}.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, D. M. Lin, V. Peng, J. Ngai, and T. P. Speed (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research}, Vol. 30, No. 4.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}

\seealso{\code{\link{maNormMain}}, \code{\link{maNorm}}.}


\examples{
# Examples use swirl dataset, for description type ? swirl
data(swirl)

# Global median normalization followed by global MAD normalization for 
# only arrays 2 and 3 in the batch swirl

mnorm1<-maNorm(swirl[,2:3], norm="m")
mnorm2<-maNormScale(mnorm1, norm="g")

}

\keyword{smooth}% at least one, from doc/KEYWORDS
\keyword{methods}% __ONLY ONE__ keyword per line
