\name{maNormMain}
\alias{maNormMain}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Main function for location and scale normalization of cDNA microarray data}
\description{
This is the main function for location and scale normalization of cDNA microarray data. Normalization is performed for a batch of arrays using location and scale normalization procedures specified by the lists of functions \code{f.loc} and \code{f.scale}. Typically, only one function is given in each list, otherwise composite normalization is performed using the weights computed by the functions \code{a.loc} and \code{a.scale}. The function operates on an object of class \code{"\link{marrayRaw}"} (or possibly \code{"\link{marrayNorm}"}, if normalization is performed in several steps) and returns an object of class \code{"\link{marrayNorm}"}. Simple wrapper functions are provided by \code{\link{maNorm}} and \code{\link{maNormScale}}.
}

\usage{
maNormMain(mbatch, f.loc=list(maNormLoess()), f.scale=NULL,
a.loc=maCompNormEq(), a.scale=maCompNormEq(), Mloc=TRUE, Mscale=TRUE, echo=FALSE)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{mbatch}{An object of class \code{"\link{marrayRaw}"}, containing intensity data for the batch of arrays to be normalized. An object of class \code{"\link{marrayNorm}"} may also be passed if normalization is performed in several steps.}
  \item{f.loc}{A list of location normalization functions, e.g., \code{\link{maNormLoess}}, \code{\link{maNormMed}}, or \code{\link{maNorm2D}}.}
  \item{f.scale}{A list of scale normalization functions, .e.g, \code{\link{maNormMAD}}.}
  \item{a.loc}{For composite normalization, a function for computing the weights used in combining several location normalization functions, e.g., \code{\link{maCompNormA}}.}
  \item{a.scale}{For composite normalization, a function for computing
  the weights used in combining several scale normalization functions.}
  \item{Mloc}{If \code{TRUE}, the location normalization values are stored in the slot \code{maMloc} of the object of class \code{"\link{marrayNorm}"} returned by the function, if \code{FALSE}, these values are not retained.}
  \item{Mscale}{If \code{TRUE}, the scale normalization values are stored in the slot \code{maMscale} of the object of class \code{"\link{marrayNorm}"} returned by the function, if \code{FALSE}, these values are not retained.}
  \item{echo}{If \code{TRUE}, the index of the array currently being
  normalized is printed.}
}

\details{
When both location and scale normalization functions (\code{f.loc} and \code{f.scale}) are passed, location normalization is performed before scale normalization. That is, scale values are computed for the location normalized log-rations. The same results could be obtained by two applications of the function \code{maNormMain}, first with only the location normalization function and \code{f.scale=NULL}, and second with only the scale normalization function and \code{f.loc=NULL}.
}

\value{
  \item{mnorm}{An object of class \code{"\link{marrayNorm}"}, containing the normalized intensity data.}
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, and T. P. Speed (2001). Normalization for cDNA microarray data. In M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), \emph{Microarrays: Optical Technologies and Informatics}, Vol. 4266 of \emph{Proceedings of SPIE}.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, D. M. Lin, V. Peng, J. Ngai, and T. P. Speed (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research}, Vol. 30, No. 4.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}


\seealso{\code{\link{maNorm}}, \code{\link{maNormScale}}, \code{\link{maNormLoess}}, \code{\link{maLoess}},
\code{\link{maNormMAD}}, \code{\link{maMAD}},
\code{\link{maNormMed}}, \code{\link{maMed}},
\code{\link{maNorm2D}}, \code{\link{ma2D}},
\code{\link{maCompNormA}}, \code{\link{maCompNormEq}}.}

\examples{
# Examples use swirl dataset, for description type ? swirl
data(swirl)

# Within-print-tip-group loess location normalization of batch swirl
# - Default normalization
swirl.norm<-maNormMain(swirl)

boxplot(swirl.norm)
boxplot(swirl.norm[,3])
plot(swirl.norm[,3])

# Global median normalization for arrays 81 and 82
swirl.norm <- maNormMain(swirl[,1:2], f.loc = list(maNormMed(x=NULL,y="maM")))

# Global loess normalization for array 81
swirl.norm <- maNormMain(swirl[,1], f.loc = list(maNormLoess(x="maA",y="maM",z=NULL)))

# Composite normalization as in Yang et al. (2002)
# No MSP controls are available here, so all spots are used for illustration
# purposes
swirl.norm <- maNormMain(swirl[,1], f.loc = list(maNormLoess(x="maA",y="maM",z=NULL),maNormLoess(x="maA",y="maM",z="maPrintTip")), a.loc=maCompNormA())
}
\keyword{smooth}% at least one, from doc/KEYWORDS
\keyword{methods}% __ONLY ONE__ keyword per line
