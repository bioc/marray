\name{maLoess}
\alias{maLoess}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Stratified univariate robust local regression}
\description{
 This function performs robust local regression of a variable \code{y} on predictor variable \code{x}, separately within values of a third variable \code{z}. It is used by \code{\link{maNormLoess}} for intensity dependent location normalization.
}
\usage{
maLoess(x, y, z, w=NULL, subset=TRUE, span=0.4, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{A numeric vector of predictor variables.}
  \item{y}{A numeric vector of responses.}
  \item{z}{Variables used to stratify the data.}
  \item{w}{An optional numeric vector of weights.}
  \item{subset}{A "logical" or "numeric" vector indicating the subset of points used to compute the fits. }
  \item{span}{The argument \code{span} which controls the degree of
    smoothing in the  \code{\link{loess}} function.}
    \item{...}{Misc arguments}
}

\details{
\code{y} is regressed on \code{x}, separately within values of \code{z} using the \code{\link{loess}} function.  
}

\value{
A numeric vector of fitted values.
}

\references{S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, and T. P. Speed (2001). Normalization for cDNA microarray data. In M. L. Bittner, Y. Chen, A. N. Dorsel, and E. R. Dougherty (eds), \emph{Microarrays: Optical Technologies and Informatics}, Vol. 4266 of \emph{Proceedings of SPIE}.\cr \cr

Y. H. Yang, S. Dudoit, P. Luu, D. M. Lin, V. Peng, J. Ngai, and T. P. Speed (2002). Normalization for cDNA microarray data: a robust composite method addressing single and multiple slide systematic variation. \emph{Nucleic Acids Research}, Vol. 30, No. 4.}}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}.}


\seealso{\code{\link{maNormMain}}, \code{\link{maNormLoess}}, \code{\link{loess}}.}


\examples{
# See examples for maNormMain.

}
\keyword{smooth}% at least one, from doc/KEYWORDS

