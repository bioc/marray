\name{maPalette}
\alias{maPalette}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Microarray color palette}
\description{
This function returns a vector of color names corresponding to a range of colors specified in the arguments.  
}
\usage{
maPalette(low = "white", high = c("green", "red"), mid=NULL, k =50)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{low}{Color for the lower end of the color palette, specified using any of the three kinds of R colors, i.e., either a color name (an element of \code{colors}), a hexadecimal string of the form \code{"#rrggbb"}, or an integer \code{i} meaning \code{palette()[i]}.}
  \item{high}{Color for the upper end of the color palette, specified
  using any of the three kinds of R colors, i.e., either a color name
  (an element of \code{colors}), a hexadecimal string of the form
  \code{"#rrggbb"}, or an integer \code{i} meaning \code{palette()[i]}.}
 \item{mid}{Color for the middle portion of the color palette, specified using any of the three kinds of R colors, i.e., either a color name (an element of \code{colors}), a hexadecimal string of the form \code{"#rrggbb"}, or an integer \code{i} meaning \code{palette()[i]}.}
  \item{k}{Number of colors in the palette.}
}

\value{
  A "character" vector of color names. This can be used to create a user-defined color palette for subsequent graphics by \code{palette}, in a \code{col=} specification in graphics functions, or in \code{\link{par}}.
}

\references{
S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}, Yee Hwa (Jean) Yang.}

\seealso{\code{\link{image}}, \code{\link{maColorBar}}, \code{\link{maImage}}, \code{\link{maImage.func}}.}

\examples{
par(mfrow=c(1,4))
pal <- maPalette(low="red", high="green")
maColorBar(seq(-2,2, 0.2), col=pal, horizontal=FALSE, k=21)
pal <- maPalette(low="red", high="green", mid="yellow")
maColorBar(seq(-2,2, 0.2), col=pal, horizontal=FALSE, k=21)
pal <- maPalette()
maColorBar(seq(-2,2, 0.2), col=pal, horizontal=FALSE, k=21)
pal <- maPalette(low="purple", high="purple",mid="white")
maColorBar(seq(-2,2, 0.2), col=pal, horizontal=FALSE, k=21)
}
\keyword{color}% at least one, from doc/KEYWORDS

