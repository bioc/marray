\name{maDiagnPlots1}
\alias{maDiagnPlots1}
\alias{maRawPlots}
\alias{maNormPlots}

%- Also NEED an `\alias' for EACH other topic documented here.
\title{Basic diagnostic plots}
\description{
  These three functions produce sets of basic diagnostic plots for pre-- and post--normalization cDNA microarray data.
}
\usage{
maDiagnPlots1(mraw, title=NULL, save=TRUE,
	fname=paste(as.character(maLabels(maTargets(mraw))[1]), ".ps", sep = ""), 
	dev=c("postscript", "jpeg"))

maRawPlots(mraw, title=NULL, save = TRUE, 
	fname = paste(as.character(maLabels(maTargets(mraw))[1]), ".ps", sep = ""), 
	dev = c("postscript", "jpeg"))


maNormPlots(mnorm, title=NULL, save = TRUE, 
	fname = paste(as.character(maLabels(maTargets(mnorm))[1]),".ps", sep = ""), 
	dev = c("postscript", "jpeg"))

}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{mraw}{Microarray object of class \code{"\link{marrayRaw}"}, for \code{maDiagnPlots1} \code{maRawPlots}. By default, plots are done for the first array in the batch.}
 \item{mnorm}{Microarray object of class \code{"\link{marrayNorm}"}, for \code{maNormPlots}. By default, plots are done for the first array in the batch.}
 \item{title}{A "character" string for the plot title. If \code{NULL} a default title is shown.}
  \item{save}{If \code{TRUE}, the figures will be saved to the file named in \code{fname}.}
  \item{fname}{A "character" string naming the output file.}
  \item{dev}{A "character" string naming the graphics device, "postscript" or "jpeg".}
}
\details{
\code{maDiagnPlots1} produces eight plots of pre-- and post--normalization cDNA microarray data: color images of Cy3 and Cy5 background intensities, and of pre-- and post--normalization log--ratios M; boxplots of pre-- and post--normalization log--ratios M by print--tip--group; MA--plots of pre-- and post--normalization log--ratios M by print--tip--group.\cr \cr

\code{maRawPlots} produces six plots of pre--normalization cDNA microarray data: color images of Cy3 and Cy5 background intensities and log--ratios M; boxplots of log--ratios M by print--tip--group and plate; MA--plot of log--ratios M by print--tip--group.\cr \cr

\code{maNormPlots} produces four plots of post--normalization cDNA microarray data: color images of location normalization values Mloc and normalized log--ratios M; boxplots of normalized log--ratios M by print--tip--group; MA--plot of normalized log--ratios M by print--tip--group.
}

\references{
S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for exploratory analysis and normalization of cDNA microarray data. In G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors, \emph{The Analysis of Gene Expression Data: Methods and Software}, Springer, New York.
}

\author{Sandrine Dudoit, \url{http://www.stat.berkeley.edu/~sandrine}, Yee Hwa (Jean) Yang.}

\note{This function maynot work properly when "save=FALSE".}

\seealso{\code{\link{maPlot}}, \code{\link{maBoxplot}}, \code{\link{maImage}}.}

\examples{
library(marrayNorm)

# Examples use swirl dataset, for description type ? swirl
data(swirl)

maPlate(swirl)<-maCompPlate(swirl,384)

# Pre-normalization
maRawPlots(swirl[,3])

# Post-normalization
mnorm<-maNorm(swirl[,3], norm="twoD")
maNormPlots(mnorm, title="Swirl array 93 -- 2D normalization")

}
\keyword{hplot}
