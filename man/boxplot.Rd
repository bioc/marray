\name{boxplot}
\docType{methods}
\alias{boxplot,marrayRaw-method}
\alias{boxplot,marrayNorm-method}



\title{Boxplots for cDNA microarray spot statistics}

\description{
  The function \code{boxplot} produces boxplots of microarray spot
  statistics for the classes  \code{"\link{marrayRaw}"},
  \code{"\link{marrayNorm}"}.
  We encourage users to use \code{boxplot} rather than \code{maBoxplot}.
  The name of the arguments have changed slightly.
}

\usage{
\S4method{boxplot}{marrayRaw}(x, xvar="maPrintTip", yvar="maM", ...)
\S4method{boxplot}{marrayNorm}(x, xvar="maPrintTip", yvar="maM", ...)
}
%- maybe also `usage' for other objects documented here.

\arguments{
  \item{x}{Microarray object of class \code{"\link{marrayRaw}"}, \code{"\link{marrayNorm}"}}
  \item{xvar}{Name of accessor method for the spot statistic used to
    stratify the data, typically a slot name for the microarray layout
    object (see \code{"\link{marrayLayout}"}) such as \code{maPlate} or
    a method such as \code{maPrintTip}. If \code{x} is NULL, the data
    are not stratified.} 
  \item{yvar}{Name of accessor method for the spot statistic of
    interest, typically a slot name for the microarray object \code{m},
    such as \code{maM}.} 
  \item{\dots}{Optional graphical parameters, see \code{\link{par}}.}
}

\details{
If there are more than one array in the batch, the function produces a
boxplot for each array in the batch. Such plots are useful when
assessing the need for between array normalization, for example, to deal
with scale differences among different arrays. Default graphical
parameters are chosen for convenience using the function
\code{\link{maDefaultPar}} (e.g. color palette,  axis labels,  plot
title) but the user has the option to overwrite these parameters at any
point. 
}

\references{
S. Dudoit and Y. H. Yang. (2002). Bioconductor R packages for
exploratory analysis and normalization of cDNA microarray data. In
G. Parmigiani, E. S. Garrett, R. A. Irizarry and S. L. Zeger, editors,
\emph{The Analysis of Gene Expression Data: Methods and Software},
Springer, New York. 
}

\author{Jean Yang and Sandrine Dudoit}

\seealso{\code{\link{maBoxplot}}, \code{\link{maDefaultPar}}.}

\examples{

# To see the demo type demo(marrayPlots)

# Examples use swirl dataset, for description type ? swirl
data(swirl)

# Boxplots of pre-normalization log-ratios M for each of the 16 
# print-tip-groups for the Swirl 93 array. 
# - Default arguments
boxplot(swirl[,3])  

# All spots 
boxplot(swirl[,3], xvar=NULL, col="green")  

# Boxplots of pre-normalization red foreground intensities for each grid row
# for the Swirl 81 array. 
boxplot(swirl[,1], xvar="maGridRow", yvar = "maRf", main = "Swirl array 81: pre-normalization red foreground intensity") 

# Boxplots of pre-normalization log-ratios for each array in swirl
boxplot(swirl, main="Swirl arrays: pre-normalization log-ratios")
}
\keyword{hplot}% at least one, from doc/KEYWORDS

