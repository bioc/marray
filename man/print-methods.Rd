\name{print-methods}
\docType{methods}

\alias{print-methods}
\alias{show-methods}
\alias{print}

\title{Printing methods for microarray objects}

\description{
Print methods were defined for the microarray classes, \code{"\link{marrayInfo}"}, \code{"\link{marrayLayout}"}, \code{"\link{marrayRaw}"}, \code{"\link{marrayNorm}"}, \code{"\link{marraySpots}"}, and \code{"\link{marrayTwo}"}. These methods produce summaries of the intensity and textual data stored in different classes of microarray objects.
}

\section{Methods}{\describe{
\item{x = ANY}{generic print method}
\item{x = marrayLayout}{for an object of class \code{"\link{marrayLayout}"}, the method prints main layout parameters such as the number of spots and the dimensions of the spot and grid matrices.}
\item{x = marrayInfo}{for an object of class \code{"\link{marrayInfo}"}, the method prints the first 10 rows of the \code{"maInfo"} and \code{"maLabels"} slots.}
\item{x = marrayRaw}{for an object of class \code{"\link{marrayRaw}"}, the method prints a short description of the microarray layout \code{"maLayout"} and the target samples hybridized to the arrays \code{"maTargets"}, and a summary of the distribution of the log-ratio statistics \code{"maM"}.}
\item{x = marrayNorm}{for an object of class \code{"\link{marrayNorm}"}, the method prints a short description of the microarray layout \code{"maLayout"} and the target samples hybridized to the arrays \code{"maTargets"}, and a summary of the distribution of the log-ratio statistics \code{"maM"}.}
\item{x = marrayTwo}{for an object of class \code{"\link{marrayTwo}"}, the method prints a short description of the microarray layout \code{"maLayout"}, and a summary of the distribution of the statistics \code{"maX"} and \code{"maY"}.}
\item{x = marraySpots}{for an object of class \code{"\link{marraySpots}"}, the method prints a short description of the microarray layout \code{"maLayout"} and the probe sequences spotted on the arrays \code{"maGnames"}.}
}}
\keyword{methods}

