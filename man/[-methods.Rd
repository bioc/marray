\name{[-methods}
\docType{methods}

\alias{[-methods}
\alias{[}

\title{Subsetting methods for microarray objects}


\description{
Subsetting methods were defined for the microarray classes,
\code{"\link{marrayInfo}"}, \code{"\link{marrayLayout}"},
\code{"\link{marrayRaw}"} and \code{"\link{marrayNorm}"}}. These
methods create instances of the given class, for a subset of spots
and/or arrays in a batch. 
}

\section{Methods}{\describe{
\item{x = ANY}{generic method.}
\item{x = marrayInfo}{\code{x[i, j]} extract object of class \code{"\link{marrayInfo}"} for spots or arrays with indices i and labels with indices j.}
\item{x = marrayLayout}{\code{x[i]} extract object of class \code{"\link{marrayLayout}"} for spots with indices i.}
\item{x = marrayRaw}{\code{x[i, j]} extract object of class \code{"\link{marrayRaw}"} for spots with indices i and arrays with indices j.}
\item{x = marrayNorm}{\code{x[i, j]} extract object of class \code{"\link{marrayNorm}"} for spots with indices i and arrays with indices j.}
}}

\keyword{methods}

