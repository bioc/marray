\name{coerce-methods}
\docType{methods}


\alias{coerce-methods}
\alias{coerce}
\alias{as}

\title{Coerce an object to belong to a given microarray class}

\description{
Coercing methods were defined to convert microarray objects of one class
into objects of another class, e.g., instances of the
\code{"\link{marrayRaw}"} class into instances of the
\code{"\link{marrayNorm}"} class. 
}

\section{Methods}{\describe{
\item{from = marrayRaw, to = marrayNorm}{convert an object of class
  \code{"\link{marrayRaw}"} into an object of class
  \code{"\link{marrayNorm}"}.} 
}}

\examples{
##  Use Package convert to convert object to other data types such as
exprSet, MAList
}

\keyword{methods}


