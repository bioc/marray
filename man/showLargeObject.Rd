\name{ShowLargeObject-class}
\docType{class}

\alias{ShowLargeObject-class}
\alias{show,ShowLargeObject-method}
%\alias{print-methods}

\title{Show Large Data Object - class}
      
\description{
A virtual class including the data classes \code{marrayRaw},
\code{marrayNorm}, \code{marrayInfo}, \code{marrayLayout},
\code{PrinterInfo}, \code{RGData} and \code{MAData}, all of which typically contain large
quantities of numerical data in vector, matrices and data.frames. 
}

\section{Methods}{
A \code{show} method is defined for objects of class
\code{ShowLargeObject} which uses \code{printHead} to print only the
leading elements or rows of components or slots which contain large
quantities of data. 
}

\author{modifid from Gordon Smyth's function}

\keyword{classes}
\keyword{data}
