\name{maBayesian}
\alias{maBayesian}

\title{Calculates an Odds Ratio for Each Gene in a Multi-slide
  Microarray Experiment. } 

\description{This function calls the function "lm.series" and "ebayes"
  in the package "limma" to calculate the log-odds of differential
     expression by empirical Bayes shrinkage of the standard errors
     towards a common value.  }

\usage{
maBayesian(mraw, ...)
}

\arguments{
  \item{mraw}{An object of class \code{matrix}, \code{dataframe},
    \code{marrayRaw}, \code{marrayNorm} or \code{exprSet} that contains
    gene expression data.} 
  \item{...}{arguments for the fuction \code{lm.series} and \code{ebayes}}
}


\author{Yee Hwa Yang \email{jean@biostat.ucsf.edu} \cr}


\examples{
## data(swirl)
## M <- sweep(maM(swirl), 2, c(1,-1,1,-1), FUN="*")
## est <- maBayesian(M)
}

\keyword{models}








