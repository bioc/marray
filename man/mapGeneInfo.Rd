\name{mapGeneInfo}
\alias{mapGeneInfo}
\alias{widget.mapGeneInfo}
\alias{URLstring}
\alias{SFGL}
\alias{UCBFGL}

\title{Creating URL strings for external database links }

\description{
  These functions are used with \code{\link{htmlPage}}.
  The function \code{mapGeneInfo}, takes all the arguments and generate
  a character matrix of two columns.  The first columns representing the
  name of the argument and the second columns represents the value of an
  argument. 
  The function \code{widget.mapGeneInfo} allows the user to enter this
  information  interactively. 
}

\usage{
mapGeneInfo(widget = FALSE, Gnames, Name = "pubmed", ID =
                 "genbank", ACC = "SMDacc", ...)
widget.mapGeneInfo(Gnames) 
}

%- maybe also `usage' for other objects documented here.
\arguments{
  \item{widget}{A logical value specifying if widgets should be used.}
  \item{Name}{The external database for spot description, E.g. "pubmed".}
  \item{ID}{The external database for spot ID, E.g. "operon", "Riken", "locuslink".}
  \item{ACC}{The external database for gene accession number, E.g. "genebank".}
  \item{Gnames}{An object of class \code{matrix}, \code{data.frame} or \code{marrayInfo}
    which contains  description of spotted probe sequences.}
  \item{\dots}{Other column names}
}
\details{
  The function \code{mapGeneInfo} generates a character matrix with the
  first column representing the column headings of "Gnames" and the
  second column representing the corresponding names in the list
  \code{URLstring}.  For example, if a particular column in "Gnames"
  with column names "ID" contains genebank accession number, then the
  function \code{mapGeneInfo} generates a row containing "ID" in the
  first column and "genbank" in the second.  Examples are \code{SFGL}
  and \code{UCBFGL}.\cr
  
  \code{URLstring} is a list contains the URL to various external
  database, E.g. operon, Riken, genbank. \cr
  The current choices are: 
  "pubmed",  "locuslink",  "riken",   "SMDclid",    "SMDacc",    "operonh2",  "operonh1" ,
  "operonm2",  "operonm1" and  "genbank" .  
  "SMDclid" and "SMDacc" are links to Stanford Microarray Databases.
 }

 \author{Jean Yee Hwa Yang}

 \examples{
 mapGeneInfo(ID="genebank", ll="locuslink")
 mapGeneInfo(ID="locuslink", Sample.ID="riken")
}

\keyword{manip}% at least one, from doc/KEYWORDS


