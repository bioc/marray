\name{htmlPage}
\alias{htmlPage}
\alias{table2html}

\title{Display gene list as a HTML page}
\description{
Given a set of index to a data.frame containing gene names information.
We create a web page with one element per genes that contains URLs links
to various external database links.
E.g Operon oligodatabase , Riken, GenBank and PubMed web sites. 
}

\usage{
htmlPage(genelist, filename = "GeneList.html", geneNames =
                 Gnames, mapURL = SFGL, othernames, title, table.head,
                 table.center = TRUE, disp = c("browser", "file")[1])

table2html(restable, filename = "GeneList.html", mapURL = SFGL,
                 title, table.head, table.center = TRUE, disp =
                 c("browser", "file")[1])
}

\arguments{
  \item{restable}{A data.frame that contains only the information you
    wish to display in the html file.  The rows corresponds to a
    different DNA spots.}
  \item{genelist}{A numeric vector of index to a data.frame}
  \item{filename}{The name of the file to store the HTML in. }
  \item{geneNames}{A data.frame containing the information related the
    each DNA spots. }
  \item{mapURL}{A matrix of characters containing the URL for various
    external database.  E.g \code{\link{SFGL}}.}
  \item{othernames}{A data.frame containing other information.}
  \item{title}{Title of the HTML page}
  \item{table.head}{ A character vector of column labels for the table}
  \item{table.center}{A logical indicating whether the table should be
    centered}
  \item{disp}{Either "File" or "Browser" (default is Browser). File will
    save the information in html file, while Browser will create an html
    files and display information in the user's browser.}
}
\details{
  This function is an extension to \code{ll.htmlpage}
}

\value{
  No value is return, the function produce a html file "filename" and
  output the results in a browser.
}
\author{Yee Hwa Yang }

\seealso{\code{ll.htmlpage}, \code{\link{URLstring}}, \code{\link{widget.mapGeneInfo}}}

\examples{
##library(annotate)
data(swirl)
Gnames <- maGeneTable(swirl)
swirlmap <- mapGeneInfo(Name = "none", ID="genbank")
## htmlPage(100:110, geneNames = Gnames, mapURL = swirlmap, title="Swirl")

moreinfo <- round(maM(swirl), 2)
swirlmap <- mapGeneInfo(Name = "pubmed", ID="genbank")
##htmlPage(100:110, geneNames = Gnames, mapURL = swirlmap, othernames=moreinfo, title="Swirl", disp="file")
}
\keyword{file}% at least one, from doc/KEYWORDS
