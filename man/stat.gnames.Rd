%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file was automatically generated by Rdoc on 2001-11-08 :05:09.
% Do NOT modify this file, instead modify the source (Rarray.R).
% Rdoc is copyright (C) 2001, Henrik Bengtsson, henrikb@braju.com.
% Download the Rdoc compiler at http://www.braju.com/R/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{stat.gnames}

\alias{stat.gnames}

\title{Sort Genes According to the Value of a Statistic}

\description{
Lists genes and corresponding statistics in decreasing order of the
statistics. This function applies to any type of statistic, including
log ratios, one and two-sample t-statistics, and F-statistics. Missing
values are ignored, as in \code{\link{sort}(..., na.last=NA)}. 
}

\usage{
stat.gnames(x, gnames, crit= 50)
}

\arguments{
 \item{x}{a numeric vector containing the statistics for each
 gene. Missing values (NAs) are allowed. }
 
\item{gnames}{a character vector containing the gene names.}

 \item{crit}{specifies the number of genes to be returned. If crit <
 1, the crit*100\% genes with the largest x values are listed. If crit
 >= 1, the crit genes with the largest x values are listed. }
}

\value{
List containing the following components 
  \item{gnames}{gene names sorted in decreasing order of the
 statistics in x.}
 \item{t}{statistics sorted in decreasing order.}
}

\author{
  Yee Hwa Yang, \email{yeehwa@stat.berkeley.edu} \cr
  Sandrine Dudoit, \email{sandrine@stat.berkeley.edu} }

\seealso{\code{\link{order}}, \code{\link{sort}}.}

\examples{
data(swirl)
aveM <- apply(maM(swirl), 1, mean.na)
Gnames <- maGeneTable(swirl)

stat.gnames(abs(aveM), Gnames, crit=10)
stat.gnames(aveM, Gnames, crit=0.01)

}

\keyword{misc}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
