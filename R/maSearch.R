###################################################################
## Examples
## library(marrayInput)
## data(swirl)
## findID("fb24a09", swirl, ID="ID")
## findID("geno1", swirl)
###################################################################

findID <-
  function(text,
           Gnames=gnames,
           ID = "Name")
{
  switch(data.class(Gnames),
         exprSet = G <- phenoData(Gnames),
         marrayRaw = G <- maGeneTable(Gnames),
         marrayNorm = G <-maGeneTable(Gnames),
         marrayInfo = G <- maInfo(Gnames),
         G <- Gnames
         )
  ind <- grep(ID, colnames(G))
  y <- as.vector(G[,ind])
  x <- grep(text, y)
  return(x)
}

