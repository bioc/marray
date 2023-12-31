###########################################################################
#
# Output functions
#
# Date : September 19, 2002
# Modify : October, 14, 2002, April 10 2005 (Agnes)
#
# Runs on R 1.6 and above
#
# This file contains wrapper functions for analysis
#
# source("~/Projects/maTools/R/maInOut.R")
#

write.marray <- function(mraw, file="maRawResults.xls", val="maM", ...)
{
  # Was tmp <- cbind(maGeneTable(mraw), eval(call(val,mraw)), row.names=NULL)
  res <- c()
  for(i in 1:length(val))
    res <- cbind(res, eval(call(val[i],mraw)))
  tmp <- cbind(maGeneTable(mraw), res, row.names=NULL)
  write.table(tmp, file, row.names=FALSE, col.names=TRUE, sep="\t", ...)
  return(NULL)
}
  

## Apr, 09 based on suggestion by
## Matthew Fero <mfero@fhcrc.org>
write.xls <- function(res, file="test.xls", ...)
  {
    ## col.names = NA and row.names = TRUE
    opt = list(...)
    defs <- list(x = res, file = file, row.names=FALSE, col.names=TRUE, sep="\t")  
    write.args <- maDotsMatch(opt, maDotsMatch(defs, formals(args("write.table"))))    
##    print(write.args)
    do.call(write.table, write.args)
    ##    write.table(res, file, row.names=FALSE, col.names=TRUE, sep="\t", ...) 
  }


write.list <- function (x,
                        filename = "data",
                        append = FALSE,
                        closefile = TRUE,
                        outfile
                        ) 
{
  if(!append)
    {
      outfile <- file(filename, "w")
      cat(file=outfile, append= append)
    }
  for(i in 1:length(x))
    {
      cat(paste(names(x)[i], ":"), file = outfile, append = TRUE)
      cat("\n", file = outfile, append=TRUE)
      if(!is.null(x[[i]]))
        {
          switch(data.class(x[[i]]),
                 matrix = write.table(x[[i]], sep="\t", file = outfile, append = TRUE),
                 table = if(!is.null(names(x[[i]])))
                 {
                   write.table(rbind(names(x[[i]]), x[[i]]),file = outfile, 
                               append = TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
                 }
                 else
                 {
                   write(x[[i]], file=outfile, append=TRUE)
                 },
                 list = write.list(x[[i]], outfile=outfile, append=TRUE, closefile=FALSE),
                 if(!is.null(names(x[[i]])))
                 {
                   write.table(rbind(names(x[[i]]), x[[i]]),file = outfile, 
                               append = TRUE, row.names=FALSE, col.names=FALSE, sep="\t")
                 }
                 else
                 {
                   write(x[[i]], ncolumns = length(x[[i]]), file = outfile, append = TRUE)
                 }
                 )
        }
      cat("\n", file = outfile, append=TRUE)
    }
  if(closefile) close(outfile)
}
