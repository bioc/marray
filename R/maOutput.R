###########################################################################
#
# Output functions
#
# Date : September 19, 2002
# Modify : October, 14, 2002
#
# Runs on R 1.6 and above
#
# This file contains wrapper functions for analysis
#
# source("~/Projects/maTools/R/maInOut.R")
#

write.marray <- function(mraw, file="maRawResults.xls", val="maM", ...)
{
  tmp <- cbind(maGeneTable(mraw), eval(call(val,mraw)))
  write.table(tmp, file, row.names=FALSE, col.names=TRUE, sep="\t", ...)
  return(NULL)
}
  
  
write.xls <- function(res, file="test.xls", ...)
  {
    write.table(res, file, row.names=FALSE, col.names=TRUE, sep="\t", ...)
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
