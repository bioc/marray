###########################################################################
# Date : September 19, 2002
# Modify : October, 18, 2002
#
# Runs on R 1.5.1 and above
#
# This file contains wrapper functions for Spot files
#
# source("~/Projects/maTools/R/spotWrap.R")
###########################################################################

###########################################################################
## This is a wrapper function specifially to generate diagnostic plots and
## quality file for every genepix files in the current working directory

spotTools <- function(fnames,
                      path=".",
                      galfile,
                      bg=TRUE,
                      plot=TRUE,
                      quality=TRUE,
                      fill=TURE,
                      raw=FALSE,
                      echo=TRUE,
                      ...)
  {
    opt <- list(...)
    normM <- normA <- NULL
    if(missing(fnames)) fnames <- dir(path, pattern="*\\.spot$")
    if(missing(galfile))
      {
        tmp <- dir(path, pattern="*\\.gal$")
        ifelse(length(tmp)==0, galfile <- fnames[1], galfile <- tmp)
      }

    if(echo) cat("Reading Gal file ...")
    args <- maDotsMatch(c(opt, galfile=galfile, path=path), formals(args("read.Galfile")))
    info <- do.call("read.Galfile", args)
    maControls(info$layout) <- maGenControls(info$gnames)
    if(echo) cat("done \n ")
    
    if(quality) Q.res <- NULL
    if(raw) rawdata <- new("marrayRaw")
    
    for(i in fnames)
      {
        defs <- list(fnames =i, path=path,  name.W="circularity",
                     layout = info$layout,
                     gnames=info$gnames,
                     fill=TRUE, quote="")
        args <- maDotsMatch(c(defs, opt), formals(args("read.Spot")))
        coredata <- do.call("read.Spot", args)
        tmp <- maW(coredata) 
        maW(coredata) <- apply(tmp >= 1, 2, as.numeric)
        if(raw)
          {
            if(length(maGf(rawdata)) == 0)
               rawdata <- coredata
            else
              {
                maGf(rawdata) <- cbind(maGf(rawdata), maGf(coredata))
                maRf(rawdata) <- cbind(maRf(rawdata), maRf(coredata))
                maGb(rawdata) <- cbind(maGb(rawdata), maGb(coredata))
                maRb(rawdata) <- cbind(maRb(rawdata), maRb(coredata))
                maW(rawdata) <- cbind(maW(rawdata), maW(coredata))
              }
          }
        if(!bg){
          nbgdata <- coredata
          slot(nbgdata, "maGb") <- matrix(0,0,0)
          slot(nbgdata, "maRb") <- matrix(0,0,0)
          data <- nbgdata; rm(nbgdata, coredata)
          fileM <- "normMnbg.xls"
          fileA <- "normAnbg.xls"
          fileMA <- "normMAnbg.xls"
          fileQ <- "qualityNbg.xls"
        }
        else
          {
            data <- coredata
            rm(coredata)
            fileM <- "normM.xls"
            fileA <- "normA.xls"
            fileMA <- "normMA.xls"
            fileQ <- "quality.xls"
          }
        gc()
        
        ## Diagnostic Plots
        if(plot){
          if(echo) cat("Generating ...");
          maDiagnPlots(data, save=TRUE)
          }
        
        ## Quality
        if(quality){
          if(echo) cat("Calculating quality info ...")
          tmp <- maQualityMain(data, path=path, output=TRUE)
          Q.res <- cbind(Q.res, unlist(tmp))
          if(echo) cat("Done \n")
        }

        ## Normalization
        if(!raw){
          defs <- list(norm="p")
          args <- maDotsMatch(maDotsDefaults(opt, defs), formals(args("maNorm")))
          normdata <- do.call("maNorm", c(list(data),args))
          normM <- cbind(normM, maM(normdata))
          normA <- cbind(normA, maA(normdata))
        }
      }


    ## Clean up plots
    if(plot){
      dir.create("DiagnPlots")
      file.copy(dir(pattern="^Plot"), "DiagnPlots", overwrite=TRUE)
      file.remove(dir(pattern="^Plot"))
    }
    
    ## Writing Quality
    if(quality){
      indtmp <- c(1, grep("Flag", rownames(Q.res)),
                  grep("Mean", rownames(Q.res)),  2:nrow(Q.res))
      write.table(cbind(rownames(Q.res[indtmp,]), Q.res[indtmp,]),
                  file=fileQ, col.names=FALSE, row.names=FALSE, quote=F, se="\t")
      if(echo) print(paste("Write to file", fileQ))
      dir.create("QualityXLS")
      file.copy(dir(pattern="^Q\\."), "QualityXLS", overwrite=TRUE)
      file.remove(dir(pattern="^Q\\."))
      assign("QualityXLS", Q.res, envir=.GlobalEnv)
    }


    ## For Normalization
    if(!raw){
      colnames(normM)<- colnames(normA) <- fnames
      normarray <- new("marrayNorm", maA=normA, maM=normM,
                       maLayout=info$layout,
                       maGnames=info$gnames)
      ## Writing M
      write.xls(cbind(maGeneTable(normarray), round(normM, 5)), file=fileM)
      if(echo) print(paste("Write to file", fileM))
      ## Writing A
      write.xls(cbind(maGeneTable(normarray), round(normA, 5)), file=fileA)
      if(echo) print(paste("Write to file", fileA))
      ## Writing M and A into one file
      ind <- as.vector(rbind(1:length(fnames), (1:length(fnames)) + length(fnames)))
      tmp <- round(cbind(normM, normA), 5)[,ind]
      write.xls(cbind(maGeneTable(normarray), tmp), file=fileMA)
      if(echo) print(paste("Write to file", fileMA))
    }
    
    if(raw)
      return(rawdata)
    if(!raw)
      return(normarray)
  }

##################################################################
## END OF FILE
##################################################################
