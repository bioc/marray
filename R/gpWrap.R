###########################################################################
# Date : September 19, 2002
# Modify : April 14, 2003
#
# Runs on R 1.5.1 and above
#
# This file contains wrapper functions for GenePix files
## source("~/Projects/maTools/R/gpWrap.R")
## source("~/Projects/maTools/R/maAnalysis.R")
###########################################################################

###########################################################################
## This is a wrapper function specifially to generate diagnostic plots and
## quality file for every genepix files in the current working directory

gpTools <- function(fnames,
                    path=".",
                    galfile,
                    bg=TRUE,
                    plot=TRUE,
                    quality=TRUE,
                    fill = TRUE,
                    raw = FALSE,
                    echo=TRUE,
                    ...)
  {
    ## SetUP QualityInfo:
    QualityInfo <- c("file","Date","Pmt",
                     "Flaginfo.-75","Flaginfo.-50","Flaginfo.0", "RS2N.Median", "GS2N.Median",
                     "CtlA.Empty.Median", "CtlA.Negative.Median","CtlA.Positive.Median",
                     "CtlM.Empty.Median", "CtlM.Negative.Median",  "CtlM.Positive.Median", 
                     "CtlA.probes.Min.","CtlA.probes.1st Qu.","CtlA.probes.Median",
                     "CtlA.probes.3rd Qu.","CtlA.probes.Max.",
                     "CtlM.probes.Min.","CtlM.probes.1st Qu.","CtlM.probes.Median",
                     "CtlM.probes.3rd Qu.","CtlM.probes.Max.",
                     "CTLNum.Empty","CTLNum.Negative","CTLNum.Positive","CTLNum.probes",
                     "Layout1","Layout.GridR","Layout.GridC","Layout.SpotR","Layout.SpotC")
    
    opt <- list(...)
    normM <- normA <- NULL
    if(missing(fnames)) fnames <- dir(path, pattern="*\\.gpr$")
    if(missing(galfile))
      {
        tmp <- dir(path, pattern="*\\.gal$")
        ifelse(length(tmp)==0, galfile <- fnames[1], galfile <- tmp)
      }

   
    if(quality) Q.res <- NULL
    if(raw) rawdata <- new("marrayRaw")
    
    for(i in fnames)
      {
        args <- maDotsMatch(maDotsMatch(opt, list(galfile=i, path=path)),
                            formals(args("read.Galfile")))
        core.info <- do.call("read.Galfile", args)
        maControls(core.info$layout) <- maGenControls(core.info$gnames)
        defs <- list(fnames=i,
                     path=path,  name.W="Flags",
                     name.Gf = "F532 Median",
                     name.Rf = "F635 Median",                               
                     layout = core.info$layout,
                     gnames=core.info$gnames,
                     comment.char="",
                     fill=TRUE, quote="\"")
        args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.GenePix")))
        coredata <- do.call("read.GenePix", args)

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

        ## Number of missing values
        missingvalues <- (sum(is.na(maM(data)[,1]))/maNspots(data))


        ## Diagnostic Plots
        if(plot){
          if(echo) cat("Generating ...");
          if(table(maW(data))["0"] / maNspots(data) > 0.45)
            {
              maDiagnPlots(data, save=TRUE); if(echo) cat("Done. \n")
            } else
          cat("Plot Failed: Percentage of good spots are", table(maW(data))["0"] / maNspots(data), "\n")
        }

        ## Quality
        if(quality){
          if(echo) cat("Calculating quality info ...")
          tmp <- maQualityMain(data, path=path, output=TRUE)
          Q.res <- cbind(Q.res, unlist(tmp)[QualityInfo])
          if(echo) cat("Done. \n")
        }

        ## Normalization
        if(!raw){
          defs <- list(norm="p")
          args <- maDotsMatch(maDotsDefaults(opt, defs), formals(args("maNorm")))
          if(missingvalues < 0.7)
            {
              normdata <- do.call("maNorm", c(list(data),args))
              normM <- cbind(normM, maM(normdata))
              normA <- cbind(normA, maA(normdata))
            } else
          {
            cat("Normalization Failed: Percentage of missing values are", missingvalues, "\n")
            cat("Use Median normalization instead. \n")
            args$norm <- "m"
            normdata <- do.call("maNorm", c(list(data),args))
            normM <- cbind(normM, maM(normdata))
            normA <- cbind(normA, maA(normdata))
          }
        } ## Normalization
      } ## loop through different files

   
    ## Clean up plots 
    if(plot){ 
      dir.create("DiagnPlots")
      if(!length(dir(pattern="^Plot"))==0)
        {
          file.copy(dir(pattern="^Plot"), "DiagnPlots", overwrite=TRUE)
          file.remove(dir(pattern="^Plot"))
        }
    }

     
    ## Writing Quality
    if(quality){
##      indtmp <- c(1, grep("Flag", rownames(Q.res)),
##                  grep("Mean", rownames(Q.res)),  2:nrow(Q.res))
##      write.table(cbind(rownames(Q.res[indtmp,]), Q.res[indtmp,]),
##                  file=fileQ, col.names=FALSE, row.names=FALSE, quote=F, sep="\t")
      write.table(Q.res, file=fileQ, col.names=FALSE, row.names=TRUE, quote=FALSE, sep="\t")
      if(echo) print(paste("Write to file", fileQ))
      dir.create("QualityXLS")
      file.copy(dir(pattern="^Q\\."), "QualityXLS", overwrite=TRUE)
      file.remove(dir(pattern="^Q\\."))
      assign("QualityXLS", Q.res, envir=.GlobalEnv)
    }

    ## Writing M (Normalization)
    if(!raw){
      if(echo) cat("Assuming all data comes from the same print-run \n")
      if(echo) cat("Reading Gal file ...")
      args <- maDotsMatch(c(opt, path=path), formals(args("read.Galfile")))
      info <- do.call("read.Galfile", c(list(galfile), list(path), args))
      maControls(info$layout) <- maGenControls(info$gnames)
      if(echo) cat("done \n ")
      
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
