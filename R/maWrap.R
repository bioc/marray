###########################################################################
# Date : September 19, 2002
# Modify : October, 14, 2002
#
# Runs on R 1.5.1 and above
#
# This file contains wrapper functions for analysis
#
# source("~/Projects/maTools/R/maWrap.R")
#
###########################################################################

widget.TwoSamples <- function(output=TRUE)
  {

    haveTkW <- require("tkWidgets", character.only=TRUE)
    if (!haveTkW)
      stop("This feature requires tkWidgets")

    wlist <- list()
    targetfile <- list(Name="Target File", Value=".txt",
                       toText=function(x) paste(x,collapse = ","),
                       fromText=NULL, canEdit=TRUE, buttonFun = fileBrowser,
                       buttonText = "Browse")
    inputdata <- list(Name="Normalized data", Value="data",
                     toText=function(x) paste(x,collapse = ","),
                     fromText=NULL, canEdit=TRUE, buttonFun = NULL,
                     buttonText = NULL)
    info <- c("Trt", "Ctl", "targetID", "slidesID", "dyesID", "NumID")
    infoValues <- c("Trt", "Ctl", "TargetName", "Slides", "Dyes", "5")

    for(i in 1:length(info))
      {
        test <- list(Name=info[i], Value=infoValues[i],
                     toText=function(x) paste(x,collapse = ","),
                     fromText=NULL, canEdit=TRUE, buttonFun = NULL,
                     buttonText = NULL)
        wlist <- c(wlist, list(test))
      }
    names(wlist) <- info
    widget1 <- list(wList = c(targetfile=list(targetfile), inputdata=list(inputdata), wlist))
    res <- widgetRender(widget1, "Two Samples Adjustment")
    resValues <- values.Widget(res)

    args <- list()
    argsName <- c()
    for(i in 1:length(resValues))
      {
        argsName <- c(argsName, resValues[[i]]$Entry)
        args <-  c(args, list(resValues[[i]]$Value))
      }
    names(args) <- argsName
    args$inputdata <- eval(as.name(args$inputdata))
    names(args)[names(args) == "inputdata"] <- "normdata"
    names(args)[names(args) == "NumID"] <- "RedID"
    outdata <- do.call(maTwoSamples, c(args, list(output=output)))
    return(outdata)
  }

maTwoSamples <- function(targetfile,
                         normdata,
                         Trt,
                         Ctl,
                         targetID="TargetName",
                         slidesID="Slides",
                         dyesID="Dyes",
                         RedID=5,
                         path=".",
                         output=TRUE)
  {
    ## normdata is an marrayNorm objects
     dyeflip <- function(x)
      {
        if(!setequal(as.vector(unique(x[,targetID])), c(Trt, Ctl)))
          res <- NULL
        else
          {
            id <- c(1,2)[as.vector(x[,targetID])==Trt]
            if((as.vector(x[id, dyesID])== RedID))
              res <- 1
            if((as.vector(x[id, dyesID])!= RedID))
              res <- -1
          }
        return(res)
      }

     if(missing(normdata)) normdata <- gpTools()

     target <- maInfo(read.marrayInfo(targetfile))

     dyetype <- unique(target[,"Dyes"])
     if(length(dyetype) !=2)
       stop("Error: More than 2 types of dyes")
     
     if(missing(Trt))    Trt <- as.vector(target[1,targetID])
     if(missing(Ctl))    Ctl <- as.vector(target[2,targetID])
     
     sorttarget <- split(target, target[,slidesID])
     dyeswitch <- unlist(lapply(sorttarget, dyeflip))

     switchM <- switchA <- newtarget <- NULL
     fnames <- colnames(maM(normdata))
     for(i in names(dyeswitch))
       {
         M <- maM(normdata)[,fnames == i]
         A <- maA(normdata)[,fnames == i]
         dyeflip <- dyeswitch[names(dyeswitch)==i]
         switchM <- cbind(switchM, M*dyeflip)
         switchA <- cbind(switchA, A)
         id <- as.vector(target[,slidesID]) == i
         newtarget <- rbind(newtarget, data.frame(t(unlist(apply(target[id,], 2, unique)))))
       }
     rownames(newtarget) <-  as.character(c(1:dim(newtarget)[1]))
     newdata <- normdata
     colnames(switchM) <- colnames(switchA) <- names(dyeswitch)
     slot(newdata, "maM") <- switchM
     slot(newdata, "maA") <- switchA
     
     maTargets(newdata) <- new("marrayInfo",
                               maLabels = as.vector(newtarget[,1]),maInfo = newtarget,
                               maNotes="Generate from maTwoSamples")
     fname <- paste(Trt,"over",Ctl, ".xls", sep="") 
     if(output) write.xls(cbind(maGeneTable(newdata), round(switchM, 5)), file=fname)
     print(paste("Write to file", fname, "\n", sep=" "))
     
     return(newdata)
   }





##################################################################
## END OF FILE
##################################################################
