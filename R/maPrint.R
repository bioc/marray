###########################################################################
# maPrint.R
#
# Print methods for microarray classes
# Change Print to summary
###########################################################################
# marrayLayout class


setMethod("summary", signature(object="marrayLayout"), function(object) {
    cat("Array layout:\t Object of class marrayLayout. \n\n")

    cat("Total number of spots:\t\t\t")
    cat(maNspots(object))
    cat("\n")

    cat("Dimensions of grid matrix:\t\t")
    cat(paste( maNgr(object), "rows by", maNgc(object), "cols"))
    cat("\n")
    
    cat("Dimensions of spot matrices:\t\t")
    cat(paste(maNsr(object), "rows by", maNsc(object), "cols"))
    cat("\n\n")

    if(length(maNspots(object))==1)
      nsub<-length((1:maNspots(object))[maSub(object)])
    else
      nsub<-NULL
    cat(paste("Currently working with a subset of ", nsub, "spots.\n\n", sep=""))

    cat("Control spots: \n")
    if(length(maControls(object))!=0){
      tmp <- table(maControls(object))
      cat(paste("There are  ",  length(tmp),  "types of controls : \n"))
      print(tmp[1:min(10, length(tmp))])
      if(length(tmp) > 10) cat ("...")
    }
    cat("\n\n")
    cat("Notes on layout: \n", maNotes(object), "\n")
  })
          

#########################
# marrayInfo class

setMethod("summary",
          signature(object="marrayInfo"),
          function(object){
  cat("Object of class marrayInfo. \n\n")
            
  nr <- length(maLabels(object))
  ds <- dim(maInfo(object))
    
  if(!is.null(ds))
    {
      if(nr!=0)
        {
          tmp<-data.frame(maLabels(object),maInfo(object))
          dimnames(tmp)[[2]][1]<-"maLabels"
          dimnames(tmp)[[2]][-1]<-dimnames(maInfo(object))[[2]]
        }
      else
	tmp<-maInfo(object)
      print(tmp[1:min(nrow(tmp), 10),])
      if(nrow(tmp) > 10)  cat("... \n")
    }

  if(is.null(ds)&(nr!=0))
    {
      tmp <- as.vector(maLabels(object))
      print(tmp[1:min(length(tmp), 10)])
      if(nr > 10)  cat("... \n")
    }
  
  cat("\nNumber of labels: ", nr, " \n")
  cat("Dimensions of maInfo matrix: ", ds[1], " rows by ", ds[2], " columns\n\n")
  cat("Notes: \n", maNotes(object), "\n")
})

#########################
# marrayRaw class
setMethod("summary", signature(object="marrayRaw"), function(object){
  if((length(maGf(object))==0) | (length(maRf(object))==0))
    cat("Input is empty \n")
  else
    {
      cat("Pre-normalization intensity data:\t Object of class marrayRaw. \n\n")
      
      cat(paste("Number of arrays: \t",ncol(maM(object)), " arrays.\n \n", sep=""))
      cat("A) Layout of spots on the array: \n")
      summary(maLayout(object))

      cat("\n")

      cat("B) Samples hybridized to the array: \n")
      summary(maTargets(object))
      cat("\n")

      cat("C) Summary statistics for log-ratio distribution: \n")
      results <- apply(maM(object),2, function(x){round(summary(x),2)})
      if(is.matrix(results))
        print(t(results))

      if(is.list(results)){
        x.min <- x.1q <- x.med <- x.mean <- x.3q <- x.max <- x.NA <- NULL
        x.names <- names(results)
        if(is.null(x.names)) x.names <- as.character(1:length(results))
        for(i in 1:length(results)){
          x.min <- c(x.min, results[[i]][1])
          x.1q <- c(x.1q, results[[i]][2])
          x.med <- c(x.med, results[[i]][3])
          x.mean <- c(x.mean, results[[i]][4])
          x.3q <- c(x.3q, results[[i]][5])
          x.max <- c(x.max, results[[i]][6])
          x.NA <- c(x.NA, results[[i]][7])
        }
        y <- data.frame(x.names, x.min, x.1q, x.med, x.mean, x.3q, x.max, as.vector(x.NA))
        dimnames(y)[[2]] <-
          c("   ", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max", "NA")
        print(y)
      }

      cat("\n")

      cat("D) Notes on intensity data: \n")
      cat(maNotes(object))
      cat("\n")
    }
})

#########################
# marrayNorm class

setMethod("summary",
          signature(object="marrayNorm"),
          function(object)
          {
            if((length(maA(object))==0) & (length(maM(object))==0))
              cat("Input is empty \n")
            else
              {
                cat("Normalized intensity data:\t Object of class marrayNorm.\n \n")
                
                cat("Call to normalization function:\n")
                if(length(maNormCall(object)) > 1)
                  print(maNormCall(object))
                
                cat(paste("\nNumber of arrays: \t",ncol(maM(object)), " arrays.\n", sep=""))

                cat("\nA) Layout of spots on the array: \n")
                summary(maLayout(object))
                cat("\n")

                cat("B) Samples hybridized to the array: \n")
                summary(maTargets(object))
                cat("\n")
                
                cat("C) Summary statistics for log-ratio distribution: \n")
                if(length(maM(object))!=0)
                  results <- apply(maM(object),2, function(x){round(summary(x),2)})
                if(is.matrix(results))
                  print(t(results))

                if(is.list(results)){
                  x.min <- x.1q <- x.med <- x.mean <- x.3q <- x.max <- x.NA <- NULL
                  x.names <- names(results)
                  if(is.null(x.names)) x.names <- as.character(1:length(results))
                  for(i in 1:length(results)){
                    x.min <- c(x.min, results[[i]][1])
                    x.1q <- c(x.1q, results[[i]][2])
                    x.med <- c(x.med, results[[i]][3])
                    x.mean <- c(x.mean, results[[i]][4])
                    x.3q <- c(x.3q, results[[i]][5])
                    x.max <- c(x.max, results[[i]][6])
                    x.NA <- c(x.NA, results[[i]][7])
                  }
                  y <- data.frame(x.names, x.min, x.1q, x.med, x.mean, x.3q, x.max, as.vector(x.NA))
        dimnames(y)[[2]] <-
          c("   ", "Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max", "NA")
                  print(y)
                }
                cat("\n")
                
                cat("D) Notes on intensity data: \n")
                cat(maNotes(object))
                cat("\n")
              }
          })

setClass("ShowLargeObject")
setIs("marrayRaw","ShowLargeObject")
setIs("marrayNorm","ShowLargeObject")
setIs("marrayInfo","ShowLargeObject")
setIs("marrayLayout","ShowLargeObject")

setMethod("show","ShowLargeObject",
#  Print and show method large data objects
#  Modified from Gordon Smyth
#  March 2004
function(object) {
  cat("An object of class \"",class(object),"\"\n",sep="")
  for (what in slotNames(object)) {
    x <- slot(object,what)
    if((length(x) > 0) | !is.null(x)) {
      cat("@",what,"\n",sep="")
      printHead(x)
      cat("\n")
    }
  }
})
