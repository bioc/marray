## Assuming same Target information
rbind.marrayInfo <- function(..., deparse.level=1)
  {
    data <- list(...)
    newx <- data[[1]]
    for(i in 2:length(data))
      {
        x <- data[[i]]
        if(length(maLabels(x))!=0)
          slot(newx,"maLabels") <- c(maLabels(newx), maLabels(x))
        if(length(maInfo(x))!=0)
          slot(newx,"maInfo")<-  rbind(maInfo(newx), maInfo(x))
        if(length(maNotes(x)) != 0)
          slot(newx,"maNotes")<- paste(maNotes(newx), maNotes(x))
      }
    return(newx)
  }

## Assuming same print-run
## Do not cbind genes
cbind.marrayRaw <- function(..., deparse.level = 1)
  {
    data <- list(...)
    newx <- data[[1]]
    for(x in data[2:length(data)])
      {
        if(length(maGf(x))!=0)
          maGf(newx) <- cbind(maGf(newx), maGf(x))
        if(length(maRf(x))!=0)
          maRf(newx) <- cbind(maRf(newx), maRf(x))
        if(length(maGb(x))!=0)
          maGb(newx) <- cbind(maGb(newx), maGb(x))
        if(length(maRb(x))!=0)
          maRb(newx) <- cbind(maRb(newx), maRb(x))
        if(length(maW(x))!=0)
          maW(newx) <- cbind(maW(newx), maW(x))
        maTargets(newx) <- rbind(maTargets(newx), maTargets(x))
        if(length(maNotes(x)) != 0)
          slot(newx,"maNotes")<- paste(maNotes(newx), maNotes(x))
      }
    return(newx)
  }

cbind.marrayNorm <-  function(..., deparse.level = 1)
  {
    data <- list(...)
    newx <- data[[1]]
    for(x in data[2:length(data)])
      {
        if(length(maM(x))!=0)
          maM(newx) <- cbind(maM(newx), maM(x))
        if(length(maA(x))!=0)
          maA(newx) <- cbind(maA(newx), maA(x))
        if(length(maMloc(x))!=0)
          maMloc(newx) <- cbind(newx@maMloc, x@maMloc)
        if(length(maMscale(x))!=0)
          maMscale(newx) <- cbind(maMscale(newx), maMscale(x))
        if(length(maW(x))!=0)
          maW(newx) <- cbind(maW(newx), maW(x))
        maTargets(newx) <- rbind(maTargets(newx), maTargets(x))
        if(length(maNotes(x)) != 0)
          slot(newx,"maNotes")<- paste(maNotes(newx), maNotes(x))
      }
    return(newx)
  }


