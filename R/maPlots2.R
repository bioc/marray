## R function for finding gene patterns
## Created Date: Sept. 18, 2002

##########################################################3
## Quality control plots
## Created from maDiagnPlots.R
## Oct, 9, 2003
## Modified by Agnes Paquet
## Things to do
## 1) Spatial plots on ranked data
## 2) Remove flagged spots from MA-plots
## 3) Other diagnostic plots?
## 4) More efficient error handling with gpTools
## 5) More controls on the parameter, to handle possible errors
## 6) Add the png format

maDiagnPlots <-  function(mraw,
                          mNorm = NULL,
                          save = TRUE,
                          fname = NULL,
                          dev = "png",  #set default to be png
                          pch,
                          col,
                          DEBUG=FALSE,
                          ...) 
{
  rm.na <- 
  function (x) 
    {
      ind <- is.na(x) | is.nan(x) | is.infinite(x)
      return(x[!ind])
    }

  if (DEBUG) print("function starting")
  
  opt <- list(...)
  data <- mraw[,1]
  if(is.null(mNorm))
    {
      ifelse(length(maW(data)) != 0, 
             defs <- list(norm="p", subset=eval(maW(data)==0)),
             defs <- list(norm="p"))
      args <-  maDotsDefaults(opt, defs)
      normdata <- do.call("maNorm", c(list(data),args))
    }
 
  else
    {
      ifelse(length(maW(data)) != 0, 
             defs <- list(norm="p", subset=eval(maW(data)==0)),
             defs <- list(norm="p"))
      args <-  maDotsDefaults(opt, defs)
      normdata <- mNorm[,1]
    }
  
  if(is.null(fname))
    f <- colnames(maGf(data))
   else f <- fname

  tmp <- unlist(strsplit(f, "\\."))
  ifelse(length(maGb(data))!=0, bg <- "Plot", bg <- "Plot.nbg")
  fstart <- paste(tmp[-length(tmp)], collapse=".")

  if(length(grep(f, dir())) != 0)
    {
      tmp <- readLines(f, n=40)
      subnames <- paste(tmp[grep("DateTime", tmp)], tmp[grep("PMT", tmp)])
      if(length(subnames) ==0) subnames <- ""
    }
  else
    subnames <- ""
  
  if(dev == "jpeg" | dev == "jpg")
    {
      dev <- "jpeg"  #R recognizes only jpeg
      fname <- paste(bg, fstart, "jpeg", sep=".")
      def <- list(dev=list(quality=100, width=1400, height=1400),
                  main= paste(fname, ": Pre- and post- normalization"))
    }

  if(dev == "postscript")
    {
      fname <- paste(bg, fstart, "ps", sep=".")
      def <- list(dev=list(paper="special", width=14, height=14),
                  main= paste(fname, ": Pre- and post- normalization"))
    }

  if(dev == "png")
    {
      fname <- paste(bg, fstart, "png", sep=".")
      def <- list(dev=list(width=1400, height=900), ## modified
                  main= paste(fname, ": Pre- and post- normalization"))
    }

  if(!is.element(dev, c("jpeg","png","postscript", "jpg")))
    {
      print("Format error, format will be set to PNG")
      fname <- paste(bg, fstart, "png", sep=".")
      def <- list(dev=list(width=1400, height=900),  ## modified
                  main= paste(fname, ": Pre- and post- normalization"))
    }

   
  if (!is.null(opt)) 
    def <- maDotsDefaults(opt, def)
  args <- c(list(fname), def$dev)
  
  if(save) 
    do.call(dev, args)
  
  ## Some controls info
  ifelse(missing(col), colcode <- unique(as.integer(maControls(data))+1),
         colcode <- col)
  ctlcode <- levels(maControls(data))
  names(colcode) <- ctlcode


  if(DEBUG) print("start layout")
  ## Layout 
  layout(matrix(c(13, 1,1,2,2,13,3,3,5,5,13,4,4,6,6,13,7,7,
                  9,10,13,8,8,9,10, 13,11,11,12,12), 5, 6),
         height=c(1, 5, 5, 5, 5), 
         width = c(5.5, 5.5, 2 ,5.5, 2, 5.5))
  
  if(DEBUG) print("start 1")
  ## 1) MA-plot (Before Normalization)
  par(mar=c(3,2,3,2))
  y <- max(maM(data), na.rm=TRUE) + 2
  x <- min(maA(normdata), na.rm=TRUE) 
  defs <- maDefaultPar(data, x="maA", y="maM", z="maPrintTip")
  flags <- as.integer(as.factor(maW(data)))
  badspot.func <- maText(maW(data) != 0,
                         labels=as.character(flags[maW(data) != 0]),
                         col=flags[maW(data) != 0], cex=0.5)
  args <- maDotsMatch(maDotsDefaults(opt, defs), formals(args("maPlot")))
  maPlot(data, ylim=c(min(maM(data), na.rm=TRUE),y), xlim=c(0,16), text.func = badspot.func,
         legend.func=NULL, main="MA-plot: raw", cex=0.6)
  legend.func <- do.call("maLegendLines", defs$def.legend)
  legend.func(x, y)


  
  if(DEBUG) print("start 2")
  ## 2) MA-plot (After Normalization)
  par(mar=c(3,2,3,2))
  ##  y <- max(maM(normdata), na.rm=TRUE) + 2
  ##  x <- min(maA(normdata), na.rm=TRUE) 
  defs <- maDefaultPar(normdata, x="maA", y="maM", z="maPrintTip")
  Cindex <- as.character(maControls(data)) != "probes"
  filternum <- -30
  if(length(c(1:length(Cindex))[Cindex]) != 0)
    {
      speccode <- rep(NA, length(Cindex))
      for(i in 1:length(ctlcode))
        speccode[maControls(data)==ctlcode[i]] <- colcode[i]
      newCindex <- Cindex[maW(normdata) > filternum]
      qualspot.func <- maText(newCindex,
                              labels=rep("*", length(newCindex))[newCindex],
                              col=speccode[maW(normdata) > filternum][newCindex], cex=1)
    }
  else
    {
      qualspot.func <- NULL
    }
  ## use the same scale as MA-plot before normalization
  maPlot(normdata[maW(normdata) > filternum,], ylim=c(min(maM(data), na.rm=TRUE),y), xlim=c(0,16),
         text.func = qualspot.func,lines.func=NULL, legend.func =NULL, main="MA-plot: Norm", cex=0.6)
  if(length(c(1:length(Cindex))[Cindex]) != 0)
    legend(x, y,  ctlcode[ctlcode!="probes"], col=colcode, pch="*")

  
  if(DEBUG) print("start 3, 4")
  ## 3 & 4) maM (Before Normalization)
  par(mar=c(2,3,5,2))
  BYcol <- maPalette(high="yellow", low="blue", mid="grey", k=50)
  #convert marrayRaw into marrayNorm, needed for assignment
  tmpNorm <- maNorm(data, norm="none")
  tmpNorm@maM <- as.matrix(rank(maM(data)))
  tmp <- maImage(tmpNorm, x="maM", main="Spatial: M-Raw", bar=FALSE, col=BYcol)
  par(mar=c(2,1,2,4))
  maColorBar(tmp$x.bar, horizontal = FALSE, col = BYcol,  main = "")

  if(DEBUG) print("start 5 & 6")
  ## 5 & 6) maM (After Normalization)
  par(mar=c(2,3,5,2))
  tmpNorm <- normdata
  tmpNorm@maM <- as.matrix(rank(maM(normdata)))
  tmp <- maImage(tmpNorm, x="maM", main="Spatial: M-Norm", bar=FALSE, col=BYcol)
  par(mar=c(2,1,2,4))
  maColorBar(tmp$x.bar, horizontal = FALSE, col = BYcol,  main = "")

  if(DEBUG) print("start 7 & 8")
  ## 7 & 8) maA 
  par(mar=c(2,3,5,2))
  Bcol <- maPalette(high="blue", low="white", k=50)
  tmp <- maImage(data, x="maA", main="Spatial: A", bar=FALSE, col=Bcol)
  par(mar=c(2,1,2,4))
  maColorBar(tmp$x.bar, horizontal = FALSE, col = Bcol,  main = "")

  if(DEBUG) print("start 9")
  ## 9
  ifelse(length(maRb(data))!=0 , RS2N <- as.vector(log(maRf(data) / maRb(data),2)),
         RS2N <- as.vector(log(maRf(data),2)))
  lab <- paste("mean:", round(mean(RS2N), 2),",", "var:", round(var(rm.na(RS2N)), 2))
  hist(rm.na(RS2N), main=lab, col="red", freq=FALSE, ylim=c(0,1.1));

  if(length(maControls(data))!=0)
    {
      tmp <- split(RS2N, maControls(data))
      tmp2 <- lapply(tmp, function(x){density(rm.na(x), sd(rm.na(RS2N))/4, na.rm=TRUE)})
      for(i in 1:length(tmp2))
        lines(tmp2[[i]], lwd=2, col=colcode[i])
      xrange <- range(rm.na(RS2N))
      xcood <- xrange[1] + (xrange[2]-xrange[1]) * 0.7
      legend(xcood, 1, names(tmp), lty=1, lwd=2, col=colcode, cex=0.8)
    }
  
  if(DEBUG) print("start 10")
  ##10
  ifelse(length(maGb(data))!=0, GS2N <- as.vector(log(maGf(data) / maGb(data),2)),
         GS2N <- as.vector(log(maGf(data),2)))
  lab <- paste("mean:", round(mean(GS2N), 2),",", "var:", round(var(rm.na(GS2N)), 2))
  hist(rm.na(GS2N), main=lab, col="green", freq=FALSE, ylim=c(0,1.1));

  if(length(maControls(data))!=0)
    {
      tmp <- split(GS2N, maControls(data))
      tmp2 <- lapply(tmp, function(x){density(rm.na(x), sd(rm.na(GS2N))/4, na.rm=TRUE)})
      for(i in 1:length(tmp2))
        lines(tmp2[[i]], lwd=2, col=colcode[i])
      xrange <- range(rm.na(RS2N))
      xcood <- xrange[1] + (xrange[2]-xrange[1]) * 0.7
      legend(xcood, 1, names(tmp), lty=1, lwd=2, col=colcode, cex=0.8)
    }

  if(DEBUG) print("start 11")
  ## 11
  if(length(maControls(data))!=0)
    maDotPlots(data, x="maM", col=colcode)

  if(DEBUG) print("start 12")
  ## 12
  if(length(maControls(data))!=0)
    maDotPlots(data, x="maA",  col=colcode)

  if(DEBUG) print("start 13")
  ## 13
  layout(1)
  par(mar=c(2,2,4,2))
  mtext(def$main, line=3)
  mtext(subnames, line=2, cex = 0.7)
  mtext(paste("Call:", maNormCall(normdata)[3]), line=1, cex = 0.7)
  

  if(DEBUG) print("Done...")
  ## Finishing
  if (save == TRUE) {
    cat(paste("save as", fname, "\n"))
    dev.off()
  }
}


maDotPlots <- function(data,
                       x=list("maA"),
                       id="ID",
                       pch,
                       col,
                       nrep=3,
                       ...)
  {
    newdata <- NULL
    for(i in x)
      newdata <- cbind(newdata, eval(call(i, data)))

    if(!is.null(newdata))
      {
        colnames(newdata) <- x
        xlim <- range(newdata, na.rm=TRUE)
      }
    else
      stop("No specified data")

    Cindex <- maControls(data) != "probes"
    ifelse(missing(pch), pchcode <- (1:ncol(newdata))+15, pchcode <- pch)
    ifelse(missing(col), colcode <- unique(as.integer(maControls(data))+1),
           colcode <- col)
    names(colcode) <- levels(maControls(data))
    Ctl <- cbind(maInfo(maGnames(data)), maControls(data))

    IDindex <- grep(id, colnames(Ctl))
    y <- split(Ctl, Ctl[,ncol(Ctl)])
    if(length(y[names(y) != "probes"]) != 0)
      {
        exty <- lapply(y[names(y) != "probes"], function(x){
          ext <- split(x, x[, IDindex])
          extid <- lapply(ext, function(xx){as.integer(row.names(xx))})
          extid[lapply(extid, length) > nrep]
        })
        exty <- exty[lapply(exty, length) != 0]
        
        ylim <- c(1, sum(unlist(lapply(exty, length))))
        
        par(mar=c(4,7,2,2))
        plot(1,1, type="n", xlim=xlim, ylim=ylim, axes=FALSE, xlab=unlist(x), ylab="")
        ii <- 1
        for(i in 1:length(exty))
          for(j in 1:length(exty[[i]]))
            {
              ind <- exty[[i]][[j]]
              for(k in 1:ncol(newdata))
                {
                  points(newdata[ind,k], rep(ii, length(newdata[ind,k])), pch=pchcode[k],
                         col=colcode[names(exty)[i]])
                  points(median(newdata[ind, k], na.rm=TRUE), ii, pch=18, col="black")
                }
              ii <- ii + 1
            }
        axis(1)
        lab <- paste(unlist(lapply(exty, names)), " (n=",
                     unlist(lapply(exty, lapply, length)), ") ", sep="")
        axis(2, at=1:ylim[2], labels=lab, las=2, cex.axis=0.6) 
        box()
      }
    else
      {
        plot(1, 1, axes=FALSE, xlab="", ylab="", type="n")
        text(1, 1, "No Control Genes")
        box()
      }
    return()
  }
