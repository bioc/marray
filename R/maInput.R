###########################################################################
##  Function: read.marrayLayout
##
##
##
###########################################################################

read.fname <-
  function(fname,
           skip,
           sep="\t",
           quote="\"",
           ...)
  {
    h<-strsplit(readLines(fname, n=skip+1),split=sep)
    h<-as.list(unlist(h[[length(h)]]))
    names(h)<- unlist(h)
    dat<-scan(fname,quiet=TRUE,what=h, sep=sep, skip=skip+1, quote=quote, ...)
    return(dat)
  }


read.marrayLayout<-
  function(fname=NULL,
           ngr,
           ngc,
           nsr,
           nsc,
           pl.col=NULL,
           ctl.col=NULL,
           sub.col=NULL,
           notes=fname,
           skip,
           sep="\t",
           quote="\"",
           ...)
{
  # Read layout file
  if(!is.null(fname)){
    h<-strsplit(readLines(fname, n=skip+1),split=sep)
    h<-as.list(unlist(h[[length(h)]]))
    names(h)<-unlist(h)
    dat<-scan(fname,quiet=TRUE,what=h, sep=sep, skip=skip+1, fill=TRUE, quote=quote, ...)
  }
  else
    notes <- "No Input File"

  # Generate spot coordinates
  nspots<-ngr * ngc * nsr * nsc
  layout <- new("marrayLayout", maNgr=ngr, maNgc = ngc,
                maNsr=nsr, maNsc=nsc, maNspots=nspots, maNotes=notes)

  if(!is.null(fname)){
  # Get Uneven layout
    if(!is.null(sub.col)){
      subvalue <- intersect(c(1:nspots), as.numeric(dat[[sub.col]]))
      maSub(layout)  <- rep(0, nspots)
      maSub(layout)[subvalue] <- 1
    }

  # Get plate id
    if(!is.null(pl.col))
      maPlate(layout) <- as.factor(dat[[pl.col]])

  # Get control
    if(!is.null(ctl.col))
      maControls(layout) <- as.factor(dat[[ctl.col]])
  }
  return(layout)
}

###########################################################################
# Function: read.marrayInfo
#
###########################################################################
read.marrayInfo <-
  function(fname,
           info.id=NULL,
           labels=NULL,
           notes=fname,
           sep="\t",
           skip=0,
           quote="\"",
           ...)
{
  h<-strsplit(readLines(fname, n=skip+1),split=sep)
  h<-as.list(unlist(h[[length(h)]]))
  names(h)<-gsub("\"", "", unlist(h))
  dat <- read.table(fname, sep=sep, skip=skip+1, fill=TRUE, quote=quote, ...)
  colnames(dat) <- h

  descript <- new("marrayInfo", maNotes=notes)

  ## Enter Description
  if(is.null(info.id)) info.id <- 1:ncol(dat)
  maInfo(descript) <- as.data.frame(dat[,info.id])

  ## Enter Predefine Labels
  if(length(labels) == nrow(dat))
    maLabels(descript) <- as.vector(labels)
  else
    {
      if(is.null(labels))
        labels <- 1
      maLabels(descript) <- as.character(as.vector(dat[,labels]))
    }
  return(descript)
}

###########################################################################
# Function: read.marrayRaw
#
###########################################################################
read.marrayRaw<-
  function(fnames,
           path=".",
           name.Gf=NULL,
           name.Gb=NULL,
           name.Rf=NULL,
           name.Rb=NULL,
           name.W=NULL,
           layout=NULL,
           gnames=NULL,
           targets=NULL,
           notes=NULL,
           skip=0,
           sep="\t",
           quote="\"",
           DEBUG=FALSE,
           ...)
{
  if(is.null(path))
    fullfnames <- fnames
  else
    fullfnames <- file.path(path, fnames)

  fname <- fullfnames[1]

  
  # Intensity data
  Gf<-Gb<-Rf<-Rb<-W<- NULL
  if(is.null(name.Rf)) Rf <- matrix(0,0,0)
  if(is.null(name.Gf)) Gf <- matrix(0,0,0)
  if(is.null(name.Gb)) Gb <- matrix(0,0,0)
  if(is.null(name.Rb)) Rb <- matrix(0,0,0)

  for(f in fullfnames)
  {
    ## Calculate Skip
    if(DEBUG) cat("Calculating skip  ... ")
    y <- readLines(f, n=100)
    skip <- grep(name.Gf, y)[1] - 1
    if(DEBUG) cat(skip, "done \n")

    cat("Reading ... ", f)
    h<-strsplit(readLines(f, n=skip+1),split=sep)
    h<-as.list(unlist(h[[length(h)]]))
    names(h)<-gsub("\"", "", unlist(h))
    
    dat<-scan(f,quiet=TRUE,what=h, sep=sep, skip=skip+1, quote=quote, ...)
    if(!is.null(name.Gf)) Gf<-cbind(Gf,as.numeric(dat[[name.Gf]]))
    if(!is.null(name.Gb)) Gb<-cbind(Gb,as.numeric(dat[[name.Gb]]))
    if(!is.null(name.Rf)) Rf<-cbind(Rf,as.numeric(dat[[name.Rf]]))
    if(!is.null(name.Rb)) Rb<-cbind(Rb,as.numeric(dat[[name.Rb]]))
    if(!is.null(name.W)) W <-cbind(W,as.numeric(dat[[name.W]]))
  }

  if(!is.null(name.W)) colnames(W) <- fnames
  if(!is.null(name.Gb)) colnames(Gb) <- fnames
  if(!is.null(name.Rb)) colnames(Rb) <- fnames
  if(!is.null(name.Rf)) colnames(Rf) <- fnames
  if(!is.null(name.Gf)) colnames(Gf) <- fnames

  ## Add Notes
  if(is.null(notes)) notes <- ""

  mraw<-new("marrayRaw", maRf=Rf, maRb=Rb, maGf=Gf, maGb=Gb, maNotes =notes)
  
  ## Add other informations ad Weights
  if(!is.null(layout)) maLayout(mraw) <- layout
  if(!is.null(gnames)) maGnames(mraw) <- gnames
  if(!is.null(targets)) maTargets(mraw) <- targets
  if(!is.null(W)) maW(mraw) <- W
  cat(" ... done \n")
  return(mraw)
}

#########################  END FIRST PART
###
### TEST FUNCTION
## read.layout(fname, 4, 4, 20, 20, pl.col="PLATE",
## ctl.col=NULL, gnames.col=c("SUID","LUID","Gene Name"), notes="BLAH, BLAH")
### END FUNCTION
###########################################################################

###########################################################################
#  Read SPOT
#
read.Spot <-  function(fnames = NULL,
                       path=".",
                       name.Gf = "Gmean",
                       name.Gb = "morphG",
                       name.Rf = "Rmean",
                       name.Rb = "morphR",
                       name.W= NULL,
                       layout = NULL,
                       gnames = NULL,
                       targets = NULL,
                       notes=NULL,
                       skip=0,
                       sep="\t",
                       quote="\"",
                       ...)
  {
    ## If fnames not specified, read everything in the dir
    if(is.null(fnames))
      fnames <- dir(path=path, pattern=paste("*", "spot", sep="\."))

    ## Calculate Skip
    if(skip == 0)
      {
        y <- readLines(file.path(path, fnames[1]), n=100)
        skip <- grep(name.Gf, y)[1] - 1
      }
    
    if(is.null(notes)) notes <- "Spot Data"

    mraw <- read.marrayRaw(fnames =fnames,
                           path=path,
                           name.Gf = name.Gf,
                           name.Gb = name.Gb,
                           name.Rf = name.Rf,
                           name.Rb = name.Rb,
                           name.W= name.W,
                           layout = layout,
                           gnames = gnames,
                           targets = targets,
                           notes = notes,
                           skip= skip,
                           sep= sep,
                           quote=quote,
                           ...)
    return(mraw)
  }



###########################################################################
#  Read SPOT
#
read.GenePix <-  function(fnames = NULL,
                          path=NULL,
                          name.Gf = "F532 Median",
                          name.Gb = "B532 Median",
                          name.Rf = "F635 Median",
                          name.Rb = "B635 Median",
                          name.W= "Flags",
                          layout = NULL,
                          gnames = NULL,
                          targets = NULL,
                          notes=NULL,
                          skip=0,
                          sep="\t",
                          quote="\"",
                          DEBUG=FALSE,
                          ...)
  {
    
    opt <- list(...)
    ## If fnames not specified, read everything in the dir
    if(is.null(fnames))
      fnames <- ifelse(is.null(path), dir(pattern="*\\.gpr$"), dir(path, pattern="*\\.gpr$"))
    else{
      if(!is.null(targets))
        fnames <- maInfo(targets)[,1]
    }
    
    if(DEBUG) print(fnames)
    if(is.null(gnames) | is.null(layout))
      {
        if(DEBUG) print("Reading Galfile ... ")
        opt <- list(...)
        defs <- list(galfile = fnames[1], path=path, info.id = c("Name", "ID"),
                     labels = "ID", sep = sep, quote=quote, fill=TRUE, check.names=FALSE,
                     as.is=TRUE, ncolumns = 4)
        gal.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.Galfile")))        
        gal <- do.call("read.Galfile", gal.args)
        if(is.null(gnames)) gnames <- gal$gnames
        if(is.null(layout)) layout <- gal$layout
        if(DEBUG) print("done \n ")
      }

    if(DEBUG) cat("Setting up controls status ... ")
    layout@maControls <- as.factor(maGenControls(gnames))
    if(DEBUG) cat("done \n ")
    
    if(is.null(notes)) notes <- "GenePix Data"

    if(DEBUG) cat("Calling read.marrayRaw ... \n")
    defs <- list(fnames = fnames, path=path,
                 name.Gf = name.Gf, name.Gb=name.Gb, name.Rf=name.Rf, name.Rb=name.Rb,
                 name.W=name.W, layout = layout, gnames=gnames, targets=targets,
                 notes = notes, skip=skip, sep=sep, quote=quote, fill=TRUE,
                 check.names=FALSE,  as.is=TRUE)
    maRaw.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.marrayRaw")))        
    mraw <- do.call("read.marrayRaw", maRaw.args)

    return(mraw)
  }


read.SMD <-  function(fnames = NULL,
                      path=".",
                      name.Gf = "CH1I_MEAN",
                      name.Gb = "CH1B_MEDIAN",
                      name.Rf = "CH2I_MEAN",
                      name.Rb = "CH2B_MEDIAN",
                      name.W= NULL,
                      layout = NULL,
                      gnames = NULL,
                      targets = NULL,
                      notes=NULL,
                      skip=0,
                      sep="\t",
                      quote="",
                      ...)
  {

    ## If fnames not specified, read everything in the dir
    if(is.null(fnames))
      fnames <- dir(path=path, pattern=paste("*", "xls", sep="\."))


    ## Calculate Skip
    if(skip == 0)
      {
        y <- readLines(file.path(path, fnames[1]), n=100)
        skip <- grep(name.Gf, y)[1] - 1
      }
    
    if(is.null(notes)) notes <- "SMD Data"

    mraw <- read.marrayRaw(fnames =fnames,
                           path=path,
                           name.Gf = name.Gf,
                           name.Gb = name.Gb,
                           name.Rf = name.Rf,
                           name.Rb = name.Rb,
                           name.W= name.W,
                           layout = layout,
                           gnames = gnames,
                           targets = targets,
                           notes = notes,
                           skip= skip,
                           sep= sep,
                           quote=quote,
                           ...)
    return(mraw)
  }


###########################################################################

read.Galfile <- function (galfile,
                          path=".",
                          info.id = c("Name", "ID"),
                          labels = "ID",
                          notes = "",
                          sep = "\t",
                          skip = NULL,
                          ncolumns = 4,
                          ...)
{
  opt <- list(...)

  if(!is.null(path))
    f <- file.path(path, galfile)
  else
    f <- galfile

  y <- readLines(f, n=100)

  skip <- intersect(grep("ID", y), grep("Name", y))[1] - 1
  defs <- list(file=f, path = path, sep=sep, skip=skip,
               fill = TRUE, quote = "\"", check.names=FALSE, as.is=TRUE, comment.char="", header=TRUE)
  read.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.table")))
  dat <- do.call("read.table", read.args)

  ## Gnames
  descript <- new("marrayInfo", maNotes = notes)
  if (is.null(info.id))
    info.id <- 1:ncol(dat)
  maInfo(descript) <- data.frame(dat[,info.id])

  if (length(labels) == nrow(dat))
    maLabels(descript) <- as.vector(labels)
  else {
    if (is.null(labels))
      labels <- 1
    maLabels(descript) <- as.character(dat[,labels])
  }
  ## Layout
##  headerinfo <- readGPRHeaders(f)
  idB <- grep("Block", colnames(dat));  ##Lblock <- dat[,id]
  idR<- grep("Row", colnames(dat));  ##Lrow <- dat[,id]
  idC <- grep("Column", colnames(dat)); ## Lcolumn <- dat[,id]
##  if(is.null(headerinfo))
##    ncolumns <- 
  mlayout <- maCompLayout(dat[,c(idB, idR, idC)], ncolumns)

  ##  ngr <- max(Lblock) / ncolumns
##   ngc <- ncolumns
##   nsr <- max(Lrow)
##   nsc <- max(Lcolumn)
##  nspots <- as.integer(ngr) * as.integer(ngc) * as.integer(nsr) * as.integer(nsc)
##  temp <- rep(FALSE, nspots)
##  ind <- (nsr * nsc) * (Lblock - 1) + (Lrow - 1) * nsc + Lcolumn
##  temp[ind] <- TRUE
##  mlayout <- new("marrayLayout", maNgr = as.integer(ngr),
##                 maNgc = as.integer(ngc),
##                 maNsr = as.integer(nsr),
##                 maNsc = as.integer(nsc),
##                 maNspots = nspots,
##                 maSub=temp)
  return(list(gnames = descript, layout=mlayout))
}

read.SMD2 <- function(fnames = NULL, path = ".", name.Gf = "CH1I_MEAN", 
    name.Gb = "CH1B_MEDIAN", name.Rf = "CH2I_MEAN", name.Rb = "CH2B_MEDIAN", 
    name.W = NULL, layout = NULL, gnames = NULL, targets = NULL, notes = NULL, 
    skip = 0, sep = "\t", quote = "", ...) {

    if (is.null(fnames)) 
        fnames <- dir(path = path, pattern = paste("*", "xls", sep = "."))
    if (is.null(path)) 
        fullfnames <- fnames
    else
        fullfnames <- file.path(path, fnames)
    y <- readLines(fullfnames[1], n = 100)
    skip <- grep(name.Gf, y)[1] - 1
    smdTable <- NULL
    
    if (is.null(layout)) {
    
        cat("Generating layout from ", fnames[1], "\n", sep="")
        smdTable <- read.table(fullfnames[1], header=TRUE, sep="\t", 
                               quote = "", skip = skip, comment.char = "")
    
        numSectors <- max(smdTable$SECTOR)
        xsize <- max(smdTable$X.COORD) - min(smdTable$X.COORD)
        ysize <- max(smdTable$Y.COORD) - min(smdTable$Y.COORD)
        
        maNgr <- round(sqrt(numSectors*ysize/xsize))
        maNgc <- round(sqrt(numSectors*xsize/ysize))
        if (is.na(maNgr)) {
            row <- grep("Exptid", y)[1]
            exptid <- strsplit(y[row], "=")[[1]][2]
            cat("Image: http://genome-www5.stanford.edu/cgi-bin/SMD/clickable.pl?exptid=",
                exptid, "\n", sep = "")
            options(warn = getOption("warn")-1)
            repeat {
                cat("Enter number of vertical sectors (", numSectors, 
                    " total sectors): ", sep = "")
                maNgr <- as.integer(readLines(n = 1))
                if (!is.na(maNgr) && maNgr > 0 && maNgr < numSectors &&
                    numSectors/maNgr == as.integer(numSectors/maNgr))
                    break
            }
            options(warn = getOption("warn")+1)
            maNgc <- numSectors / maNgr
        }
        
        row <- grep("Rows per Sector", y)[1]
        maNsr <- as.integer(strsplit(y[row], "=")[[1]][2])
        row <- grep("Columns per Sector", y)[1]
        maNsc <- as.integer(strsplit(y[row], "=")[[1]][2])
        
        maNspots <- maNgr * maNgc * maNsr * maNsc
        
        maSub <- rep(FALSE, maNspots)
        maSub[smdTable$SPOT] <- TRUE
        
        row <- grep("Printname", y)[1]
        printname <- strsplit(y[row], "=")[[1]][2]
        row <- grep("Tip Configuration", y)[1]
        tipconfig <- strsplit(y[row], "=")[[1]][2]
        maNotes <- paste("Print Name: ", printname, 
                         "\nTip Configuration: ", tipconfig, sep = "")
        
        layout <- new("marrayLayout", maNgr = maNgr, maNgc = maNgc, 
                      maNsr = maNsr, maNsc = maNsc, maNspots = maNspots, 
                      maSub = maSub, maNotes = maNotes)
    }
    
    if (is.null(gnames)) {
    
        cat("Generating probe sequence info from ", fnames[1], "\n", sep="")
        if (is.null(smdTable))
            smdTable <- read.table(fullfnames[1], header=TRUE, sep="\t", 
                                   quote = "", skip = skip, comment.char = "")
        maLabels <- as.character(smdTable$SUID)
        cols <- 2:(match("SUID", colnames(smdTable))-1)
        maInfo <- smdTable[,cols]
        
        gnames <- new("marrayInfo", maLabels = maLabels, maInfo = maInfo)
    }
    
    if (is.null(targets)) {
    
        cat("Generating target sample info from all files\n")
        maLabels <- character(0)
        maInfo <- data.frame()
        for (i in 1:length(fnames)) {
             z <- readLines(fullfnames[i], n = skip)
             row <- grep("Exptid", z)[1]
             maLabels <- c(maLabels, strsplit(z[row], "=")[[1]][2])
             row <- grep("Experiment Name", z)[1]
             Experiment <- strsplit(z[row], "=")[[1]][2]
             row <- grep("Channel 1 Description", z)[1]
             Cy3 <- strsplit(z[row], "=")[[1]][2]
             row <- grep("Channel 2 Description", z)[1]
             Cy5 <- strsplit(z[row], "=")[[1]][2]
             row <- grep("SlideName", z)[1]
             SlideName <- strsplit(z[row], "=")[[1]][2]
             maInfo <- rbind(maInfo, data.frame(Experiment = Experiment,
                                                Cy3 = Cy3, Cy5 = Cy5, 
                                                SlideName = SlideName))
        }
        rownames(maInfo) <- 1:dim(maInfo)[1]
        
        targets <- new("marrayInfo", maLabels = maLabels, maInfo = maInfo)
    }
    
    if (is.null(notes)) {
        
        cat("Generating notes from ", fnames[1], "\n", sep="")
        row <- grep("Organism", y)[1]
        organism <- strsplit(y[row], "=")[[1]][2]
        row <- grep("Category", y)[1]
        category <- strsplit(y[row], "=")[[1]][2]
        row <- grep("Subcategory", y)[1]
        subcategory <- strsplit(y[row], "=")[[1]][2]
        row <- grep("Description", y)[1]
        description <- strsplit(y[row], "=")[[1]][2]
        row <- grep("Experimenter", y)[1]
        experimenter <- strsplit(y[row], "=")[[1]][2]
        row <- grep("Contact email", y)[1]
        email <- strsplit(y[row], "=")[[1]][2]
        row <- grep("Scanning Software", y)[1]
        software <- strsplit(y[row], "=")[[1]][2]
        row <- grep("Software version", y)[1]
        version <- strsplit(y[row], "=")[[1]][2]
        row <- grep("Scanning parameters", y)[1]
        parameters <- strsplit(y[row], "=")[[1]]
        if (length(parameters) > 1)
            parameters <- paste(parameters[2:length(parameters)], collapse = ", ")
        else
            parameters <- NA
        
        notes <- paste("Organism: ", organism, 
                       "\nCategory: ", category, 
                       "\nSubcategory: ", subcategory, 
                       "\nDescription: ", description, 
                       "\nExperimenter: ", experimenter, 
                       "\nE-Mail: ", email,
                       "\nScanning Software: ", software, " ", version,
                       "\nScanning Parameters: ", parameters, sep = "")
    }
    
    mraw <- read.marrayRaw(fnames = fnames, path = path, name.Gf = name.Gf, 
        name.Gb = name.Gb, name.Rf = name.Rf, name.Rb = name.Rb, 
        name.W = name.W, layout = layout, gnames = gnames, targets = targets, 
        notes = notes, skip = skip, sep = sep, quote = quote, 
        ...)
    
    return(mraw)
}
