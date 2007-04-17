###########################################################################
##  Function: read.marrayLayout
##
##  source("~/Projects/madman/Rpacks/marray/R/maInput.R")
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
           skip=NULL,
           sep="\t",
           quote="\"",
           DEBUG=FALSE,
           ...)
{
  if(DEBUG) print("in read.marrayraw")
  if(DEBUG) cat("Path", path, "\n")
  if(is.null(path))
    fullfnames <- fnames
  else
    fullfnames <- file.path(path, fnames)

  fname <- fullfnames[1]

  if(DEBUG) print(skip)

  # Intensity data
  if(is.null(skip))
    {
      if (DEBUG) print("In is.null(skip) ")
      y <- readLines(fullfnames[1], n=100)
      #modify name.Gf for special character
      regName.Gf <- gsub("\\(", "\\\\\\(", name.Gf)
      regName.Gf <- gsub("\\)", "\\\\\\)", regName.Gf)
      skip2 <- grep(regName.Gf, y)[1] - 1
      if(DEBUG) cat(skip2, "done \n")
    }
  else {
    if(DEBUG) print("skip != NULL")
    skip2 <- skip
    if(DEBUG) print(paste("skip2 =", skip2, sep=" "))
  }

  if(is.null(layout))
    {
      if (DEBUG) print("in is.null(layout)")
      nspots <- length(readLines(fullfnames[1])) - skip2 - 1
      if(DEBUG) cat(nspots, "of rows \n")
    }
  else
    {
      if(sum(layout@maSub) == 1)
        nspots  <- layout@maNgr * layout@maNgc * layout@maNsr * layout@maNsc
      else
        nspots <- sum(layout@maSub)
    }
  ##     nspots <- sum(layout@maSub)  ## Aug 3, maSub is default set to 1.

  if (DEBUG) print(paste("nspots = ", nspots, sep=""))

  Y <- matrix(0, nspots, length(fullfnames))
  colnames(Y) <- fullfnames
  if(!is.null(name.Rf)) Rf <- Y
  if(!is.null(name.Gf)) Gf <- Y
  if(!is.null(name.Gb)) Gb <- Y
  if(!is.null(name.Rb)) Rb <- Y
  if(!is.null(name.W)) W <- Y

  for(f in fullfnames)
  {
    cat("Reading ... ", f, "\n")

    ## Calculate Skip
    if(is.null(skip))
      {
        if(DEBUG) print("in is.null(skip), part 2")
        if(DEBUG) cat("Calculating skip  ... ")
        y <- readLines(f, n=100)
        #modify name.Gf for special character
        regName.Gf <- gsub("\\(", "\\\\\\(", name.Gf)
        regName.Gf <- gsub("\\)", "\\\\\\)", regName.Gf)

        skip2 <- grep(regName.Gf, y)[1] - 1
        if(DEBUG) cat(skip2, "done \n")
      }
    else
      {
        if (DEBUG) print("in !is.null(skip) 2")
        skip2 <- skip
      }

    dat <- read.table(f, skip = skip2, header = TRUE,
                      sep = sep, quote = quote, check.names = FALSE,
                      as.is = TRUE, comment.char = "", nrows = nspots, ...)

    if(!is.null(name.Gf)) Gf[,f]<- as.matrix(dat[, name.Gf])
    if(!is.null(name.Gb)) Gb[,f]<- as.matrix(dat[, name.Gb])
    if(!is.null(name.Rf)) Rf[,f]<-as.matrix(dat[, name.Rf])
    if(!is.null(name.Rb)) Rb[,f]<-as.matrix(dat[, name.Rb])
    if(!is.null(name.W)) W[,f] <-as.matrix(dat[, name.W])
  }

  ## Add Notes

  if(is.null(notes)) notes <- ""

  mraw <- new("marrayRaw", maNotes=notes)
  if(!is.null(name.Gf)) mraw@maGf <- Gf
  if(!is.null(name.Gb)) mraw@maGb <- Gb
  if(!is.null(name.Rf)) mraw@maRf <- Rf
  if(!is.null(name.Rb)) mraw@maRb <- Rb
  if(!is.null(layout)) maLayout(mraw) <- layout
  if(!is.null(gnames)) maGnames(mraw) <- gnames
  if(!is.null(targets)) maTargets(mraw) <- targets
  if(!is.null(name.W)) maW(mraw) <- W
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
                       skip=NULL,
                       sep="\t",
                       quote="\"",
                       DEBUG=FALSE,
                       ...)
  {
    ## If fnames not specified, read everything in the dir
    opt <- list(...)
    if(is.null(fnames))
      {
        if(!is.null(targets))
          fullnames <- as.vector(maInfo(targets)[,1]) ## modified
        else
          {
            if(is.null(path))
              fullnames <- dir(pattern="\\.spot$")
            else
              fullnames <- dir(path, pattern="\\.spot$")
          }
      }
    else
      fullnames <- fnames

    setgal <- FALSE
    if(DEBUG) print(fullnames)
    if(is.null(gnames) | is.null(layout))
      {
        if(DEBUG) print("Reading Galfile ... ")
        galfile <- dir(pattern=".\\gal", path=path)
        defs <- list(galfile = galfile, path=path, info.id = c("ID", "Name"),
                     labels = "ID", sep = sep, quote=quote, fill=TRUE, check.names=FALSE,
                     as.is=TRUE, ncolumns = 4)
        gal.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.Galfile")))
        gal <- do.call("read.Galfile", gal.args)
        if(is.null(gnames)) gnames <- gal$gnames
        if(is.null(layout)) layout <- gal$layout
        if(DEBUG) print("done \n ")
        setgal <- TRUE
      }

    if(DEBUG) print(table(layout@maControls))
    if(DEBUG) cat("Setting up controls status ... ")
    if(length(layout@maControls) == 0)
      {
        if(DEBUG) cat("Generate controls \n")
        GenControls.args <- maDotsMatch(maDotsMatch(opt, list(Gnames=gnames)), formals(args("maGenControls")))
        layout@maControls <- as.factor(do.call("maGenControls", GenControls.args))
        if(DEBUG) print(table(layout@maControls))
      }
    if(DEBUG) cat("done \n ")
    if(is.null(notes)) notes <- "Spot Data"

    if(DEBUG) cat("Calling read.marrayRaw ... \n")
    defs <- list(fnames = fullnames, path=path,
                 name.Gf = name.Gf, name.Gb=name.Gb, name.Rf=name.Rf, name.Rb=name.Rb,
                 name.W=name.W, layout = layout, gnames=gnames, targets=targets,
                 notes = notes, skip=skip, sep=sep, quote=quote, fill=TRUE)
    maRaw.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.marrayRaw")))
    maRaw.args <- c(maRaw.args, defs[names(defs[setdiff(names(defs), names(maRaw.args))])])
    mraw <- do.call("read.marrayRaw", maRaw.args)

    if(setgal)
      {
        ## Checking orders
        mraw <- mraw[gal$neworder,]
        ## make sure the reordering doesn't affect the subset fucntion.
        mraw@maLayout@maSub <- layout@maSub
      }
    return(mraw)
  }



###########################################################################
#  Read GENEPIX
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
                          skip=NULL,
                          sep="\t",
                          quote="\"",
                          DEBUG=FALSE,
                          ...)
  {

    opt <- list(...)
    ## If fnames not specified, read everything in the dir
    if(is.null(fnames))
      {
        if(!is.null(targets))
          fullnames <- as.vector(maInfo(targets)[,1]) ## modified
        else
          {
            if(is.null(path))
              fullnames <- dir(pattern=".*\\.gpr$")
            else
              fullnames <- dir(path, pattern=".*\\.gpr$")
          }
      }
    else
      fullnames <- fnames

    setgal <- FALSE
    if(DEBUG) print(fullnames)
    if(is.null(gnames) | is.null(layout))
      {
        if(DEBUG) print("Reading Galfile ... ")
        defs <- list(galfile = fullnames[1], path=path, info.id = c("ID", "Name"),
                     labels = "ID", sep = sep, quote=quote, fill=TRUE, check.names=FALSE,
                     as.is=TRUE, ncolumns = 4, skip=skip)
        gal.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.Galfile")))
        gal <- do.call("read.Galfile", gal.args)
        if(is.null(gnames)) gnames <- gal$gnames
        if(is.null(layout)) layout <- gal$layout
        if(DEBUG) print("done \n ")
        setgal <- TRUE
      }

    if(DEBUG) print(table(layout@maControls))
    if(DEBUG) cat("Setting up controls status ... ")
    if(length(layout@maControls) == 0)
      {
        if(DEBUG) cat("Generate controls \n")
        GenControls.args <- maDotsMatch(maDotsMatch(opt, list(Gnames=gnames)), formals(args("maGenControls")))
        layout@maControls <- as.factor(do.call("maGenControls", GenControls.args))
        if(DEBUG) print(table(layout@maControls))
      }
    if(DEBUG) cat("done \n ")

    if(is.null(notes)) notes <- "GenePix Data"

    if(DEBUG) cat("Calling read.marrayRaw ... \n")
    defs <- list(fnames = fullnames, path=path,
                 name.Gf = name.Gf, name.Gb=name.Gb, name.Rf=name.Rf, name.Rb=name.Rb,
                 name.W=name.W, layout = layout, gnames=gnames, targets=targets,
                 notes = notes, skip=skip, sep=sep, quote=quote, fill=TRUE)
    maRaw.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.marrayRaw")))
    maRaw.args <- c(maRaw.args, defs[names(defs[setdiff(names(defs), names(maRaw.args))])])
    mraw <- do.call("read.marrayRaw", maRaw.args)

    if(setgal)
      {
        ## Checking orders
        mraw <- mraw[gal$neworder,]
        ## make sure the reordering doesn't affect the subset fucntion.
        mraw@maLayout@maSub <- layout@maSub
      }

    return(mraw)
  }




###########################################################################
#  Read Agilent
#
read.Agilent <-  function(fnames = NULL,
                          path=NULL,
                          name.Gf = "gMedianSignal",
                          name.Gb = "gBGMedianSignal",
                          name.Rf = "rMedianSignal",
                          name.Rb = "rBGMedianSignal",
                          name.W= NULL,
                          layout = NULL,
                          gnames = NULL,
                          targets = NULL,
                          notes=NULL,
                          skip=NULL,
                          sep="\t",
                          quote="\"",
                          DEBUG=FALSE,
                          info.id=NULL,
                          ...)
  {

    opt <- list(...)
    ## If fnames not specified, read everything in the dir
    if(is.null(fnames))
      {
        if(!is.null(targets))
          fullnames <- as.vector(maInfo(targets)[,1]) ## modified
        else
          {
            if(is.null(path))
              fullnames <- dir(pattern=".*\\.txt$")
            else
              fullnames <- dir(path, pattern=".*\\.txt$")
          }
      }
    else
      fullnames <- fnames

    if(DEBUG) print(fullnames)
    if(is.null(gnames) | is.null(layout))
      {
        if(DEBUG) cat("Setting up Gene Annotation information  ... ")
        ## Gnames : reading data in (the first one)
        y <- readLines(fullnames[1], n=100)
        if(class(info.id) == "character")
          skip <- intersect(grep(info.id[1], y), grep("Row", y))[1] - 1
        else
          skip <- intersect(grep("Row", y), grep("ProbeName", y))[1] - 1

        defs <- list(file=fullnames[1], path = path, sep=sep, skip=skip,
                     fill = TRUE, quote = "\"", check.names=FALSE,
                     as.is=TRUE, comment.char="", header=TRUE)
        read.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.table")))
        dat <- do.call("read.table", read.args)

        descript <- new("marrayInfo", maNotes = "Agilent")
        if(is.null(info.id))
          {
            info.id <- intersect(c("ProbeName", "ProbeUID", "SystematicName", "Description", "ControlType", "SwissProt", "GenBank", "Primate", "Sequence", "GenPept"), colnames(dat))
          }
        maInfo(descript) <- data.frame(dat[,info.id])
        rownames(descript@maInfo) <- make.names(as.vector(dat[,info.id[1]]), unique=TRUE)
        descript@maLabels <- as.character(dat[,info.id[1]])
        if(DEBUG) cat("done \n")

        ## Layout
        if(DEBUG) cat("Setting up Layout information  ... ")
##        idR<- grep(layout.id["Row"], colnames(dat));  ##Lrow <- dat[,id]
##        idC <- grep(layout.id["Column"], colnames(dat)); ## Lcolumn <- dat[,id]
        layout <- maCompLayout(cbind(1, 1, dat[,c("Row", "Col")]), ncolumns=1)
        tmp <- rep(TRUE, max(as.integer(dat[,"FeatureNum"])))
        tmp[-dat[,"FeatureNum"]] <- FALSE
        maSub(layout) <- tmp
        if(DEBUG) cat("done \n")
      }

    if(DEBUG) cat("Setting up controls status ... ")
    if(length(layout@maControls) == 0) {
      tmp <- as.character(dat[,"ControlType"])
      tmp[tmp == "0"] <- "probes"
      tmp[tmp == "1"] <- "Positive"
      tmp[tmp == "-1"] <- "Negative"
      layout@maControls <- as.factor(tmp)
      rm(tmp)
    }
    if(DEBUG) cat("done \n ")

    if(is.null(notes)) notes <- "Agilent Data"

    if(DEBUG) cat("Calling read.marrayRaw ... \n")
    if(DEBUG) cat("To read", fullnames, "... \n")

    defs <- list(fnames = fullnames, path=path,
                 name.Gf = name.Gf, name.Gb=name.Gb, name.Rf=name.Rf, name.Rb=name.Rb,
                 name.W=name.W, layout = layout, gnames=descript, targets=targets,
                 notes = notes, skip=skip, sep=sep, quote=quote, fill=TRUE,
                 check.names=FALSE,  as.is=TRUE)
    maRaw.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.marrayRaw")))
    mraw <- do.call("read.marrayRaw", maRaw.args)

    return(mraw)
  }




###########################################################################

read.Galfile <- function (galfile,
                          path=".",
                          info.id = c("ID", "Name"),
                          layout.id =c(Block="Block", Row="Row", Column="Column"),
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

  ## Fix no check whether skip is null or not ##
  ## Pointed out by Dustin Potter,  Feb 9, 2006 ##

  if(is.null(skip)){
    if(class(info.id) == "character")
      skip <- intersect(grep(info.id[1], y), grep(layout.id[1], y))[1] - 1
    else
      skip <- intersect(grep("ID", y), grep("Name", y))[1] - 1
  }
  
  defs <- list(file=f, path = path, sep=sep, skip=skip,
               fill = TRUE, quote = "\"", check.names=FALSE,
               as.is=TRUE, comment.char="", header=TRUE)
  read.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.table")))
  dat <- do.call("read.table", read.args)

  ## Gnames
  descript <- new("marrayInfo", maNotes = notes)
  if (is.null(info.id))
    info.id <- 1:ncol(dat)
  maInfo(descript) <- data.frame(dat[,info.id])

  ## Fix for R2.0 because data.frame needed unique rownames
  rownames(descript@maInfo) <- make.names(as.vector(dat[,info.id[1]]), unique=TRUE)

  if (length(labels) == nrow(dat))
    maLabels(descript) <- as.vector(labels)
  else {
    if (is.null(labels))
      labels <- 1
    maLabels(descript) <- as.character(dat[,labels])
  }

  ## Layout
  idB <- grep(layout.id["Block"], colnames(dat));  ##Lblock <- dat[,id]
  idR<- grep(layout.id["Row"], colnames(dat));  ##Lrow <- dat[,id]
  idC <- grep(layout.id["Column"], colnames(dat)); ## Lcolumn <- dat[,id]
  mlayout <- maCompLayout(dat[,c(idB, idR, idC)], ncolumns)
  newmat <- cbind(gr = ((dat[,idB] - 1)%/%ncolumns) + 1,
                  gc = ((dat[,idB] - 1)%%ncolumns) + 1, sr = dat[,idR], sc = dat[,idC])
  ngr<-maNgr(mlayout); ngc<-maNgc(mlayout)
  nsr<-maNsr(mlayout); nsc<-maNsc(mlayout)
  ind<-(nsr * nsc)* ((newmat[,1] - 1) * ngc + (newmat[,2] - 1)) + (newmat[,3] - 1) * nsc + newmat[,4]
  ## ind provides the proper location of the entry
  ## Assume no NAs and no missing values
  return(list(gnames = descript, layout=mlayout, neworder=order(ind)))
}



##############################################################
## Read SMD, default arguments are for Hs genome
## For other genomes, you will need to modify the info.id field

read.SMD <- function(fnames = NULL, path = NULL,
                     name.Gf = "Ch1 Intensity (Median)",
                     name.Gb = "Ch1 Background (Median)",
                     name.Rf = "Ch2 Intensity (Median)",
                     name.Rb = "Ch2 Background (Median)",
                     name.W = NULL,
                     info.id = c("Name", "Clone ID"),
                     layout = NULL, gnames = NULL, targets = NULL, notes = NULL,
                     skip = NULL, sep = "\t", quote = "\"", DEBUG=FALSE, ...)
{
  ## If fnames is NULL check target, if target is also NULL then read from dir()
  if(is.null(fnames))
    {
      if(!is.null(targets))
        fullnames <- as.vector(maInfo(targets)[,1]) ## modified
      else
        {
          if(is.null(path))
            fullnames <- dir(pattern="\\.xls$")
          else
            fullnames <- dir(path, pattern="\\.xls$")
        }
    }
  else
    fullnames <- fnames

  if(DEBUG) print(fnames)

  ## Jean Yang, April 18, 2007
  ## Add opt = list(...)
  ## Bug pointed out by Michael Gormley <mpg33@drexel.edu>
  opt <- list(...)
  if(is.null(gnames) | is.null(layout))
      {
        if(DEBUG) cat("Generating layout from ", fullnames[1], "\n", sep="")
        opt <- list(...)
        defs <- list(galfile = fullnames[1], path=path, info.id = c("Name", "Clone ID"),
                     layout.id=c(Block="Sector", Row="X Grid Coordinate \\(within sector\\)",
                       Column="Y Grid Coordinate \\(within sector\\)"),
                     labels = "Spot", sep = sep, quote=quote, fill=TRUE, check.names=FALSE,
                     as.is=TRUE, ncolumns = 4)
        gal.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.Galfile")))
        if (DEBUG) print("Reading Gal file")
        gal <- do.call("read.Galfile", gal.args)
        if(is.null(gnames)) gnames <- gal$gnames
        if(is.null(layout)) layout <- gal$layout
        if(DEBUG) print("done \n ")
      }

    if(is.null(targets))
    {
      cat("Generating target sample info from all files\n")
      maLabels <- character(0)
      maInfo <- data.frame()
      for (i in 1:length(fullnames)) {
        #modify name.Gf for special character
        regName.Gf <- gsub("\\(", "\\\\\\(", name.Gf)
        regName.Gf <- gsub("\\)", "\\\\\\)", regName.Gf)
        if (DEBUG) print(regName.Gf)
        comment <- grep(regName.Gf, readLines(fullnames[i], n=100))
        if (DEBUG) print(paste("comment = ", comment, sep=""))
        z <- readLines(fullnames[i], n = comment)
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
    }

  rownames(maInfo) <- 1:dim(maInfo)[1]
  targets <- new("marrayInfo", maLabels = maLabels, maInfo = maInfo)

  defs <- list(fnames = fullnames, path=path,
               name.Gf = name.Gf, name.Gb=name.Gb, name.Rf=name.Rf, name.Rb=name.Rb,
               name.W=name.W, layout = layout, gnames=gnames, targets=targets,
               notes = notes, skip=skip, sep=sep, quote=quote, fill=TRUE,
               check.names=FALSE,  as.is=TRUE)
  if(DEBUG) cat("Calling read.marrayRaw ... \n")
  maRaw.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.marrayRaw")))
  mraw <- do.call("read.marrayRaw", maRaw.args)
  return(mraw)
}

###########################################################################
# Function: checkTargetInfo
# ADD:  check that the foreground and backgruond intensities are
# stored in the same order as provided in the first column of target file.
# Date:  Jan 05, 2005
###########################################################################
checkTargetInfo <- function(mraw)
  {
    if(length(mraw@maTargets@maInfo) == 0)
      stop("Missing Target Information \n")
    targetFnames <- as.vector(mraw@maTargets@maInfo[,1])

    if(class(mraw) == "marrayRaw")
      {
        if((length(mraw@maGf) == 0 ) & length(mraw@maRf) == 0)
          stop("Missing intensities information in both channels\n")
        if((length(mraw@maGf) != 0 ))
          colFnames <- colnames(mraw@maGf)
        else
          colFnames <- colnames(mraw@maRf)
      }
    if(class(mraw) == "marrayNorm")
      {
        if((length(mraw@maM) == 0 ))
          stop("Missing Log-ratio information \n")
        colFnames <- colnames(mraw@maM)
      }

    res <- sapply(targetFnames, grep, colFnames) == 1:length(colFnames)
    return(sum(res) == length(colFnames))

  }

############################################
## END OF FILE
############################################
