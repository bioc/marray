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
  if(is.null(path))
    fullfnames <- fnames
  else
    fullfnames <- file.path(path, fnames)

  fname <- fullfnames[1]

  
  # Intensity data
  if(is.null(skip))
    {
      y <- readLines(fullfnames[1], n=100)
      skip <- grep(name.Gf, y)[1] - 1
      if(DEBUG) cat(skip, "done \n")
    }
    
  if(is.null(layout))
    {
      nspots <- length(readLines(fullfnames[1])) - skip - 1
      if(DEBUG) cat(nspots, "of rows \n")
    }
  else
    nspots <- maNspots(layout)
  
  Y <- matrix(0, nspots, length(fullfnames))
  colnames(Y) <- fullfnames
  if(!is.null(name.Rf)) Rf <- Y
  if(!is.null(name.Gf)) Gf <- Y
  if(!is.null(name.Gb)) Gb <- Y
  if(!is.null(name.Rb)) Rb <- Y
  if(!is.null(name.W)) W <- Y
  
  for(f in fullfnames)
  {
    cat("Reading ... ", f)

    ## Calculate Skip
    if(is.null(skip))
      {
        if(DEBUG) cat("Calculating skip  ... ")
        y <- readLines(f, n=100)
        skip <- grep(name.Gf, y)[1] - 1
        if(DEBUG) cat(skip, "done \n")
      }
    dat <- read.table(f, skip = skip, header = TRUE, 
                      sep = sep, as.is = TRUE, quote = quote, check.names = FALSE, 
                      comment.char = "", nrows = nspots, ...)

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
                       ...)
  {
    ## If fnames not specified, read everything in the dir

    if(is.null(fnames))
      {
        if(!is.null(targets))
          fullnames <- maInfo(targets)[,1]  
        else
          {
            if(is.null(path))
              fullnames <- dir(pattern="*\\.spot$")
            else
              fullnames <- dir(path, pattern="*\\.spot$")
          }
      }
    else
      fullnames <- fnames

    if(is.null(notes)) notes <- "Spot Data"

    mraw <- read.marrayRaw(fnames =fullnames,
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
          fullnames <- maInfo(targets)[,1]  
        else
          {
            if(is.null(path))
              fullnames <- dir(pattern="*\\.gpr$")
            else
              fullnames <- dir(path, pattern="*\\.gpr$")
          }
      }
    else
      fullnames <- fnames
    
    if(DEBUG) print(fullnames)
    if(is.null(gnames) | is.null(layout))
      {
        if(DEBUG) print("Reading Galfile ... ")
        opt <- list(...)
        defs <- list(galfile = fullnames[1], path=path, info.id = c("Name", "ID"),
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
    defs <- list(fnames = fullnames, path=path,
                 name.Gf = name.Gf, name.Gb=name.Gb, name.Rf=name.Rf, name.Rb=name.Rb,
                 name.W=name.W, layout = layout, gnames=gnames, targets=targets,
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

  skip <- intersect(grep(info.id[1], y), grep(layout.id[1], y))[1] - 1
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
  rownames(descript@maInfo) <- as.vector(dat[,info.id[1]])
  
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
  return(list(gnames = descript, layout=mlayout))
}

read.SMD <- function(fnames = NULL, path = NULL, name.Gf = "CH1I_MEDIAN", 
                     name.Gb = "CH1B_MEDIAN", name.Rf = "CH2I_MEDIAN", name.Rb = "CH2B_MEDIAN", 
                     name.W = NULL, layout = NULL, gnames = NULL, targets = NULL, notes = NULL, 
                     skip = NULL, sep = "\t", quote = "\"", DEBUG=FALSE, ...)
{
  ## If fnames is NULL check target, if target is also NULL then read from dir()
  if(is.null(fnames))
    {
      if(!is.null(targets))
        fullnames <- maInfo(targets)[,1]  
      else
        {
          if(is.null(path))
            fullnames <- dir(pattern="*\\.xls$")
          else
            fullnames <- dir(path, pattern="*\\.xls$")
        }
    }
  else
    fullnames <- fnames
  
  if(DEBUG) print(fnames)
  if(is.null(gnames) | is.null(layout))
      {
        if(DEBUG) cat("Generating layout from ", fullnames[1], "\n", sep="")
        opt <- list(...)
        defs <- list(galfile = fullnames[1], path=path, info.id = c("NAME", "Clone ID"),
                     layout.id=c(Block="SECTOR", Row="SECTORROW", Column="SECTORCOL"),
                     labels = "SPOT", sep = sep, quote=quote, fill=TRUE, check.names=FALSE,
                     as.is=TRUE, ncolumns = 4)
        gal.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.Galfile")))        
        gal <- do.call("read.Galfile", gal.args)
        if(is.null(gnames)) gnames <- gal$gnames
        if(is.null(layout)) layout <- gal$layout
        if(DEBUG) print("done \n ")
      }

  if(is.null(targets))
    {
      skip <- grep(name.Gf, readLines(fullnames[1], n=100))
      cat("Generating target sample info from all files\n")
      maLabels <- character(0)
      maInfo <- data.frame()
      for (i in 1:length(fullnames)) {
        z <- readLines(fullnames[i], n = skip)
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

#    if (is.null(notes)) {
#      cat("Generating notes from ", fnames[1], "\n", sep="")
#      row <- grep("Organism", y)[1]
#      organism <- strsplit(y[row], "=")[[1]][2]
#      row <- grep("Category", y)[1]
#      category <- strsplit(y[row], "=")[[1]][2]
#      row <- grep("Subcategory", y)[1]
#      subcategory <- strsplit(y[row], "=")[[1]][2]
#      row <- grep("Description", y)[1]
#      description <- strsplit(y[row], "=")[[1]][2]
#      row <- grep("Experimenter", y)[1]
#      experimenter <- strsplit(y[row], "=")[[1]][2]
#      row <- grep("Contact email", y)[1]
#      email <- strsplit(y[row], "=")[[1]][2]
#      row <- grep("Scanning Software", y)[1]
#      software <- strsplit(y[row], "=")[[1]][2]
#      row <- grep("Software version", y)[1]
#      version <- strsplit(y[row], "=")[[1]][2]
#      row <- grep("Scanning parameters", y)[1]
#      parameters <- strsplit(y[row], "=")[[1]]
#      if (length(parameters) > 1)
#        parameters <- paste(parameters[2:length(parameters)], collapse = ", ")
#      else
#        parameters <- NA
#      
#      notes <- paste("Organism: ", organism, 
#                     "\nCategory: ", category, 
#                     "\nSubcategory: ", subcategory, 
#                     "\nDescription: ", description, 
#                     "\nExperimenter: ", experimenter, 
#                     "\nE-Mail: ", email,
#                     "\nScanning Software: ", software, " ", version,
#                     "\nScanning Parameters: ", parameters, sep = "")
#    }

  if(DEBUG) cat("Calling read.marrayRaw ... \n")
  defs <- list(fnames = fullnames, path=path,
               name.Gf = name.Gf, name.Gb=name.Gb, name.Rf=name.Rf, name.Rb=name.Rb,
               name.W=name.W, layout = layout, gnames=gnames, targets=targets,
               notes = notes, skip=skip, sep=sep, quote=quote, fill=TRUE,
               check.names=FALSE,  as.is=TRUE)
  maRaw.args <- maDotsMatch(maDotsMatch(opt, defs), formals(args("read.marrayRaw")))        
  mraw <- do.call("read.marrayRaw", maRaw.args)
  return(mraw)
}

