###########################################################################
## WIDGET for marrayLayout
## setwd("C:\\MyDoc\\Projects\\BioC\\marrayInput\\R")
###########################################################################

widget.marrayLayout<-function(path="",
                              skip=0,
                              sep="\t",
                              quote="",
                              ...)
{
  require(tcltk)
  require(tkWidgets)
  
  ## Functions:
  ## Calling read.marrayLayout
  inputLayout <- function(name.plate,
                          name.controls,
                          skip=skip,
                          sep=sep,
                          quote=quote,
                          ...)
    {
      fname <- tclvalue("galfname")

      ## Name of the new marrayLayout objects in R
      if(tclvalue("layoutName")=="")
        newname <- paste(unlist(strsplit(tclvalue("galfname"), "\\."))[1],
                         "Layout", sep="")
      else
        newname <- tclvalue("layoutName")

      ## Values of notes
      if(tclvalue("notes") =="")
        notes <- fname
      else
        notes <- tclvalue("notes")

      ## automatically look for skip for gal file.
      if(fname != ""){
        x <- unlist(strsplit(fname, "\\."))
        xx <- x[length(x)]
        if(xx == "gal"){
          y <- readLines(fname, n=100)
          skip <- grep("Row", y)[1] - 1
        }
      }
      else
        fname <- NULL

      ## Calling read.marrayLayout
      newLayout <- read.marrayLayout(fname=fname,
                                     ngr=as.numeric(tclvalue("gr")),
                                     ngc=as.numeric(tclvalue("gc")),
                                     nsr=as.numeric(tclvalue("sr")),
                                     nsc=as.numeric(tclvalue("sc")),
                                     pl.col=name.plate,
                                     ctl.col=name.controls,
                                     notes=notes,
                                     skip=skip,
                                     sep=sep,
                                     quote=quote,
                                     ...
                                     )
      ## Assign name to newname)
      assign(newname, newLayout, env = .GlobalEnv)
      cat(paste("\n Finish creating a new marrayLayout: ", newname, "\n", sep=""))
    } ## end inputLayout


  ## List column names from file
  listColnames <- function(fname,
                           skip=skip,
                           sep=sep,
                           quote=quote,
                           ...)
    {
      x <- unlist(strsplit(fname, "\\."))
      xx <- x[length(x)]
      if(xx == "gal"){
        y <- readLines(fname, n=100)
        skip <- grep("Row", y)[1] - 1
      }

      h<-strsplit(readLines(fname, n=skip+1),split=sep)
      h<-as.list(unlist(h[[length(h)]]))
      columnHeadings<-gsub("\"", "", unlist(h))

      base <- tktoplevel()
      tkwm.title(base, "Select Column Information")
      mainfrm <- tkframe(base, borderwidth=2)

      ## Plate
      filevarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
      tkpack(filevarfr, tklabel(filevarfr, text="Plate"), side='top')
      for(ch in columnHeadings)
        tkpack(filevarfr, tkcheckbutton(filevarfr, text=ch,
                                        variable=paste('pl', ch, sep=".")))

      ## Controls
      file2varfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
      tkpack(file2varfr, tklabel(file2varfr, text="Controls"), side='top')
      for(ch in columnHeadings)
      tkpack(filevarfr, tkcheckbutton(file2varfr, text=ch,
                                      variable=paste('ct', ch, sep=".")))

      tkpack(mainfrm, filevarfr, side="left", fill="both")
      tkpack(mainfrm, filevarfr, side="right", fill="both")

      getResults <- function()
        {
          for(ch in columnHeadings){
            check <- eval(parse(text=tclvalue(paste("pl", ch, sep="."))))
            if(check == "1") name.plate <<- ch
            check <- eval(parse(text=tclvalue(paste("ct", ch, sep="."))))
            if(check == "1") name.controls <<- ch
            tkdestroy(base)
          }
        }
      butfrm <- tkframe(mainfrm, borderwidth= 1, relief="groove")
      b.but <- tkbutton(butfrm, command=getResults, text="Go Back")
      tkpack(butfrm, b.but, side="bottom")
      tkpack(mainfrm, butfrm, side="bottom")
    }

  # Jianhua added this function in
  fileBrow <- function()
    {
      tkdelete(fileEntry, 0, "end")
      temp <- fileBrowser(path)
      temp <- paste(temp, sep = "", collapse = "," )
      tkinsert(fileEntry, 0, temp)
    }

  newname <- newLayout <- name.plate <- name.controls <- NULL

  base <- tktoplevel()
  tkwm.title(base, "marrayLayout Builder")
  mainfrm <- tkframe(base, borderwidth=2)

  # Input Notes info
  notesvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  tkpack(notesvarfr, tklabel(notesvarfr,
                             text="Name of marrayLayout object:"), side='top')
  tkpack(notesvarfr, tkentry(notesvarfr,
                             width=30, textvariable="layoutName"), side='left')
  tkpack(mainfrm, notesvarfr, side="top")

  # Select spot coord
  coordvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  heading <- tklabel(coordvarfr, text="Spot coordinates")

  gr.label <- tklabel(coordvarfr, text="Grid row")
  gc.label <- tklabel(coordvarfr, text="Grid col")
  sr.label <- tklabel(coordvarfr, text="Spot row")
  sc.label <- tklabel(coordvarfr, text="Spot col")

  gr.entry <- tkentry(coordvarfr, width=7)
  gc.entry <- tkentry(coordvarfr, width=7)
  sr.entry <- tkentry(coordvarfr, width=7)
  sc.entry <- tkentry(coordvarfr, width=7)

  tkgrid(heading, columnspan = 4)
  tkgrid(gr.label, gr.entry, gc.label, gc.entry)
  tkgrid(sr.label, sr.entry, sc.label, sc.entry)
  tkgrid.configure(gr.label, gr.entry, gc.label, gc.entry)
  tkgrid.configure(sr.label, sr.entry, sc.label, sc.entry)

  tkconfigure(gr.entry, textvariable = "gr")
  tkconfigure(gc.entry, textvariable = "gc")
  tkconfigure(sr.entry, textvariable = "sr")
  tkconfigure(sc.entry, textvariable = "sc")

  tkpack(mainfrm, coordvarfr, side="top")

  # Select .gal file
  galvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  fileFrame <- tkframe(galvarfr)
  fileLabel <- tklabel(fileFrame, text="Info file:")
  fileEntry <- tkentry(fileFrame, width=30, textvariable="galfname")
  fileBrow <- tkbutton(fileFrame, text = "Browse", command = fileBrow)
  tkgrid(fileLabel, columnspan = 2)
  tkgrid(fileEntry, fileBrow)
  tkpack(fileFrame, side = "top")
  tkpack(mainfrm, galvarfr, side="top")

  # Input Notes info
  notesvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  tkpack(notesvarfr, tklabel(notesvarfr, text="Notes:"), side='top')
  tkpack(notesvarfr, tkentry(notesvarfr, width=30, textvariable="notes"), side='left')
  tkpack(mainfrm, notesvarfr, side="top")

  # Input print-tip and controls column
  butfrm <- tkframe(mainfrm, borderwidth= 1, relief="groove")
  a.but <- tkbutton(butfrm, command=function()
                    {
                      if(length(tclvalue("galfname")) != 0)
                        listColnames(tclvalue("galfname"),
                                     skip=skip, sep=sep, quote=quote, ...)
                    },
                    text="GetInfo")
  b.but <- tkbutton(butfrm, command=function()tkdestroy(base), text="Quit")
  c.but <- tkbutton(butfrm, command=function()
                    inputLayout(name.plate, name.controls,
                                skip=skip, sep=sep, quote=quote, ...),
                    text="Build")
  tkpack(butfrm, a.but, c.but, side="left")
  tkpack(butfrm, b.but, side="right")
  tkpack(mainfrm, butfrm,anchor='w')
  return(invisible())
}


###########################################################################
## WIDGET : marrayInfo
###########################################################################
widget.marrayInfo <- function(path=".",
                              skip=0,
                              sep="\t",
                              quote="",
                              ...)
{
  require(tcltk)
  require(tkWidgets)
   
  ## Functions:
  ## Calling read.marrayInfo
  inputInfo <- function(skip=skip,
                        sep=sep,
                        quote=quote,
                        ...)
    {
      fname <- tclvalue("Fname")


      ## Name of the new marrayInfo objects in R
      if(tclvalue("infoName")=="")
        newname <- paste(unlist(strsplit(tclvalue("Fname"), "\\."))[1], "Info", sep="")
      else
        newname <- tclvalue("infoName")

      ## Name of the new marrayNotes objects in R
      if(tclvalue("notes") =="")
        notes <- fname
      else
        notes <- tclvalue("notes")

      ## Setting some defaults for GenePix Gal file
      x <- unlist(strsplit(fname, "\\."))
      xx <- x[length(x)]
      if(xx == "gal"){
        y <- readLines(fname, n=100)
        skip <- grep("Row", y)[1] - 1
      }

      info.id <- eval(parse(text=tclvalue("infoID")))
      labels.id <- eval(parse(text=tclvalue("labels")))
      newInfo <- read.marrayInfo(fname=fname,
                                 info.id = info.id,
                                 labels = labels.id,
                                 notes=notes,
                                 skip=skip,
                                 sep=sep,
                                 quote=quote,
                                 ...
                                 )
      cat(paste("\n Finish creating a new marrayInfo: ", newname, "\n", sep=""))
      assign(newname, newInfo, env = .GlobalEnv)
    }

  # Jianhua put this function in
  fileBrow <- function()
    {
      tkdelete(fileEntry, 0, "end")
      temp <- fileBrowser(path)
      temp <- paste(temp, sep = "", collapse = ",")
      tkinsert(fileEntry, 0, temp)
    }

  newName <- NULL
  base <- tktoplevel()
  tkwm.title(base, "MarrayInfo Builder")
  mainfrm <- tkframe(base, borderwidth=2)

  # Input Names info
  notesvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  tkpack(notesvarfr, tklabel(notesvarfr, text="Name of marrayInfo object:"), side='top')
  tkpack(notesvarfr, tkentry(notesvarfr, width=30, textvariable="infoName"), side='left')
  tkpack(mainfrm, notesvarfr, side="top")

  # Select Fname file
  galvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  fileFrame <- tkframe(galvarfr)
  fileLbl <- tklabel(fileFrame, text="File Name:")
  fileEntry <- tkentry(fileFrame, width=30, textvariable="Fname")
  fileBrow <- tkbutton(fileFrame, text = "Browse",
                       command = fileBrow)
  tkgrid(fileLbl, columnspan = 2)
  tkgrid(fileEntry, fileBrow)
  tkpack(fileFrame, side = "top")
  tkpack(mainfrm, galvarfr, side="top")

  # Input ID info
  idvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  tkgrid(idvarfr, columnspan = 2)
  info.lab <- tklabel(idvarfr, text="Info Index")
  info.entry <- tkentry(idvarfr, width=15, textvariable="infoID")
  labels.lab <-  tklabel(idvarfr, text="Labels")
  labels.entry <-  tkentry(idvarfr, width=15, textvariable="labels")
  tkgrid(info.lab, labels.lab)
  tkgrid(info.entry, labels.entry)
  tkgrid.configure(info.lab, labels.lab)
  tkgrid.configure(info.entry, labels.entry)
  tkpack(mainfrm, idvarfr, side="top")

  # Input Notes info
  notesvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  tkpack(notesvarfr, tklabel(notesvarfr, text="Notes:"), side='top')
  tkpack(notesvarfr, tkentry(notesvarfr, width=30, textvariable="notes"), side='left')
  tkpack(mainfrm, notesvarfr, side="top")

  # Input print-tip and controls column
  butfrm <- tkframe(mainfrm, borderwidth= 1, relief="groove")
  a.but <- tkbutton(butfrm, command=function()
                    inputInfo(sep=sep,
                              skip=skip,
                              quote=quote,
                              ...), text="Build")
  b.but <- tkbutton(butfrm, command=function()tkdestroy(base), text="End")

  tkpack(butfrm, a.but, side="left")
  tkpack(butfrm, b.but, side="right")
  tkpack(mainfrm, butfrm,anchor='w')
  return(invisible())
}

###########################################################################
## WIDGET for marrayRaw
##
###########################################################################

widget.marrayRaw<-function(ext=c("spot", "xls", "gpr"),
                           skip=0,
                           sep="\t",
                           quote="\"",
                           ...)

{
  require(tcltk)
  require(tkWidgets)
  
  ## Jianhua added this in to catch the "..."
  args <- list(...)

  ## Calling read.marrayLayout
  inputRaw <- function(fnames,
                       skip=skip,
                       sep=sep,
                       quote=quote,
                       ...)
    {
      if(tclvalue("rawName")=="")
        newname <- paste(ext, "data", sep="")
      else
        newname <- tclvalue("rawName")

      if(tclvalue("notes") =="")
        notes <- paste("File from", ext)
      else
        notes <- tclvalue("notes")

      ## Deal with Names
      if(tclvalue("nameW") == "") name.W <- NULL else name.W <- tclvalue("nameW")
      if(tclvalue("nameRf") == "")
        stop("Input Foreground intensities") else name.Rf <- tclvalue("nameRf")
      if(tclvalue("nameGf") == "")
        stop("Input Foreground intensities") else name.Gf <- tclvalue("nameGf")
      if(tclvalue("nameRb") == "") name.Rb <- NULL else name.Rb <- tclvalue("nameRb")
      if(tclvalue("nameGb") == "") name.Gb <- NULL else name.Gb <- tclvalue("nameGb")

      ## Calculate Skip
      y <- readLines(fnames[1], n=100)
      skip <- grep(name.Gf, y)[1] - 1

      ## Call read.marrayRaw
      newRaw <- read.marrayRaw(fnames=fnames,
                               path=NULL,
                               name.Gf=name.Gf,
                               name.Gb=name.Gb,
                               name.Rf=name.Rf,
                               name.Rb=name.Rb,
                               name.W=name.W,
                               layout=eval(parse(text=tclvalue("layout"))),
                               gnames=eval(parse(text=tclvalue("genenames"))),
                               targets=eval(parse(text=tclvalue("targets"))),
                               notes=tclvalue("notes"),
                               skip=skip,
                               sep=sep,
                               quote=quote,
                                ...)

      cat(paste("\n Finish creating a new marrayInfo: ", newname, "\n", sep=""))
      assign(newname, newRaw, env = .GlobalEnv)
    }

  ### Jianhua
  # This is the name I put in to store the files when "Files" is
  # clicked. Change it to whatever you would like to
  # Jean: Name it "fnames"
  fnames <- NULL

  # Jianhua added the following functions in
  filesClick <- function(){
      if(!is.null(args$path))
          fnames <<- fileBrowser(args$path)
      else
          fnames <<- fileBrowser()
  }

  layBrowser <- function(){
      browserToEntry(layoutEntry)
  }

  samBrowser <- function(){
      browserToEntry(targetEntry)
  }

  genBrowser <- function(){
      browserToEntry(geneEntry)
  }

  browserToEntry <- function(anEntry){
      first = TRUE
      temp <- objectBrowser()
      tkdelete(anEntry, 0, "end")
      for(i in 1:length(temp)){
          if(first){
              tkinsert(anEntry, "end", temp[[i]]$name)
              first <- FALSE
          }else{
              toPut <- paste(",", temp[[i]]$name, sep = "")
              tkinsert(anEntry, "end", toPut)
          }
      }
  }

  layoutPressed <- function(){
      if(!is.null(args$path))
          widget.marrayLayout(path = args$path)
      else
          widget.marrayLayout(path="")
  }

  genePressed <- function(){
      if(!is.null(args$path))
          widget.marrayInfo(path = args$path)
      else
          widget.marrayInfo(path="")
  }
  # end of Jianhua's functions

  newName <- NULL

  base <- tktoplevel()
  tkwm.title(base, "MarrayRaw builder")
  mainfrm <- tkframe(base, borderwidth=2)

  # Select files
  filevarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  tkpack(filevarfr, tkbutton(filevarfr, text="Files",
                             command = filesClick), side='top')

  # Input Names
  notesvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  tkpack(notesvarfr, tklabel(notesvarfr, text="Name of the marrayRaw object:"),
         side='top')
  tkpack(notesvarfr, tkentry(notesvarfr, width=50, textvariable="rawName"), side='left')
  tkpack(mainfrm, notesvarfr, side="top")

  ## Files Names
  tkpack(mainfrm, filevarfr, side="top")

  # Select fg and bg
  coordvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  heading <- tklabel(coordvarfr, text="Foreground and background intensities")

  gf.label <- tklabel(coordvarfr, text="Green Foreground")
  gb.label <- tklabel(coordvarfr, text="Green Background")
  rf.label <- tklabel(coordvarfr, text="Red Foreground")
  rb.label <- tklabel(coordvarfr, text="Red Background")
  W.label <- tklabel(coordvarfr, text="Weights")

  gf.entry <- tkentry(coordvarfr, width=10, textvariable = "nameGf" )
  gb.entry <- tkentry(coordvarfr, width=10, textvariable = "nameGb" )
  rf.entry <- tkentry(coordvarfr, width=10, textvariable = "nameRf" )
  rb.entry <- tkentry(coordvarfr, width=10, textvariable = "nameRb" )
  W.entry <-  tkentry(coordvarfr, width=10, textvariable = "nameW" )


  tkgrid(heading, columnspan = 4)
  tkgrid(gf.label, gf.entry, gb.label, gb.entry)
  tkgrid(rf.label, rf.entry, rb.label, rb.entry)
  tkgrid(W.label, W.entry)
  tkgrid.configure(gf.label, gf.entry, gb.label, gb.entry)
  tkgrid.configure(rf.label, rf.entry, rb.label, rb.entry)
  tkgrid.configure(W.label, W.entry)

  tkpack(mainfrm, coordvarfr, side="top")

##  Jean: Jianhua, this is not quite working... what should I do?
##  if(ext == "spot")
##    {
##      tkinsert(gf.entry, "0.0", text="Gmean")
##      tkinsert(gb.entry, "0.0", text="morphG")
##      tkinsert(rf.entry, "0.0", text="Rmean")
##      tkinsert(rb.entry, "0.0", text="morphR")
##    }

  # Input Targets, Layout, GeneNames info
  notesvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)

  # I modified the following to put the browsers in
  ## Layout
  layoutFrame <- tkframe(notesvarfr)
  layoutLbl <- tklabel(layoutFrame, text="Layout:")
  layoutEntry <- tkentry(layoutFrame, width=43, textvariable="layout")
  layoutBrowser <- tkbutton(layoutFrame, text = "Browse",
                            command = layBrowser)
  tkgrid(layoutLbl, columnspan = 2)
  tkgrid(layoutEntry, layoutBrowser)
  tkpack(layoutFrame, side = "top")
  tkpack(mainfrm, notesvarfr, side="top")

  ## Target
  notesvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  targetFrame <- tkframe(notesvarfr)
  targetLbl <- tklabel(targetFrame, text="Target Information:")
  targetEntry <- tkentry(targetFrame, width=43, textvariable="targets")
  targetBrowser <- tkbutton(targetFrame, text = "Browse",
                            command = samBrowser)
  tkgrid(targetLbl, columnspan = 2)
  tkgrid(targetEntry, targetBrowser)
  tkpack(targetFrame, side = "top")
  tkpack(mainfrm, notesvarfr, side="top")

  ## Gene / Probe
  notesvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  geneFrame <- tkframe(notesvarfr)
  geneLbl <- tklabel(geneFrame, text="Gene Information:")
  geneEntry <- tkentry(geneFrame, width=43, textvariable="genenames")
  geneBrowser <- tkbutton(geneFrame, text = "Browse",
                          command = genBrowser)
  tkgrid(geneLbl, columnspan = 2)
  tkgrid(geneEntry, geneBrowser)
  tkpack(geneFrame, side = "top")
  tkpack(mainfrm, notesvarfr, side="top")

  ## Notes
  notesvarfr <- tkframe(mainfrm, relief="groove", borderwidth=2)
  tkpack(notesvarfr, tklabel(notesvarfr, text="Notes:"), side='top')
  tkpack(notesvarfr, tkentry(notesvarfr, width=50, textvariable="notes"), side='left')
  tkpack(mainfrm, notesvarfr, side="top")

  ## Button

  butfrm <- tkframe(mainfrm, borderwidth= 1, relief="groove")
#  a.but <- tkbutton(butfrm, command=function()
#                    {widget.marrayLayout()}, text="Layout")
  a.but <- tkbutton(butfrm, command=layoutPressed, text="Layout")
  b.but <- tkbutton(butfrm, command=function()
                    {widget.marrayInfo()}, text="Target")
#  c.but <- tkbutton(butfrm, command=function()
#                    {widget.marrayInfo()}, text="Genes")
  c.but <- tkbutton(butfrm, command=genePressed, text="Genes")
  d.but <- tkbutton(butfrm, command=function()
                    inputRaw(fnames,
                             sep=sep,
                             skip=skip,
                             quote=quote,
                             ...), text="Build")
  e.but <- tkbutton(butfrm, command=function()tkdestroy(base), text="Quit")
  tkpack(butfrm, a.but, b.but, c.but, d.but, e.but, side="left")
  tkpack(mainfrm, butfrm,anchor='w')

  # Jianhua add this in so that the widget is modal
  ##  tkwait.window(base)
  return(invisible())
}


















