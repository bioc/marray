###################################################################
##
## Date: October 11, 2002
## 
## source("~/Projects/maTools/R/maAnnotate.R")
## 
###################################################################

##########################################################################
## Widget for html page

mapGeneInfo <- function(widget=FALSE, Gnames, Name="pubmed", ID="genbank", ACC="SMDacc",  ...)
  {
    if(widget)
      {
        res <- widget.mapGeneInfo(Gnames)
        return(res)
      }
    else
      {
        opt <- list(...)
        base <- matrix(c("Grid", "Spot", "Row", "Column", "Block",
                         "cood", "cood", "cood", "cood", "cood"), ncol=2)
        rownames(base) <- c("Grid", "Spot", "Row", "Column", "Block")
        newinfo <- rbind(c("Name", Name),
                         c("ID", ID),
                         c("ACC", ACC))
        rownames(newinfo) <- c("Name", "ID", "ACC")
        return(rbind(newinfo, cbind(names(opt), unlist(opt)),base))
      }
  }

widget.mapGeneInfo <- function(Gnames)
  {
    print("widget")
    startfun <- function()
      {
        print("The URL choices are:")
        print(names(URLstring))
      }
    
    require(tcltk)
    require(tkWidgets)
    switch(data.class(Gnames),
           marrayNorm = headings <- colnames(maInfo(maGnames(Gnames))),
           marrayRaw= headings <- colnames(maInfo(maGnames(Gnames))),
           data.frame = headings <- colnames(Gnames),
           headings <- colnames(Gnames)
           )
    
    headings <- headings[-unique(c(grep("Grid", headings),
                                   grep("Spot", headings),
                                   grep("Row", headings),
                                   grep("Column", headings),
                                   grep("Block", headings)))]

    wlist <- list()
    for(hvalue in headings)
      {
        test <- list(Name=hvalue, Value=hvalue,
                     toText=function(x) paste(x,collapse = ","),
                     fromText=NULL, canEdit=TRUE, buttonFun = NULL,
                     buttonText = "Choices")
        wlist <- c(wlist, list(test))
      }
    names(wlist) <- headings
    widget1 <- list(wList = wlist,
                    preFun = startfun)
    res <- widgetRender(widget1, "Map Gene Names")

    resValues <- values.Widget(res)
    base <- matrix(c("Grid", "Spot", "Row", "Column", "Block",
                     "cood", "cood", "cood", "cood", "cood"), ncol=2)
    for(i in 1:length(resValues))
      base <- rbind(base, c(resValues[[i]]$Entry, resValues[[i]]$Value))
    return(base)
  }


##########################################################################
htmlPage <- function(genelist,
                     filename="GeneList.html",
                     geneNames=Gnames,
                     mapURL=SFGL,
                     othernames,
                     title,
                     table.head,
                     table.center=TRUE,
                     disp=c("browser", "file")[1])
{
  switch(class(geneNames),
         data.frame= data <- geneNames,
         marrayRaw = data <- maGeneTable(geneNames),
         marrayNorm = data <- maGeneTable(geneNames),
         marrayInfo = data <- maInfo(geneNames),
         matrix = data <- data.frame(geneNames),
         data <- geneNames)

  if(missing(othernames))
    restable <- data[genelist,] else
  restable <- cbind(data, othernames)[genelist,]
  
  args <- list(filename = filename, mapURL = mapURL,
               table.center = table.center,disp = disp)
  if(!missing(title)) args <- c(args, list(title=title))
  if(!missing(table.head)) args <- c(args, list(table.head=table.head))
  do.call("table2html", c(list(restable), args))
  return()
}

##ll.htmlpage <- function(genelist, filename, title, othernames, table.head, 
##                        table.center = TRUE) 
##  {
##    restable <- cbind(LocusLink = genelist, othernames)
##    args <- list(table.center=table.center)
##    if(!missing(filename)) args <- c(args, list(filename="GeneList.html"))
##    if(!missing(title)) args <- c(args, list(title=title))
##    if(!missing(table.head)) args <- c(args, list(table.head=table.head))
##    do.call("table2html", c(list(restable), args))
##    return()
##  }
                        

#####################################################
## Base Function
##
tablegen <-  function(input)
  {
    HTwrap <-   function (x, tag = "TD") {
      paste("<", tag, ">", x, "</", tag, ">", sep = "")}
    
    HTwrap.matrix <- function(input)
      {
        output <- ""
        for (nm in 1:ncol(input))
          output <- paste(output, HTwrap(input[,nm]), sep = "")
        return(output)
      }

    HTwrap.list <- function(input)
      {
        output <- ""
        for (nm in 1:length(input))
          output <- paste(output, HTwrap(input[[nm]]), sep = "")
        return(output)
      }
    
    switch(data.class(input),
           vector = output <- HTwrap(input),
           matrix = output <- HTwrap.matrix(input),
           list = output <- HTwrap.list(input),
           output <- HTwrap(input)
           )
    return(output)
  }


opVersionID <- function(opID)
  {
    code <- unlist(lapply(strsplit(as.vector(opID), split=""),
                          function(x){paste(x[1:2], collapse="")}))
    tmp <- table(code)
    code2 <- names(tmp)[tmp==max(tmp)]
    switch(code2,
           M2 = res <- "operonm2",
           M0 = res <- "operonm1",
           H2 = res <- "operonh2",
           H0 = res <- "operonh1"
           )
    return(res)
  }


gsubAnchor <-function (id, urlString) 
{
  test <-  function(x){
    if(!is.na(x))
      res <- gsub(pattern="UNIQID", replacement=x, urlString)
    else
      res <- x
    return(res)
  }
  paste("<A HREF=", sapply(as.character(id), test), ">", id, "</A>", sep = "")
}
#####################################################
## Table 2 HTML
## Date: Feb 16, 2003
table2html <- function (restable, filename = "GeneList.html",
                        mapURL = SFGL, title, table.head, table.center = TRUE, 
                        disp = c("browser", "file")[1]) 
{

  HTwrap <- function(x, tag = "TD") {
    paste("<", tag, ">", x, "</", tag, ">", sep = "")
  }
  
  outfile <- file(filename, "w")
  cat("<html>", file = outfile)
  cat(HTwrap(HTwrap("BioConductor Gene Listing", tag = "TITLE"), 
             tag = "head"), file = outfile)
  cat("<body bgcolor=\"#FFFFEF\">", "<H1 ALIGN=CENTER > BioConductor Gene Listing </H1>", 
      file = outfile, sep = "\n")
  if (!missing(title)) 
    cat("<CENTER><H1 ALIGN=\"CENTER\">", title, " </H1></CENTER>\n", 
        file = outfile, sep = "\n")
  if (table.center) 
    cat("<CENTER> \n", file = outfile)
  cat("<TABLE BORDER=4>", file = outfile, sep = "\n")
  if (!missing(table.head)) {
    headout <- paste("<TH>", table.head, "</TH>")
    cat("<TR>", headout, "</TR>", file = outfile, sep = "\n")
  }
  if (is.null(mapURL)) 
    mapURL <- widget.mapGeneInfo(restable)

  oldGnamesID <- colnames(restable)
  GnamesID <- rep("none", length(oldGnamesID))
  for (i in 1:nrow(mapURL))
    GnamesID[grep(mapURL[i, 1], oldGnamesID)] <- mapURL[i,2]
  
  if (sum(GnamesID == "operon") != 0) 
    GnamesID[grep("operon", GnamesID)] <- opVersionID(restable[1:100, grep("operon", GnamesID)])
  mainTable <- Headings <- NULL
  
  for (i in 1:length(GnamesID)) {
    info <- GnamesID[i]
    x <- as.vector(restable[, i])
    if(!is.null(class(x))) if(class(x) == "numeric") x <- round(x, 2)
    if ((info != "") | is.null(info)) {
      switch(info, cood = mainTable <- paste(mainTable, HTwrap(x), sep = ""),
             none = mainTable <- paste(mainTable,  HTwrap(x), sep = ""),
             mainTable <- paste(mainTable,
                                HTwrap(gsubAnchor(x, urlString = URLstring[[info]])), sep = ""))
    }
    Headings <- c(Headings, colnames(restable)[i])
  }
  
  cat(paste(HTwrap(Headings), collapse = ""), file = outfile)
  cat("\n", file = outfile)
  cat(HTwrap(mainTable, tag = "TR"), file = outfile, sep = "\n")
  cat("</TABLE>", "</body>", "</html>", sep = "\n", file = outfile)
  if (table.center) 
    cat("</CENTER> \n", file = outfile)
  close(outfile)

  if (disp == "browser") 
    browseURL(paste("file://", getwd(), filename, sep = "/"))
  ##  openBrowser(paste("file://", getwd(), filename, sep = "/"))
  return()
}



###################################################################
## predefine info
URLstring <- list(
 pubmed = "http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&db=PubMed&term=UNIQID",
 locuslink = "http://www.ncbi.nlm.nih.gov/LocusLink/LocRpt.cgi?l=UNIQID",
 riken = "http://read.gsc.riken.go.jp/chipinfo.php?defkey=&chiprearrayid=UNIQID",
 SMDclid = "http://genome-www4.stanford.edu/cgi-bin/SMD/source/sourceResult?option=CloneID&criteria1=IMAGE:UNIQID&choice=cDNA",
 SMDacc = "http://genome-www4.stanford.edu/cgi-bin/SMD/source/sourceResult?option=Number&criteria=UNIQID&choice=Gene",
 operonh2 = "http://oparray.operon.com/human2/index.php?single_query=UNIQID",
 operonh1 = "http://oparray.operon.com/~operon/human/index.php?single_query=UNIQID",
 operonm2 = "http://oparray.operon.com/mouse2/index.php?single_query=UNIQID",
 operonm1 = "http://oparray.operon.com/~operon/mouse/index.php?single_query=UNIQID",
 operonST="http://sandlertest.ucsf.edu/NOMAD/nomad-cgi/query_annot.pl?UNIQID",
 genbank ="http://www.ncbi.nlm.nih.gov/entrez/viewer.fcgi?DB=nucleotide&val=UNIQID",
 unigeneMm="http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Mm&CID=UNIQID",
 unigeneHS="http://www.ncbi.nlm.nih.gov/UniGene/clust.cgi?ORG=Hs&CID=UNIQID")
                  

SFGL <- mapGeneInfo(ID="operonST",
                    ACC="SMDacc",
                    LocusLink="locuslink",
                    Cluster="unigeneMm",
                    LOCUSLINK="locuslink",
                    GenBank="genbank",
                    Name="none")

UCBFGL <- mapGeneInfo(ID="riken",
                      ACC="SMDacc")

###################################################################
