##########################################################################
# maComp.R
#
# Simple computations for cDNA microarray objects
#
###########################################################################
# Convert integer vector of indices to logical vector

maNum2Logic<-function(n=length(subset), subset=TRUE)
{
  if(is.logical(subset))
    return(subset)
  if(is.numeric(subset))
  {
    which<-rep(FALSE,n)
    which[subset]<-TRUE 
    return(which)
  }
}

###########################################################################
# Produce a table of with spot coordinates and gene names

maGeneTable <- function(object)
{
  tmp <- data.frame(
  GR=maGridRow(object),
  GC=maGridCol(object),
  SR=maSpotRow(object),
  SC=maSpotCol(object),
  maInfo(maGnames(object)))
  colnames(tmp) <- c("Grid.R", "Grid.C", "Spot.R", "Spot.C", names(maInfo(maGnames(object))))
  return(tmp)
}

###########################################################################
# Generate plate IDs from dimensions of grid and spot matrices
## Modified Feb 27, 2004

maCompPlate <- function (x, n = 384)
{
    is.int <- function(z) trunc(z)==z
    totalPlate <- maNspots(x)/n
    if (is.int(totalPlate))
      {
        tmp <- n/(maNgr(x) * maNgc(x))
        factor(rep(rep(1:totalPlate, rep(tmp, totalPlate)), (maNgr(x) *
                                                             maNgc(x))))[maSub(x)]
      }
    else
      {
        totalPlate.c <- ceiling(totalPlate)
        totalPlate.f <- floor(totalPlate)
        tmp1 <- n/(maNgr(x) * maNgc(x))
        tmp2 <- ((totalPlate-totalPlate.f) * n ) / (maNgr(x) * maNgc(x))
        factor(rep(rep(1:totalPlate.c, c(rep(tmp1, totalPlate.f),tmp2)),(maNgr(x)* maNgc(x))))[maSub(x)]
      }
  }

##maCompPlate<-function(x, n=384){
##  totalPlate <- maNspots(x) /n
##  tmp <- n / (maNgr(x) * maNgc(x))
##  factor(rep(rep(1:totalPlate, rep(tmp, totalPlate)), (maNgr(x) * maNgc(x))))
##}

###########################################################################
# Convert spot index to grid and spot matrix coordinates (4 coords)

maInd2Coord <- function (x, L)
{
    coord<-cbind(maGridRow(L), maGridCol(L), maSpotRow(L), maSpotCol(L))[x,]
    colnames(coord) <- c("Grid.R", "Grid.C", "Spot.R", "Spot.C")
    coord
}

# Convert grid and spot matrix coordinates (4 coords) to spot index.
# Works for subsetted arrays

maCoord2Ind <- function (x, L)
{
  ngr<-maNgr(L)
  ngc<-maNgc(L)
  nsr<-maNsr(L)
  nsc<-maNsc(L)
  n<-maNspots(L)
  ind<-(nsr * nsc)* ((x[,1] - 1) * ngc + (x[,2] - 1)) + (x[,3] - 1) *
        nsc + x[,4]
  ord<-order(ind)
  ind<-ind[ord]
  coord<-x[ord,]
  sub<-(1:n)[maSub(L)]
  ind<-intersect(ind, sub)
  ind<-order(sub)[sub %in% ind]
  ind
}

###########################################################################
# Generate grid and spot matrix coordinates from ranges of rows and
# columns for the grid and spot matrices

maCompCoord<-function(grows, gcols, srows, scols)
{
  ngr <- length(grows)
  ngc <- length(gcols)
  nsr <- length(srows)
  nsc <- length(scols)
  t1 <- rep(grows, rep(nsr * nsc * ngc, ngr))
  t2 <- rep(rep(gcols, rep(nsr * nsc, ngc)), ngr)
  t3 <- rep(rep(srows, rep(nsc, nsr)), ngc * ngr)
  t4 <- rep(scols, nsr * ngc * ngr)
  coord<-cbind(t1, t2, t3, t4)
  colnames(coord) <- c("Grid.R", "Grid.C", "Spot.R", "Spot.C")
  coord
}

# Generate spot index from ranges of rows and columns for the grid 
# and spot matrices

maCompInd<-function(grows, gcols, srows, scols, L)
{
  coord<-maCompCoord(grows, gcols, srows, scols)
  maCoord2Ind(coord, L)
}


###########################################################################
##
## Added on April 5, 2004
##
###########################################################################

maCompLayout <- function(mat, ncolumns=4)
  {
    if (dim(mat)[2]==3)
      {
        Blocks <- mat[,1]
        gr <- mat[,1] / ncolumns
        newmat <- cbind(gr=(mat[,1] - 1) %/% ncolumns, gc=((mat[,1] -1) %% 4) + 1,
                        sr=mat[,2], sc=mat[,3])
      }  else
    newmat <- mat
    
    ngr <- max(newmat[,1]);  ngc <- max(newmat[,2])
    nsr <- max(newmat[,3]);  nsc <- max(newmat[,4])
    nspots <- as.integer(ngr) * as.integer(ngc) * as.integer(nsr) *   as.integer(nsc)
    
    mlayout <- new("marrayLayout", maNgr = as.integer(ngr),
                   maNgc = as.integer(ngc),
                   maNsr = as.integer(nsr),
                   maNsc = as.integer(nsc),
                   maNspots = nspots)
    maSub(mlayout) <- maCoord2Ind(newmat, mlayout)
    return(mlayout)
  }


###########################################################################
##
##  Move controlCode and maGenControls from marrayTools to marrayClasses
##  May 7, 2003
###########################################################################
controlCode <-
structure(c("Buffer", "Empty", "EMPTY", "AT", "NC", "M200009348", "M200012700", 
"M200016219", "M200016205", "M200013499", "M200003425", "M200006376", 
"M200001318", "M200004477", "M200001732", "M200006590", "M200000829", 
"H200000553", "H200000680", "H200001719", "H200001847", "H200007830", 
"H200008181", "H200008484", "H200008489", "H200009216", "H200009498", 
"H200011103", "H200019704", "18S", "SSC", "mCNR", "NLG", "GABA", 
"pNICE", "M13", "Cot-1", "orward", "everse", "Genomic", "p133", 
"Blam", "T7/SP6", "ephrin", "Buffer", "Empty", "Empty", "Negative", "Negative", 
"Positive", "Positive", "Positive", "Positive", "Positive", "Positive", 
"Positive", "Positive", "Positive", "Positive", "Positive", "Positive", 
"Positive", "Positive", "Positive", "Positive", "Positive", "Positive", 
"Positive", "Positive", "Positive", "Positive", "Positive", "Positive", 
"Positive", "Buffer", "con", "con", "con", "con", "con", "Positive", 
"Negative", "Negative", "Negative", "Negative", "Negative", "con", 
"con"), .Dim = c(44, 2), .Dimnames = list(c("1", "2", "3", "4", 
"5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", 
"16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", 
"27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", 
"38", "39", "40", "41", "42", "43", "44"), c("Pattern", "Name")))

## But Fix March 03, 2003
## Replace controlCode by controlcode (1.0.6)
## Add controlcode = controlCode (1.0.9), add defaults for gpTools and spotTools to work

maGenControls <- function(Gnames, controlcode=controlCode, id="ID")
{
  if(class(Gnames) == "marrayInfo")
    {
      ifelse(is.numeric(id), tmp <- id, tmp <- grep(id, colnames(maInfo(Gnames))))
      if(length(tmp) == 0)
        {
          tmp <- 1
          print("Controls are generated from an arbitaray columns\n")
        }
      ID <- as.vector(maInfo(Gnames[,tmp])[[1]])
    }
  else
    {
      if(is.null(dim(Gnames)))
        ID <- Gnames
      else
        {
          ifelse(is.numeric(id), tmp <- id, tmp <- grep(id, colnames(Gnames)))
          ID <- as.vector(Gnames[,tmp])
        }
    }
  Control <- rep("probes", length(ID))
  for(i in 1:nrow(controlcode))
    {
      position <- grep(as.vector(controlcode[i,"Pattern"]), ID)
      if(length(position) > 0)
        Control[position] <- as.vector(controlcode[i, "Name"])
    }
  return(Control)
}
