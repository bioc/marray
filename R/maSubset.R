###########################################################################
# maSubset.R
#
# Subsetting methods for microarray classes
#
###########################################################################
# marrayInfo class


setMethod("[", "marrayInfo", function(x, i, j, ..., drop=FALSE)
{
  newx<-x
  if(missing(j))
   j<-TRUE
  if(missing(i))
   i<-TRUE
  if(length(maLabels(x))!=0)
    slot(newx,"maLabels")<-maLabels(x)[i]
  if(length(maInfo(x))!=0)
    slot(newx,"maInfo")<-maInfo(x)[i, j, drop=FALSE]
  return(newx)
})

#########################
# marrayLayout class

setMethod("[", "marrayLayout", function(x, i, j, ..., drop=FALSE)
{
  newx<-x
  if(!missing(i))
  {
    if(length(maPlate(x))!=0)
      slot(newx,"maPlate")<-maPlate(x)[i]
    if(length(maControls(x))!=0)
    slot(newx,"maControls")<-maControls(x)[i]
    if(length(maNspots(x))!=0)
    {
      sub<-rep(FALSE,maNspots(x))
      sub[maSub(x)][i]<-TRUE
      slot(newx, "maSub")<-sub
    }
  }
  return(newx)
})

#########################
# marrayRaw class

setMethod("[", "marrayRaw", function(x, i, j, ..., drop=FALSE)
{
  newx<-x
  if(missing(j))
   j<-TRUE
  if(missing(i))
   i<-TRUE

  if(length(maRf(x))!=0)
    slot(newx,"maRf")<-maRf(x)[i,j,drop=FALSE]
  if(length(maGf(x))!=0)
    slot(newx,"maGf")<-maGf(x)[i,j,drop=FALSE]
  if(length(maRb(x))!=0)
    slot(newx,"maRb")<-maRb(x)[i,j,drop=FALSE]
  if(length(maGb(x))!=0)
    slot(newx,"maGb")<-maGb(x)[i,j,drop=FALSE]
 if(length(maW(x))!=0)
    slot(newx,"maW")<-maW(x)[i,j,drop=FALSE]
  slot(newx,"maLayout")<-maLayout(x)[i]
  slot(newx,"maGnames")<-maGnames(x)[i]
  slot(newx,"maTargets")<-maTargets(x)[j]
  return(newx)
})

#########################
# marrayNorm class

setMethod("[", "marrayNorm", function(x, i, j, ..., drop=FALSE)
{
  newx<-x
  if(missing(j))
   j<-TRUE
  if(missing(i))
   i<-TRUE


  if(length(maA(x))!=0)
    slot(newx,"maA")<-maA(x)[i,j,drop=FALSE]
  if(length(maM(x))!=0)
    slot(newx,"maM")<-maM(x)[i,j,drop=FALSE]
  if(length(maMloc(x))!=0)
    slot(newx,"maMloc")<-maMloc(x)[i,j,drop=FALSE]
  if(length(maMscale(x))!=0)
    slot(newx,"maMscale")<-maMscale(x)[i,j,drop=FALSE]
  if(length(maW(x))!=0)
    slot(newx,"maW")<-maW(x)[i,j,drop=FALSE]
  slot(newx,"maLayout")<-maLayout(x)[i]
  slot(newx,"maGnames")<-maGnames(x)[i]
  slot(newx,"maTargets")<-maTargets(x)[j]
  slot(newx,"maNormCall")<-maNormCall(x)
  return(newx)
})

###########################################################################
