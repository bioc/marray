############################################################################
# maBasicPlots.R
# TO BE REMOVE IN THE NEXT RELEASE [MARCH 15, 2004]
# Wrapper for diagnostic plots for two-color cDNA microarrays
#
###########################################################################
# Pre- and post-normalization plots

maDiagnPlots1<-function(mraw, title=NULL, save=TRUE,
	fname=paste(as.character(maLabels(maTargets(mraw))[1]),".ps",sep=""),
	dev=c("postscript","jpeg"))
{

  .Deprecated("maQualityPlots in package arrayQuality")

   mraw<-mraw[,1]
   # Default loess normalization within print-tip-group
   mnorm<-maNorm(mraw,norm="p")
  
  if(save==TRUE)
    do.call(dev,list(fname))
  
  layout(matrix(c(1:4,9,10,5:8,11,12),2,6,byrow=TRUE),width=c(5.5,3,5.5,3,6,6))

  # maImage
  stats<-c("maGb", "maRb", "maM","maM")
  titl<-c("Gbg","Rbg","Unnormalized M","Normalized M")
  Gcoltmp <- maPalette(low="white", high="green", k=50)
  Rcoltmp <- maPalette(low="white", high="red", k=50)
  RGcoltmp <- maPalette(low="green", high="red", k=50)
  cols<-list(Gcoltmp,Rcoltmp,RGcoltmp,RGcoltmp)
  m<-c("mraw","mraw","mraw","mnorm")
  for(i in 1:4)
  {
    x.bar<-do.call("maImage",list(m=eval(as.symbol(m[i])), x=stats[i], subset=TRUE, col=cols[[i]], contours=FALSE, bar=FALSE,main=titl[i]))$x.bar
     maColorBar(x.bar,horizontal=FALSE,col=cols[[i]],main="")
  }

  # maBoxplot
  maBoxplot(mraw,"maPrintTip","maM",main="Unnormalized M")
  maBoxplot(mnorm,"maPrintTip","maM",main="Normalized M")

  # maPlot
  defs<-maDefaultPar(mraw[,1],"maA","maM","maPrintTip")
  legend.func<-do.call("maLegendLines",defs$def.legend)
  args.lines<-c(list(TRUE,f=0.3),defs$def.lines)
  lines.func<-do.call("maLowessLines",args.lines)
  
  maPlot(mraw,"maA","maM","maPrintTip",lines.func,text.func=maText(),legend.func,main="Unnormalized MA-plot")
  maPlot(mnorm,"maA","maM","maPrintTip",lines.func,text.func=maText(),legend.func,main="Normalized MA-plot")

  # Back to defaults
  layout(1)

  if(!is.null(title))
    mtext(title, line=3)
  else
    mtext(paste(as.character(maLabels(maTargets(mraw))[1]),": Pre- and post- print-tip-group loess normalization",sep=""),line=3)

  if(save==TRUE)
    dev.off()
}


########################################################################### 
# Pre-normalization plots

maRawPlots<-function(mraw, title=NULL, save=TRUE, 
	fname=paste(as.character(maLabels(maTargets(mraw))[1]),".ps",sep=""), 
	dev=c("postscript","jpeg"))
{

  .Deprecated("maQualityPlots in package arrayQuality")
  mraw<-mraw[,1]
  if(save==TRUE)
    do.call(dev,list(fname))

  layout(matrix(c(1:6,7,7,8,8,9,9),2,6,byrow=TRUE),width=c(6,2.5,6,2.5,6,2.5))

  # maImage
  stats<-c("maGb", "maRb", "maM")
  titl<-c("Gbg","Rbg","Unnormalized M")
  m<-c("mraw","mraw","mraw")
  Gcoltmp <- maPalette(low="white", high="green", k=50)
  Rcoltmp <- maPalette(low="white", high="red", k=50)
  RGcoltmp <- maPalette(low="green", high="red", k=50)
  cols<-list(Gcoltmp,Rcoltmp,RGcoltmp)
  for(i in 1:3)
  {
    x.bar<-do.call("maImage",list(m=eval(as.symbol(m[i])), x=stats[i], subset=TRUE, col=cols[[i]], contours=FALSE, bar=FALSE,main=titl[i]))$x.bar
    maColorBar(x.bar,horizontal=FALSE,col=cols[[i]],main="")
  }

  # maBoxplot
  maBoxplot(mraw,x="maPrintTip", y="maM",main="Unnormalized M")
  maBoxplot(mraw,x="maPlate", y="maM",main="Unnormalized M")
  
  # maPlot
  defs<-maDefaultPar(mraw,"maA","maM","maPrintTip")
  legend.func<-do.call("maLegendLines",defs$def.legend)
  args.lines<-c(list(TRUE,f=0.3),defs$def.lines)
  lines.func<-do.call("maLowessLines",args.lines)
  
  maPlot(mraw,"maA","maM","maPrintTip",lines.func,text.func=maText(),legend.func,main="Unnormalized MA-plot")
  
  # Back to defaults
  layout(1)

  if(!is.null(title))
    mtext(title, line=3)
  else
    mtext(paste(as.character(maLabels(maTargets(mraw))[1]),": Pre-normalization",sep=""),line=3)

  
  if(save==TRUE)
    dev.off()
}


########################################################################## 
# Post-normalization plots

maNormPlots<-function(mnorm, title=NULL, save=TRUE, 
	fname=paste(as.character(maLabels(maTargets(mnorm))[1]),".ps",sep=""),
	dev=c("postscript","jpeg"))
{

  .Deprecated("maQualityPlots in package arrayQuality")
  mnorm<-mnorm[,1]
  if(save==TRUE)
    do.call(dev,list(fname))

  layout(matrix(c(1:4,5,5,6,6),2,4,byrow=TRUE),width=c(7.5,2.5,7.5,2.5))

  # maImage
  stats<-c("maMloc", "maM")
  titl<-c("Loc. normalization", "Normalized M")
  RGcoltmp <- maPalette(low="green", high="red", k=50)
  cols<-list(RGcoltmp,RGcoltmp)
  for(i in 1:2)
  {
    x.bar<-do.call("maImage",list(m=mnorm, x=stats[i], subset=TRUE, col=cols[[i]], contours=FALSE, bar=FALSE,main=titl[i]))$x.bar
    maColorBar(x.bar,horizontal=FALSE,col=cols[[i]],main="")
  }

  # maBoxplot
  maBoxplot(mnorm,x="maPrintTip",y="maM",main="Normalized M")
  
  # maPlot
  defs<-maDefaultPar(mnorm,"maA","maM","maPrintTip")
  legend.func<-do.call("maLegendLines",defs$def.legend)
  args.lines<-c(list(TRUE,f=0.3),defs$def.lines)
  lines.func<-do.call("maLowessLines",args.lines)
  
  maPlot(mnorm,"maA","maM","maPrintTip",lines.func,text.func=maText(),legend.func,main="Normalized MA-plot")
  

  # Back to defaults
  layout(1)

  if(!is.null(title))
    mtext(title, line=3)
  else
    mtext(paste(as.character(maLabels(maTargets(mnorm))[1]), ": ", as.character(list(maNormCall(mnorm)))," normalization",sep=""),line=3)

  if(save==TRUE)
    dev.off()
}

########################################################################### 
