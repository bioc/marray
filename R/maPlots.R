############################################################################
# maPlot.R
# maPlot3.R
#
# S4 methods
# Diagnostic plots for two-color cDNA microarrays
#
###########################################################################

## S3 + S4 settings
## March 15, 2004
plot.marrayRaw <-
  function (x, xvar = "maA", yvar = "maM", zvar="maPrintTip", lines.func,text.func,legend.func,...)
  {maPlot(m=x, x=xvar, y=yvar, lines.func=lines.func, text.func=text.func, legend.func=legend.func,...)}

plot.marrayNorm <-
  function (x, xvar = "maA", yvar = "maM", zvar="maPrintTip", lines.func,text.func,legend.func,...)
  {maPlot(m=x, x=xvar, y=yvar, lines.func=lines.func, text.func=text.func, legend.func=legend.func,...)}

setMethod("boxplot", signature(x="marrayRaw"), function (x, xvar = "maPrintTip", yvar = "maM", ...){
            maBoxplot(m=x, x=xvar, y=yvar, ...)})
setMethod("boxplot", signature(x="marrayNorm"), function (x, xvar = "maPrintTip", yvar = "maM", ...){
              maBoxplot(m=x, x=xvar, y=yvar, ...)})
setMethod("image", signature(x="marrayRaw"),
          function (x, xvar = "maM", subset = TRUE, col, contours = FALSE,  bar = TRUE, overlay=NULL, ...){
            maImage(m=x, x=xvar, subset=subset, col=col, contours=contours, bar=bar, overlay=overlay, ... )})
setMethod("image", signature(x="marrayNorm"),
          function (x, xvar = "maM", subset = TRUE, col, contours = FALSE,  bar = TRUE, overlay=NULL, ...){
              maImage(m=x, x=xvar, subset=subset, col=col, contours=contours, bar=bar, overlay=overlay, ... )})

## points, lines, text
setMethod("points", signature(x="marrayRaw"), function (x, xvar = "maA", yvar = "maM", ...)
          {addPoints(object=x, xvar=xvar, yvar=yvar, ...)})
setMethod("points", signature(x="marrayNorm"), function (x, xvar = "maA", yvar = "maM", ...)
          {addPoints(object=x, xvar=xvar, yvar=yvar, ...)})
setMethod("text", signature(x="marrayRaw"), function (x, xvar = "maA", yvar = "maM",...)
          {addText(object=x, xvar=xvar, yvar=yvar, ...)})
setMethod("text", signature(x="marrayNorm"), function (x, xvar = "maA", yvar = "maM", ...)
          {addText(object=x, xvar=xvar, yvar=yvar, ...)})
setMethod("lines", signature(x="marrayRaw"), function (x, xvar = "maA", yvar = "maM", zvar="maPrintTip",...)
          {addLines(object=x, xvar=xvar, yvar=yvar, ...)})
setMethod("lines", signature(x="marrayNorm"), function (x, xvar = "maA", yvar = "maM", zvar="maPrintTip",...)
          {addLines(object=x, xvar=xvar, yvar=yvar, ...)})


addText <- function(object, xvar="maA", yvar="maM", subset=NULL, labels = as.character(1:length(subset)), ...)
  {
    text.func <- maText(subset=subset, labels=labels, ...)
    text.func(as.numeric(eval(call(xvar,object))), as.numeric(eval(call(y,m))))
  }


addPoints <- function(object, xvar="maA", yvar="maM", subset=TRUE, ...)
  {
    xx <- as.numeric(eval(call(xvar,object)))
    yy <- as.numeric(eval(call(yvar,object)))
    points(xx[subset], yy[subset], ...)
  }

addLines <- function(object, xvar="maA", yvar="maM", zvar="maPrintTip", subset=TRUE, ...)
  {
    defs <- maDefaultPar(m=object,x=xvar,y=yvar,z=zvar)
    lines.func<-do.call("maLoessLines",
                        c(list(subset=subset, loess.args=list(span=0.4, degree=1, family="symmetric",
                                              control=loess.control(trace.hat="approximate",
                                                iterations=5,surface="direct"))),defs$def.lines))
    xx <- as.numeric(eval(call(xvar,object)))
    yy <- as.numeric(eval(call(yvar,object)))
    zz <- as.numeric(eval(call(zvar,object)))
    lines.func(xx, yy, zz)
  }
 

###########################################################################
# Default plotting parameters for microarray objects
###########################################################################
# Compare default and ... plotting parameters and let ... overwrite defaults

maDotsDefaults<-function(dots, defaults)
{
  args<-c(dots,defaults[setdiff(names(defaults),names(dots))])
  return(args)
}

#########################
# Default parameters for microarray objects of class marrayRaw and marrayNorm
maDefaultPar<-function(m,x,y,z)
{
    m<-m[,1]
    main<-as.character(maLabels(maTargets(m)))

    xlab<-ylab<-zlab<-""
    col<-2
    lty<-1
    lwd<-2.5
    las<-1
    names<-""
    def.legend<-list()

    ylab<-strsplit(y,"ma")[[1]]
    if(ylab[1]=="")
     ylab<-ylab[2]

    if(!is.null(x))
    {
      xlab<-strsplit(x,"ma")[[1]]
      if(xlab[1]=="")
        xlab<-xlab[2]
    }

    if(!is.null(z))
    {
      zz<-eval(call(z,m))
      zlab<-strsplit(z,"ma")[[1]]
      if(zlab[1]=="")
        zlab<-zlab[2]

      if(z!="maPrintTip")
      {
        ## names<-paste(zlab,unique(zz),sep=" ")  ## BUG [modified Feb 27, 2004]
        ifelse(is.factor(zz), tmp <- levels(zz), tmp<- unique(zz))
        names<-paste(zlab, tmp, sep=" ")  
        col<-(1:length(names))
        lty<-rep(1,length(names))
        las<-1
        ncol<-1
        ord<-order(lty,col)
        def.legend<-list(legend=names[ord],col=col[ord], lty=lty[ord], lwd=lwd, ncol=ncol)
      }

      if(z=="maPrintTip")
      {
        which<-unique(zz)
        Ig<-maNgr(m)
        Jg<-maNgc(m)
        lg.names<- paste("(",sort(rep(1:Ig,Jg)),",",rep(1:Jg,Ig),")",sep="")
        lg.col<-sort(rep(2:(Ig+1),Jg))
        lg.lty<-rep(1:Jg,Ig)
        names<-lg.names[which]
        col<-lg.col[which]
        lty<-lg.lty[which]
        las<-3
        ncol<-Jg
        ord<-order(lg.lty,lg.col)
        def.legend<-list(legend=lg.names[ord],col=lg.col[ord], lty=lg.lty[ord], lwd=lwd, ncol=ncol)
      }
   }

    def.box<-list(xlab=xlab,ylab=ylab,names=names,col=col,las=las,main=main)
    def.plot<-list(xlab=xlab,ylab=ylab,pch=20,col=1,main=main)
    def.lines<-list(col=col,lty=lty,lwd=lwd)
    def.text<-list(pch=16,col="purple")
    return(list(def.box=def.box,def.plot=def.plot,def.lines=def.lines,
		def.legend=def.legend,def.text=def.text))
}

###########################################################################
# maBoxplot: Boxplot methods
###########################################################################
# Boxplots for a single and multiple arrays

maBoxplot<- function(m, x="maPrintTip", y="maM", ...) {
    opt<-list(...)

    if(maNsamples(m)==1)
    {

      yy<-as.numeric(eval(call(y,m)))

      if(is.null(x))
        xx<-rep(1,length(yy))
      if(!is.null(x))
        xx<-eval(call(x,m))

      # Optional graphical parameter defaults
      def<-maDefaultPar(m,x,y,x)$def.box
      if(!is.null(opt))
        def<-maDotsDefaults(opt,def)

      args<-c(list(yy~xx),def)
      do.call("boxplot",args)
   }
   if(maNsamples(m)>1)
   {
     yy<-as.data.frame(eval(call(y,m)))

     # Optional graphical parameter defaults
     if(length(maLabels(maTargets(m))) != 0)
       def <- list(names=maLabels(maTargets(m)),ylab=strsplit(y,"ma")[[1]][2],col=2)
     else
       def <- list(names=dimnames(yy)[[2]], ylab=strsplit(y,"ma")[[1]][2],col=2)

     if(!is.null(opt))
        def<-maDotsDefaults(opt,def)

      args<-c(list(yy),def)
      do.call("boxplot",args)
    }
##    if(y=="maM") abline(h=0,col="gray",lwd=2.5)
}

###########################################################################
##
# maPlot: Scatter-plot methods with fitted lines and points highlighted
##
###########################################################################
# General function for scatter-plot of y vs. x with fitted lines within
# values of z and subset of points highlighted

maPlot.func<-function(x, y, z,
	lines.func=maLowessLines(subset=TRUE,f=0.3,col=1:length(unique(z)),
		lty=1,lwd=2.5),
	text.func=maText(),
	legend.func=maLegendLines(legend=as.character(unique(z)),
		col=1:length(unique(z)), lty=1,lwd=2.5,ncol=1),
	...)
{
  plot(x,y,...)

  # Plot fitted curves
  if(!is.null(lines.func))
    lines.func(x,y,z)

  # Legend
  if(!is.null(legend.func))
    legend.func(x,y)

  # Label a subset of points
  if(!is.null(text.func))
    text.func(x,y)

}

#########################
# Label a subset of points
# Jean: Modify (Oct 21, 2002)  line "tmp <- length(c(1:length(subset))[subset])"

maText <-function (subset = NULL, labels = as.character(1:length(subset)),
    ...)
{
    function(x, y) {
      tmp <- length(c(1:length(subset))[subset])
      if (tmp > 0) {
            if (length(subset) < length(labels))
                text(x[subset], y[subset], labels[subset], ...)
            if (length(subset) > length(labels))
                text(x[subset], y[subset], labels, ...)
            if ((length(subset) == length(labels)) & is.logical(subset))
                text(x[subset], y[subset], labels[subset], ...)
            if ((length(subset) == length(labels)) & is.numeric(subset))
                text(x[subset], y[subset], labels, ...)
        }
    }
}


#########################
# Plot fitted lines

# Lowess
maLowessLines<-function(subset=TRUE, f=0.3, col=2, lty=1, lwd=2.5,...)
{
  function(x,y,z)
  {
    subset<-maNum2Logic(length(x), subset)
    g<-unique(z[subset])
    if(length(col)<length(g))
      col<-rep(col[1],length(g))
    if(length(lty)<length(g))
      lty<-rep(lty[1],length(g))

    for(i in (1:length(g)))
    {
      which<-z[subset]==g[i]
      xx<-x[subset][which]
      yy<-y[subset][which]
      ind <- is.na(xx) | is.na(yy) | is.infinite(xx) | is.infinite(yy)
      fit<- lowess(xx[!ind], yy[!ind], f=f)
      lines(fit,col=col[i],lty=lty[i],lwd=lwd,...)
    }
  }
}

# Loess

maLoessLines<-function(subset=TRUE, weights=NULL,
                        loess.args=list(span=0.4, degree=1, family="symmetric",
                          control=loess.control(trace.hat="approximate",
                            iterations=5,surface="direct")),col=2, lty=1, lwd=2.5, ...)
{
  function(x,y,z)
  {
    subset<-maNum2Logic(length(x), subset)
    g<-unique(z[subset])
    if(length(col)<length(g))
      col<-rep(col[1],length(g))
    if(length(lty)<length(g))
      lty<-rep(lty[1],length(g))

    for(i in (1:length(g)))
    {
      which<-z[subset]==g[i]
      xx<-x[subset][which]
      yy<-y[subset][which]
      ww<-weights[subset][which]
      args<-c(list(yy ~ xx, weights=ww),loess.args)
      fit<-do.call("loess",args)
      xf<-seq(quantile(xx,0.005,na.rm=TRUE),quantile(xx,0.995,na.rm=TRUE),length=100)
      yf<-predict(fit,data.frame(xx=xf))
      lines(xf,yf,col=col[i],lty=lty[i],lwd=lwd,...)
    }
  }
}

#########################
# Add legend to existing plot

maLegendLines<-function(legend="", col=2, lty=1, lwd=2.5, ncol=1, ...)
{
  function(x,y)
  {
    a<-min(x[!(is.na(x)|is.infinite(x))])
    b<-max(y[!(is.na(y)|is.infinite(y))])
    legend(a,b,legend=as.character(legend),col=col,lty=lty,lwd=lwd,ncol=ncol,...)
  }
}

###########################################################################
# Methods for microarray objects: wrapper around maPlot.func
##
## Jean April 9,2003 modified default to maLoessLines
##

maPlot <- function(m, x="maA", y="maM", z="maPrintTip",lines.func,text.func,legend.func, ...)
{

  m<-m[,1]
  # Default plotting arguments
  defs<-maDefaultPar(m,x,y,z)

  if(missing(lines.func))
    lines.func<-do.call("maLoessLines", c(list(subset=TRUE,
                                               loess.args=list(span=0.4, degree=1, family="symmetric",
                                                 control=loess.control(trace.hat="approximate",
                                                   iterations=5,surface="direct"))),
                                               defs$def.lines))
  if(missing(text.func))
    text.func<-maText()
  if(missing(legend.func))
    legend.func<-do.call("maLegendLines",defs$def.legend)

  xx<-as.numeric(eval(call(x,m)))
  yy<-as.numeric(eval(call(y,m)))
  if(is.null(z))
     zz<-rep(1,length(xx))
  if(!is.null(z))
    zz<-eval(call(z,m))

  opt<-list(...)
  if(!is.null(opt))
    def.plot<-maDotsDefaults(opt,defs$def.plot)

  do.call("maPlot.func", c(list(x=xx,y=yy,z=zz,lines.func=lines.func,text.func=text.func,legend.func=legend.func),def.plot))
  if(y=="maM") abline(h=0,col="gray",lwd=2.5)
}

###########################################################################
# maImage: image methods for microarray objects
# add overlay, NULL which means no overlay or a subset of information
# Modified March 15, 2004 to add overlay option (contribution from katie)
###########################################################################

maImage.func<-function(x, L, subset=TRUE, col=heat.colors(12), contours=FALSE, overlay=NULL, ...)
{
  
  ## x is a matrix
  if(missing(L))
    L <- new(marrayLayout, maNgr=1, maNgc=1, maNsr=nrow(x), maNsc=ncol(x))

  ## When only a subset of spots are stored in marray object, pad with NA
  subset<-maNum2Logic(maNspots(L), subset)
  z<-rep(NA,maNspots(L))
  z[maSub(L)][subset]<-x[subset]
 
  # Create a "full layout"
  Ig<-maNgr(L);  Jg<-maNgc(L);  Is<-maNsr(L);  Js<-maNsc(L)
  L0<-read.marrayLayout(ngr=Ig, ngc=Jg, nsr=Is, nsc=Js)
  nr<-Is*Ig;  nc<-Js*Jg
  row.ind<-(maGridRow(L0)-1)*Is+maSpotRow(L0)
  col.ind<-(maGridCol(L0)-1)*Js+maSpotCol(L0)
  ord<-order(row.ind,col.ind)

  z<-matrix(z[ord],nrow=nr,ncol=nc,byrow=TRUE)
  z<-t(z)[,nr:1]

  # Image of spot statistics
  image(1:nc, 1:nr, z, axes=FALSE, col=col, ...)
  axis(3, at = (0:(Jg-1))*Js + Js/2, labels = 1:Jg)
  axis(2, at = (0:(Ig-1))*Is + Is/2, labels = Ig:1,las=1)
  if(contours)
    contour(1:nc,1:nr,z,add=TRUE, ...)
  box(lwd=4)
  abline(v=((1:Jg-1)*Js + 0.5),lwd=3)
  abline(h=((1:Ig-1)*Is + 0.5),lwd=3)

 # Only for overlay
  if(!is.null(overlay))
    {
      zsub<-rep(NA,maNspots(L))
      subset.over<-maNum2Logic(maNspots(L), overlay)
      zsub[maSub(L)][subset] <- subset.over[subset]
      ## This is subset of subset so we assume that subset.over
      ## (or overlay) is a subset of the subset data a subset on the original

      ## Plotting
      zsub<-matrix(zsub[ord],nrow=nr,ncol=nc,byrow=TRUE)
      zsub<-t(zsub)[,nr:1]
      ## we hard code this part for the moment
      points(row(zsub)[zsub],col(zsub)[zsub], pch=22, col="black", cex=0.9)  
    }
 
}

#########################
# Methods for microarray objects: wrapper around maImage.func

# Methods for microarray objects: wrapper around maImage.func
# Modified by Jean : Sept 15, 2002 to include centering of color
# Modified March 15, 2004 to add overlay option (contribution from katie)

maImage <- function(m, x="maM", subset=TRUE, col, contours=FALSE, bar=TRUE, overlay=NULL, ...)
{
  subset<-maNum2Logic(maNspots(m), subset)
  m<-m[,1]
  if(missing(col))
  {
    col<-rainbow(50)
    if(is.element(x,c("maGb","maGf","maLG")))
      col<-maPalette(low="white", high="green", k=50)
    if(is.element(x,c("maRb","maRf","maLR")))
      col<-maPalette(low="white", high="red", k=50)
    if(is.element(x,c("maM","maMloc","maMscale")))
      col<-maPalette(low="blue", mid="gray", high="yellow", k=50)
    if(is.element(x,c("maA")))
      col<-maPalette(low="white", high="blue", k=50)
  }

  xx<-as.numeric(eval(call(x,m)))


  # Set color range [ensure it's centered]
  tmp<-xx[subset]
  tmp<-tmp[!(is.na(tmp)|is.infinite(tmp))]

  zmax <- ceiling(max(tmp))
  zmin <- floor(min(tmp))
  if(zmin < 0){
    ztmp <- max(abs(zmin), zmax)
    zrange <- c(-ztmp, ztmp)
  }
  else
    zrange <- c(zmin, zmax)

  # Optional graphical parameter defaults
  def<-list(xlab="",ylab="",main=paste(maLabels(maTargets(m)), ": image of ",
                              strsplit(x,"ma")[[1]][2], sep=""), zlim=zrange)
  opt<-list(...)
  if(!is.null(opt))
    def<-maDotsDefaults(opt,def)
  args<-c(list(x=xx, L=maLayout(m), subset=subset, col=col, contours=contours, overlay=overlay),def)
  x.bar <- seq(args$zlim[1], args$zlim[2], length=41)

  if(!bar)
    do.call("maImage.func",args)

  if(bar)
  {
    layout(matrix(c(1,2),1,2),width=c(8,2))
    par(mar=c(4,4,5,3))
    do.call("maImage.func",args)
    par(mar=c(3,0,3,1))
    maColorBar(x.bar,horizontal=FALSE,col=col,main="")
    layout(1)
    par(mar=c(5,4,4,2) + 0.1)
  }

  return(list(x.col=col[1:length(x.bar)], x.bar=x.bar,
              summary=summary(xx[subset])))
}

###########################################################################
# Color bar for calibration
###########################################################################

maPalette <- function(low = "white", high = c("green", "red"), mid=NULL, k =50)
{
    low <- col2rgb(low)/255
    high <- col2rgb(high)/255

    if(is.null(mid)){
        r <- seq(low[1], high[1], len = k)
        g <- seq(low[2], high[2], len = k)
        b <- seq(low[3], high[3], len = k)
      }
    if(!is.null(mid)){
      k2 <- round(k/2)
      mid <- col2rgb(mid)/255
      r <- c(seq(low[1], mid[1], len = k2),
             seq(mid[1], high[1], len = k2))
      g <- c(seq(low[2], mid[2], len = k2),
             seq(mid[2], high[2], len = k2))
      b <- c(seq(low[3], mid[3], len = k2),
             seq(mid[3], high[3], len = k2))
    }
    rgb(r, g, b)
}


maColorBar<-function(x, horizontal = TRUE, col=heat.colors(50),
scale=1:length(x), k=10,  ...)
{
  if(is.numeric(x))
  {
    x <- x
    colmap <- col
  }
  else
  {
    colmap <- x
    low<-range(scale)[1]
    high<-range(scale)[2]
    x <- seq(low, high, length=length(x))
  }

  if(length(x)>k)
    x.small<-seq(x[1], x[length(x)],length=k)
  else
    x.small<-x

  if(horizontal)
  {
    image(x, 1, matrix(x,length(x),1), axes=FALSE, xlab="", ylab="", col=colmap, ...)
    axis(1, at=rev(x.small), labels=signif(rev(x.small),2), srt=270)
  }
  if(!horizontal)
  {
    image(1, x, matrix(x,1,length(x)), axes=FALSE, xlab="", ylab="", col=colmap, ...)
    par(las=1)
    axis(4, at=rev(x.small), labels=signif(rev(x.small), 2))
    par(las=0) # Back to default
  }
  box()
}

###########################################################################
# Functions to filter genes, return a logical vector

maTop<-function(x, h=1, l=1)
{
  x<-as.vector(x)
  return((x>=quantile(x, 1-h, na.rm=TRUE)) | (x<=quantile(x, l, na.rm=TRUE)))
}

###########################################################################
