###########################################################################
# maNorm.R
#
# Functions for location and scale normalization for two-color cDNA microarrays
#  source("~/Projects/madman/Rpacks/marray/R/maNorm.R")
###########################################################################

############################################################################
# maNormMain
# Main within-slide location and scale normalization function

maNormMain<-function(mbatch,
                     f.loc=list(maNormLoess()),
                     f.scale=NULL,
                     a.loc=maCompNormEq(),
                     a.scale=maCompNormEq(),
                     Mloc=TRUE, Mscale=TRUE, echo=FALSE)
{

    if (is(mbatch, "marrayRaw"))
    {
        mnorm<-as(mbatch, "marrayNorm")
        M<-Ml<-Ms<-NULL
    }
    if (is(mbatch, "marrayNorm"))
        mnorm<-mbatch

  slot(mnorm, "maNormCall") <- match.call()

  if(length(f.loc)>0)
    M<-Ml<-NULL
  if(length(f.scale)>0)
    M<-Ms<-NULL

  for(i in 1:ncol(maM(mbatch)))
  {

    if(echo)
      cat(paste("Normalizing array ", i, ".\n", sep=""))
    m<-mbatch[,i]
    M1<-M2<-NULL
    Mnorm<-maM(m)

    # Location
    if(length(f.loc)>0)
    {
      for(func in f.loc)
        M1<-cbind(M1, func(m))

      if(length(f.loc)>1)
          M1 <- rowSums(M1*a.loc(maA(m),length(f.loc)))
      Ml<-cbind(Ml,M1)
      Mnorm<-(Mnorm - M1)
    }


    # Scale
    if(length(f.scale)>0)
    {
      # Scale computed for location normalized data
      m<-mnorm[,i]
      slot(m, "maM")<-Mnorm
      for(func in f.scale)
        M2<-cbind(M2, func(m))

      if(length(f.scale)>1)
          M2 <- rowSums(M2*a.scale(maA(m),length(f.scale)))
      Ms<-cbind(Ms,M2)
      Mnorm<-Mnorm/M2
    }

    M<-cbind(M, Mnorm)
  }

  slot(mnorm, "maM")<-M
  if(length(f.loc)>0 & Mloc)
    slot(mnorm, "maMloc")<-Ml
  if(length(f.scale)>0 & Mscale)
    slot(mnorm, "maMscale")<-Ms

  return(mnorm)
}

############################################################################
# maNorm
# Wrapper function around maNormMain for simple normalization procedures
#
# From Gordon's Mar 17 E-mail
# loess(y~x,span=0.4,degree=1,family="symmetric",control=loess.control(trace.hat="ap
# proximate",iterations=5,surface="direct"))
#

maNorm<-function(mbatch,
	norm=c("printTipLoess", "none", "median", "loess", "twoD", "scalePrintTipMAD"),
	subset=TRUE, span=0.4, Mloc=TRUE, Mscale=TRUE, echo=FALSE, ...)
{
  ## Set normalization defaults
  opt <- list(...)

  ## norm <- unlist(strsplit(norm, split=""))[1]
  norm.method <- match.arg(norm)
  if(echo)
    cat(paste("Normalization method: ", norm.method, ".\n", sep=""))

  switch(norm.method,
         none = maNormMain(mbatch, f.loc=NULL,
           Mloc=Mloc, Mscale=Mscale, echo=echo),

         median = maNormMain(mbatch,
           f.loc=list(maNormMed(x=NULL, y="maM", subset=subset)),
           Mloc=Mloc, Mscale=Mscale, echo=echo),

         loess = maNormMain(mbatch,
           f.loc=list(maNormLoess(x="maA", y="maM", z=NULL, w=NULL,
             subset=subset, span=span, ...)),
           Mloc=Mloc, Mscale=Mscale, echo=echo),

         twoD = maNormMain(mbatch,
           f.loc=list(maNorm2D(x = "maSpotRow", y = "maSpotCol",
             z = "maM", g = "maPrintTip",
             w = NULL, subset = subset, span = span, ...)),
           Mloc=Mloc, Mscale=Mscale, echo=echo),

         printTipLoess = maNormMain(mbatch,
           f.loc=list(maNormLoess(x="maA", y="maM", z="maPrintTip",
             w=NULL, subset=subset, span=span, ...)),
           Mloc=Mloc, Mscale=Mscale, echo=echo),

         scalePrintTipMAD = maNormMain(mbatch,
           f.loc=list(maNormLoess(x="maA", y="maM", z="maPrintTip",
             w=NULL, subset=subset, span=span, ...)),
           f.scale=list(maNormMAD(x="maPrintTip", y="maM", geo=TRUE,
             subset=subset)),
           Mloc=Mloc, Mscale=Mscale, echo=echo))
}

############################################################################
# maNormScale
# Wrapper function around maNormMain for simple scale normalization procedures

maNormScale<-function(mbatch,
		norm=c("globalMAD", "printTipMAD"),
		subset=TRUE, geo=TRUE, Mscale=TRUE, echo=FALSE)
{
  ## norm <- unlist(strsplit(norm, split=""))[1]

  norm.method <- match.arg(norm)
  if(echo)
    cat(paste("Normalization method: ", norm.method, ".\n", sep=""))

  if(norm.method == "globalMAD")
  {
     # Not a very nice way to do between slide norm
     if(!geo)
       mnorm<-maNormMain(mbatch,
          f.loc=NULL,
	  f.scale=list(maNormMAD(x=NULL, y="maM", geo=geo, subset=subset)),
	  Mscale=Mscale, echo=echo)
     if(geo)
     {
       mnorm<-maNormMain(mbatch,
         f.loc=NULL,
	 f.scale=list(maNormMAD(x=NULL, y="maM", geo=FALSE, subset=subset)),
	 Mscale=TRUE, echo=echo)
         #mgeo<-apply(maMscale(mnorm), 1, function(z) exp(mean(log(z))))
         mgeo<- exp(mean(log(maMscale(mnorm)[1,])))

       maM(mnorm)<-maM(mnorm)*mgeo
       if(Mscale)
         maMscale(mnorm)<-maMscale(mnorm)/mgeo
       if(!Mscale)
         maMscale(mnorm)<-NULL
     }
  }

  if(norm.method == "printTipMAD")
     mnorm<-maNormMain(mbatch,
        f.loc=NULL,
	f.scale=list(maNormMAD(x="maPrintTip", y="maM", geo=geo,
		subset=subset)),
	Mscale=Mscale, echo=echo)

  return(mnorm)
}

###########################################################################
# Median location normalization
###########################################################################
# maMed
# General function: median of y within values of x

maMed<-function(x, y, subset=TRUE)
{
  subset<-maNum2Logic(length(y), subset)
  yfit<-rep(NA,length(y))
  for(i in unique(x))
    yfit[x==i]<-median(y[(x==i)&subset], na.rm=TRUE)
  return(yfit)
}

#########################
# maNormMed
# Function for objects of class marrayRaw

maNormMed<-function(x=NULL, y="maM", subset=TRUE)
{
  function(m)
  {
    if(is.character(x))
      maMed(x=eval(call(x,m)), y=eval(call(y,m)), subset=subset)
    else
      maMed(x=TRUE, y=eval(call(y,m)), subset=subset)
  }
}

###########################################################################
# Loess location normalization
###########################################################################
# maLoess
# General function: regress y on x within values of z
# In our case, x=NA iff y=NA
## Jean modify Nov 5th, April 9, 2003
## Very ugly fix
## Jean modify Sep 14, 2004
## Sample Loess if the number is too large

maLoess <-
function (x, y, z, w = NULL, subset = TRUE, span = 0.4, ...)
{
  opt <- list(...)
  subset <- maNum2Logic(length(y), subset)
  yfit <- rep(NA, length(y))
  good <- !(is.infinite(x) | is.infinite(y) | is.na(x) | is.na(y))
  if(length(z) == 1) info <- z
  else
    info <- z[subset]
  for (i in unique(info))
    {
      which <- z == i
      ##      cat(i, " :: ", sum(which), "::", unique(z), "\n")
      if((sum(which) == 1) | (length(unique(z)) == 1)) {
        samplesub <- rep(TRUE, length(x))
        samplesub[-sample(1:length(x), min(length(x), 2000))] <- FALSE
        ## cat("in if ", sum(samplesub), " :: ", length(samplesub), "\n")
      }
      else
        samplesub <- TRUE

      defs <- list(degree=1,
                   family="symmetric",
                   control=loess.control(trace.hat="approximate",iterations=5, surface="direct"),
                   subset = samplesub & which & subset & good,
                   span=span,
                   na.action = na.omit,
                   weights=w)
      args <- maDotsMatch(c(defs, opt), formals(args("loess")))
      fit <- loess(y ~ x,  weights = args$w, subset = args$subset,
                   span = args$span, na.action = args$na.action,
                   degree = args$degree, family=args$family,
                   control=args$control)
      ##      gc(TRUE)
      ##      fit <- loess(y ~ x, weights = w, subset = which & subset &
      ##                   good, span = span, na.action = na.omit)
      yfit[which & good] <- predict(fit, x[which & good])
    }
    return(yfit)
}


#########################
# maNormLoess
# Function for objects of class marrayRaw

maNormLoess<-function(x="maA", y="maM", z="maPrintTip", w=NULL, subset=TRUE, span=0.4, ...)
{
  function(m)
  {
    if(is.character(z))
      maLoess(x=eval(call(x,m)), y=eval(call(y,m)), z=eval(call(z,m)), w=w, subset=subset, span=span, ...)
    else
      maLoess(x=eval(call(x,m)), y=eval(call(y,m)), z=TRUE, w=w, subset=subset, span=span, ...)
  }
}

###########################################################################
# 2D spatial location normalization
###########################################################################
# ma2D
# General function: regress z on x*y within values of g

ma2D<-function(x, y, z, g, w=NULL, subset=TRUE, span=0.4, ...)
{
  opt <- list(...)
  subset<-maNum2Logic(length(y), subset)
  zfit<-rep(NA,length(z))
  # Deal with NA, NaN, Inf
  good <- !(is.infinite(z) | is.na(z))

  for(i in unique(g))
  {
    which<-g==i
    defs <- list(degree=1,
                 family="symmetric",
                 control=loess.control(trace.hat="approximate",iterations=5, surface="direct"),
                 subset = which & subset & good,
                 span=span,
                 na.action = na.omit,
                 weights=w)
    args <- maDotsMatch(c(defs, opt), formals(args("loess")))
    fit <- loess(z ~ x * y,  weights = args$w, subset = args$subset,
                 span = args$span, na.action = args$na.action,
                 degree = args$degree, family=args$family,
                 control=args$control)
    ##    gc(TRUE)
    ##    fit<-loess(z ~ x * y, weights=w, subset=subset&which&good, span=span, na.action=na.omit, ...)
    zfit[which]<-predict(fit, data.frame(x=x[which],y=y[which]))
   }
   return(zfit)
}

#########################
# maNorm2D
# Function for objects of class marrayRaw

maNorm2D<-function(x="maSpotRow", y="maSpotCol", z="maM", g="maPrintTip", w=NULL, subset=TRUE, span=0.4, ...)
{
  function(m)
    ma2D(x=eval(call(x,m)), y=eval(call(y,m)), z=eval(call(z,m)), g=eval(call(g,m)), w=w, subset=subset, span=span, ...)
}

###########################################################################
# MAD scale normalization
###########################################################################
# maMAD
# General function: MAD of y within values of x

maMAD<-function(x, y, geo=TRUE, subset=TRUE)
{
  subset<-maNum2Logic(length(y), subset)
  yfit<-rep(NA,length(y))
  m<-rep(NA,length(unique(x)))
  for(i in unique(x))
  {
    m[i]<-mad(y[(x==i)&subset], na.rm=TRUE)
    yfit[x==i]<-m[i]
  }
  if(geo)
   yfit<-yfit/exp(mean(log(m)))
  return(yfit)
}

#########################
# maNormMAD
# Function for objects of class marrayRaw

maNormMAD<-function(x=NULL, y="maM", geo=TRUE, subset=TRUE)
{
  function(m)
  {
    if(is.character(x))
      maMAD(x=eval(call(x,m)), y=eval(call(y,m)), geo=geo, subset=subset)
    else
      maMAD(x=TRUE, y=eval(call(y,m)), geo=geo, subset=subset)
  }
}

###########################################################################
# Weights for composite normalization
###########################################################################
# maCompNormEq
# Equal weights for composite normalization

maCompNormEq<-function()
{
  function(x,n)
  {
    return(matrix(1/n,length(x),n))
  }
}

#########################
# maCompNormA
# Weights for composite intensity dependent normalization

maCompNormA<-function()
{
  function(x,n)
  {
    if(n!=2)
      stop("Can only do composite location normalization for 2 functions in f.loc")
    f<-ecdf(x)
    a<-f(x)
    return(as.matrix(cbind(a,1-a)))
  }
}

###########################################################################
