##########################################################################
# maGet.R
#
# Accessor methods for two-color cDNA microarrays classes
#
###########################################################################

###########################################################################
# Accessor methods for marrayInfo class

if(!isGeneric("maLabels"))
   setGeneric("maLabels", function(object)
   standardGeneric("maLabels"))

setMethod("maLabels", "marrayInfo", function(object) slot(object, "maLabels"))

if(!isGeneric("maInfo"))
   setGeneric("maInfo", function(object) standardGeneric("maInfo"))
setMethod("maInfo", "marrayInfo", function(object) slot(object, "maInfo"))

if(!isGeneric("maNotes"))
   setGeneric("maNotes", function(object) standardGeneric("maNotes"))
setMethod("maNotes", "marrayInfo", function(object) slot(object,
   "maNotes") )

##########################################################################
# Accessor methods for marrayLayout class

if(!isGeneric("maNspots"))
   setGeneric("maNspots", function(object)
   standardGeneric("maNspots"))

setMethod("maNspots", "marrayLayout",
	function(object) {
		n<-slot(object, "maNspots")
		if(length(n)==0)
		  n<-maNgr(object)*maNgc(object)*maNsr(object)*maNsc(object)
		n
	})

if(!isGeneric("maNgr"))
   setGeneric("maNgr", function(object) standardGeneric("maNgr"))
setMethod("maNgr", "marrayLayout", function(object) slot(object, "maNgr"))

if(!isGeneric("maNgc"))
   setGeneric("maNgc", function(object) standardGeneric("maNgc"))
setMethod("maNgc", "marrayLayout", function(object) slot(object, "maNgc"))

if(!isGeneric("maNsr"))
   setGeneric("maNsr", function(object) standardGeneric("maNsr"))
setMethod("maNsr", "marrayLayout", function(object) slot(object, "maNsr"))

if(!isGeneric("maNsc"))
   setGeneric("maNsc", function(object) standardGeneric("maNsc"))
setMethod("maNsc", "marrayLayout", function(object) slot(object, "maNsc"))

if(!isGeneric("maSub"))
   setGeneric("maSub", function(object) standardGeneric("maSub"))
setMethod("maSub", "marrayLayout", function(object) slot(object, "maSub"))

if(!isGeneric("maPlate"))
   setGeneric("maPlate", function(object) standardGeneric("maPlate"))
setMethod("maPlate", "marrayLayout", function(object)
          {if(length(object@maPlate) == 0)
             maCompPlate(object)
          else
            slot(object, "maPlate")})

if(!isGeneric("maControls"))
   setGeneric("maControls", function(object) standardGeneric("maControls"))
setMethod("maControls", "marrayLayout",
          function(object) slot(object, "maControls"))

if(!isGeneric("maNotes"))
   setGeneric("maNotes", function(object) standardGeneric("maNotes"))
setMethod("maNotes", "marrayLayout", function(object) slot(object, "maNotes"))

#########################
# Methods for quantities that are not slots of marrayLayout

if(!isGeneric("maGridRow"))
   setGeneric("maGridRow", function(object) standardGeneric("maGridRow"))
setMethod("maGridRow", "marrayLayout",
	function(object) {
		ngr<-maNgr(object)
		ngc<-maNgc(object)
		nsr<-maNsr(object)
		nsc<-maNsc(object)
                gr<-numeric(0)
 		if(length(ngr*ngc*nsr*nsc)!=0)
		{
		  gr <- rep(1:ngr, rep(nsr * nsc * ngc, ngr))
		  gr <- gr[as.logical(maSub(object))]
                }
		gr
	})

if(!isGeneric("maGridCol"))
   setGeneric("maGridCol", function(object) standardGeneric("maGridCol"))
setMethod("maGridCol", "marrayLayout",
	function(object) {
		ngr<-maNgr(object)
		ngc<-maNgc(object)
		nsr<-maNsr(object)
		nsc<-maNsc(object)
		gc<-numeric(0)
		if(length(ngr*ngc*nsr*nsc)!=0)
		{
		  gc <- rep(rep(1:ngc, rep(nsr * nsc, ngc)), ngr)
		  gc <- gc[as.logical(maSub(object))]
		}
		gc
   	})


if(!isGeneric("maSpotRow"))
   setGeneric("maSpotRow", function(object) standardGeneric("maSpotRow"))
setMethod("maSpotRow", "marrayLayout",
	function(object) {
		ngr<-maNgr(object)
		ngc<-maNgc(object)
		nsr<-maNsr(object)
		nsc<-maNsc(object)
		sr<-numeric(0)
		if(length(ngr*ngc*nsr*nsc)!=0)
		{
		  sr <- rep(rep(1:nsr, rep(nsc, nsr)), ngc * ngr)
		  sr <- sr[as.logical(maSub(object))]
		}
		sr
	})

if(!isGeneric("maSpotCol"))
   setGeneric("maSpotCol", function(object) standardGeneric("maSpotCol"))
setMethod("maSpotCol", "marrayLayout",
	function(object) {
		ngr<-maNgr(object)
		ngc<-maNgc(object)
		nsr<-maNsr(object)
		nsc<-maNsc(object)
		sc<-numeric(0)
		if(length(ngr*ngc*nsr*nsc)!=0)
		{
  		  sc <- rep(1:nsc, nsr * ngc * ngr)
		  sc <- sc[as.logical(maSub(object))]
		}
		sc
	})


if(!isGeneric("maPrintTip"))
   setGeneric("maPrintTip", function(object) standardGeneric("maPrintTip"))
setMethod("maPrintTip", "marrayLayout",
	function(object) {
		ngr<-maNgr(object)
		ngc<-maNgc(object)
		nsr<-maNsr(object)
		nsc<-maNsc(object)
		pt<-numeric(0)
		if(length(ngr*ngc*nsr*nsc)!=0)
		{
		  gr<- rep(1:ngr, rep(nsr * nsc * ngc, ngr))
		  gc <- rep(rep(1:ngc, rep(nsr * nsc, ngc)), ngr)
		  pt<-(gr-1) * ngc + gc
		  pt<-pt[as.logical(maSub(object))]
		}
		pt
	})


###########################################################################
# Accessor methods for marrayRaw class

if(!isGeneric("maRf"))
   setGeneric("maRf", function(object) standardGeneric("maRf"))
setMethod("maRf", "marrayRaw", function(object) slot(object, "maRf"))

if(!isGeneric("maGf"))
   setGeneric("maGf", function(object) standardGeneric("maGf"))
setMethod("maGf", "marrayRaw", function(object) slot(object, "maGf"))

if(!isGeneric("maRb"))
   setGeneric("maRb", function(object) standardGeneric("maRb"))
setMethod("maRb", "marrayRaw", function(object) slot(object, "maRb"))

if(!isGeneric("maGb"))
   setGeneric("maGb", function(object) standardGeneric("maGb"))
setMethod("maGb", "marrayRaw", function(object) slot(object, "maGb"))

if(!isGeneric("maW"))
   setGeneric("maW", function(object) standardGeneric("maW"))
setMethod("maW", "marrayRaw", function(object) slot(object, "maW"))

#########################
# Methods for quantities that are not slots of marrayRaw

if(!isGeneric("maLR"))
   setGeneric("maLR", function(object) standardGeneric("maLR"))
setMethod("maLR", "marrayRaw",
          function(object)
          {
            r<-maRf(object)
            if(length(maRf(object))!=0)
	    {
              if(length(maRb(object))!=0)
                r<-maRf(object) - maRb(object)
              else
                r<-maRf(object)
              r<-log(ifelse(r>0, r, NA),2)
	    }
      	   r
          }
          )



if(!isGeneric("maLG"))
   setGeneric("maLG", function(object) standardGeneric("maLG"))
setMethod("maLG", "marrayRaw",
          function(object)
          {
            g<-maGf(object)
            if(length(maGf(object))!=0)
	    {
              if(length(maGb(object))!=0)
                g<-maGf(object) - maGb(object)
              else
                g<-maGf(object)
              g<-log(ifelse(g>0, g, NA),2)
	    }
      	   g
          }
          )


if(!isGeneric("maM"))
   setGeneric("maM", function(object) standardGeneric("maM"))
setMethod("maM", "marrayRaw",
	function(object)
	{
           M<-matrix(nr=0,nc=0)
	   if((length(maLR(object))!=0) & (length(maLG(object))!=0))
 	      M<-maLR(object)-maLG(object)
	   M
        }
	)


if(!isGeneric("maA"))
   setGeneric("maA", function(object) standardGeneric("maA"))
setMethod("maA", "marrayRaw",
	function(object)
	{
           A<-matrix(nr=0,nc=0)
           if((length(maLR(object))!=0) & (length(maLG(object))!=0))
             A<-(maLR(object)+maLG(object))/2
	   A
        }
	)


#########################
# maLayout slots

if(!isGeneric("maLayout"))
   setGeneric("maLayout", function(object) standardGeneric("maLayout"))
setMethod("maLayout", "marrayRaw", function(object) slot(object, "maLayout"))

if(!isGeneric("maNspots"))
   setGeneric("maNspots", function(object) standardGeneric("maNspots"))
setMethod("maNspots", "marrayRaw", function(object) maNspots(maLayout(object)))

if(!isGeneric("maNgr"))
   setGeneric("maNgr", function(object) standardGeneric("maNgr"))
setMethod("maNgr", "marrayRaw", function(object) maNgr(maLayout(object)))

if(!isGeneric("maNgc"))
   setGeneric("maNgc", function(object) standardGeneric("maNgc"))
setMethod("maNgc", "marrayRaw", function(object) maNgc(maLayout(object)))

if(!isGeneric("maNsr"))
   setGeneric("maNsr", function(object) standardGeneric("maNsr"))
setMethod("maNsr", "marrayRaw", function(object) maNsr(maLayout(object)))

if(!isGeneric("maNsc"))
   setGeneric("maNsc", function(object) standardGeneric("maNsc"))
setMethod("maNsc", "marrayRaw", function(object) maNsc(maLayout(object)))

if(!isGeneric("maGridRow"))
   setGeneric("maGridRow", function(object) standardGeneric("maGridRow"))
setMethod("maGridRow", "marrayRaw", function(object) maGridRow(maLayout(object)))

if(!isGeneric("maGridCol"))
   setGeneric("maGridCol", function(object) standardGeneric("maGridCol"))
setMethod("maGridCol", "marrayRaw", function(object) maGridCol(maLayout(object)))

if(!isGeneric("maSpotRow"))
   setGeneric("maSpotRow", function(object) standardGeneric("maSpotRow"))
setMethod("maSpotRow", "marrayRaw", function(object) maSpotRow(maLayout(object)))

if(!isGeneric("maSpotCol"))
   setGeneric("maSpotCol", function(object) standardGeneric("maSpotCol"))
setMethod("maSpotCol", "marrayRaw", function(object) maSpotCol(maLayout(object)))

if(!isGeneric("maPrintTip"))
   setGeneric("maPrintTip", function(object) standardGeneric("maPrintTip"))
setMethod("maPrintTip", "marrayRaw", function(object) maPrintTip(maLayout(object)))

if(!isGeneric("maSub"))
   setGeneric("maSub", function(object) standardGeneric("maSub"))
setMethod("maSub", "marrayRaw", function(object) maSub(maLayout(object)))

if(!isGeneric("maPlate"))
   setGeneric("maPlate", function(object) standardGeneric("maPlate"))
setMethod("maPlate", "marrayRaw", function(object) maPlate(maLayout(object)))

if(!isGeneric("maControls"))
   setGeneric("maControls", function(object) standardGeneric("maControls"))
setMethod("maControls", "marrayRaw", function(object) maControls(maLayout(object)))

if(!isGeneric("maNsamples"))
   setGeneric("maNsamples", function(object) standardGeneric("maNsamples"))
setMethod("maNsamples", "marrayRaw", function(object) max(ncol(maRf(object)), ncol(maGf(object))))

#########################
if(!isGeneric("maGnames"))
   setGeneric("maGnames", function(object) standardGeneric("maGnames"))
setMethod("maGnames", "marrayRaw", function(object) slot(object, "maGnames"))

if(!isGeneric("maTargets"))
   setGeneric("maTargets", function(object) standardGeneric("maTargets"))
setMethod("maTargets", "marrayRaw", function(object) slot(object, "maTargets"))

if(!isGeneric("maNotes"))
   setGeneric("maNotes", function(object) standardGeneric("maNotes"))
setMethod("maNotes", "marrayRaw", function(object) slot(object, "maNotes"))

###########################################################################
# Accessor methods for marrayNorm class

if(!isGeneric("maA"))
   setGeneric("maA", function(object) standardGeneric("maA"))
setMethod("maA", "marrayNorm", function(object) slot(object, "maA"))

if(!isGeneric("maM"))
   setGeneric("maM", function(object) standardGeneric("maM"))
setMethod("maM", "marrayNorm", function(object) slot(object, "maM"))

if(!isGeneric("maMloc"))
   setGeneric("maMloc", function(object) standardGeneric("maMloc"))
setMethod("maMloc", "marrayNorm", function(object) slot(object, "maMloc"))

if(!isGeneric("maMscale"))
   setGeneric("maMscale", function(object) standardGeneric("maMscale"))
setMethod("maMscale", "marrayNorm", function(object) slot(object, "maMscale"))

if(!isGeneric("maW"))
   setGeneric("maW", function(object) standardGeneric("maW"))
setMethod("maW", "marrayNorm", function(object) slot(object, "maW"))

##############
# Methods for quantities that are not slots of marrayRaw

if(!isGeneric("maLR"))
   setGeneric("maLR", function(object) standardGeneric("maLR"))
setMethod("maLR", "marrayNorm",
          function(object)
          {
            M <- maM(object)
            A <- maA(object)
            r <- M/2 + A
            r
          }
          )

if(!isGeneric("maLG"))
  setGeneric("maLG", function(object) standardGeneric("maLG"))
setMethod("maLG", "marrayNorm",
          function(object)
          {
            M <- maM(object)
            A <- maA(object)
            g <- A - M/2
            g
          }
          )

#########################
# maLayout slots

if(!isGeneric("maLayout"))
   setGeneric("maLayout", function(object) standardGeneric("maLayout"))
setMethod("maLayout", "marrayNorm", function(object) slot(object, "maLayout"))

if(!isGeneric("maNspots"))
   setGeneric("maNspots", function(object) standardGeneric("maNspots"))
setMethod("maNspots", "marrayNorm", function(object) maNspots(maLayout(object)))

if(!isGeneric("maNgr"))
   setGeneric("maNgr", function(object) standardGeneric("maNgr"))
setMethod("maNgr", "marrayNorm", function(object) maNgr(maLayout(object)))

if(!isGeneric("maNgc"))
   setGeneric("maNgc", function(object) standardGeneric("maNgc"))
setMethod("maNgc", "marrayNorm", function(object) maNgc(maLayout(object)))

if(!isGeneric("maNsr"))
   setGeneric("maNsr", function(object) standardGeneric("maNsr"))
setMethod("maNsr", "marrayNorm", function(object) maNsr(maLayout(object)))

if(!isGeneric("maNsc"))
   setGeneric("maNsc", function(object) standardGeneric("maNsc"))
setMethod("maNsc", "marrayNorm", function(object) maNsc(maLayout(object)))

if(!isGeneric("maGridRow"))
   setGeneric("maGridRow", function(object) standardGeneric("maGridRow"))
setMethod("maGridRow", "marrayNorm", function(object) maGridRow(maLayout(object)))

if(!isGeneric("maGridCol"))
   setGeneric("maGridCol", function(object) standardGeneric("maGridCol"))
setMethod("maGridCol", "marrayNorm", function(object) maGridCol(maLayout(object)))

if(!isGeneric("maSpotRow"))
   setGeneric("maSpotRow", function(object) standardGeneric("maSpotRow"))
setMethod("maSpotRow", "marrayNorm", function(object) maSpotRow(maLayout(object)))

if(!isGeneric("maSpotCol"))
   setGeneric("maSpotCol", function(object) standardGeneric("maSpotCol"))
setMethod("maSpotCol", "marrayNorm", function(object) maSpotCol(maLayout(object)))

if(!isGeneric("maPrintTip"))
   setGeneric("maPrintTip", function(object) standardGeneric("maPrintTip"))
setMethod("maPrintTip", "marrayNorm", function(object) maPrintTip(maLayout(object)))

if(!isGeneric("maSub"))
   setGeneric("maSub", function(object) standardGeneric("maSub"))
setMethod("maSub", "marrayNorm", function(object) maSub(maLayout(object)))

if(!isGeneric("maPlate"))
   setGeneric("maPlate", function(object) standardGeneric("maPlate"))
setMethod("maPlate", "marrayNorm", function(object) maPlate(maLayout(object)))

if(!isGeneric("maControls"))
   setGeneric("maControls", function(object) standardGeneric("maControls"))
setMethod("maControls", "marrayNorm", function(object) maControls(maLayout(object)))

if(!isGeneric("maNsamples"))
   setGeneric("maNsamples", function(object) standardGeneric("maNsamples"))
setMethod("maNsamples", "marrayNorm", function(object) ncol(maM(object)))

#########################
if(!isGeneric("maGnames"))
   setGeneric("maGnames", function(object) standardGeneric("maGnames"))
setMethod("maGnames", "marrayNorm", function(object) slot(object, "maGnames"))

if(!isGeneric("maTargets"))
   setGeneric("maTargets", function(object) standardGeneric("maTargets"))
setMethod("maTargets", "marrayNorm", function(object) slot(object, "maTargets"))

if(!isGeneric("maNotes"))
   setGeneric("maNotes", function(object) standardGeneric("maNotes"))
setMethod("maNotes", "marrayNorm", function(object) slot(object, "maNotes"))

if(!isGeneric("maNormCall"))
   setGeneric("maNormCall", function(object) standardGeneric("maNormCall"))
setMethod("maNormCall", "marrayNorm", function(object) slot(object, "maNormCall"))

