##########################################################################
# maSet.R
#
# Assignment methods for two-color cDNA microarrays classes
#
###########################################################################

###########################################################################
# Assignment methods for maLabels slot of marrayInfo class

if( !isGeneric("maLabels<-") )
      setGeneric("maLabels<-", function(object, value)
               standardGeneric("maLabels<-"))

setReplaceMethod("maLabels", signature("marrayInfo","character"),
             function(object, value){
               slot(object, "maLabels")<-value
               object
             })

setReplaceMethod("maLabels", signature("marrayInfo", "numeric"),
             function(object, value){
               info <- maInfo(object)
               if(length(value)==1)
		 labels<-info[,value]
	       else
                 labels <- apply(info[,value], 1, paste, collapse=" ")
               slot(object, "maLabels") <- labels
               object
             })

###########################################################################
# Assignment methods for maInfo slot of marrayInfo class

if( !isGeneric("maInfo<-") )
      setGeneric("maInfo<-", function(object, value)
               standardGeneric("maInfo<-"))

setReplaceMethod("maInfo", signature("marrayInfo","data.frame"),
             function(object, value) {
               slot(object, "maInfo")<-value
               object
             })

###########################################################################
# Assignment methods for maNgr slot of marrayLayout class

if( !isGeneric("maNgr<-") )
      setGeneric("maNgr<-", function(object, value)
               standardGeneric("maNgr<-"))

setReplaceMethod("maNgr", signature("marrayLayout","numeric"),
  	function(object, value) {
     	  slot(object,"maNgr")<- value
	  object
	})

setReplaceMethod("maNgr", signature("marrayRaw","numeric"),
	function(object, value) {
	  maNgr(slot(object,"maLayout"))<- value
	  object
	})

setReplaceMethod("maNgr", signature("marrayNorm","numeric"),
  	function(object, value) {
	  maNgr(slot(object,"maLayout"))<- value
	  object
	})
#########################################################################
# Assignment methods for maNgc slot of marrayLayout class

if( !isGeneric("maNgc<-") )
      setGeneric("maNgc<-", function(object, value)
               standardGeneric("maNgc<-"))

setReplaceMethod("maNgc", signature("marrayLayout","numeric"),
  	function(object, value) {
     	  slot(object,"maNgc")<- value
	  object
	})

setReplaceMethod("maNgc", signature("marrayRaw","numeric"),
	function(object, value) {
	  maNgc(slot(object,"maLayout"))<- value
	  object
	})

setReplaceMethod("maNgc", signature("marrayNorm","numeric"),
  	function(object, value) {
	  maNgc(slot(object,"maLayout"))<- value
	  object
	})

###########################################################################
# Assignment methods for maNsr slot of marrayLayout class

if( !isGeneric("maNsr<-") )
      setGeneric("maNsr<-", function(object, value)
               standardGeneric("maNsr<-"))

setReplaceMethod("maNsr", signature("marrayLayout","numeric"),
  	function(object, value) {
     	  slot(object,"maNsr")<- value
	  object
	})

setReplaceMethod("maNsr", signature("marrayRaw","numeric"),
	function(object, value) {
	  maNsr(slot(object,"maLayout"))<- value
	  object
	})

setReplaceMethod("maNsr", signature("marrayNorm","numeric"),
  	function(object, value) {
	  maNsr(slot(object,"maLayout"))<- value
	  object
	})

###########################################################################
# Assignment methods for maNsc slot of marrayLayout class

if( !isGeneric("maNsc<-") )
      setGeneric("maNsc<-", function(object, value)
               standardGeneric("maNsc<-"))

setReplaceMethod("maNsc", signature("marrayLayout","numeric"),
  	function(object, value) {
     	  slot(object,"maNsc")<- value
	  object
	})

setReplaceMethod("maNsc", signature("marrayRaw","numeric"),
	function(object, value) {
	  maNsc(slot(object,"maLayout"))<- value
	  object
	})

setReplaceMethod("maNsc", signature("marrayNorm","numeric"),
  	function(object, value) {
	  maNsc(slot(object,"maLayout"))<- value
	  object
	})

###########################################################################
# Assignment methods for maNspots slot of marrayLayout class

if( !isGeneric("maNspots<-") )
      setGeneric("maNspots<-", function(object, value)
               standardGeneric("maNspots<-"))

setReplaceMethod("maNspots", signature("marrayLayout","numeric"),
  	function(object, value) {
     	  slot(object,"maNspots")<- value
	  object
	})

setReplaceMethod("maNspots", signature("marrayRaw","numeric"),
	function(object, value) {
	  maNspots(slot(object,"maLayout"))<- value
	  object
	})

setReplaceMethod("maNspots", signature("marrayNorm","numeric"),
  	function(object, value) {
	  maNspots(slot(object,"maLayout"))<- value
	  object
	})

###########################################################################
# Assignment methods for maSub slot of marrayLayout class

if( !isGeneric("maSub<-") )
      setGeneric("maSub<-", function(object, value)
               standardGeneric("maSub<-"))

setReplaceMethod("maSub", signature("marrayLayout","logical"),
             function(object, value) {
               slot(object, "maSub")<-value
               object
             })

setReplaceMethod("maSub", signature("marrayLayout","numeric"),
             function(object, value){
               slot(object, "maSub")<-maNum2Logic(maNspots(object),value)
               object
             })

setReplaceMethod("maSub", signature("marrayRaw"),
             function(object,value) {
                maSub(slot(object, "maLayout")) <- value
                object
             })

setReplaceMethod("maSub", signature("marrayNorm"),
             function(object,value){
                maSub(slot(object, "maLayout")) <- value
                object
             })

###########################################################################
# Assignment methods for maPlate slot of marrayLayout class

if( !isGeneric("maPlate<-") )
      setGeneric("maPlate<-", function(object, value)
               standardGeneric("maPlate<-"))

setReplaceMethod("maPlate", "marrayLayout",
  	function(object, value) {
     	  slot(object,"maPlate")<- factor(value)
	  object
	})

setReplaceMethod("maPlate", "marrayRaw",
	function(object, value) {
	  maPlate(slot(object,"maLayout"))<- factor(value)
	  object
	})

setReplaceMethod("maPlate", "marrayNorm",
  	function(object, value) {
	  maPlate(slot(object,"maLayout"))<- factor(value)
	  object
	})


###########################################################################
# Assignment methods for maControls slot of marrayLayout class

if( !isGeneric("maControls<-") )
      setGeneric("maControls<-", function(object, value)
               standardGeneric("maControls<-"))

setReplaceMethod("maControls", "marrayLayout",
  function(object, value) {
     slot(object,"maControls")<- factor(value)
     object
  })

setReplaceMethod("maControls", "marrayRaw",
  function(object, value) {
     maControls(slot(object,"maLayout"))<- factor(value)
     object
  })

setReplaceMethod("maControls", "marrayNorm",
  function(object, value) {
     maControls(slot(object,"maLayout")) <- factor(value)
     object
  })


###########################################################################
# Assignment methods for maLayout slot of marrayRaw, marrayNorm etc. classes

if( !isGeneric("maLayout<-") )
      setGeneric("maLayout<-", function(object, value)
               standardGeneric("maLayout<-"))

setReplaceMethod("maLayout", signature("marrayRaw", "marrayLayout"),
  function(object, value) {
     slot(object,"maLayout")<- value
     object
  })

setReplaceMethod("maLayout", signature("marrayNorm", "marrayLayout"),
  function(object, value) {
     slot(object,"maLayout")<- value
     object
  })

###########################################################################
# Assignment methods for maGnames slot of marrayRaw, marrayNorm

if( !isGeneric("maGnames<-") )
      setGeneric("maGnames<-", function(object, value)
               standardGeneric("maGnames<-"))

setReplaceMethod("maGnames", signature("marrayRaw", "marrayInfo"),
  function(object, value) {
     slot(object,"maGnames")<- value
     object
  })

setReplaceMethod("maGnames", signature("marrayNorm", "marrayInfo"),
  function(object, value) {
     slot(object,"maGnames")<- value
     object
  })

###########################################################################
# Assignment methods for maTargets slot of marrayRaw, marrayNorm

if( !isGeneric("maTargets<-") )
      setGeneric("maTargets<-", function(object, value)
               standardGeneric("maTargets<-"))

setReplaceMethod("maTargets", signature("marrayRaw", "marrayInfo"),
  function(object, value) {
     slot(object,"maTargets")<- value
     object
  })

setReplaceMethod("maTargets", signature("marrayNorm", "marrayInfo"),
  function(object, value) {
     slot(object,"maTargets")<- value
     object
  })


###########################################################################
# Assignment methods for maW slot of marrayRaw and marrayNorm classes

if( !isGeneric("maW<-") )
      setGeneric("maW<-", function(object, value)
               standardGeneric("maW<-"))

setReplaceMethod("maW", signature("marrayRaw", "matrix"),
  function(object, value) {
     slot(object,"maW")<- value
     object
  })

setReplaceMethod("maW", signature("marrayNorm", "matrix"),
  function(object, value) {
     slot(object,"maW")<- value
     object
  })

###########################################################################
# Assignment methods for maGf slot of marrayRaw class

if( !isGeneric("maGf<-") )
      setGeneric("maGf<-", function(object, value)
               standardGeneric("maGf<-"))

setReplaceMethod("maGf", signature("marrayRaw", "matrix"),
  function(object, value) {
     slot(object,"maGf")<- value
     object
  })

###########################################################################
# Assignment methods for maRf slot of marrayRaw class

if( !isGeneric("maRf<-") )
      setGeneric("maRf<-", function(object, value)
               standardGeneric("maRf<-"))

setReplaceMethod("maRf", signature("marrayRaw", "matrix"),
  function(object, value) {
     slot(object,"maRf")<- value
     object
  })

###########################################################################
# Assignment methods for maGb slot of marrayRaw class

if( !isGeneric("maGb<-") )
      setGeneric("maGb<-", function(object, value)
               standardGeneric("maGb<-"))

setReplaceMethod("maGb", signature("marrayRaw", "matrix"),
  function(object, value) {
     slot(object,"maGb")<- value
     object
  })

setReplaceMethod("maGb", signature("marrayRaw", "NULL"),
  function(object, value) {
     slot(object,"maGb")<- value
     object
  })

###########################################################################
# Assignment methods for maRb slot of marrayRaw class

if( !isGeneric("maRb<-") )
      setGeneric("maRb<-", function(object, value)
               standardGeneric("maRb<-"))

setReplaceMethod("maRb", signature("marrayRaw", "matrix"),
  function(object, value) {
     slot(object,"maRb")<- value
     object
  })

setReplaceMethod("maRb", signature("marrayRaw", "NULL"),
  function(object, value) {
     slot(object,"maRb")<- value
     object
  })

###########################################################################
# Assignment methods for maA slot of marrayNorm class

if( !isGeneric("maA<-") )
      setGeneric("maA<-", function(object, value)
               standardGeneric("maA<-"))

setReplaceMethod("maA", signature("marrayNorm", "matrix"),
  function(object, value) {
     slot(object,"maA")<- value
     object
  })
###########################################################################
# Assignment methods for maM slot of marrayNorm class

if( !isGeneric("maM<-") )
      setGeneric("maM<-", function(object, value)
               standardGeneric("maM<-"))

setReplaceMethod("maM", signature("marrayNorm", "matrix"),
  function(object, value) {
     slot(object,"maM")<- value
     object
  })

###########################################################################
# Assignment methods for maMloc slot of marrayNorm class

if( !isGeneric("maMloc<-") )
      setGeneric("maMloc<-", function(object, value)
               standardGeneric("maMloc<-"))

setReplaceMethod("maMloc", signature("marrayNorm", "matrix"),
  function(object, value) {
     slot(object,"maMloc")<- value
     object
  })

###########################################################################
# Assignment methods for maMscale slot of marrayNorm class

if( !isGeneric("maMscale<-") )
      setGeneric("maMscale<-", function(object, value)
               standardGeneric("maMscale<-"))

setReplaceMethod("maMscale", signature("marrayNorm", "matrix"),
  function(object, value) {
     slot(object,"maMscale")<- value
     object
  })

###########################################################################
# Assignment methods for maNotes slot of various classes

if( !isGeneric("maNotes<-") )
      setGeneric("maNotes<-", function(object, value)
               standardGeneric("maNotes<-"))

setReplaceMethod("maNotes", signature("marrayLayout", "character"),
  function(object, value) {
     slot(object,"maNotes")<- value
     object
  })

setReplaceMethod("maNotes", signature("marrayRaw", "character"),
  function(object, value) {
     slot(object,"maNotes")<- value
     object
  })

setReplaceMethod("maNotes", signature("marrayNorm", "character"),
  function(object, value) {
     slot(object,"maNotes")<- value
     object
  })

setReplaceMethod("maNotes", signature("marrayInfo", "character"),
  function(object, value) {
     slot(object,"maNotes")<- value
     object
  })
