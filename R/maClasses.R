##########################################################################
# maClasses.R
#
# Basic class definitions for two-color cDNA microarrays
#
###########################################################################

###########################################################################
# marrayInfo
# Class for describing either samples hybed to array, spots, or other objects
## Jean: Add prototype June-19,03

setClass("marrayInfo",
         representation(maLabels="character",
			maInfo="data.frame",
			maNotes="character"),
         prototype=list(maInfo=data.frame()))

##########################################################################
# marrayLayout
# Class for the microarray layout:  spot coords, grid coords, print-tip group,
# plate, controls (e.g. with expected log-ratio M=0), etc.

setClass("marrayLayout",
         representation(maNgr="numeric",
			maNgc="numeric",
			maNsr="numeric",
                        maNsc="numeric",
			maNspots="numeric",
			maSub="logical",
                        maPlate="factor",
			maControls="factor",
			maNotes="character"),
         prototype=list(maSub=TRUE, maPlate=factor(numeric(0)),
           maControls=factor(numeric(0))))

###########################################################################
# marrayRaw
# Class for pre-normalization intensity data and layout for several arrays

setClass("marrayRaw",
         representation(maRf="matrix",
 			maGf="matrix",
			maRb="matrix",
			maGb="matrix",
                        maW="matrix",
			maLayout="marrayLayout",
                        maGnames="marrayInfo",
			maTargets="marrayInfo",
			maNotes="character"))

###########################################################################
# marrayNorm
# Class used for normalization of intensity data (before and after)
# Use match.call() inside normalization function to return the call

setClass("marrayNorm",
         representation(maA="matrix",
			maM="matrix",
			maMloc="matrix",
			maMscale="matrix",
                        maW="matrix",
			maLayout="marrayLayout",
                        maGnames="marrayInfo",
			maTargets="marrayInfo",
			maNotes="character",
			maNormCall="character")
         )


dim.marrayRaw <- function(x) dim(x@maRf)
dim.marrayNorm <- function(x) dim(x@maM)
dim.marrayInfo <- function(x) dim(x$maInfo)

