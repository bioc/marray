###########################################################################
## Date : October 25, 2002
##
## source("~/Projects/maTools/R/maRankGenes.R")
##
###########################################################################

maSelectGnames <- function(statdata,
                           crit1=50,
                           crit2=crit1,
                           sub=TRUE,
                           selectstat,
                           operate=c("intersect", "union"))
  {
    operate.list <- function(x, operate)
      {
        res <- x[[1]]
        for(i in 2:length(x))
          res <- do.call(operate, list(res, x[[i]]))
        return(res)
      }
    
    gene.ID <- function(x)
      {
        x$gnames
      }

    if(is.integer(sub))
      {
        tmp <- rep(FALSE, nrow(statdata))
        tmp[sub] <- TRUE
        sub <- tmp
      }

    if(missing(selectstat)) selectstat <- 1:ncol(statdata)
    Gnames <- 1:nrow(statdata)
    
    list.id <- list()
    for(i in selectstat)
      {
        switch(data.class(statdata),
                           matrix = headings <- colnames(statdata)[i],
                           data.frame = headings <-  dimnames(statdata)[[2]][i],
                           headings <- colnames(statdata)[i]
                           )
        if(headings == "bayesFun"){
          list.id <- c(list.id,
                       list(gene.ID(stat.gnames(statdata[sub,i],
                                                Gnames[sub], crit=crit1))))
        }
        else
          {
            tmp1 <- gene.ID(stat.gnames(statdata[sub,i], Gnames[sub], crit=crit1)) 
            tmp2 <- gene.ID(stat.gnames(-(statdata[sub,i]), Gnames[sub], crit=crit2)) 
            list.id <- c(list.id, list(c(tmp1, tmp2)))
          }
      }
    finalid <- operate.list(list.id, operate[1])
    return(finalid)
  }



###################################################################
## Select values based on intensities binning
###################################################################
stat.confband.text <-
function (M, A, crit1 = 0.025, crit2 = crit1, nclass = 5)
{
    if (crit1 >= 1)
        crit1 <- crit1/length.na(M)
    if (crit2 >= 1)
        crit2 <- crit2/length.na(M)
    txtA <- (rep(FALSE, length(A)))
    Abin <- quantile.na(A, probs = seq(0, nclass, 1)/nclass)
    for (i in 1:nclass) {
        tmpind <- (Abin[i] <= A) & (A < Abin[i + 1])
        xtmp <- M
        xtmp[!tmpind] <- NA
        n1 <- sum.na(tmpind)
        cutoff <- quantile.na(xtmp, probs = c(crit1, (1 - crit2)))
        vals <- ((xtmp < cutoff[1]) | (xtmp > cutoff[2]))
        txtA[vals] <- TRUE
    }
    res <- c(1:length(txtA))[txtA]
    tmp <- res[rev(order(M[res]))]
    return(tmp)
}


###########################################################################
# Statistics for Microarray Analysis
# Exploratory analysis - Mainly preprocessing.
#
# Date : August 9, 2000
# Last update : May 17, 2001
#
# History:
#   May 17, 2001: Fix to norm.scale.func
#   March, 19: Splitting Rarray in to smaller files.  
#              Including Comments at the start of each function.
#   Nov, 20: Change the argument on plot.mva...it's not usable otherwise.
#            Bug fix ma.func
#   Nov, 13: Ben's Bug fix on stat.ma
#   Nov, 10: Change data structure from matrix to list of matrix.  
#   Sept, 28: Bug fix: ma.func
#
# Authors: Sandrine Dudoit and Yee Hwa (Jean) Yang.
##########################################################################


##########################################################################
#  stat.gnames
#  History:  
#     March 19, 2001:  remove infinite values from the ordering.
#
##########################################################################

stat.gnames<-function(x, gnames, crit=50)
{
    ind <- is.infinite(x)
    x <- x[!ind]
    if (crit < 1) {
        which <- rev(order.na(x, na.last = FALSE))[1:(round(length(x) * 
            crit))]
        if (sum(is.na(x)) > (length(x) - round(length(x) * crit))) 
            warning("NA exists under your selection criteria")
    }
    if (crit >= 1) {
        which <- rev(order.na(x, na.last = FALSE))[1:crit]
        if (sum(is.na(x)) > (length(x) - crit)) 
            warning("NA exists under your selection criteria")
    }
    if (is.matrix(gnames) | is.data.frame(gnames)) 
      {
	gnames <- gnames[!ind, ]
        res <- list(gnames = gnames[which, ], t = x[which])
      }
    if (is.vector(gnames)) 
      {
	gnames <- gnames[!ind]
        res <- list(gnames = gnames[which], t = x[which])
      }
    res
}


##########################################################################
#                                End of file
##########################################################################
