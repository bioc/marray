###########################################################################
## Date : October 11, 2002
##
## 
###########################################################################

maDotsMatch <- function(dots, defaults)
  {
    ind <- intersect(names(dots), setdiff(names(defaults), "..."))
    for(i in ind)
      defaults[[i]] <- dots[[i]]
    return(defaults[!names(defaults)=="..."])
  }
