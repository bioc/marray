############################################################################
# maPlot3.R
#
# S4 methods
# Diagnostic plots for two-color cDNA microarrays
#
###########################################################################


  if(!isGeneric("boxplot"))  setGeneric("boxplot")
  setMethod("boxplot", signature(x="marrayRaw"), function (x, xvar = "maPrintTip", yvar = "maM", ...)
            {
              maBoxplot(m=x, x=xvar, y=yvar, ...)
            }
            )
  setMethod("boxplot", signature(x="marrayNorm"), function (x, xvar = "maPrintTip", yvar = "maM", ...)
            {
              maBoxplot(m=x, x=xvar, y=yvar, ...)
            }
            )


  if(!isGeneric("image"))  setGeneric("image")
  setMethod("image", signature(x="marrayRaw"), function (x, xvar = "maM", subset = TRUE, col, contours = FALSE,  bar = TRUE, ...)
            {
              maImage(m=x, x=xvar, subset=subset, col=col, contours=contours, bar=bar, ... )
            }
            )

  setMethod("image", signature(x="marrayNorm"), function (x, xvar = "maM", subset = TRUE, col, contours = FALSE,  bar = TRUE, ...)
            {
              maImage(m=x, x=xvar, subset=subset, col=col, contours=contours, bar=bar, ... )
            }
            )

