###########################################################################
## TO BE REMOVE
## Date : October 11, 2002
##
## Modified from Sandrine's Code
##
## source("~/Projects/maTools/R/maStat.R")
##
## Sandrine Test this  
## X <- matrix(rnorm(1000, 10), nc=10)
## Y <- sample(1:2,ncol(X),replace=TRUE)
## maStat(X, funNames=c("bayesFun", "meanFun"))
## maStat(X, funNames=c("bayesFun", "meanFun"), y=Y)
###########################################################################


##################################################################
## Widget Wrapper
##

widget.Stat <- function(expr, outputName="statres", funNames,... )
  {
    LABELFONT <- "Helvetica 12"
    BUTWIDTH <- 10
    BUTTONLAST <- NULL
    CANCEL <- FALSE
    END <- FALSE
    if(missing(funNames))
      FunctionLists <- c("bayesFun", "meanFun", "ttestFun", "numNAFun")
    else
      FunctionLists <- funNames
    
    require(tcltk) || stop("tcltk support is absent")
    
    cancel <- function() {
      CANCEL <<- TRUE
      tkdestroy(base)
    }

    calculate <- function(...)
      {
        newname <- outputName
        funNames <- c()
        for(i in FunctionLists)
          {
            check <- eval(parse(text=tclvalue(i)))
            if(check == '1') funNames <- c(funNames, i)
          }
        res <<- maStat(expr, funNames = funNames, ...)
        write.xls(res, paste(newname, "xls", sep="."))
        assign(newname, res, envir = .GlobalEnv)
        cat(paste("\n Finish calculation, results:", newname, "\n", sep=""))
        tkdestroy(base)
        return()
      }
    
    base <- tktoplevel()
    tkwm.title(base, "Calculation")
    mainfrm <- tkframe(base, borderwidth=2)
    
    ## Buttons
    buttonfr <- tkframe(base)
    for(n in FunctionLists)
      tkpack(buttonfr,
             tkcheckbutton(buttonfr, text=n, variable=n),
             anchor='w')
    tkpack(buttonfr)
        
    butFrame <- tkframe(base)
    cancelBut <- tkbutton(butFrame, text = "Cancel", width = BUTWIDTH, 
                          command = cancel)
    calBut <- tkbutton(butFrame, text = "Calculate", width = BUTWIDTH, 
                       command = calculate)
    tkgrid(calBut, cancelBut)
    tkpack(butFrame)
    tkwait.window(base)
    return(invisible())
  }  

###########################################################################
## Wrapper function
##
##

maStat <- function(expr, funNames, ...)
  {

    switch(data.class(expr),
           exprSet = M <- exprs(expr),
           marrayRaw = M <- maM(expr),
           marrayNorm = M <- maM(expr),
           M <- expr
           )

    opt <- list(...)
    res <- resNames <- c()
    
    if(!class(funNames) == "list")
      {
        print(funNames)
          for(fun in funNames)
          {
            argsfun <- maDotsMatch(opt, formals(args(fun)))
            if(is.null(argsfun))  test <- eval(call(fun))
            if(!is.null(argsfun))  test <- eval(call(fun), argsfun)
            tmp <- test(M)
            ifelse(is.null(colnames(tmp)), tmp2 <- fun, tmp2 <- colnames(tmp))
            resNames <- c(resNames, tmp2)
            res <- cbind(res, tmp)
            colnames(res) <- resNames
          }
      }

    if(class(funNames) == "list")
      for(i in 1:length(funNames))
        {
          tmp <- eval(funNames[[i]])(M)
          ifelse(is.null(colnames(tmp)),
                 tmp2 <- ifelse(is.null(names(funNames)[i]), i, names(funNames)[i]),
                 tmp2 <- colnames(tmp))
          resNames <- c(resNames, tmp2)
          res <- cbind(res, tmp)
          colnames(res) <- resNames
        }
      
    return(res)
  }

###########################################################################
##
## this is filterfunc in genefilter library
## maStatFun <- function (...)
## {
##    flist <- list(...)
##    if (length(flist) == 1 && is.list(flist[[1]]))
##        flist <- flist[[1]]
##    f <- function(x) {
##	fval <- NULL
##        for (fun in flist) {
##            fval <- cbind(fval,fun(x))
##      }
##    }
##    class(f) <- "filterfun"
##    return(f)
## }

###########################################################################
## Functions that calculates various statistics
##
bayesFun <- function(...)
{
  ## take only matrix
  function(M) {
    res <- maBayesian(M, ...)
    return(res$lods)
  }
}

meanFun <- function(y=NULL,
                    na.rm = TRUE)
{
  function(M) {
    meansub <- function(x, y=NULL, na.rm=TRUE)
      {
        if (na.rm) {
          ind <- is.na(x) | is.nan(x) | is.infinite(x)
          x <- x[!ind]
          y <- y[!ind]
        }
        ifelse(is.null(y),
               res <- mean(x),
               res <- diff(unlist(lapply(split(x,y), mean))))
        names(res)<- "Mean"
        return(res)
      }
    switch(data.class(M),
           matrix = apply(M, 1, meansub, y=y, na.rm=na.rm),
           list = unlist(lapply(M, meansub, y=y, na.rm=na.rm)),
           meansub(M, y=y, na.rm=na.rm)
           )
  }
}

ttestFun <- function(y=NULL,
                      var.equal = FALSE,
                      alternative =c("two.sided", "less", "greater"),
                      na.rm = TRUE
                      )
{
  function(M) {
    ttestsub <- function(x, y=NULL,
                      var.equal = FALSE,
                      alternative =c("two.sided", "less", "greater"),
                      na.rm=TRUE)
      {
        if (na.rm) {
          ind <- is.na(x) | is.nan(x) | is.infinite(x)
          x <- x[!ind]
          y <- y[!ind]
        }
        if(length(x) == 0)
          {
            res <- rep(NA, 4)
            names(res) <- c("statistic","estimate","parameter","p.value")
          }
        else
          {
            ifelse(is.null(y),
                   res <- t.test(x, var.equal=var.equal, alternative=alternative),
                   res <- t.test(x ~ y, var.equal=var.equal, alternative=alternative))
            res<-unlist(res[c("statistic","estimate","parameter","p.value")])
            if(is.null(y))
              names(res) <- c("statistic","estimate","parameter","p.value")
            if(!is.null(y))
              names(res) <- c("statistic","estimate group1","estimate group2","parameter","p.value")
          }
	return(res)
      }
    switch(data.class(M),
           matrix = t(apply(M, 1, ttestsub, y=y, var.equal=var.equal,
             alternative=alternative, na.rm==na.rm)),
           list = unlist(lapply(M, ttestsub,  y=y, var.equal=var.equal,
             alternative=alternative, na.rm==na.rm)),
           ttestsub(M,  y=y, var.equal=var.equal,
                    alternative=alternative, na.rm==na.rm)
           )
  }
}
         
numNAFun <- function()
  {
    function(M)
      {
        switch(data.class(M),
               matrix = apply(M, 1, function(x)sum(is.na(x))),
               list = unlist(lapply(M, function(x)sum(is.na(x)))),
               sum(is.na(M))
               )
      }
  }
