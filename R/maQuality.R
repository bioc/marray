maQuality <- function(mraw, path=".")
  {
    res <- NULL
    for(i in 1:ncol(maM(mraw)))
      res <- cbind(res, unlist(maQualityMain(mraw[,i], path=path)))
    return(res)
  }

maQualityMain <- function(mraw, path=".", fname, output=FALSE)
  {
    data <- mraw[,1]

    if(missing(fname))
      f <- colnames(maGf(data)) 
    tmp <- unlist(strsplit(f, "\\."))
    ifelse(length(maGb(data))!=0, bg <- "Q", bg <- "Q.nbg")
    fstart <- paste(tmp[-length(tmp)], collapse=".")
    fname <- paste(bg, fstart, "xls", sep=".")

    ## info
    Date <- PMT <- NULL
    tmp <- readLines(file.path(path,f), n=40)    
    
    if(length(grep("DateTime", tmp)) != 0)
      Date <- gsub("\"", "",strsplit(tmp[grep("DateTime", tmp)], split="=")[[1]][2])

    if(length(grep("PMT", tmp)) != 0)
      {
        PMT1 <- gsub("\"", "",strsplit(tmp[grep("PMT", tmp)], split="=")[[1]][2])
        PMT <- gsub("\t", ",", PMT1)
      }
    
    ## Flag Info
    if(length(maW(data)) != 0)
      {
        tmp <- table(maW(data))
        Flaginfo <- round(100*tmp/maNspots(data), 2)
      }
    else
      {
        Flaginfo <- NULL
      }

    ## S 2 N
    ifelse(length(maRb(data))!=0 , RS2N <- as.vector(log(maRf(data) / maRb(data),2)),
           RS2N <- as.vector(log(maRf(data),2)))
    RS2Ninfo <- summary(RS2N)[1:6]

    ifelse(length(maGb(data))!=0 , GS2N <- as.vector(log(maGf(data) / maGb(data),2)),
           GS2N <- as.vector(log(maGf(data),2)))
    GS2Ninfo <- summary(GS2N)[1:6]

    ## Control Spots
    A <- split(maA(data),  maControls(data))
    CtlA <- lapply(A, function(x){summary(x)[1:6]})
    CTLNum <- table(maControls(data))
    
    M <- split(maM(data),  maControls(data))
    CtlM <- lapply(M, function(x){summary(x)[1:6]})

    ## Spatial Effects

    ## Reuslts
    res <- list(file = f,
                Date = Date,
                Pmt =  PMT,
                Layout = c(maNspots(data),
                  GridR = maNgr(data),
                  GridC = maNgc(data),
                  SpotR = maNsr(data),
                  SpotC = maNsc(data)),
                Flaginfo = Flaginfo,
                RS2N = RS2Ninfo,
                GS2N = GS2Ninfo,
                CTLNum = CTLNum,
                CtlA = CtlA,
                CtlM = CtlM)
    if(output) write.list(res, file=fname)
    return(res)
  }
