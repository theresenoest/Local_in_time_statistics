## ============================================================================================
##
## analyseData
##
## ============================================================================================

analyseData <- function(path, startStopTimes, n, caseStr, data, followUpTime,
                        data0=c(), followUpTime0=c(), weights=c(), weights0=c()) {

  print(paste(path, n))
  set.seed(100)

  funcList <- computeFuncList(caseStr)
  commonNullDist <- funcList$commonNullDist
  sampleData <- funcList$sampleData
  computeSortedStat <- funcList$computeSortedStat
  computeNofExtreme <- funcList$computeNofExtreme
  xlab <- funcList$xlab
  ylab <- funcList$ylab
  
  nofPairsVec <- c()
  nullDistList <- list()
  nullDist <- c()
  nofPeriods <- length(startStopTimes)-n+1

  nofPairs0vec <- rep(0, nofPeriods)
  nofPairs1vec <- rep(0, nofPeriods)
  pVals <- matrix(NA, nrow=nofPeriods, ncol=nrow(data))
  for (p in 1:nofPeriods) {
    print(paste("  ", p, "of", nofPeriods))
    curData <- data[, followUpTime>=startStopTimes[p] &
                    followUpTime<=startStopTimes[p+n-1]]   
    curWeights <- c()
    if (length(weights)>0)
      curWeights <- weights[, followUpTime>=startStopTimes[p] &
                            followUpTime<=startStopTimes[p+n-1]]
    nofPairs <- nofPairs1vec[p] <- ncol(curData)
    nofPairs0 <- 0
    compStat <- TRUE
    if (caseStr=="Weight") {
      if ((sum(followUpTime0>=startStopTimes[p] &
               followUpTime0<=startStopTimes[p+n-1])<3) | nofPairs<3) {
        compStat <- FALSE
        nofPairs0vec[p] <- sum(followUpTime0>=startStopTimes[p] &
                               followUpTime0<=startStopTimes[p+n-1])
      } else {
        curData <- cbind(curData, data0[, followUpTime0>=startStopTimes[p] &
                                        followUpTime0<=startStopTimes[p+n-1]])
        if (length(weights)>0)
          curWeights <- cbind(curWeights, weights0[, followUpTime0>=startStopTimes[p] &
                                                   followUpTime0<=startStopTimes[p+n-1]])
        nofPairs0 <- nofPairs0vec[p] <- ncol(curData)-nofPairs
      }
    }

    if (!compStat) {
      pVals[p, ] <- 1
    } else {
      expSd <- computeSortedStat(curData, nofPairs, nofPairs0, curWeights)
      if (!commonNullDist || length(nofPairsVec)==0 || !(nofPairsVec %in% nofPairs)) {
        ## Compute null distribution
        runif(1)
        mc.reset.stream()
        workerFunc <- function(s) {
          curDlist <- sampleData(data, curData, nofPairs,  weights, curWeights)
          curD <- curDlist$data
          curW <- curDlist$weights
          res <- computeSortedStat(curD, nofPairs, nofPairs0, curW)
          rm(curDlist, curD, curW)
          res
        }
        
        res <- mclapply(values, workerFunc, mc.cores = numWorkers)
        nullDist <- matrix(NA,  nrow=nSim, ncol=nrow(curData))
        for (s in 1:nSim)
          nullDist[s, ] <- res[[s]]
        
        if (commonNullDist) {
          nofPairsVec <- c(nofPairsVec, nofPairs)
          nullDistList[[length(nofPairsVec)]] <- nullDist
        }
      } else { ## commonNullDist && (nofPairs %in% nofPairsVec)
        nullDist <- nullDistList[[which(nofPairsVec %in% nofPairs)]]
      }
      
      nofExtreme <-  computeNofExtreme(expSd, nullDist)
      pVals[p, ] <- nofExtreme/nSim
    }
  }

  write.table(pVals, paste(path, "pVals", n, ".txt", sep=""), col.names=F, row.names=F, quote=F)

  pdf(paste(path, "plotPvals", n, ".pdf", sep=""))
  for (p in 1:nofPeriods) {
    main <- ""
    if (caseStr=="Weight")
      main <- paste("Nof with spread:", nofPairs1vec[p], "   Nof without spread:", nofPairs0vec[p]) 
    plot(pVals[p,], type="l", ylim=c(0,1), log="x", main=main)
    abline(h=0.05, col=2)
  }
  dev.off()
    
  pdf(paste(path, "pVals", n, "AllGenes.pdf", sep=""))
  ind <- 1:ncol(pVals)
  mm <- 4000
  nn <- 100    
  xVals <- 6:(length(pVals[1,ind])-5)   
  plot(xVals, runmed(pVals[1,ind],11)[xVals], type="l", ylim=c(0,1), ylab="p-value",
       xlab=xlab)
  for (p in 1:nofPeriods) 
    lines(xVals, runmed(pVals[p,ind],11)[xVals])
  for (p in (nofPeriods-5):nofPeriods) {
    pV <- runmed(pVals[p,ind],11)[xVals]
    lines(pV,col=2)
    number <- nofPeriods-p+1
    points(xVals[mm+nn*number], pV[mm+nn*number], pch=paste(number))
  }
  for (p in 1:6) {
    pV <- runmed(pVals[p,ind],11)[xVals]
    lines(pV,col=4)
    number <- p
    points(xVals[mm+nn*number], pV[mm+nn*number], pch=paste(number))
  }
  abline(h=0.05, col=3)
  dev.off()
    
  pdf(paste(path, "pVals", n, "2000Genes.pdf", sep=""))
  ind <- 1:2000
  mm <- 300
  nn <- 30    
  xVals <- 6:(length(pVals[1,ind])-5)   
  plot(xVals, runmed(pVals[1,ind],11)[xVals], type="l", ylim=c(0,1), ylab="p-value",
       xlab=xlab)  
  for (p in 1:nofPeriods) 
    lines(xVals, runmed(pVals[p,ind],11)[xVals])
  for (p in (nofPeriods-5):nofPeriods) {
    pV <- runmed(pVals[p,ind],11)[xVals]
    lines(pV,col=2)
    number <- nofPeriods-p+1
    points(xVals[mm+nn*number], pV[mm+nn*number], pch=paste(number))
  }
  for (p in 1:6) {
    pV <- runmed(pVals[p,ind],11)[xVals]
    lines(pV,col=4)
    number <- p
    points(xVals[mm+nn*number], pV[mm+nn*number], pch=paste(number))
  }
  abline(h=0.05, col=3)
  dev.off()

  ## Make plots of p-values over time
  
  ## smooth pVals - median-filter over genes - window l
  pValsSmooth <- pVals
  for (p in 1:nofPeriods)
    pValsSmooth[p,] <- runmed(pVals[p,],99)
  
  if (length(followUpTime0)>length(followUpTime)) {
    meanTimeForPeriod <- rep(NA, nofPeriods)
    for (p in 1:nofPeriods)
      meanTimeForPeriod[p] <- (followUpTime0[p+n-1]+followUpTime0[p])/2
  } else {
    meanTimeForPeriod <- rep(NA, nofPeriods)
    for (p in 1:nofPeriods)
      meanTimeForPeriod[p] <- (followUpTime[p+n-1]+followUpTime[p])/2
  }
  
  days <- min(dataList$followUpTime):max(dataList$followUpTime)
  periodForDay <- rep(NA, length(days))
  periodForDay[days < mean(meanTimeForPeriod[1:2])] <- 1
  periodForDay[days > mean(meanTimeForPeriod[(length(meanTimeForPeriod)-1):
                                             (length(meanTimeForPeriod))])] <- nofPeriods
  for (p in 2:(nofPeriods-1))
    periodForDay[days >= mean(meanTimeForPeriod[(p-1):p]) &
                 days <= mean(meanTimeForPeriod[(p):(p+1)])] <- p
  
  library(caTools)
  
  for (winLen in c(1, 30, 90, 365)) {
    pdf(paste(path, "pValsTime", caseStr, winLen, ".pdf", sep=""), height=6, width=7) 
    
    pValVec <- (pValsSmooth[,50])[periodForDay]
    pValVec <- caTools::runmean(pValVec, winLen)
    plot(days, pValVec, type="l", xlim=rev(range(dataList$followUpTime)), ylim =c(0,1),
         ylab="p-value", xlab="Follow up time (days)", cex=1.0, cex.axis=1.0, cex.lab=1.0)
    text(days[1], pValVec[1], labels=paste(50), cex=0.6)
    col <- 1
    for (r in c(200, 500, 1000, 2000)) {
      col <- col+1
      pValVec <- (pValsSmooth[,r])[periodForDay]
      pValVec <- caTools::runmean(pValVec, winLen)
      lines(days, pValVec, col=col)
      text(days[1], pValVec[1], labels=paste(r), col=col, cex=0.6)
    }
    for (y in 0:7)
      abline(v=y*365)
    abline(h=0.05, lty=2)
    
    dev.off()
  }
}
