## ============================================================================================
##
## predictDiagnosis
##
## ============================================================================================

computePredResMat <- function(predRes, predRes0, followUpTime, followUpTime0) {
  selPeriodStrs <- c("Year 1", "Year 1-2", "Year 3-4", "Year 5-6")
  selPeriodInterval <- matrix(NA, nrow=4, ncol=2)
  selPeriodInterval[1,] <- c(1,365)
  selPeriodInterval[2,] <- c(1,2*365)
  selPeriodInterval[3,] <- c(2*365+1,4*365)
  selPeriodInterval[4,] <- c(4*365+1,6*365)
  
  resMat <- matrix(NA, nrow=5, ncol=5)
  rownames(resMat) <- c("All years", selPeriodStrs)
  colnames(resMat) <- c("FN", "TP", "FP", "TN", "P-value")
  
  for (i in 0:length(selPeriodStrs)) {
    fn <- sum(!predRes>0)
    tp <- sum(predRes>0)
    fp <- sum(!predRes0<0)
    tn <- sum(predRes0<0)
    if (i>0) {
      sel <- followUpTime>selPeriodInterval[i,1] & followUpTime<=selPeriodInterval[i,2]
      sel0 <- followUpTime0>selPeriodInterval[i,1] & followUpTime0<=selPeriodInterval[i,2]
      fn <- sum(!predRes[sel]>0)
      tp <- sum(predRes[sel]>0)
      fp <- sum(!predRes0[sel0]<0)
      tn <- sum(predRes0[sel0]<0)
    }
    
    resMat[i+1,] <-
      c(fn, tp, fp, tn,
        fisher.test(matrix(c(tp, fn, fp, tn), ncol=2, byrow=F),
                    alternative="g")$p.value)
  }

  resMat
}


predictDiagnosis  <- function(dataList, n, nofGenes,  strFollowUpTime, strPlotFracCorr, strResMat,
                              m, m2, selPeriodStr, selPeriodInterval, selPeriodFromWith,
                              hasInsitu=FALSE) {

  ## Predict diagnosis using a leave-one-out approach, n=periodLength
  caseStr <- "Weight"
  caseN <- "N1"
    
  curDataList <- selectDataset(dataList, caseN, dataList$withSpread)
  data <- curDataList$data
  followUpTime <- curDataList$followUpTime
  weights <- curDataList$weights
  curDataList0 <- selectDataset(dataList, caseN, dataList$withoutSpread)
  data0 <- curDataList0$data
  followUpTime0 <- curDataList0$followUpTime
  weights0 <- curDataList0$weights
  if (hasInsitu) {
    curDataListInsitu <- selectDataset(dataList, caseN, dataList$insitu)
    dataInsitu <- curDataListInsitu$data
    followUpTimeInsitu <- curDataListInsitu$followUpTime
    weightsInsitu <- curDataListInsitu$weights
  }
  selFollowUpTime <- followUpTime0
  notSelFollowUpTime <- followUpTime
  if (selPeriodFromWith) {
    selFollowUpTime <- followUpTime
    notSelFollowUpTime <- followUpTime0
  }
  nofPeriods <- length(selFollowUpTime)-n+1
    
  ## Select genes for predictor
  ## Find sorting of genes from analyses for "Weights" in predictor
  ## Weights with sign for genes selected for predictor:  Number of genes selected?
    
  predRes <- matrix(NA, ncol=length(nofGenes), nrow=ncol(data))
  predRes0 <- matrix(NA, ncol=length(nofGenes), nrow=ncol(data0))
  rownames(predRes) <- colnames(data)
  rownames(predRes0)  <- colnames(data0)
  colnames(predRes) <- colnames(predRes0) <- nofGenes
  if (hasInsitu) {
    predResInsitu <- matrix(NA, ncol=length(nofGenes), nrow=ncol(dataInsitu))
    colnames(predResInsitu) <- nofGenes
  }
  
  for (p in 1:nofPeriods)
    print(c(selFollowUpTime[p+n-1]-selFollowUpTime[p],
            sum(notSelFollowUpTime %in% selFollowUpTime[p]:selFollowUpTime[p+n-1])))
  
  meanTimeForPeriod <- rep(NA, nofPeriods)
  for (p in 1:nofPeriods)
    meanTimeForPeriod[p] <- (selFollowUpTime[p+n-1]+selFollowUpTime[p])/2
  
  bestPeriod <- rep(NA, length(followUpTime))
  bestPeriod0 <- rep(NA, length(followUpTime0))
  bestPeriod[followUpTime < mean(meanTimeForPeriod[1:2])] <- 1
  bestPeriod0[followUpTime0 < mean(meanTimeForPeriod[1:2])] <- 1
  bestPeriod[followUpTime > mean(meanTimeForPeriod[(length(meanTimeForPeriod)-1):
                                                   (length(meanTimeForPeriod))])] <- nofPeriods
  bestPeriod0[followUpTime0 > mean(meanTimeForPeriod[(length(meanTimeForPeriod)-1):
                                                     (length(meanTimeForPeriod))])] <- nofPeriods
  if (hasInsitu) {
    bestPeriodInsitu <- rep(NA, length(followUpTimeInsitu))
    bestPeriodInsitu[followUpTimeInsitu < mean(meanTimeForPeriod[1:2])] <- 1
    bestPeriodInsitu[followUpTimeInsitu >
                     mean(meanTimeForPeriod[(length(meanTimeForPeriod)-1):
                                            (length(meanTimeForPeriod))])] <- nofPeriods
  } else {
    bestPeriodInsitu <- c()
  }

  for (p in 2:(nofPeriods-1)) {
    bestPeriod[followUpTime >= mean(meanTimeForPeriod[(p-1):p]) &
               followUpTime <= mean(meanTimeForPeriod[(p):(p+1)])] <- p
    bestPeriod0[followUpTime0 >= mean(meanTimeForPeriod[(p-1):p]) &
                followUpTime0 <= mean(meanTimeForPeriod[(p):(p+1)])] <- p
    if (hasInsitu) {
      bestPeriodInsitu[followUpTimeInsitu >= mean(meanTimeForPeriod[(p-1):p]) &
                       followUpTimeInsitu <= mean(meanTimeForPeriod[(p):(p+1)])] <- p
    }
  }
  
  for (p in unique(sort(c(bestPeriod, bestPeriod0, bestPeriodInsitu)))) {
    selInd <-  (1:length(followUpTime))[followUpTime>=selFollowUpTime[p] &
                                        followUpTime<=selFollowUpTime[p+n-1]]
    selInd0 <- (1:length(followUpTime0))[followUpTime0>=selFollowUpTime[p] &
                                         followUpTime0<=selFollowUpTime[p+n-1]]
    
    leaveOut <- which(bestPeriod==p)
    leaveOut0 <- which(bestPeriod0==p)
      if (hasInsitu) {
        leaveOutInsitu <- which(bestPeriodInsitu==p)
      }
    
    ## Predict for group1
    if (length(leaveOut)>0)
      for (i in leaveOut) { 
        curData <- data[, selInd[selInd!=i]]
        curWeights <- weights[, selInd[selInd!=i]]
        curFollowUpTime <- followUpTime[selInd[selInd!=i]]
        curData0 <- data0[, selInd0]
        curWeights0 <- weights0[, selInd0]
        curFollowUpTime0 <- followUpTime0[selInd0]
        
        mu1 <- computeMean(curData, curWeights)
        mu0 <- computeMean(curData0, curWeights0)
        var1 <- computeVar(curData, curWeights)
        var0 <- computeVar(curData0, curWeights0)
        
        weight <- as.numeric((mu1-mu0)/sqrt(var1 + var0))
        sortList <- sort.list(abs(weight), decreasing=TRUE)
        
        for (j in 1:length(nofGenes)) {
          selGenes <- sortList[1:nofGenes[j]]
          predRes[i, j] <- sum(weight[selGenes]*data[selGenes,i])
        }
      }
    
    ## Predict for group0
    if (length(leaveOut0)>0)
      for (i in leaveOut0) { 
        curData <- data[, selInd] 
        curWeights <- weights[, selInd] 
        curFollowUpTime <- followUpTime[selInd]
        curData0 <- data0[, selInd0[selInd0!=i]]
        curWeights0 <- weights0[, selInd0[selInd0!=i]]
        curFollowUpTime0 <- followUpTime0[selInd0[selInd0!=i]]
        
        mu1 <- computeMean(curData, curWeights)
        mu0 <- computeMean(curData0, curWeights0)
        var1 <- computeVar(curData, curWeights)
        var0 <- computeVar(curData0, curWeights0)
        
        weight <- as.numeric((mu1-mu0)/sqrt(var1+var0))
        sortList <- sort.list(abs(weight), decreasing=TRUE)
        for (j in 1:length(nofGenes)) {
          selGenes <- sortList[1:nofGenes[j]]
          predRes0[i, j] <- sum(weight[selGenes]*data0[selGenes,i])
        }
      }
    
    ## Predict for groupInsitu
    if (hasInsitu) {
      if (length(leaveOutInsitu)>0) {
        curData <- data[, selInd] 
        curWeights <- weights[, selInd] 
        curFollowUpTime <- followUpTime[selInd]
        curData0 <- data0[, selInd0]
        curWeights0 <- weights0[, selInd0]
        curFollowUpTime0 <- followUpTime0[selInd0]
        
        mu1 <- computeMean(curData, curWeights)
        mu0 <- computeMean(curData0, curWeights0)
        var1 <- computeVar(curData, curWeights)
        var0 <- computeVar(curData0, curWeights0)
        
        weight <- as.numeric((mu1-mu0)/sqrt(var1+var0))
        sortList <- sort.list(abs(weight), decreasing=TRUE)
        
        for (i in leaveOutInsitu) { 
          for (j in 1:length(nofGenes)) {
            selGenes <- sortList[1:nofGenes[j]]
            predResInsitu[i, j] <- sum(weight[selGenes]*dataInsitu[selGenes,i])
          }
        }
      }
    }
  }

  for (nGenes in nofGenes) {

    resMat <- computePredResMat(predRes[,paste(nGenes)], predRes0[,paste(nGenes)],
                                followUpTime, followUpTime0)
    write.table(resMat, paste(strResMat, nGenes, ".txt", sep=""), sep="\t")
                
    ## === Print correctly and wrongly classified and p-value ===
    print(paste("Number of genes selected:", nGenes))
    ##
    fn <- sum(!predRes[,paste(nGenes)]>0)
    tp <- sum(predRes[,paste(nGenes)]>0)
    fp <- sum(!predRes0[,paste(nGenes)]<0)
    tn <- sum(predRes0[,paste(nGenes)]<0)
    print("FN TP FP TN:")
    print(c(fn, tp, fp, tn))
    print("p-value in Fisher test, all years:")
    print(fisher.test(matrix(c(tp, fn, fp, tn), ncol=2, byrow=F), alternative="g")$p.value)
    ##      
    sel <- followUpTime>selPeriodInterval[1] & followUpTime<=selPeriodInterval[2]
    sel0 <- followUpTime0>selPeriodInterval[1] & followUpTime0<=selPeriodInterval[2]
    fn <- sum(!predRes[,paste(nGenes)][sel]>0)
    tp <- sum(predRes[,paste(nGenes)][sel]>0)
    fp <- sum(!predRes0[,paste(nGenes)][sel0]<0)
    tn <- sum(predRes0[,paste(nGenes)][sel0]<0)
    print("FN TP FP TN:")
    print(c(fn, tp, fp, tn))
    print(paste("p-value in Fisher test, ", selPeriodStr, ":", sep=""))          
    print(fisher.test(matrix(c(tp, fn, fp, tn), ncol=2, byrow=F), alternative="g")$p.value)
    ##
    print("")
    
    ## === Plot correctly and wrongly classified ===
    plotDataList <- dataList
    plotCorrectPredWithS <- predRes[,paste(nGenes)]>0
    plotCorrectPredWithoutS <- predRes0[,paste(nGenes)]<=0
    ##save(plotDataList, plotCorrectPredWithS, plotCorrectPredWithoutS,
    ##     file=paste(strFollowUpTime, nGenes, ".pdf.R", sep=""))
    plotFollowUpTimesWithPredRes(dataList, predRes[,paste(nGenes)]>0,
                                 predRes0[,paste(nGenes)]<=0, 
                                 paste(strFollowUpTime, nGenes, ".pdf", sep=""))
    
    ## === Plot fraction of correctly classified ===
    
    pdf(paste(strPlotFracCorr, nGenes, ".pdf", sep=""), height=7, width=15)
    
    tVecStart <- min(followUpTime):(max(followUpTime)-m2)
    tVecStop <- tVecStart+m2
    tVecMean <- tVecStart+m
    isCorr <- predRes[,paste(nGenes)]>0
    fracCorr <- rep(NA, length(tVecStart))
    for (i in 1:length(tVecStart))
      fracCorr[i] <- mean(isCorr[followUpTime>=tVecStart[i] & followUpTime<=tVecStop[i]])
    keep <- !is.na(fracCorr)
    ##
    plotX <- tVecMean[keep]
    plotY <- runmed(fracCorr[keep],m2)
    plotXlim <- rev(range(dataList$followUpTime))
    ##save(plotX, plotY, plotXlim, file=paste(strPlotFracCorr, nGenes, ".pdf.red.R", sep=""))
    ##
    plot(tVecMean[keep], runmed(fracCorr[keep],m2), xlim=rev(range(dataList$followUpTime)),
         type="l", ylim =c(0,1), xlab="Follow up time (days)", ylab="Fraction correctly classified",
         cex=1.4, cex.axis=1.4, cex.lab=1.4, col=2)
    
    tVecStart <- min(followUpTime0):(max(followUpTime0)-m2)
    tVecStop <- tVecStart+m2
    tVecMean <- tVecStart+m
    isCorr <- predRes0[,paste(nGenes)]<=0
    fracCorr <- rep(NA, length(tVecStart))
    for (i in 1:length(tVecStart))
      fracCorr[i] <- mean(isCorr[followUpTime0>=tVecStart[i] & followUpTime0<=tVecStop[i]])
    keep <- !is.na(fracCorr)
    ##
    plotX <- tVecMean[keep]
    plotY <- runmed(fracCorr[keep],m2)
    ##save(plotX, plotY, plotXlim, file=paste(strPlotFracCorr, nGenes, ".pdf.black.R", sep=""))
    ##
    lines(tVecMean[keep], runmed(fracCorr[keep], m2),  cex=1.4, cex.axis=1.4, cex.lab=1.4, col=1)
    
    abline(h=0.5, lty=2)
    for (i in 0:7)
      abline(v=365*i)
    
    dev.off()
  }
}

