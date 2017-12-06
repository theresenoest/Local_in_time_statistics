require(lumi)
require(lumiHumanIDMapping)
library(multtest)
library(limma)
library(splines)
library(sva)
library(glmnet)
library(glmnetcr)
library(GSVA)
library(combinat)
library(limma)

source("readData.r")
source("analyseData.r")
source("predictDiagnosis.r")

## ============================================================================================
##
## getDate
## - Convert from str on format ddMMMyyyy  to ISOdate
##
## ============================================================================================

getDate <- function(str) {
  year <- as.numeric(substr(str,6,9))
    day <- as.numeric(substr(str,1,2))

  month <- substr(str,3,5)
  if (month=="JAN")
    month <- 1
  else if (month=="FEB")
    month <- 2
  else if (month=="MAR")
    month <- 3
  else if (month=="APR")
    month <- 4
  else if (month=="MAY")
    month <- 5
  else if (month=="JUN")
    month <- 6
  else if (month=="JUL")
    month <- 7
  else if (month=="AUG")
    month <- 8
  else if (month=="SEP")
    month <- 9
  else if (month=="OCT")
    month <- 10
  else if (month=="NOV")
    month <- 11
  else if (month=="DEC")
    month <- 12
  
  ISOdate(year, month, day)
}

## ============================================================================================
##
## selectDataset
##
## ============================================================================================

selectDataset <- function(dataList, caseN, selVec) {

  data <- dataList$dataN1
  if (caseN=="N2")
    data <- dataList$dataN2
  else if (caseN=="N3")
    data <- dataList$dataN3

  data <- data[, selVec]
  followUpTime <- dataList$followUpTime[selVec]
  weights <- c()
  if (length(dataList$runNo)>0) {
    weights <- t(dataList$weight[dataList$runNo[selVec],])
    rownames(weights) <- rownames(data)
  }

  list(data=data,
       followUpTime=followUpTime,
       weights=weights)
}

## ============================================================================================
##
## computeFuncList
##
## ============================================================================================

computeMean <- function(data, weights=c()) {
  if (length(weights)==0) {    
    (1/ncol(data))  * rowSums(data)
  } else {
    (1/rowSums(weights))  * rowSums(data*weights)
  }
}

computeVar <- function(data, weights=c()) {
  if (length(weights)==0) {
    mu <- matrix(rep(computeMean(data),ncol(data)), ncol=ncol(data))
    expr <- (data-mu)
    (1/(ncol(expr)-1)) * rowSums(expr*expr)
  } else {
    mu <- matrix(rep(computeMean(data, weights),ncol(data)), ncol=ncol(data))
    expr <- (data-mu) * weights
    (1/(ncol(expr)-1)) * rowSums(expr*expr)
  }
}

computeSd <- function(data, weights=c()) {
    sqrt(computeVar(data, weights))
}

###

computeFuncList <- function(caseStr) {
  
  ## caseStr=="Sd"
  commonNullDist <- TRUE
  sampleData <- function(data, curData, n, weights=c(), curWeights=c()) {
    ind <- sample(1:ncol(data), n)
    if (length(weights)==0) {
      list(data=data[, ind], weights=c())
    } else {
      list(data=data[, ind], weights=weights[, ind])
    }
  }
  computeSortedStat <- function(data, n, n0, weights=c()) { 
    sort(as.numeric(computeSd(data, weights)))    
  }
  computeNofExtreme <- function(expSd, nullDist) {
    colSums(nullDist<=matrix(rep(expSd,each=nSim), nrow=nSim))
  }
  xlab="genes ordered with respect to standard deviation"
  ylab="standard deviation"
  if (caseStr=="Exp1") {
    commonNullDist <- TRUE
    sampleData <- function(data, curData, n, weights=c(), curWeights=c()) {
      ind <- sample(1:ncol(data), n)
      if (length(weights)==0) { 
        list(data=data[, ind], c())
      } else {
        list(data=data[, ind], weights=weights[, ind])
      }
    }
    computeSortedStat <- function(data, n, n0, weights=c()) {
      sort(abs(as.numeric(computeMean(data, weights))), decreasing=TRUE)
    }
    computeNofExtreme <- function(expSd, nullDist) {
      colSums(nullDist>=matrix(rep(expSd,each=nSim), nrow=nSim))
    }
    xlab="genes ordered with respect to absolute value of expectation"
    ylab="absolute value of expectation"
  } else if (caseStr=="Exp2") {
    commonNullDist <- FALSE
    sampleData <- function(data, curData, n, weights=c(), curWeights=c()) {
      mat <- matrix(rep(sample(c(-1,1), ncol(curData), replace=T), each=nrow(curData)),
                    nrow=nrow(curData))
      list(data=curData * mat, weights=curWeights)
    }
    computeSortedStat <- function(data, n, n0, weights=c()) {
      sort(abs(as.numeric(computeMean(data, weights))), decreasing=TRUE)
    }
    computeNofExtreme <- function(expSd, nullDist) {
      colSums(nullDist>=matrix(rep(expSd,each=nSim), nrow=nSim))
    }
    xlab="genes ordered with respect to absolute value of expectation"
    ylab="absolute value of expectation"
  } else if (caseStr=="ExpSd") {
    commonNullDist <- FALSE
    sampleData <- function(data, curData, n, weights=c(), curWeights=c()) {
      mat <- matrix(rep(sample(c(-1,1), ncol(curData), replace=T), each=nrow(curData)),
                    nrow=nrow(curData))
      list(data=curData * mat, weights=curWeights)
    }
    computeSortedStat <- function(data, n, n0, weights=c()) {
       sort(abs(as.numeric(computeMean(data, weights)/computeSd(data, weights))), decreasing=TRUE)
    }
    computeNofExtreme <- function(expSd, nullDist) {
      colSums(nullDist>=matrix(rep(expSd,each=nSim), nrow=nSim))
    }
    xlab="genes ordered with respect to absolute value of expectation/standard deviation"
    ylab="absolute value of expectation/standard deviation"
  } else if (caseStr=="Weight") {
    commonNullDist <- FALSE
    sampleData <- function(data, curData, n, weights=c(), curWeights=c()) {
      ind <- sample(1:ncol(curData))
      if (length(weights)==0) {
        list(data=curData[, ind], weights=c())
      } else {
        list(data=curData[, ind], weights=curWeights[, ind])
      }
    }
    computeSortedStat <- function(data, n, n0, weights=c()) {
      sort(abs(as.numeric((computeMean(data[,1:n], weights[,1:n])-
                           computeMean(data[,(n+1):(n+n0)], weights[,(n+1):(n+n0)]))/
                          sqrt(computeVar(data[,1:n], weights[,1:n])+
                               computeVar(data[,(n+1):(n+n0)], weights[,(n+1):(n+n0)])))),
           decreasing=TRUE)
    }
    computeNofExtreme <- function(expSd, nullDist) {
      colSums(nullDist>=matrix(rep(expSd,each=nSim), nrow=nSim))
    }
    xlab="genes ordered with respect to absolute value of weight"
    ylab="absolute value of weight"
  }

  list(commonNullDist=commonNullDist,
       sampleData=sampleData,
       computeSortedStat=computeSortedStat,
       computeNofExtreme=computeNofExtreme,
       xlab=xlab,
       ylab=ylab)
}

## ============================================================================================
##
## makeGeneList
##
## ============================================================================================

makeGeneList <- function(nofPeriods, periodLength, startStopTimes, selYear, filename) {

  caseN <- "N1"
  meanTimeForPeriod <- rep(NA, nofPeriods)
  for (p in 1:nofPeriods)
    meanTimeForPeriod[p] <- (startStopTimes[p+periodLength -1]+startStopTimes[p])/2
    
  days <- sort(c((1:7)*365, (0:7)*365+182))
  periodForDay <- rep(NA, length(days))
  periodForDay[days < mean(meanTimeForPeriod[1:2])] <- 1
  periodForDay[days > mean(meanTimeForPeriod[(length(meanTimeForPeriod)-1):(length(meanTimeForPeriod))])] <- nofPeriods
  for (p in 2:(nofPeriods-1))
    periodForDay[days >= mean(meanTimeForPeriod[(p-1):p]) & days <= mean(meanTimeForPeriod[(p):(p+1)])] <- p
  keep <- abs(days-meanTimeForPeriod[periodForDay])<60 ## Ikke avvik på mer enn 2 mnd.
  days <- days[keep]
  periodForYear <- periodForDay[keep]
  years <- round(days/365,1)

  curDataList <- selectDataset(dataList, caseN, dataList$withSpread)
  data <- curDataList$data
  followUpTime <- curDataList$followUpTime
  weights <- curDataList$weights
  
  curDataList0 <- selectDataset(dataList, caseN, dataList$withoutSpread)
  data0 <- curDataList0$data
  followUpTime0 <- curDataList0$followUpTime
  weights0 <- curDataList0$weights

  p <- periodForYear[years==selYear]
  curData <- data[, followUpTime>=startStopTimes[p] & followUpTime<=startStopTimes[p+periodLength-1]]   
  curData0 <- data0[, followUpTime0>=startStopTimes[p] & followUpTime0<=startStopTimes[p+periodLength-1]]
  curWeights <- weights[, followUpTime>=startStopTimes[p] & followUpTime<=startStopTimes[p+periodLength-1]]   
  curWeights0 <- weights0[, followUpTime0>=startStopTimes[p] & followUpTime0<=startStopTimes[p+periodLength-1]]

  genes <- rownames(curData)
  logFC <- as.numeric(computeMean(curData, curWeights)-computeMean(curData0, curWeights0))     
  weightsForScore <-
    logFC/as.numeric(sqrt(computeVar(curData, curWeights) + computeVar(curData0, curWeights0)))

  res <- cbind(genes, round(weightsForScore,2), round(logFC,2))
  sortList <- sort.list(abs(weightsForScore), decreasing=T)
  res <- res[sortList,]
  colnames(res) <- c("gene", "weight", "logFC") 

  write.table(res, filename, row.names=F, quote=F, sep="\t")
}

## ============================================================================================
##
## computeScore
##
## ============================================================================================

computeScore <- function(data, weights) {
  if (length(weights)!=1)
    colSums((matrix(rep(weights, ncol(data)), ncol=ncol(data)))*data, na.rm=TRUE)
  else
    weights*data
}


plotFollowUpTimesWithPredRes <- function(dataList, correctPredWithS, correctPredWithoutS, filename) {
  f <- dataList$followUpTime
  withS <- dataList$withSpread
  withoutS <- dataList$withoutSpread

  pdf(filename, height=3.25, width=15) 
  plot(f[withS][!correctPredWithS], rep(0, length(f[withS][!correctPredWithS])), pch="|",
       xlim=rev(range(dataList$followUpTime)), ylim =c(-2,1), ylab="", yaxt="n",
       xlab="Follow up time (days)", cex=1.4, cex.axis=1.4, cex.lab=1.4)
  abline(v=-0)
  abline(v=365)
  abline(v=365*2)
  abline(v=365*3)
  abline(v=365*4)
  abline(v=365*5)
  abline(v=365*6)
  abline(v=365*7)

  ## Wrong result: red
  points(f[withS][!correctPredWithS], rep(0, length(f[withS][!correctPredWithS])), pch="|", cex=1.4, col=2)
  points(f[withoutS][!correctPredWithoutS], rep(-1.5, length(f[withoutS][!correctPredWithoutS])), pch="|", cex=1.4, col=2)

  ## Correct result: green
  points((f[withS][correctPredWithS]), rep(-0.1, length((f[withS])[correctPredWithS])), pch="|", cex=1.4, col=3)
  points((f[withoutS][correctPredWithoutS]), rep(-1.6, length((f[withoutS])[correctPredWithoutS])), pch="|", cex=1.4, col=3)

  f5 <- sort(f[withS])[(1:1000)*5]
  points(f5, rep(0.5, length(f5)), cex=1.4)

  f5 <- sort(f[withoutS])[(1:1000)*5]
  points(f5, rep(-1, length(f5)), cex=1.4)

  axis(2, at=c(-1.5,0), labels=c("without", "with"), xlab="", cex.axis=1.2)

  dev.off()
}
