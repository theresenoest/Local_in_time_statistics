source("tools.r", chdir=T)

## Setup for parallell simulations
nSim <- 1000
values <- 1:nSim
numWorkers <- 8

## Read data 
dataList <- readCancerData()

## Select start and stop times that will be used for selecting time periods
selStartStopTimes <- dataList$followUpTime[dataList$withSpread] 

## === Analyses that compares cases with and without spread ===
caseStr <- "Weight"
caseN <- "N1"
periodLength <- 30
curDataList <- selectDataset(dataList, caseN, dataList$withSpread)
curDataList0 <- selectDataset(dataList, caseN, dataList$withoutSpread)
analyseData("Res/Weight/", selStartStopTimes, periodLength, caseStr,
            curDataList$data, curDataList$followUpTime,
            data0=curDataList0$data, followUpTime0=curDataList0$followUpTime,
            weights=curDataList$weights, weights0=curDataList0$weights)

## === Analyses for case-ctrl pairs with spread ===
caseN <- "N1"
for (caseStr in c("Sd", "Exp1", "Exp2")) {
  periodLength <- 30
  if (caseStr=="Sd")
    periodLength <- 35
  curDataList <- selectDataset(dataList, caseN, dataList$withSpread)
  analyseData(paste("Res/", caseStr, "/", sep=""), selStartStopTimes, periodLength,
              caseStr, curDataList$data, curDataList$followUpTime,
              weights=curDataList$weights)      
}  
  
## ===== Predict diagnosis, make gene list and plots =====
nofGenes <- c(50, 100, 200, 500, 1000, 2000)
periodLength <- 30
nofPeriods <- length(selStartStopTimes)-periodLength+1  
selPeriodStr <- "year 4"
selPeriodInterval <- c(3*365,4*365)
selPeriodFromWith <- TRUE
hasInsitu <- FALSE
m <- 182
m2 <- 365
predictDiagnosis(dataList, periodLength, nofGenes, 
                 "Res/followUpTimesWithPredRes", "Res/plotFracCorr", "Res/predResMat", 
                 m, m2, selPeriodStr, selPeriodInterval, selPeriodFromWith, hasInsitu)

## Make gene list
selYear <- 3.5
makeGeneList(nofPeriods, periodLength, selStartStopTimes, selYear, "Res/genelist.txt")

