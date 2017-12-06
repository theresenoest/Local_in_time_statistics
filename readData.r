## ============================================================================================
##
## readCancerData
##   - Read preprocessed data (log2  normalized data + background information),
##     compute log2 differences in gene expression and extract information from background
##
## Auxilliary functions:
##   mapToGenes
##   makeListWithData
##
## ============================================================================================

mapToGenes <- function(data) {
  ## Map to gene symbol
  nuIDs <- rownames(data)
  mappingInfo <- nuID2RefSeqID(nuIDs, lib.mapping='lumiHumanIDMapping', returnAllInfo=TRUE)
  geneName <- as.character(mappingInfo[,3])
  geneName <- geneName[geneName!=""] 

  ## Compute mean of probes for each gene
  exprs0 <- data
  uGeneNames <- unique(geneName)
  exprs <- matrix(NA, ncol=ncol(exprs0), nrow=length(uGeneNames))
  colnames(exprs) <- colnames(exprs0)
  rownames(exprs) <- uGeneNames
  for (j in 1:length(uGeneNames))
    exprs[j,] <- colMeans(exprs0[mappingInfo[,3]==uGeneNames[j], ,drop=F])
  
  exprs
}

makeListWithData <- function(data, withoutSpread, withSpread, followUpTime, isCase) {

  data <- mapToGenes(data)

  ## Compute log2 differences and standardize data (three different standardizations)
  withoutSpread <- withoutSpread[isCase]
  withSpread <- withSpread[isCase]
  followUpTime <- followUpTime[isCase]
  ## 
  dataN1 <- data[,isCase] - (data[,!isCase]) ## case - ctrl, no standardization
  sd <- apply(dataN1, 1, sd)
  mu <- rowMeans(dataN1)
  dataN2 <- apply(dataN1, 2, function(x){(x-mu)/sd}) ## standardize to mean=0, sd=1
  dataN3 <- apply(dataN1, 2, function(x){x/sd}) ## Assume mean 0, standardize to sd=1
  
  ## Sort data with respect to followUpTime
  sortList <- sort.list(followUpTime)
  withoutSpread <- withoutSpread[sortList]
  withSpread <- withSpread[sortList]
  followUpTime <- followUpTime[sortList]
  dataN1 <- dataN1[,sortList]
  dataN2 <- dataN2[,sortList]
  dataN3 <- dataN3[,sortList]
  
  list(withoutSpread=withoutSpread,
       withSpread=withSpread,
       followUpTime=followUpTime,
       weight=matrix(1,nrow=3, ncol=nrow(dataN1)),
       dataN1=dataN1,
       dataN2=dataN2,
       dataN3=dataN3)
}

readCancerData <- function() {

  ## === Read data and background ===
  load(file="../Data/prospectiveCancerData.R")
  data <- normdata  # normdata  is the gene expression matrix
  ## The backgound matrix:
  ## Rearrange columns in data and rows in background so that case and ctrl of the same pair are adjacent  
  # ind = individuals to include
  data <- data[, ind]
  background <- background[ind, ]
  ## Now the case and control from a case-control pair are found in neighboring columns, i.e. log-differences
  ## are computed as data[, background[,"Case_ctrl"]=="case"] - data[, background[,"Case_ctrl"]=="ctrl"]:

  ## === Find followUpTime ===
  n <- length(background[,"diagnosisdate"])
  followUpTime <- rep(NA, n)
  for (i in 1:n)
    if (background[i,"diagnosisdate"]!="") {
      date1 <- getDate(background[i,"freezing_date"])
      date2 <- getDate(background[i,"diagnosisdate"])
      followUpTime[i] <- round(difftime(date2,date1,units="days") )
    }

  ## === Define strata - withoutSpread and withSpread ===
  withSpread <- background[,"Case_ctrl"]=="case" & background[,"metastasis"]>0
  withoutSpread <- background[,"Case_ctrl"]=="case" & background[,"metastasis"]==0

  ## === Define isCase ===
  isCase <- background[,"Case_ctrl"]=="case"

  makeListWithData(data, withoutSpread, withSpread, followUpTime, isCase)
}
