
library(randomForest)
library(plyr)
library(doParallel)
library(vegan)
library(ggmap)
library(reshape2)
library(RSQLite)
library(stringr)

load("data/base_map.rdata")
lazyLoad("data/forestsdb")
source("r/model.predict.RanFor.4.2.r")

setOldClass("randomForest.formula")
setOldClass("idf")

###Class Definition###
setClass("bugs", representation(bugdata="data.frame",
                                predictors="data.frame",
                                dbconn="SQLiteConnection"),
         prototype=list(dbconn = dbConnect("SQLite", "data/bug_metadata.db")))
setClass("mmi", representation(subsample = "data.frame",
                               metrics = "data.frame",
                               modelprediction = "data.frame",
                               result = "data.frame",
                               finalscore = "data.frame",
                               datalength = "numeric",
                               summary ="data.frame"),
         contains="bugs",
         prototype = list(subsample = data.frame(),
                          metrics = data.frame(),
                          residuals = data.frame(),
                          result = data.frame(),
                          finalscore = data.frame(),
                          datalength = numeric(),
                          summary = data.frame()
         )
)
setClass("oe", representation(ambiguous="data.frame",
                              oesubsample="data.frame",
                              iterations="matrix",
                              fulliterations="list",
                              datalength="numeric",
                              oeresults="data.frame"), 
         contains="bugs",
         prototype=list(ambiguous=data.frame(),
                        oesubsample=data.frame(),
                        datalength=numeric(),
                        iterations=matrix(),
                        fulliterations=list(),
                        datalength=numeric(),
                        oeresults=data.frame()
                        ))

setClass("metric.mean", representation(mean.metric="data.frame"),
         contains=c("mmi", "oe"))
###Validity Checks###
setValidity("bugs", function(object){
  bugcolumns <- c("StationCode", "SampleID", "FinalID", "BAResult", "DistinctCode")
  predictorcolumns <- c("StationCode", "New_Lat",     "New_Long",    "ELEV_RANGE",  "BDH_AVE",     "PPT_00_09",  
                  "LPREM_mean",  "KFCT_AVE",    "TEMP_00_09",  "P_MEAN",      "N_MEAN",      "PRMH_AVE",   
                  "SITE_ELEV",   "MgO_Mean",    "S_Mean",      "SumAve_P",   
                  "CaO_Mean")
  for(i in 1:5){
    if(!(bugcolumns[i] %in% names(object@bugdata)))
      return(paste("'bugdata' missing column:", bugcolumns[i]))
    if(sum(is.na(object@bugdata[, bugcolumns[i]])) != 0 & i != 5)
      return(paste("Missing data in column", bugcolumns[i], "of 'bugdata'"))
  }
  for(i in 1:17){
    if(!(predictorcolumns[i] %in% names(object@predictors)))
      return(paste("'predictors' missing column:", predictorcolumns[i]))
    if(sum(is.na(object@predictors[, predictorcolumns[i]])) != 0)
      return(paste("Missing data in column", predictorcolumns[i], "of 'bugdata'"))
    }
  if(!("AREA_SQKM" %in% names(object@predictors)) & !("LogWSA" %in% names(object@predictors)))
    return("Predictors must include a column AREA_SQKM or LogWSA")
  if(!setequal(object@bugdata$StationCode, object@predictors$StationCode))
    return("All StationCode IDs must be represented in both bug and predictor data")
  if(length(unique(object@bugdata$SampleID)) != nrow(unique(object@bugdata[, c("StationCode", "SampleID")])))
    return("SampleIDs must be unique to one StationCode")
  TRUE
})
setValidity("mmi", function(object){
  bugcolumns <- c("StationCode", "SampleID", "FinalID", "BAResult", "DistinctCode")
  predictorcolumns <- c("StationCode", "New_Lat",     "New_Long",    "ELEV_RANGE",  "BDH_AVE",     "PPT_00_09",  
                        "LPREM_mean",  "KFCT_AVE",    "TEMP_00_09",  "P_MEAN",      "N_MEAN",      "PRMH_AVE",   
                        "SITE_ELEV",   "MgO_Mean",    "S_Mean",      "SumAve_P",   
                        "CaO_Mean")
  for(i in 1:5){
    if(!(bugcolumns[i] %in% names(object@bugdata)))
      return(paste("'bugdata' missing column:", bugcolumns[i]))
    if(sum(is.na(object@bugdata[, bugcolumns[i]])) != 0 & i != 5)
      return(paste("Missing data in column", bugcolumns[i], "of 'bugdata'"))
  }
  for(i in 1:17){
    if(!(predictorcolumns[i] %in% names(object@predictors)))
      return(paste("'predictors' missing column:", predictorcolumns[i]))
    if(sum(is.na(object@predictors[, predictorcolumns[i]])) != 0)
      return(paste("Missing data in column", predictorcolumns[i], "of 'bugdata'"))
  }
  if(!("AREA_SQKM" %in% names(object@predictors)) & !("LogWSA" %in% names(object@predictors)))
    return("Predictors must include a column AREA_SQKM or LogWSA")
  if(!setequal(object@bugdata$StationCode, object@predictors$StationCode))
    return("All StationCode IDs must be represented in both bug and predictor data")
  if(length(unique(object@bugdata$SampleID)) != nrow(unique(object@bugdata[, c("StationCode", "SampleID")])))
    return("SampleIDs must be unique to one StationCode")
  TRUE
})
setValidity("oe", function(object){
  bugcolumns <- c("StationCode", "SampleID", "FinalID", "BAResult")
  predictorcolumns <- c("StationCode", "New_Lat",     "New_Long",    "ELEV_RANGE",  "BDH_AVE",     "PPT_00_09",  
                        "LPREM_mean",  "KFCT_AVE",    "TEMP_00_09",  "P_MEAN",      "N_MEAN",      "PRMH_AVE",   
                        "SITE_ELEV",   "MgO_Mean",    "S_Mean",      "SumAve_P",   
                        "CaO_Mean")
  for(i in 1:4){
    if(!(bugcolumns[i] %in% names(object@bugdata)))
      return(paste("'bugdata' missing column:", bugcolumns[i]))
    if(sum(is.na(object@bugdata[, bugcolumns[i]])) != 0 & i != 5)
      return(paste("Missing data in column", bugcolumns[i], "of 'bugdata'"))
  }
  for(i in 1:17){
    if(!(predictorcolumns[i] %in% names(object@predictors)))
      return(paste("'predictors' missing column:", predictorcolumns[i]))
    if(sum(is.na(object@predictors[, predictorcolumns[i]])) != 0)
      return(paste("Missing data in column", predictorcolumns[i], "of 'bugdata'"))
  }
  if(!("AREA_SQKM" %in% names(object@predictors)) & !("LogWSA" %in% names(object@predictors)))
    return("Predictors must include a column AREA_SQKM or LogWSA")
  if(!setequal(object@bugdata$StationCode, object@predictors$StationCode))
    return("All StationCode IDs must be represented in both bug and predictor data")
  if(length(unique(object@bugdata$SampleID)) != nrow(unique(object@bugdata[, c("StationCode", "SampleID")])))
    return("SampleIDs must be unique to one StationCode")
  TRUE
})
###MMI Methods###

setMethod("show", "bugs", function(object){
  print(head(object@bugdata))
  cat("\n")
  print(head(object@predictors))
  })

setGeneric("nameMatch", function(object, effort = "SAFIT1")
           standardGeneric("nameMatch"))
setMethod("nameMatch", "mmi", function(object, effort = "SAFIT1"){

  colnames(object@bugdata)[which(colnames(object@bugdata) == "FinalID")] <- "Taxa"
  colnames(object@bugdata)[which(colnames(object@bugdata) == "BAResult")] <- "Result"
  
  ###Clean data###
  object@bugdata$Taxa <- str_trim(object@bugdata$Taxa) 
  ###Aggregate taxa###
  object@bugdata <- ddply(object@bugdata, "SampleID", function(df){
    ddply(df, "Taxa", function(sdf){
      id <- unique(sdf[, !(colnames(sdf) %in% "Result")])
      Result <- sum(sdf$Result)
      cbind(id, Result)
    })
  })
  
  ###Match to STE###
  ibi <- dbGetQuery(object@dbconn, "SELECT * FROM ibitable2")
  object@bugdata$STE <- rep(NA, length(object@bugdata$Taxa))
  object@bugdata$STE <- ibi[match(object@bugdata$Taxa, ibi$FinalID), as.character(effort)]
  object@bugdata$STE <- as.character(object@bugdata$STE)
  if(length(object@bugdata$STE[which(is.na(object@bugdata$STE))] > 0)){
    warning(paste("The following taxa were not found in the database:",
                  paste(as.character(object@bugdata$Taxa[which(is.na(object@bugdata$STE))]),
                        collapse=" ", sep=" , ")))
  }
  object@bugdata$STE[which(is.na(object@bugdata$STE))] <- "Missing"
  
  
  ###Determine Distinctiveness###
  distinctsorter <- function(taxon, data){
    level <- ibi$TaxonomicLevelCode[match(taxon, ibi$FinalID)] 
    levelname <- as.character(ibi$TaxonomicLevelName[match(taxon, ibi$FinalID)])
    print(levelname)
    if(is.na(levelname))return(T)
    samelevel <- ibi$FinalID[which(ibi[, levelname] == taxon)]
    matchedlevel <- object@bugdata$STE %in% ibi$CustomSTE[match(samelevel, ibi$FinalID)]
    result <- ibi$TaxonomicLevelCode[match(data$STE[matchedlevel], ibi$FinalID)] > level
    length(which(result)) != 0
  }

  distinctlist <- dlply(object@bugdata, "SampleID", function(df){
    sapply(1:nrow(df), function(i){
      ifelse(distinctsorter(df$STE[i], df), "Non-Distinct", "Distinct")
    })})
  object@bugdata$distinct <- unlist(distinctlist)
  ###Override##
  object@bugdata$DistinctCode <- as.character(object@bugdata$DistinctCode)
  object@bugdata$distinct[which(object@bugdata$distinct == "Non-distinct" & object@bugdata$DistinctCode == "Yes")] <- "Distinct" 
  return(object)
})
setMethod("nameMatch", "oe", function(object, effort = "SAFIT1__OTU_a"){
  
  colnames(object@bugdata)[which(colnames(object@bugdata) == "FinalID")] <- "Taxa"
  colnames(object@bugdata)[which(colnames(object@bugdata) == "BAResult")] <- "Result"
  
  ###Clean data###
  object@bugdata$Taxa <- str_trim(object@bugdata$Taxa)
  ###Aggregate taxa###
  object@bugdata <- ddply(object@bugdata, "SampleID", function(df){
    ddply(df, "Taxa", function(sdf){
      id <- unique(sdf[, !(colnames(sdf) %in% "Result")])
      Result <- sum(sdf$Result)
      cbind(id, Result)
    })
  })
  
  ###Match to STE###
  otu_crosswalk <- dbGetQuery(object@dbconn, "SELECT * FROM ibitable2")
  object@bugdata$STE <- rep(NA, length(object@bugdata$Taxa))
  object@bugdata$STE <- otu_crosswalk[match(object@bugdata$Taxa, otu_crosswalk$FinalID), as.character(effort)]
  object@bugdata$STE <- as.character(object@bugdata$STE)
  object@bugdata$STE <- object@bugdata$STE[which(object@bugdata$STE != "Exclude")]
  object@bugdata$STE[which(is.na(object@bugdata$STE))] <- "Missing"
  
  ###Calculate ambiguous###
  percent.ambiguous <- ddply(object@bugdata, "SampleID", function(df){
    100*sum(df$Result[df$STE == "Ambiguous"])/sum(df$Result)
  })
  taxa.ambiguous <- ddply(object@bugdata, "SampleID", function(df){
    100*length(df$Taxa[df$STE == "Ambiguous"])/length(df$Taxa)
  })
  object@ambiguous <- merge(percent.ambiguous, taxa.ambiguous, by="SampleID")
  names(object@ambiguous)[2:3] <- c("individuals", "taxa")
  object@bugdata <- object@bugdata[object@bugdata$STE != "Ambiguous",]
  return(object)
})

setGeneric("subsample", function(object)
  standardGeneric("subsample"))
setMethod("subsample", "mmi", function(object){
  if(is.null(object@bugdata$distinct)){object <- nameMatch(object)}
  object@datalength <- length(object@bugdata)
  object@bugdata$SampleID <- as.character(object@bugdata$SampleID)
  rarifydown <- function(data){unlist(sapply(unique(data$SampleID), function(sample){
    v <- data[data$SampleID==sample, "Result"]
    if(sum(v)>=500){rrarefy(v, 500)} else
    {v}
  }))}
  registerDoParallel()
  rarificationresult <- foreach(i=1:20, .combine=cbind, .packages="vegan") %dopar% {
    rarifydown(object@bugdata)
  }

  closeAllConnections()
  object@subsample <- as.data.frame(cbind(object@bugdata, rarificationresult))
  colnames(object@subsample)[(object@datalength + 1):(object@datalength + 20)]<- paste("Replicate", 1:20)
  return(object)
})
setMethod("subsample", "oe", function(object){
  if(nrow(object@ambiguous)==0){object <- nameMatch(object)}
  object@datalength <- length(object@bugdata)
  object@bugdata$SampleID <- as.character(object@bugdata$SampleID)
  rarifydown <- function(data){unlist(sapply(unique(data$SampleID), function(sample){
    v <- data[data$SampleID==sample, "Result"]
    if(sum(v)>=400){rrarefy(v, 400)} else
    {v}
  }))}
  registerDoParallel()
  rarificationresult <- foreach(i=1:20, .combine=cbind, .packages="vegan") %dopar% {
    rarifydown(object@bugdata)
  }
  closeAllConnections()
  object@oesubsample <- as.data.frame(cbind(object@bugdata, rarificationresult))
  colnames(object@oesubsample)[(object@datalength + 1):(object@datalength + 20)]<- paste("Replicate", 1:20)
  return(object)
})

setGeneric("metrics", function(object)
  standardGeneric("metrics"))
setMethod("metrics", "mmi", function(object){
  if(nrow(object@subsample) == 0){object <- subsample(object)}
  ibi <- dbGetQuery(object@dbconn, "SELECT * FROM ibi")
  object@subsample$MaxTol <- ibi$MaxTol[match(object@subsample$Taxa, ibi$FinalID)]
  object@subsample$MaxTol <- as.numeric(object@subsample$MaxTol)
  object@subsample$Class <- ibi$Class[match(object@subsample$Taxa, ibi$FinalID)]
  object@subsample$Order <- as.character(ibi$Order[match(object@subsample$Taxa, ibi$FinalID)])
  object@subsample$FunctionalFeedingGroup <- as.character(ibi$FunctionalFeedingGroup[match(object@subsample$Taxa, ibi$FinalID)])
  
  object@metrics <- as.data.frame(matrix(NA, nrow = length(unique(object@bugdata$SampleID)), ncol = 180))
  for(i in 1:20){
    ###Number of Coleoptera taxa###
    object@metrics[[i]] <- ddply(object@subsample[object@subsample$distinct == "Distinct" & object@subsample[[object@datalength + i]]>0, ], "SampleID",
                          function(d)length(unique(d$STE[d$Order == "Coleoptera"])))[, 2]  		   
    
    ###Numer of Diptera taxa
    object@metrics[[i+20]] <- ddply(object@subsample[object@subsample$distinct == "Distinct" & object@subsample[[object@datalength + i]]>0, ], "SampleID",
                             function(d)length(unique(d$STE[d$Order == "Diptera"])))[, 2]
    
    ###Number of Scraper taxa###
    object@metrics[[i+40]] <- ddply(object@subsample[object@subsample$distinct == "Distinct" & object@subsample[[object@datalength + i]]>0, ], "SampleID",
                             function(d)length(unique(d$STE[which(d$FunctionalFeedingGroup == "SC")])))[, 2]
    
    ###Number of Predator taxa###
    object@metrics[[i+60]] <- ddply(object@subsample[object@subsample$distinct == "Distinct" & object@subsample[[object@datalength + i]]>0, ], "SampleID",
                             function(d)length(unique(d$STE[which(d$FunctionalFeedingGroup == "P")])))[, 2]	
    
    ###Number of CF/CG taxa###
    object@metrics[[i+80]] <- ddply(object@subsample[object@subsample$distinct == "Distinct" & object@subsample[[object@datalength + i]]>0, ], "SampleID",
                                    function(d)length(unique(d$STE[which(d$FunctionalFeedingGroup %in% c("CF", "CG"))])))[, 2]
    
    ###Number of Shredder taxa###
    object@metrics[[i+100]] <- ddply(object@subsample[object@subsample$distinct == "Distinct" & object@subsample[[object@datalength + i]]>0, ], "SampleID",
                                    function(d)length(unique(d$STE[which(d$FunctionalFeedingGroup == "SH")])))[, 2]
    
    
    ###Percent Intolerant###
    object@metrics[[i+120]] <- ddply(object@subsample, "SampleID",
                              function(d){
                                100*sum(d[which(d$MaxTol <= 2), object@datalength + 1])/sum(d[, object@datalength + 1])
                              })[, 2]
    
    ###Percent NonInsect###
    object@metrics[[i+140]] <- ddply(object@subsample, "SampleID",
                                     function(d){
                                       100*sum(d[which(d$Class != "Insecta"), object@datalength + 1])/sum(d[, object@datalength + 1])
                                     })[, 2]
    
    ###Weighted Tolerance###
    object@metrics[[i+160]] <- ddply(object@subsample, "SampleID",
                              function(d){
                                p <- which(!is.na(d$MaxTol))
                                sum((d[p, object@datalength + i] * d$MaxTol[p]))/sum(d[p, object@datalength + i])
                              })[, 2]
    
  }
  
  ###Means###
  object@metrics$coleoptera <- apply(object@metrics[, 1:20], 1, mean)
  object@metrics$diptera <- apply(object@metrics[, 21:40], 1, mean)
  object@metrics$scraper <- apply(object@metrics[, 41:60], 1, mean)
  object@metrics$predator <- apply(object@metrics[, 61:80], 1, mean)
  object@metrics$cfcg <- apply(object@metrics[, 81:100], 1, mean)
  object@metrics$shredder <- apply(object@metrics[, 101:120], 1, mean)
  object@metrics$intolerant <- apply(object@metrics[, 121:140], 1, mean)
  object@metrics$noninsect <- apply(object@metrics[, 141:160], 1, mean)
  object@metrics$tolerance <- apply(object@metrics[, 161:180], 1, mean)
  return(object)
})

setGeneric("randomForest", function(object)
  standardGeneric("randomForest"))
setMethod("randomForest", "mmi", function(object){
  if(nrow(object@metrics) == 0){object <- metrics(object)}
  object@predictors <- merge(unique(object@bugdata[, c("StationCode", "SampleID")]), object@predictors, by="StationCode")
  object@modelprediction <- as.data.frame(matrix(NA, nrow = nrow(object@predictors)))
  
  if(is.null(object@predictors$LogWSA))
    object@predictors$LogWSA <- log10(object@predictors$AREA_SQKM)
    
  object@modelprediction$coleoptera <- predict(Number_Coleoptera_Taxa.FinalForest, object@predictors)
  object@modelprediction$diptera <- object@metrics$diptera
  object@modelprediction$scraper <- predict(Number_Scraper_Taxa.FinalForest, object@predictors)
  object@modelprediction$predator <- predict(Number_Predator_Taxa.FinalForest, object@predictors)
  object@modelprediction$cfcg <- predict(Number_CF___CG_Taxa.FinalForest, object@predictors)
  object@modelprediction$shredder <- predict(Number_ShredderTaxa.FinalForest, object@predictors)
  object@modelprediction$intolerant <- predict(Percent_Intolerant_Taxa__0_2_.FinalForest, object@predictors)
  object@modelprediction$noninsect <- predict(Percent_Non_Insecta_Taxa.FinalForest, object@predictors)
  object@modelprediction$tolerance <- predict(Tolerance_Value.FinalForest, object@predictors)
  object@modelprediction$V1 <- unique(object@predictors$SampleID)
  return(object)
})
setMethod("randomForest", "oe", function(object){
  if(nrow(object@oesubsample)==0){object <- subsample(object)}
  
  if(!("LogWSA" %in% names(object@predictors)))
    object@predictors$LogWSA <- log10(object@predictors$AREA_SQKM)
  
  names(object@predictors)[which(names(object@predictors) == "TEMP_00_09")] <- "AvgTemp00_09"
  names(object@predictors)[which(names(object@predictors) == "PPT_00_09")] <- "AvgPPT00_09"
  names(object@predictors)[which(names(object@predictors) == "LogWSA")] <- "Log_Area"
  names(object@predictors)[which(names(object@predictors) == "SITE_ELEV")] <-  "AvgOfElevation"
  


  bugspa <- dbGetQuery(object@dbconn, "Select * FROM bugscal_pa2")
  row.names(bugspa) <- bugspa[, 1]
  bugspa <- bugspa[, 2:ncol(bugspa)]
  
  object@predictors <- merge(unique(object@oesubsample[, c("StationCode", "SampleID")]), object@predictors,
                             by="StationCode")
  row.names(object@predictors) <- paste(object@predictors$StationCode, "%", object@predictors$SampleID,
                                       sep="")
  iterate <- function(rep){
    patable <- dcast(data=object@oesubsample[, c("StationCode", "SampleID", "STE", rep)],
                     StationCode + SampleID ~ STE,
                     value.var=rep,
                     fun.aggregate=function(x)sum(x)/length(x))
    patable[is.na(patable)] <- 0
    row.names(patable) <- paste(patable$StationCode, "%", patable$SampleID, sep="")
    
    iresult <- model.predict.RanFor.4.2(bugcal.pa=bugspa,
                                        grps.final=dbGetQuery(object@dbconn, "Select * FROM grps_final")$grps_final,
                                        preds.final=c("AvgTemp00_09", "Log_Area", "AvgPPT00_09", "AvgOfElevation"),
                                        ranfor.mod=rf.mod,
                                        prednew=object@predictors,
                                        bugnew=patable,
                                        Pc=0.5,
                                        Cal.OOB=FALSE)
    iresult$SampleID <- unique(patable$SampleID)
    return(iresult)
  }
  object@fulliterations <- lapply(paste("Replicate", 1:20), function(i)iterate(i))
  labels <- strsplit(row.names(object@fulliterations[[1]]), "%")
  labels <- as.data.frame(matrix(unlist(labels), nrow=length(labels), byrow=T))
  object@fulliterations <- lapply(object@fulliterations, function(l){
    row.names(l)<-labels[, 2]
    l
    })
  object@iterations <- do.call(cbind, lapply(object@fulliterations, function(l)l$OoverE))
  object@oeresults <- data.frame(labels, apply(object@iterations, 1, mean))
  names(object@oeresults) <- c("StationCode", "SampleID", "OoverE")
  object
})
  
setGeneric("score", function(object, object2)
  standardGeneric("score"))
setMethod("score", "mmi", function(object){
  if(nrow(object@modelprediction) == 0){object <- randomForest(object)}
  maxmin <- dbGetQuery(object@dbconn, "SELECT * FROM maxmin")
  object@result <- as.data.frame(matrix(NA, nrow = length(unique(object@modelprediction$V1)), ncol = 18))
  colnames(object@result) <- c("coleoptera", "diptera", "scraper", "predator", "cfcg", "shredder",
                               "intolerant", "noninsect", "tolerance", "coleoptera_score", "diptera_score", 
                               "scraper_score", "predator_score", "cfcg_score", "shredder_score", "intolerant_score", 
                               "noninsect_score", "tolerance_score")
  for(i in c(1, 3:7)){
    object@result[, i] <- (object@metrics[, 180 + i] - object@modelprediction[, i+1] - maxmin[i, 4]) /
      (maxmin[i, 3] - maxmin[i, 4])
  } 
  object@result$diptera <- (object@modelprediction$diptera -  maxmin[2, 4]) /
    (maxmin[2, 3] - maxmin[2, 4])
  object@result$noninsect <- (object@metrics$noninsect - object@modelprediction$noninsect - maxmin[8, 3]) /
    (maxmin[8, 4] - maxmin[8, 3])
  object@result$tolerance <- (object@metrics$tolerance - object@modelprediction$tolerance - maxmin[9, 3]) /
    (maxmin[9, 4] - maxmin[9, 3])
  
  for(i in 1:9)
    object@result[, i + 9] <- ifelse(object@result[, i] <= 0, 0, ifelse(
      object@result[, i] >= 1, 1, object@result[, i]))
  
  object@finalscore <- data.frame(unique(object@modelprediction$V1), 
                                  apply(object@result[, 10:18], 1, mean)/0.5974159)
  d <- data.frame(object@finalscore, object@metrics[, 181:189], 
                  object@modelprediction, object@result[, 9:18])
  d <- merge(unique(object@bugdata[, c("StationCode", "SampleID")]), d, by.x="SampleID", by.y="unique.object.modelprediction.V1.")
  colnames(d)[1:3] <- c("SampleID", "StationCode", "MMI Score")
  object@summary <- d

  return(object)
  }         
)
setMethod("score", "oe", function(object)randomForest(object))

setMethod("summary", "mmi", function(object = "mmi"){
  if(nrow(object@result) != 0){
    object@summary
  } else
    show(object)
})
setMethod("summary", "oe", function(object = "oe"){
  if(nrow(object@oeresults) != 0){
    object@oeresults
  } else
    show(object)
})
setMethod("plot", "mmi", function(x = "mmi"){
  load("data/base_map.rdata")
  x@result$MMIScore <- cut(x@finalscore[, 2], breaks=c(0, .3, .8, 1.5), labels=c("low", "medium", "high"))
  x@result <- cbind(x@result, x@predictors[, 1:4])
  ggmap(base_map) + 
    geom_point(data=x@result, aes(x=New_Long, y=New_Lat, colour=MMIScore), size=4, alpha=.6)
})
setMethod("plot", "oe", function(x = "oe"){
  load("data/base_map.rdata")
  x@result$MMIScore <- cut(x@oeresults[, 2], breaks=c(0, .3, .8, 1.5), labels=c("low", "medium", "high"))
  x@result <- cbind(x@result, x@predictors[, 1:4])
  ggmap(base_map) + 
    geom_point(data=x@result, aes(x=New_Long, y=New_Lat, colour=MMIScore), size=4, alpha=.6)
})

setMethod("initialize", "metric.mean", function(.Object="metric.mean", x="mmi", y="oe"){
  for(i in names(getSlots("mmi"))){
    slot(.Object, i) <- slot(x, i)
  }
  for(i in names(getSlots("oe"))){
    slot(.Object, i) <- slot(y, i)
  }
  .Object@mean.metric <- merge(y@oeresults, x@summary[, 1:3])
  .Object@mean.metric$Hybrid <- apply(.Object@mean.metric[, c("OoverE", "MMI Score")], 1, mean)
  .Object
})

setMethod("summary", "metric.mean", function(object = "metric.mean", report="basic"){
  if(report %in% c("basic", "standard", "detailed", "complete")){
    object@mean.metric$Count <- ddply(object@subsample, "SampleID", function(df)sum(df$Result))[, 2]
    object@mean.metric$Pcnt_Ambiguous_Individuals <- object@ambiguous$individuals
    object@mean.metric$Pcnt_Ambiguous_Taxa <- object@ambiguous$taxa
    object@mean.metric$overall_flag <- ifelse(object@mean.metric$Count >=450 & object@mean.metric$Pcnt_Ambiguous_Individuals < 20,
                                      "Adequate", "Inadequate")
    object@mean.metric <- object@mean.metric[, c(1:2, 6:9, 3:5)]
  }
  if(report %in% c("standard", "detailed", "complete")){
    object@mean.metric$mmi_count_flag <- ifelse(object@mean.metric$Count >=450, "Adequate", "Inadequate")
    object@mean.metric$ambig_count_flag <- ifelse(object@mean.metric$Pcnt_Ambiguous_Individuals < 20, "Adequate", "Inadequate")
    #object@mean.metric$ambig_taxa_flag <- ifelse(object@mean.metric$Pcnt_Ambiguous_Taxa < 25, "Adequate", "Inadequate")
    object@mean.metric <- object@mean.metric[, c(1:6, 10:11, 7:9)]
  }
  if(report %in% c("detailed", "complete")){
    names(object@metrics)[181:189] <- paste(names(object@metrics)[181:189], "_", "metric", sep="")
    names(object@modelprediction) <- paste(names(object@modelprediction), "_", "predicted_metric", sep="")
    object@mean.metric <- cbind(object@mean.metric, object@metrics[, 181:189], object@result[, 10:18], object@modelprediction[2:9])
    object@mean.metric$Expected_capture <- apply(do.call(cbind, lapply(object@fulliterations, function(l)l$E)), 1, mean)
    object@mean.metric$Observed_taxa <- apply(do.call(cbind, lapply(object@fulliterations, function(l)l$O)), 1, mean)
  }
  if(report == "complete"){
    list(object@mean.metric, object@subsample, object@oesubsample, do.call(cbind, object@fulliterations))
  } else
  if(!(report %in% c("basic", "standard", "detailed"))){
    stop("Allowed values for report: 'basic', 'standard', 'detailed', 'complete'")
  } else
    object@mean.metric
})

setMethod("plot", "metric.mean", function(x="metric.mean"){
  load("data/base_map.rdata")
  x@mean.metric$HybridScore <- cut(x@mean.metric$Hybrid, breaks=c(0, .4, .8, 1.5), labels=c("low", "medium", "high"))
  hscore <- cbind(x@mean.metric$HybridScore, x@predictors[, 1:4])
  names(hscore)[1] <- "HybridScore"
  ggmap(base_map) + 
    geom_point(data=hscore, aes(x=New_Long, y=New_Lat, colour=HybridScore), size=4, alpha=.6)
})

