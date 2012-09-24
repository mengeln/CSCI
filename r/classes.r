
library(randomForest)
library(plyr)
library(doParallel)
library(vegan)
library(ggmap)
library(reshape2)

load("data/FinalForests.Rdata")
load("data/ibiv4.RData")
load("data/taxonomy_v5.RData")
load("data/Metrics.Max.Min.Rdata")
load("data/otu_crosswalk.rdata")
load("data/OE.RF.Model.Rdata")
source("r/model.predict.RanFor.4.2.r")

###Setup###
min <- c(Number_Coleoptera_Taxa.min.met, Number_Diptera_Taxa.min.met, Number_Scraper_Taxa.min.met,
         Number_Predator_Taxa.min.met, Number_CF___CG_Taxa.min.met, Number_ShredderTaxa.min.met,
         Percent_Intolerant_Taxa__0_2_.min.met, Percent_Non_Insecta_Taxa.min.met,
         Tolerance_Value.min.met
         )
max <- c(Number_Coleoptera_Taxa.max.met, Number_Diptera_Taxa.max.met, Number_Scraper_Taxa.max.met,
         Number_Predator_Taxa.max.met, Number_CF___CG_Taxa.max.met, Number_ShredderTaxa.max.met,
         Percent_Intolerant_Taxa__0_2_.max.met, Percent_Non_Insecta_Taxa.max.met,
         Tolerance_Value.max.met
         )
maxminlabels <- c("coleoptera", "diptera", "scraper", "predator", "cfcg", "shredder",
                  "intolerant", "noninsect", "tolerance")
maxmin <- data.frame(maxminlabels, max, min)

taxonomy <- idata.frame(taxonomy_v5)
ibi <- idata.frame(ibiv4)

setOldClass("randomForest.formula")
setOldClass("idf")

###Class Definition###
setClass("bugs", representation(bugdata="data.frame",
                                predictors="data.frame"))
setClass("mmi", representation(taxonomy = "idf",
                               ibi = "idf",
                               cf.randomForest = "randomForest.formula",
                               coleoptera.randomForest = "randomForest.formula",
                               predator.randomForest = "randomForest.formula",
                               scraper.randomForest = "randomForest.formula",
                               shredder.randomForest = "randomForest.formula",
                               intolerant.randomForest = "randomForest.formula",
                               noninsect.randomForest = "randomForest.formula",
                               tolerance.randomForest = "randomForest.formula",
                               subsample = "data.frame",
                               metrics = "data.frame",
                               modelprediction = "data.frame",
                               maxmin = "data.frame",
                               result = "data.frame",
                               finalscore = "data.frame",
                               datalength = "numeric"),
         contains="bugs",
         prototype = list(ibi = ibi,
                          taxonomy = taxonomy,
                          cf.randomForest = Number_CF___CG_Taxa.FinalForest,
                          coleoptera.randomForest = Number_Coleoptera_Taxa.FinalForest,
                          predator.randomForest = Number_Predator_Taxa.FinalForest,
                          scraper.randomForest = Number_Scraper_Taxa.FinalForest,
                          shredder.randomForest = Number_ShredderTaxa.FinalForest,
                          intolerant.randomForest = Percent_Intolerant_Taxa__0_2_.FinalForest,
                          noninsect.randomForest = Percent_Non_Insecta_Taxa.FinalForest,
                          tolerance.randomForest = Tolerance_Value.FinalForest,
                          subsample = data.frame(),
                          metrics = data.frame(),
                          residuals = data.frame(),
                          maxmin = maxmin,
                          result = data.frame(),
                          finalscore = data.frame(),
                          datalength = numeric()
         )
)
setClass("oe", representation(otu_crosswalk="data.frame",
                              ambiguous="data.frame",
                              oesubsample="data.frame",
                              iterations="matrix",
                              fulliterations="list",
                              datalength="numeric",
                              rf.mod="randomForest.formula",
                              preds.final="character",
                              grps.final="integer",
                              bugcal.pa="data.frame",
                              oeresults="data.frame"), 
         contains="bugs",
         prototype=list(otu_crosswalk=otu_crosswalk,
                        ambiguous=data.frame(),
                        oesubsample=data.frame(),
                        datalength=numeric(),
                        iterations=matrix(),
                        fulliterations=list(),
                        datalength=numeric(),
                        rf.mod=rf.mod,
                        preds.final=preds.final,
                        grps.final=grps.final,
                        bugcal.pa=bugcal.pa,
                        oeresults=data.frame()
                        ))
setClass("metric.mean", representation(mean.metric="data.frame"),
         contains=c("mmi", "oe"))
###Validity Checks###
setValidity("bugs", function(object){
  bugcolumns <- c("StationCode", "SampleID", "FinalID", "BAResult", "DistinctCode")
  predictorcolumns <- c("StationCode", "New_Lat",     "New_Long",    "ELEV_RANGE",  "BDH_AVE",     "PPT_00_09",  
                  "LPREM_mean",  "KFCT_AVE",    "TEMP_00_09",  "P_MEAN",      "N_MEAN",      "PRMH_AVE",   
                  "AREA_SQKM",   "LogWSA",      "SITE_ELEV",   "MgO_Mean",    "S_Mean",      "SumAve_P",   
                  "CaO_Mean")
  for(i in 1:5){
    if(!(bugcolumns[i] %in% names(object@bugdata)))
      return(paste("'bugdata' missing column:", bugcolumns[i]))
    if(sum(is.na(object@bugdata[, bugcolumns[i]])) != 0 & i != 5)
      return(paste("Missing data in column", bugcolumns[i], "of 'bugdata'"))
  }
  for(i in 1:19){
    if(!(predictorcolumns[i] %in% names(object@predictors)))
      return(paste("'predictors' missing column:", predictorcolumns[i]))
    if(sum(is.na(object@predictors[, predictorcolumns[i]])) != 0)
      return(paste("Missing data in column", predictorcolumns[i], "of 'bugdata'"))
    }
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
  
  ###Aggregate taxa###
  object@bugdata <- ddply(object@bugdata, "SampleID", function(df){
    ddply(df, "Taxa", function(sdf){
      id <- unique(sdf[, !(colnames(sdf) %in% "Result")])
      Result <- sum(sdf$Result)
      cbind(id, Result)
    })
  })
  
  ###Match to STE###
  object@bugdata$STE <- rep(NA, length(object@bugdata$Taxa))
  object@bugdata$STE <- object@ibi[match(object@bugdata$Taxa, object@ibi$FinalID), as.character(effort)]
  object@bugdata$STE <- as.character(object@bugdata$STE)
  object@bugdata$STE[which(is.na(object@bugdata$STE))] <- "Missing"
  
  ###Determine Distinctiveness###
  distinctsorter <- function(taxon, data){
    level <- object@taxonomy$TaxonomicLevelCode[match(taxon, object@taxonomy$FinalID)] 
    levelname <- as.character(object@taxonomy$TaxonomicLevelName[match(taxon, object@taxonomy$FinalID)])
    samelevel <- object@taxonomy$FinalID[which(taxonomy[, levelname] == taxon)]
    matchedlevel <- object@bugdata$STE %in% object@ibi$CustomSTE[match(samelevel, object@ibi$FinalID)]
    result <- object@taxonomy$TaxonomicLevelCode[match(data$STE[matchedlevel], object@taxonomy$FinalID)] > level
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
  
  ###Aggregate taxa###
  object@bugdata <- ddply(object@bugdata, "SampleID", function(df){
    ddply(df, "Taxa", function(sdf){
      id <- unique(sdf[, !(colnames(sdf) %in% "Result")])
      Result <- sum(sdf$Result)
      cbind(id, Result)
    })
  })
  
  ###Match to STE###
  object@bugdata$STE <- rep(NA, length(object@bugdata$Taxa))
  object@bugdata$STE <- object@otu_crosswalk[match(object@bugdata$Taxa, object@otu_crosswalk$FinalID), as.character(effort)]
  object@bugdata$STE <- as.character(object@bugdata$STE)
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
  object@oesubsample <- as.data.frame(cbind(object@bugdata, rarificationresult))
  colnames(object@oesubsample)[(object@datalength + 1):(object@datalength + 20)]<- paste("Replicate", 1:20)
  return(object)
})

setGeneric("metrics", function(object)
  standardGeneric("metrics"))
setMethod("metrics", "mmi", function(object){
  if(nrow(object@subsample) == 0){object <- subsample(object)}
  object@subsample$MaxTol <- ibi$MaxTol[match(object@subsample$Taxa, object@ibi$FinalID)]
  object@subsample$MaxTol <- as.numeric(object@subsample$MaxTol)
  object@subsample$Class <- ibi$Class[match(object@subsample$Taxa, object@ibi$FinalID)]
  object@subsample$Order <- as.character(object@ibi$Order[match(object@subsample$Taxa, object@ibi$FinalID)])
  object@subsample$FunctionalFeedingGroup <- as.character(object@ibi$FunctionalFeedingGroup[match(object@subsample$Taxa, object@ibi$FinalID)])
  
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
  
  object@modelprediction$coleoptera <- predict(object@coleoptera.randomForest, object@predictors)
  object@modelprediction$diptera <- object@metrics$diptera
  object@modelprediction$scraper <- predict(object@scraper.randomForest, object@predictors)
  object@modelprediction$predator <- predict(object@predator.randomForest, object@predictors)
  object@modelprediction$cfcg <- predict(object@cf.randomForest, object@predictors)
  object@modelprediction$shredder <- predict(object@shredder.randomForest, object@predictors)
  object@modelprediction$intolerant <- predict(object@intolerant.randomForest, object@predictors)
  object@modelprediction$noninsect <- predict(object@noninsect.randomForest, object@predictors)
  object@modelprediction$tolerance <- predict(object@tolerance.randomForest, object@predictors)
  object@modelprediction$V1 <- unique(object@predictors$SampleID)
  return(object)
})
setMethod("randomForest", "oe", function(object){
  if(nrow(object@oesubsample)==0){object <- subsample(object)}
  names(object@predictors)[which(names(object@predictors) == "TEMP_00_09")] <- "AvgTemp00_09"
  names(object@predictors)[which(names(object@predictors) == "PPT_00_09")] <- "AvgPPT00_09"
  names(object@predictors)[which(names(object@predictors) == "LogWSA")] <- "Log_Area"
  names(object@predictors)[which(names(object@predictors) == "SITE_ELEV")] <-  "AvgOfElevation"
  Pc=0.5
  object@predictors <- merge(unique(object@oesubsample[, c("StationCode", "SampleID")]), object@predictors,
                             by="StationCode")
  row.names(object@predictors) <- paste(object@predictors$StationCode, "%", object@predictors$SampleID,
                                       sep="")
  
  iterate <- function(rep){
    patable <- dcast(object@oesubsample[, c("StationCode", "SampleID", "STE", rep)],
                     StationCode + SampleID ~ STE, value.var=rep, mean)
    patable[is.na(patable)] <- 0
    row.names(patable) <- paste(patable$StationCode, "%", patable$SampleID, sep="")
    iresult <- model.predict.RanFor.4.2(object@bugcal.pa, object@grps.final, object@preds.final, object@rf.mod, prednew=object@predictors,
                                        bugnew=patable,Pc,Cal.OOB=FALSE)
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
  
setGeneric("score", function(object)
  standardGeneric("score"))
setMethod("score", "mmi", function(object){
  if(nrow(object@modelprediction) == 0){object <- randomForest(object)}
  object@result <- as.data.frame(matrix(NA, nrow = length(unique(object@modelprediction$V1)), ncol = 18))
  colnames(object@result) <- c("coleoptera", "diptera", "scraper", "predator", "cfcg", "shredder",
                               "intolerant", "noninsect", "tolerance", "coleoptera_score", "diptera_score", 
                               "scraper_score", "predator_score", "cfcg_score", "shredder_score", "intolerant_score", 
                               "noninsect_score", "tolerance_score")
  for(i in c(1, 3:7)){
    object@result[, i] <- (object@metrics[, 180 + i] - object@modelprediction[, i+1] - maxmin[i, 3]) /
      (maxmin[i, 2] - maxmin[i, 3])
  } 
  object@result$diptera <- (object@modelprediction$diptera -  maxmin[2, 3]) /
    (maxmin[2, 2] - maxmin[2, 3])
  object@result$noninsect <- (object@metrics$noninsect - object@modelprediction$noninsect - maxmin[i, 2]) /
    (maxmin[i, 3] - maxmin[i, 2])
  object@result$tolerance <- (object@metrics$tolerance - object@modelprediction$tolerance - maxmin[i, 2]) /
    (maxmin[i, 3] - maxmin[i, 2])
  
  for(i in 1:9)
    object@result[, i + 9] <- ifelse(object@result[, i] <= 0, 0, ifelse(
      object@result[, i] >= 1, 1, object@result[, i]))
  
  object@finalscore <- data.frame(unique(object@modelprediction$V1), 
                                  apply(object@result[, 10:18], 1, mean)/0.5974159)
  return(object)
  }         
)
setMethod("score", "oe", function(object)randomForest(object))

setMethod("summary", "mmi", function(object = "mmi"){
  if(nrow(object@result) != 0){
    d <- data.frame(object@finalscore, object@metrics[, 181:189], 
                    object@modelprediction, object@result[, 9:18])
    d <- merge(unique(object@bugdata[, c("StationCode", "SampleID")]), d, by.x="SampleID", by.y="unique.object.modelprediction.V1.")
    colnames(d)[1:3] <- c("SampleID", "StationCode", "MMI Score")
    d
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
setMethod("initialize", "metric.mean", function(mmi, oe){
  
  
})
setMethod("mean", c("mmi", "oe"))
  