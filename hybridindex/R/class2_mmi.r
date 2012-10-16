lazyLoad("data/forestsdb")

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

setValidity("mmi", function(object){
  if(validity(object) & !("LifeStageCode" %in% names(object@bugdata))){return("Column LifeStageCode missing")} else
    TRUE
})

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

setMethod("metrics", "mmi", function(object){
  if(nrow(object@subsample) == 0){object <- subsample(object)}
  ibi <- dbGetQuery(object@dbconn, "SELECT * FROM ibitable2")
  object@subsample <- merge(object@subsample, ibi[, c("FinalID", "LifeStageCode",
                                                      "ToleranceValue", "Class", "Order", "FunctionalFeedingGroup")],
                            by.x=c("Taxa", "LifeStageCode"), by.y=c("FinalID", "LifeStageCode"))
#   object@subsample$MaxTol <- ibi$MaxTol[match(object@subsample$Taxa, ibi$FinalID)]
#   object@subsample$MaxTol <- as.numeric(object@subsample$MaxTol)
#   object@subsample$Class <- ibi$Class[match(object@subsample$Taxa, ibi$FinalID)]
#   object@subsample$Order <- as.character(ibi$Order[match(object@subsample$Taxa, ibi$FinalID)])
#   object@subsample$FunctionalFeedingGroup <- as.character(ibi$FunctionalFeedingGroup[match(object@subsample$Taxa, ibi$FinalID)])
  
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
                                       100*sum(d[which(d$ToleranceValue <= 2), object@datalength + 1])/sum(d[, object@datalength + 1])
                                     })[, 2]
    
    ###Percent NonInsect###
    object@metrics[[i+140]] <- ddply(object@subsample, "SampleID",
                                     function(d){
                                       100*sum(d[which(d$Class != "Insecta"), object@datalength + 1])/sum(d[, object@datalength + 1])
                                     })[, 2]
    
    ###Weighted Tolerance###
    object@metrics[[i+160]] <- ddply(object@subsample, "SampleID",
                                     function(d){
                                       p <- which(!is.na(d$ToleranceValue))
                                       sum((d[p, object@datalength + i] * d$ToleranceValue[p]))/sum(d[p, object@datalength + i])
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

setMethod("rForest", "mmi", function(object){
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

setMethod("score", "mmi", function(object){
  if(nrow(object@modelprediction) == 0){object <- rForest(object)}
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
})

setMethod("summary", "mmi", function(object = "mmi"){
  if(nrow(object@result) != 0){
    object@summary
  } else
    show(object)
})

setMethod("plot", "mmi", function(x = "mmi"){
  load("data/base_map.rdata")
  x@result$MMIScore <- cut(x@finalscore[, 2], breaks=c(0, .3, .8, 1.5), labels=c("low", "medium", "high"))
  x@result <- cbind(x@result, x@predictors[, c("StationCode", "SampleID", "New_Lat", "New_Long")])
  ggmap(base_map) + 
    geom_point(data=x@result, aes(x=New_Long, y=New_Lat, colour=MMIScore), size=4, alpha=.6)
})
            
            