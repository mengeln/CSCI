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

setMethod("nameMatch", "mmi", function(object, effort = 1){
  agg <- aggregate(BMI(object@bugdata))[[effort]]
  class(agg) <- "data.frame"
  agg$Result <- agg$BAResult
  agg$Taxa <- as.character(agg$FinalID)
  agg$distinct <- agg[, paste0("distinct_SAFIT", effort)]
  if(effort==1)x <- rename(agg, c("distinct_SAFIT1" = "distinct_SAFIT2", "SAFIT2" = "iggSAFIT2", "SAFIT1" = "SAFIT2"))
  object@bugdata <- agg
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
  BMIall_hybrid <- function(x){
    x$Habit <- as.character(x$Habit)
    replicate_list <- lapply(1:20, function(i){
      data <- x[, 1:28]
      data$BAResult <- x[, paste("Replicate", i)]
      data.table(data)
    })
    result <- lapply(replicate_list, function(x){
      x <- subset(x, BAResult > 0)
      x[, list(
        Shannon_Diversity = diversity(tapply(BAResult, as.character(SAFIT1), sum), index= "shannon"),
        Intolerant_PercentTaxa = nrow(.SD[distinct=="Distinct" & ToleranceValue <= 2])/nrow(.SD[distinct=="Distinct"]),
        ToleranceValue = sum(BAResult * ToleranceValue, na.rm=T)/sum(BAResult),
        CFCG_Taxa = nrow(.SD[distinct=="Distinct" & (FunctionalFeedingGroup == "CF" | FunctionalFeedingGroup == "CG")]),
        Shredder_Taxa = nrow(.SD[distinct=="Distinct" & (FunctionalFeedingGroup == "SH")]),
        Clinger_Taxa = nrow(.SD[distinct=="Distinct" & (Habit == "CN")]),
        Coleoptera_Taxa = nrow(.SD[distinct=="Distinct" & (Order == "Coleoptera")]),
        Noninsect_PercentTaxa = nrow(.SD[distinct=="Distinct" & (Class != "Insecta")])/nrow(.SD[distinct=="Distinct"])  
      ),
        by=SampleID]})
    result.df <- lapply(result, data.frame)
    result.df <- lapply(1:20, function(i){
      names(result.df[[i]]) <- c("SampleID", paste0(names(result.df[[i]])[2:9], i))
      result.df[[i]]
    })
    result.reduce <- Reduce(function(x,y)merge(x,y, by="SampleID"), result.df)
    names <- c("Shannon_Diversity", "Intolerant_PercentTaxa", "ToleranceValue",
               "CFCG_Taxa", "Shredder_Taxa", "Clinger_Taxa", "Coleoptera_Taxa", "Noninsect_PercentTaxa")
    means <- sapply(names, function(names)apply(result.reduce[, grep(names, names(result.reduce))], 1, mean))
    if(class(means) != "matrix")means <- t(means)
    result_final <- cbind(result.reduce, means)
    result_final
  }
  object@metrics <- BMIall_hybrid(object@subsample)
  return(object)
})

setMethod("rForest", "mmi", function(object){
  if(nrow(object@metrics) == 0){object <- metrics(object)}
  load(system.file("data", "Metrics.RFModels.RData",  package="CSCI"))
  object@predictors <- merge(unique(object@bugdata[, c("StationCode", "SampleID")]), object@predictors, by="StationCode", all.x=TRUE)
  object@modelprediction <- as.data.frame(matrix(NA, nrow = nrow(object@predictors)))
  
  object@predictors$LogWSA<-log10(object@predictors$AREA_SQKM)
  object@predictors$Log_P_MEAN<-  log10(object@predictors$P_MEAN + 0.0001)
  object@predictors$Log_N_MEAN<-  log10(object@predictors$N_MEAN + 0.00001)
  
  res <- sapply(final.forests, function(rf)predict(rf, object@predictors))
  if(class(res)!="matrix")res <- data.frame(t(res[1:8]))
  
  object@modelprediction <- as.data.frame(res)
  names(object@modelprediction) <- c("Shannon_Diversity", "Intolerant_PercentTaxa", "ToleranceValue",
                                     "CFCG_Taxa", "Shredder_Taxa", "Clinger_Taxa", "Coleoptera_Taxa", "Noninsect_PercentTaxa")
  object@modelprediction$V1 <- unique(object@predictors$SampleID)
  return(object)
})

setMethod("score", "mmi", function(object){
  if(nrow(object@modelprediction) == 0){object <- rForest(object)}
  load(system.file("data", "maxmin.rdata",  package="CSCI"))
  col_names <- c("Shannon_Diversity", "Intolerant_PercentTaxa", "ToleranceValue",
                 "CFCG_Taxa", "Shredder_Taxa", "Clinger_Taxa", "Coleoptera_Taxa", "Noninsect_PercentTaxa")
  object@metrics <- object@metrics[order(object@metrics$SampleID), ]
  object@modelprediction <- object@modelprediction[order(object@modelprediction$V1), ]
  
  object_result <- sapply(col_names, function(col){
    result <- (object@metrics[, col] - object@modelprediction[, col] - maxmin[1, col])/(maxmin[2, col] - maxmin[1, col])
    result <- ifelse(result > 1, 1, ifelse(
      result < 0, 0, result))
    result
  })
  if(class(object_result) != "matrix")object_result <- t(object_result)
  object@result <- data.frame(object_result)
  names(object@result) <- paste0(col_names, "_score")
  object@finalscore <- data.frame(unique(object@modelprediction$V1), 
                                  apply(object@result, 1, mean)/0.808484)
  d <- data.frame(object@finalscore, object@metrics[, col_names], 
                  object@modelprediction, object@result)
  d <- merge(unique(object@bugdata[, c("StationCode", "SampleID")]), d, by.x="SampleID", by.y="unique.object.modelprediction.V1.")
  colnames(d)[1:3] <- c("SampleID", "StationCode", "MMI_Score")
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
  load(system.file("data", "base_map.rdata", package="CSCI"))
  if(nrow(x@finalscore) == 0)stop("No data to plot; compute MMI scores first")
  x@result$MMIScore <- x@finalscore[, 2]
  x@result <- cbind(x@result, x@predictors[, c("StationCode", "SampleID", "New_Lat", "New_Long")])
  ggmap(base_map) + 
    geom_point(data=x@result, aes(x=New_Long, y=New_Lat, colour=MMIScore), size=4, alpha=.6) + 
    scale_color_continuous(low="red", high="green", name="MMI Score") + labs(x="", y="")
})
            
            