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

setMethod("subsample", "mmi", function(object, rand = sample.int(10000, 1)){
  if(is.null(object@bugdata$distinct)){object <- nameMatch(object)}
  
  subsample <- as.data.frame(sapply(seq(1 + rand, 20 + rand), function(i){
    commMatrix <- acast(object@bugdata, SampleID ~ Taxa + LifeStageCode + distinct, value.var="BAResult", fill=0,
                        fun.aggregate = sum, na.rm=TRUE)
    samp <- rep.int(500, nrow(commMatrix))
    samp[rowSums(commMatrix, na.rm=TRUE) < 500] <- 
      rowSums(commMatrix, na.rm=TRUE)[rowSums(commMatrix, na.rm=TRUE) < 500]
    set.seed(i)
    
    commMatrix <- rrarefy(commMatrix, samp)
    
    if(i == 1+ rand)
      melt(commMatrix)
    else
      melt(commMatrix)$value
  }
  ))
  colnames(subsample)[3:22] <- paste("Replicate", 1:20)
  subsample <- subsample[rowSums(subsample[, 3:22]) != 0, ]
  break_ids <- strsplit(as.character(subsample$Var2), "_")
  subsample$Taxa <- sapply(break_ids, `[`, 1)
  subsample$LifeStageCode <- sapply(break_ids, `[`, 2)
  subsample$distinct <- sapply(break_ids, `[`, 3)
  subsample$SampleID <- subsample$Var1
  subsample <- subsample[, c(-1, -2)]
  object@subsample <- merge(arrange(object@bugdata, Taxa, LifeStageCode), subsample, all.x=TRUE, 
                            by.x=c("SampleID", "SAFIT2", "LifeStageCode", "distinct"),
                            by.y=c("SampleID", "Taxa", "LifeStageCode", "distinct"))
  return(object)
})

setMethod("metrics", "mmi", function(object){
  if(nrow(object@subsample) == 0){object <- subsample(object, rand = sample.int(10000, 1))}
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
        Clinger_PercentTaxa = nrow(.SD[distinct=="Distinct" & (Habit == "CN")])/nrow(.SD[distinct=="Distinct"]),
        Coleoptera_PercentTaxa = nrow(.SD[distinct=="Distinct" & (Order == "Coleoptera")])/nrow(.SD[distinct=="Distinct"]),
        Taxonomic_Richness = length(unique(.SD[distinct=="Distinct", FinalID])),
        EPT_PercentTaxa = nrow(.SD[distinct=="Distinct" & (Order %in% c("Ephemeroptera", "Plecoptera", "Trichoptera"))])/nrow(.SD[distinct=="Distinct"]),
        Shredder_Taxa = nrow(.SD[distinct=="Distinct" & (FunctionalFeedingGroup == "SH")]),
        Intolerant_Percent = sum(.SD[ToleranceValue >= 8, BAResult])/sum(BAResult)
      ), by=SampleID]})
    result.df <- lapply(result, data.frame)
    endpoint <- length(csci_metrics) + 1
    result.df <- lapply(1:20, function(i){
      names(result.df[[i]]) <- c("SampleID", paste0(names(result.df[[i]])[2:endpoint], i))
      result.df[[i]]
    })
    result.reduce <- Reduce(function(x,y)merge(x,y, by="SampleID"), result.df)
    names <- csci_metrics
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
  load(system.file("data", "Metrics.RFModels_v2.RData",  package="CSCI"))
  object@predictors <- merge(unique(object@bugdata[, c("StationCode", "SampleID")]), object@predictors, by="StationCode", all.x=TRUE)
  object@modelprediction <- as.data.frame(matrix(NA, nrow = nrow(object@predictors)))
  
  if(is.null(object@predictors$LogWSA))
    object@predictors$LogWSA <-log10(object@predictors$AREA_SQKM)
  object@predictors$Log_P_MEAN <-  log10(object@predictors$P_MEAN + 0.00001)
  
  res <- sapply(final.forests, function(rf)predict(rf, object@predictors))
  if(class(res)!="matrix")res <- data.frame(t(res[1:8]))
  
  object@modelprediction <- as.data.frame(res)
  names(object@modelprediction) <- csci_metrics
  object@modelprediction$V1 <- unique(object@predictors$SampleID)
  return(object)
})

setMethod("score", "mmi", function(object){
  if(nrow(object@modelprediction) == 0){object <- rForest(object)}
  load(system.file("data", "maxmin_v2.rdata",  package="CSCI"))
  col_names <- csci_metrics
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
                                  apply(object@result, 1, mean)/0.628016448)
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

            
            