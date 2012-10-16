setClass("metricMean", representation(mean.metric="data.frame"),
         contains=c("mmi", "oe"))

setValidity("metricMean", function(object){
  if(!(setequal(object@subsample$SampleID, object@oesubsample$SampleID))){return("Incompatible objects")} else
    T
})

setMethod("initialize", "metricMean", function(.Object="metricMean", x="mmi", y="oe"){
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

setMethod("summary", "metricMean", function(object = "metricMean", report="basic"){
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

setMethod("plot", "metricMean", function(x="metricMean"){
  load("data/base_map.rdata")
  x@mean.metric$HybridScore <- cut(x@mean.metric$Hybrid, breaks=c(0, .4, .8, 1.5), labels=c("low", "medium", "high"))
  hscore <- cbind(x@mean.metric$HybridScore, x@predictors[, c("StationCode", "SampleID", "New_Lat", "New_Long")])
  names(hscore)[1] <- "HybridScore"
  ggmap(base_map) + 
    geom_point(data=hscore, aes(x=New_Long, y=New_Lat, colour=HybridScore), size=4, alpha=.6)
})