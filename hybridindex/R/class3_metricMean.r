setClass("metricMean", representation(mean.metric="data.frame"),
         contains=c("mmi", "oe"))

setValidity("metricMean", function(object){
  if(!(setequal(object@subsample$SampleID, object@oesubsample$SampleID))){return("Incompatible objects")} else
    TRUE
})

setMethod("initialize", "metricMean", function(.Object="metricMean", x="mmi", y="oe"){
  for(i in names(getSlots("mmi"))){
    slot(.Object, i) <- slot(x, i)
  }
  for(i in names(getSlots("oe"))){
    slot(.Object, i) <- slot(y, i)
  }
  .Object@mean.metric <- merge(y@oeresults, x@summary[, 1:3])
  .Object@mean.metric$CSCI <- apply(.Object@mean.metric[, c("OoverE", "MMI_Score")], 1, mean)
  .Object
})

setMethod("summary", "metricMean", function(object = "metricMean", report="core"){
  
  object@mean.metric$Count <- ddply(object@subsample, "SampleID", function(df)sum(df$Result))[, 2]
  object@mean.metric$Number_of_Iterations <- ifelse(object@mean.metric$Count >= 500, 20, 1)
  object@mean.metric$Pcnt_Ambiguous_Individuals <- object@ambiguous$individuals
  object@mean.metric$Pcnt_Ambiguous_Taxa <- object@ambiguous$taxa
  object@mean.metric$overall_flag <- ifelse(object@mean.metric$Count >=450 & object@mean.metric$Pcnt_Ambiguous_Individuals < 20,
                                            "Adequate", "Inadequate")
  object@mean.metric <- object@mean.metric[, c(1:2, 6:9, 3:5)]
  object@mean.metric$mmi_count_flag <- ifelse(object@mean.metric$Count >=450, "Adequate", "Inadequate")
  object@mean.metric$ambig_count_flag <- ifelse(object@mean.metric$Pcnt_Ambiguous_Individuals < 20, "Adequate", "Inadequate")
  object@mean.metric <- object@mean.metric[, c(1:6, 10:11, 7:9)]
  if(report == "core")return(object@mean.metric)
  
  
  names <- c("Shannon_Diversity", "Intolerant_PercentTaxa",
             "ToleranceValue", "Shredder_Taxa", "Clinger_Taxa", "Coleoptera_Taxa", 
             "Noninsect_PercentTaxa")
  if(report == "Supp1_mmi"){
    model <- object@modelprediction[, names]
    names(model) <- paste0(names(model), "_predicted")
    return(cbind(object@mean.metric[, c("SampleID", "StationCode", "MMI_Score")], object@metrics[, names], model, object@result))
  }

  predict <- predict(oe_stuff[[1]],newdata=object@predictors[,oe_stuff[[4]]],type='prob')
  colnames(predict) <- paste0("pGroup", 1:11)
  if(report == "Suppl1_grps"){
    return(cbind(object@mean.metric[, c("SampleID", "StationCode")], predict))
  }
  
  if(report == "Suppl1_OE.E"){
    return(cbind(object@mean.metric[, c("SampleID", "StationCode")], 
                 predict %*% apply(bugcal.pa,2,function(x){tapply(x,grps.final,function(y){sum(y)/length(y)})})))
  }
  
  if(report == "Suppl1_OE.O"){
    object@oesubsample$Replicate_mean <- apply(object@oesubsample[, paste("Replicate", 1:20)], 1, mean)
    dcast(object@oesubsample, SampleID + StationCode ~ STE, value.var="Replicate_mean", sum, na.rm=TRUE)
  }
})
  

setMethod("plot", "metricMean", function(x="metricMean"){
  load(system.file("data", "base_map.rdata", package="CSCI"))
  x@mean.metric$CSCIScore <- cut(x@mean.metric$CSCI, breaks=c(0, .4, .8, 1.5), labels=c("low", "medium", "high"))
  hscore <- cbind(x@mean.metric$CSCIScore, x@mean.metric$CSCI, x@predictors[, c("StationCode", "SampleID", "New_Lat", "New_Long")])
  names(hscore)[1] <- "CSCIScore"
  names(hscore)[2] <- "CSCI"
  ggmap(base_map) + 
    geom_point(data=hscore, aes(x=New_Long, y=New_Lat, colour=CSCI), size=4, alpha=.8) + 
    scale_color_continuous(low="red", high="green", name="CSCI Index Score") + labs(x="", y="")
})

