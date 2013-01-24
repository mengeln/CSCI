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
  load(system.file("data", "oe_stuff.rdata", package="CSCI"))
  arglist <- c("core", "Suppl1_mmi", "Suppl1_grps", "Suppl1_OE", "Suppl2_OE", "Suppl2_mmi")
  report <- match.arg(report, c(arglist, "all"), several.ok=TRUE)
  if(report == "all")report <- arglist
  reportlist <- list()
  add <-  function(obj){
    reportlist[[length(reportlist)+1]] <- obj
    reportlist
  }
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
  
  if("core" %in% report){
    object@mean.metric$E <- object@fulliterations[[1]]$E
    object@mean.metric$Mean_O <- apply(
      Reduce(function(x,y)cbind(x, y$O), object@fulliterations, init=rep(NA, nrow(object@fulliterations[[1]]))),
      1, mean, na.rm=TRUE)
    reportlist <- add(object@mean.metric[, c(names(object@mean.metric)[1:8], "E", "Mean_O", "OoverE", "MMI_Score",
                                             "CSCI")])
    names(reportlist) <- "core"
  }
  
  names <- c("Shannon_Diversity", "Intolerant_PercentTaxa",
             "ToleranceValue", "Shredder_Taxa", "Clinger_Taxa", "Coleoptera_Taxa", 
             "Noninsect_PercentTaxa", "CFCG_Taxa")
  if("Suppl1_mmi" %in% report){
    model <- object@modelprediction[, names]
    names(model) <- paste0(names(model), "_predicted")
    reportlist <- add(cbind(object@mean.metric[, c("SampleID", "StationCode", "MMI_Score")],
                            object@metrics[, names], model, object@result))
    names(reportlist)[length(reportlist)] <- "Supp1_mmi"
  }

  predict <- predict(oe_stuff[[1]],newdata=object@predictors[,oe_stuff[[4]]],type='prob')
  colnames(predict) <- paste0("pGroup", 1:11)
  if("Suppl1_grps" %in% report){
    reportlist <- add(cbind(object@mean.metric[, c("SampleID", "StationCode")], predict))
    names(reportlist)[length(reportlist)] <- "Suppl1_grps"
  }
  
  if("Suppl1_OE" %in% report){
    E <- cbind(object@mean.metric[, c("SampleID", "StationCode")], 
                 predict %*% apply(oe_stuff[[2]],2,function(x){tapply(x,oe_stuff[[3]],function(y){sum(y)/length(y)})}))
    object@oesubsample$Replicate_mean <- apply(object@oesubsample[, paste("Replicate", 1:20)], 1, mean)
    O <- dcast(object@oesubsample, SampleID + StationCode ~ STE, value.var="Replicate_mean", sum, na.rm=TRUE)
    ids <- c("SampleID", "StationCode")
    result <- merge(melt(E, id.vars=ids), melt(O, id.vars=ids), by=c("variable", ids), all.x=TRUE)
    names(result) <- c("OTU", "SampleID", "StationCode", "CaptureProb", "Mean Observed")
    result$Observed[is.na(result$"Mean Observed")] <- 0
    reportlist <- add(result[, c("SampleID", "StationCode", "OTU", "CaptureProb", "Mean Observed")])
    names(reportlist)[length(reportlist)] <- "Suppl1_OE"
  }
  if("Suppl2_OE" %in% report){
    E <- cbind(object@mean.metric[, c("SampleID", "StationCode")], 
               predict %*% apply(oe_stuff[[2]],2,function(x){tapply(x,oe_stuff[[3]],function(y){sum(y)/length(y)})}))
    E <- melt(E, id.vars=c("SampleID", "StationCode"))
    names(E)[3:4] <- c("OTU", "CaptureProb") 
    O <- object@oesubsample[, c("SampleID", "StationCode", "STE", paste("Replicate", 1:20))]
    O <- dcast(melt(O, id.vars=c("SampleID", "StationCode", "STE")), SampleID + StationCode + STE ~ variable,
                    value.var="value", fun.aggregate=sum)
    names(O)[3] <- c("OTU")
    result <- merge(E, O, by=c("SampleID", "StationCode", "OTU"))
    reportlist <- add(result)
    names(reportlist)[length(reportlist)] <- "Suppl2_OE"
  }
  if("Suppl2_mmi" %in% report){
    cmmi <- melt(object@metrics, id.vars="SampleID")
    cmmi$variable <- as.character(cmmi$variable)
    cmmi$replicate <- substr(cmmi$variable, nchar(match.arg(cmmi$variable, names, TRUE))+1, nchar(cmmi$variable))
    cmmi$metric <- match.arg(cmmi$variable, names, TRUE)
    cmmi$replicate[cmmi$replicate == ""] <- "Mean"
    reportlist <- add(cmmi[, c("SampleID", "metric", "replicate", "value")])
    names(reportlist)[length(reportlist)] <- "Suppl2_mmi"
  }
  if(length(reportlist)==1)transform(reportlist) else reportlist
})
  

setMethod("plot", "metricMean", function(x="metricMean"){
  load(system.file("data", "base_map.rdata", package="CSCI"))
  x@mean.metric$CSCIScore <- cut(x@mean.metric$CSCI, breaks=c(0, .4, .8, 1.5), labels=c("low", "medium", "high"))
  hscore <- cbind(x@mean.metric$CSCIScore, x@mean.metric$CSCI, x@predictors[, c("StationCode", "SampleID", "New_Lat", "New_Long")])
  names(hscore)[1] <- "CSCIScore"
  names(hscore)[2] <- "CSCI"
  ggmap(base_map) + 
    geom_point(data=hscore, aes(x=New_Long, y=New_Lat, colour=CSCI), size=4, alpha=.8) + 
    scale_color_continuous(low="red", high="green", name="CSCI Score") + labs(x="", y="")
})

