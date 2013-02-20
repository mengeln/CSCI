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

setMethod("summary", "metricMean", function(object = "metricMean", report="all"){
  load(system.file("data", "oe_stuff.rdata", package="CSCI"))
  load(system.file("data", "extent.rdata", package="CSCI"))
  arglist <- c("core", "Suppl1_mmi", "Suppl1_grps", "Suppl1_OE", "Suppl2_OE", "Suppl2_mmi")
  report <- match.arg(report, c(arglist, "all"), several.ok=TRUE)
  if(report == "all")report <- arglist
  reportlist <- list()
  add <-  function(obj){
    reportlist[[length(reportlist)+1]] <- obj
    reportlist
  }
  object@mean.metric$Count <- ddply(object@subsample, "SampleID", function(df)sum(df$Result))[, 2]
  object@mean.metric$Number_of_MMI_Iterations <- ifelse(object@mean.metric$Count >= 500, 20, 1)
  object@mean.metric$Number_of_OE_Iterations <- 
    ifelse(object@mean.metric$Count - object@mean.metric$Count*(object@ambiguous$individuals/100) >= 400, 20, 1)
  object@mean.metric$Pcnt_Ambiguous_Individuals <- object@ambiguous$individuals
  object@mean.metric$Pcnt_Ambiguous_Taxa <- object@ambiguous$taxa
  object@mean.metric$overall_flag <- ifelse(object@mean.metric$Count >=450 & object@mean.metric$Pcnt_Ambiguous_Individuals < 20,
                                            "Adequate", "Inadequate")
  object@mean.metric <- object@mean.metric[, c(1:2, 6:9, 3:5)]
  object@mean.metric$mmi_count_flag <- ifelse(object@mean.metric$Count >=450, "Adequate", "Inadequate")
  object@mean.metric$ambig_count_flag <- ifelse(object@mean.metric$Pcnt_Ambiguous_Individuals < 20, "Adequate", "Inadequate")
  
  predcheck <- function(data){
    res <- apply(
      sapply(names(extent), function(col){
        sapply(data[, col], function(x){extent[1, col] > x | extent[2, col] < x})
      }), 1, function(x)paste(names(which(x)), collapse=", "))
    res[res == ""] <- "All within range"
    res
  }
  object@mean.metric$Model_Experience_Flag <- predcheck(object@predictors)
    
  if("core" %in% report){
    object@mean.metric$E <- object@fulliterations[[1]]$E
    object@mean.metric$Mean_O <- object@mean.metric$E * object@mean.metric$OoverE
    reportlist <- add(object@mean.metric[, c("StationCode", "SampleID", "Count", "Number_of_MMI_Iterations", 
                                             "Number_of_OE_Iterations", "Pcnt_Ambiguous_Individuals", "mmi_count_flag",
                                             "ambig_count_flag", "Model_Experience_Flag", 
                                             "E", "Mean_O", "OoverE", 
                                             "MMI_Score", "CSCI")])
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
    names(reportlist)[length(reportlist)] <- "Suppl1_mmi"
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
    result$"Mean Observed"[is.na(result$"Mean Observed")] <- 0
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
    
    x <- result
    x[, 5:24] <- colwise(function(x)ifelse(x > 0, 1, 0))(x[, 5:24])
    
    
    test <- sapply(paste("Replicate", 1:20), function(rep){
      daply(x, "SampleID", function(df){
        captable <- reportlist$Suppl1_OE[reportlist$Suppl1_OE$SampleID == unique(df$SampleID), ]
        ingroup <- as.character(captable$OTU[captable$CaptureProb > 0.5])
        sum(df[df$OTU %in% ingroup, rep] > 0)/
          sum(captable$CaptureProb[captable$OTU %in% ingroup])
      })
    })
    test <- data.frame("SampleID" = row.names(test), "StationCode" = 
                         reportlist$Suppl1_OE$StationCode[match(row.names(test), reportlist$Suppl1_OE$SampleID)],
                       "OTU" = "OoverE", CaptureProb = NA, test, row.names=NULL)
    names(test)[5:24] <- paste("Replicate", 1:20)
    reportlist <- add(rbind(result, test))
    names(reportlist)[length(reportlist)] <- "Suppl2_OE"
  }
  if(all(c("Suppl2_mmi", "Suppl1_mmi") %in% report)){
    load(system.file("data", "maxmin.rdata",  package="CSCI"))
    
    cmmi <- melt(object@metrics, id.vars="SampleID")
    cmmi$variable <- as.character(cmmi$variable)
    cmmi$replicate <- substr(cmmi$variable, nchar(match.arg(cmmi$variable, names, TRUE))+1, nchar(cmmi$variable))
    cmmi$metric <- match.arg(cmmi$variable, names, TRUE)
    cmmi$replicate[cmmi$replicate == ""] <- "Mean"
    cmmi <- cmmi[cmmi$replicate != "Mean", c("SampleID", "metric", "replicate", "value")]
    
    
    prediction <- melt(reportlist$Suppl1_mmi, id.vars=c("StationCode", "SampleID"))
    prediction <- subset(prediction, grepl("predicted", variable))
    prediction$variable <- as.character(prediction$variable)
    prediction$variable <- substr(prediction$variable, 1, nchar(prediction$variable) - 10)
    names(prediction)[3:4] <- c("metric", "predicted_value")
    x <- merge(cmmi, prediction)
    x$score <- mapply(function(value, predict, metric){
      result <- (value - predict - maxmin[1, metric])/(maxmin[2, metric] - maxmin[1, metric])
      result <- ifelse(result > 1, 1, ifelse(
        result < 0, 0, result))
      result
    }, x$value, x$predicted_value, x$metric)
    x$replicate <- as.numeric(x$replicate)
    x <- arrange(x[, c("StationCode", "SampleID", "metric", "replicate", "value", "predicted_value", "score")],
                 SampleID, metric, replicate)
    
    reportlist <- add(x)
    names(reportlist)[length(reportlist)] <- "Suppl2_mmi"
  }
  if(length(reportlist)==1)transform(reportlist) else reportlist
})
  

