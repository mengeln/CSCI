cleanData <- function(data, purge=TRUE){
  meta <- loadMetaData()
  nomatch <- !(data$FinalID %in% meta$FinalID)
  if(!purge)
    bad <- data[nomatch, ]
  
  data <- data[!nomatch, ]

  lsc <- with(data, paste(FinalID, LifeStageCode)) %in%
    with(meta, paste(FinalID, LifeStageCode))
  
  data$LifeStageCode[!lsc] <- meta$DefaultLifeStage[match(data$FinalID[!lsc], 
                                                       meta$FinalID)]
  
  if(!purge)
    data <- rbind(data, bad)
  
  data
}