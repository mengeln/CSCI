library(car)
library(plyr)
library(randomForest)

data(Metrics.RFModels, package = "CSCI")
stations <- read.csv("stations.csv")#modify as needed
stations$LogWSA<-log10(stations$AREA_SQKM)
stations$Log_P_MEAN<-  log10(stations$P_MEAN + 0.0001)
stations$Log_N_MEAN<-  log10(stations$N_MEAN + 0.00001)
refcal <- stations[stations$SiteSet == "RefCal", ]
predictors <- c("LogWSA", "New_Long", "New_Lat", "SITE_ELEV", "ELEV_RANGE", 
                "TEMP_00_09", "PPT_00_09", "SumAve_P", "KFCT_AVE", "BDH_AVE", 
                "MgO_Mean", "Log_P_MEAN", "CaO_Mean", "PRMH_AVE", "S_Mean", "PCT_SEDIM", 
                "LPREM_mean", "Log_N_MEAN")

###rank importance from models###
r <- lapply(lapply(final.forests, importance), function(df){
  df <- as.data.frame(df)
  data.frame(Variable = row.names(df), rank = 18:(18-nrow(df) +1))
})
ranks <- suppressWarnings(Reduce(function(x, y)merge(x, y, by = "Variable", all=TRUE), r))
averageImportance <- data.frame(variable = ranks$Variable, 
                                Importance = apply(ranks[, -1], 1, mean, na.rm=TRUE))

msel <- lapply(lapply(final.forests, importance), function(df){
  df <- as.data.frame(df)
  data.frame(Variable = row.names(df), mse = df[, 1])
})
mse <- suppressWarnings(Reduce(function(x, y)merge(x, y, by = "Variable", all=TRUE), msel))
mse[is.na(mse)] <- 0

averageMSE <- data.frame(variable = mse$Variable,
                         Importance = apply(mse[, -1], 1, mean))


#function to check if new data are within ellipse
isInside <- function(newdata, refcaldata = refcal, averageImp = averageImportance, confidence = 0.99){
  predictors <- c("LogWSA", "New_Long", "New_Lat", "SITE_ELEV", "ELEV_RANGE", 
                  "TEMP_00_09", "PPT_00_09", "SumAve_P", "KFCT_AVE", "BDH_AVE", 
                  "MgO_Mean", "Log_P_MEAN", "CaO_Mean", "PRMH_AVE", "S_Mean", "PCT_SEDIM", 
                  "LPREM_mean", "Log_N_MEAN")
  refcal <- scale(refcaldata[, predictors])
  
  refcal <- refcal[, averageImp$variable] * averageImp$Importance
  pcrefcal <- princomp(refcal)
  pcdata <- data.frame(axis1 = pcrefcal$scores[, 1], axis2 = pcrefcal$scores[, 2],
                       axis3 = pcrefcal$scores[, 3])
  suppressWarnings(.Call("R_GD_nullDevice", PACKAGE = "grDevices"))
  test <- dataEllipse(as.matrix(pcdata[, 1:2]), levels=c(confidence))
  test2 <- dataEllipse(as.matrix(pcdata[, 2:3]), levels=c(confidence))
  dev.off()
  #Calculate the ellipse radii
  axes <- sort(unlist(union(colwise(max)(as.data.frame(test2)), colwise(max)(as.data.frame(test)))))
  xradius <-  axes[4]
  yradius <- axes[2]
  zradius <- axes[1]
  
  prep <- scale(newdata[, predictors])
  prep <- prep[, averageImp$variable] * averageImp$Importance
  pdata <- predict(pcrefcal, prep)
  x <- pdata[, 1]
  y <- pdata[, 2]
  z <- pdata[, 3]  
  data.frame(outlier = ifelse(
    (x^2)/(xradius^2) + (y^2)/(yradius^2) + (z^2)/(zradius^2) <= 1,
    0, 1), row.names = newdata$StationCode)

}

outliercheck <- isInside(stations[stations$SiteSet %in% c("RefVal", "StressVal"), ], 
                         averageImp = averageMSE, confidence =0.99)
sum(outliercheck$outlier) #number of outliers in the non-RefCal data

psaregion <- stations$PSA9c_1987[stations$StationCode %in% row.names(outliercheck)[outliercheck$outlier == 1]]
table(psaregion)/length(psaregion)

table(psaregion)/length(psaregion) - table(stations$PSA9c_1987[stations$SiteSet  %in% c("RefVal", "StressVal")])/
  length(stations$PSA9c_1987[stations$SiteSet  %in% c("RefVal", "StressVal")]) 
  


