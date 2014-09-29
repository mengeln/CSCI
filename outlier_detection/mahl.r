
stations<-read.csv("../Desktop/stations.csv")
stations$LogWSA<-log10(stations$AREA_SQKM)
stations$Log_P_MEAN<-  log10(stations$P_MEAN + 0.0001)
stations$Log_N_MEAN<-  log10(stations$N_MEAN + 0.00001)

predictors <- c("LogWSA", "New_Long", "New_Lat", "SITE_ELEV", "ELEV_RANGE", 
                "TEMP_00_09", "PPT_00_09", "SumAve_P", "KFCT_AVE", "BDH_AVE", 
                "MgO_Mean", "Log_P_MEAN", "CaO_Mean", "PRMH_AVE", "S_Mean", "PCT_SEDIM", 
                "LPREM_mean", "Log_N_MEAN")

mahalDist <- function (stations, predictors) {
  load("../Desktop/mahalData.rdata")
  
  covlist <- lapply(split(datmat, grps.final), cov)
  covlist.df <- lapply(1:length(grpsize), function(i) (grpsize[i]-1)*covlist[[i]])
  covpool <- Reduce(`+`, covlist.df)
  covpool <- covpool/(sum(grpsize) - ngrps)
  
  covinvMatrix <- solve(covpool)
  
  grpmns <- apply(datmat,2,function(x)tapply(x,grps.final,mean))
  
  res <- t(sapply(split(stations[, predictors], 1:nrow(stations)), function(row){
    dist <- mahalanobis(grpmns, as.matrix(row), covinvMatrix, inverted=T)
  }))
  row.names(res) <- stations$StationCode
  res
}

outlier.flag <- function(dist){
  dff <- ncol(dist) - 1
  
  crit.01 <- qchisq(1 - 0.01, df = dff)
  crit.05 <- qchisq(1 - 0.5, df = dff)
  
  min <- apply(dist, 1, min)
  
  data.frame(crit_05 = ifelse(min > crit.05, 1, 0),
             crit_01 = ifelse(min > crit.01, 1, 0),
             row.names = row.names(dist))
}

distanceM <- mahalDist(stations, predictors)

outlier <- outlier.flag(distanceM)
colMeans(outlier[row.names(outlier) %in% as.character(stations$StationCode[stations$SiteSet == "RefCal"]),])            
                     