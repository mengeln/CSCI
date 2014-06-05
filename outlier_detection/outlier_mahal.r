load("L:/CSCI_ME/FamilyIndex/data/station.rdata")

csci_predictors <- c("BDH_AVE", "CaO_Mean", "ELEV_RANGE", "KFCT_AVE",
                     "P_MEAN", "LogWSA", "LPREM_mean", "New_Lat",
                     "New_Long", "PPT_00_09", "PRMH_AVE", "S_Mean", 
                     "SITE_ELEV", "SumAve_P", "TEMP_00_09")

stations <- na.omit(stations[, c("StationCode", "SiteStatus", csci_predictors)])


refsites <- as.matrix(stations[stations$SiteStatus == "Reference",# really should be just refCal
                               csci_predictors])
nonrefsites <- as.matrix(stations[stations$SiteStatus != "Reference", csci_predictors])
rownames(nonrefsites) <- stations$StationCode[stations$SiteStatus != "Reference"]

# Calculate distance reference site distance from their center
ref <- mahalanobis(refsites , colMeans(refsites), cov(refsites))^0.5
thres <- mean(ref) + sd(ref)*4 # cutoff for outlier detection

outlier <- mahalanobis(nonrefsites, # new data
                       colMeans(refsites), # reference center
                       cov(refsites) # ref covariance
                       )^0.5 > thres


sum(outlier) # out of nrow(nonrefsites)
names(outlier[outlier])  # outliers by StationCode

