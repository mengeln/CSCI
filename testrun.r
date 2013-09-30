library(devtools)
library(BMIMetrics)
library(plyr)
library(randomForest)
library(vegan)
library(stringr)
library(reshape2)
library(data.table)
library(ggplot2)
library(foreach)
library(doParallel)
files <- list.files("hybridindex/R/", full.names=TRUE)
invisible(sapply(files, source))



bugs <- read.csv("L:/CSCI_ME/bug_mmi/validation/tahoe.bugs.csv")
stations <- read.csv("L:/CSCI_ME/bug_mmi/validation/tahoe.stations.csv")

test <- CSCI(bugs, stations, purge=TRUE)

# load("L:/CSCI_ME/FamilyIndex/data/bugs.sub400.rdata")
# load("L:/CSCI_ME/FamilyIndex/data/station.rdata")
# bugs <- bugs.sub400 
# 
preds <- c("BDH_AVE", "CaO_Mean", "ELEV_RANGE", "KFCT_AVE",
           "P_MEAN", "LogWSA", "LPREM_mean", "New_Lat",
           "New_Long", "PPT_00_09", "PRMH_AVE", "S_Mean", 
           "SITE_ELEV", "SumAve_P", "TEMP_00_09")
# group <- cbind(as.character(unique(bugs$SampleID)),
#                rep(1:71, each = 50)[1:length(unique(bugs$SampleID))])
# bugs$group <- group[match(bugs$SampleID, group[, 1]), 2]
# bugssplit <- split(bugs, bugs$group)
# 
# test <- lapply(bugssplit, function(x){
#   res <- try(CSCI(x, stations[, c("StationCode", preds)], purge=TRUE))
#   gc()
#   res
# })


result <- test$core
stations$SiteStatus <- factor(stations$SiteStatus, levels(stations$SiteStatus)[c(2, 1, 3)])
result <- join(result, stations, by="StationCode")


means <- ddply(result, .(SiteStatus), function(x)mean(x$CSCI))
thres <- means$V1[means$SiteStatus=="Reference"] - 
  sd(result$CSCI[result$SiteStatus == "Reference"])

ggplot(result, aes(SiteStatus, CSCI)) + geom_boxplot() +
  geom_hline(yintercept=thres, colour="red") +
  facet_wrap(~PSA9c_1987)
  theme_bw()

ggplot(result, aes(CSCI, MaxOfCOND)) + geom_point()

