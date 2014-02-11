# BUILD

library(devtools)

build("hybridindex/", path="../CSCI_bin/", binary=TRUE)

install_github("bug_mmi", "mengeln", subdir="hybridindex")


stats <- as.matrix(stations[stations$SiteStatus == "Reference", csci_predictors])
nonrefstats <- as.matrix(stations[stations$SiteStatus != "Reference", csci_predictors])

ref <- mahalanobis(stats, colMeans(stats), cov(stats))^0.5
mean(ref) + sd(ref)*4

mahalanobis(nonrefstats, colMeans(stats), cov(stats))^0.5 < mean(ref) + sd(ref)*4