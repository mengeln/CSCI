library(CSCI)

stations <- read.csv("C:/Documents and Settings/gisuser/Desktop/demo.stations.csv")
bugs <- read.csv("C:/Documents and Settings/gisuser/Desktop/demo.bugs.csv")
bugs$SampleID <- with(bugs, paste0(StationCode, SampleDate, BugCode))


CSCI <- function (bugs, stations) {
  mmi <- new("mmi", bugs, stations)
  mmi_s <- score(mmi)
  
  oe <- new("oe", bugs, stations)
  oe_s <- score(oe)
  
  res <- new("metricMean", mmi_s, oe_s)
  summary(res, report="all")
}

report <- CSCI(bugs, stations)
View(report$Suppl2_OE)
View(report$Suppl2_mmi)




