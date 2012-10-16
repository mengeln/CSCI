library(randomForest)
library(plyr)
library(doParallel)
library(vegan)
library(ggmap)
library(reshape2)
library(RSQLite)
library(stringr)
        
lapply(list.files("r"), function(file)source(paste("r/", file, sep="")))

install.packages("P:/PartTimers/MarkEngeln/bug_mmi/hybridindex_1.3.tar.gz", repos = NULL, type="source")
library(hybridindex)

val.bugs2 <- read.csv("P:/PartTimers/MarkEngeln/bug_mmi/validation/val.bugs2.csv", stringsAsFactors=F)
val.pred2 <- read.csv("P:/PartTimers/MarkEngeln/bug_mmi/validation/val.pred2.csv", stringsAsFactors=F)

mmidata <- new("mmi", bugdata=val.bugs2, predictors=val.pred2)
mmiresults <- score(mmidata)

oedata <- new("oe", bugdata=val.bugs2, predictors=val.pred2)
oeresults <- score(oedata)

hybrid <- new("metricMean", mmiresults, oeresults)
View(summary(hybrid, report="detailed"))
