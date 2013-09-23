library(CSCI)
library(foreach)
library(doParallel)


bugs <- read.csv("P:/MarkEngeln/from Rafi/AllData/AllBugs.csv")
stations <- read.csv("P:/MarkEngeln/from Rafi/AllData/AllStations.csv")

s <- proc.time()
cl <- makeCluster(8)
registerDoParallel(cl)

groups <- 1:length(unique(bugs$SampleID))
groups <- split(groups, c(rep(1:50, each=70), 
                          rep(51, 18)))


listresult <- foreach(i = groups, .packages="CSCI") %dopar% {
  CSCI(bugs[bugs$SampleID %in% unique(bugs$SampleID)[i], ], stations)
}
stopCluster(cl)

result <- do.call(function(...)mapply(rbind, ...), listresult)
print(proc.time() - s)
