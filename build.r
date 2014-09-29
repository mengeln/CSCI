# BUILD

library(devtools)
library(roxygen2)

roxygenize("hybridindex")

build("hybridindex/", path="../CSCI_bin/", binary=TRUE)
install("hybridindex/", quick=FALSE)
# install_github("bug_mmi", "mengeln", subdir="hybridindex")

