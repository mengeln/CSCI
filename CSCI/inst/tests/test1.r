library(testthat)
library(BMIMetrics)
library(reshape2)
library(randomForest)

lapply(paste0("L:/CSCI_ME/bug_mmi/hybridindex/R/", list.files("L:/CSCI_ME/bug_mmi/hybridindex/R/")), source)
data(bugs_stations, package="CSCI")

context("High level run through")

test_that("Runs through example data", {
  expect_that(length(CSCI(bugs_stations[[1]], bugs_stations[[2]])) == 6 , is_true())
})


test_that("Works with just AREA_SQKM or LogWSA", {
  expect_that(local({bugs_stations[[2]] <- bugs_stations[[2]][, c(-2)]
                     length(CSCI(bugs_stations[[1]], bugs_stations[[2]])) == 6}), is_true())
  expect_that(local({bugs_stations[[2]] <- bugs_stations[[2]][, c(-20)]
                     length(CSCI(bugs_stations[[1]], bugs_stations[[2]])) == 6}), is_true())
})
