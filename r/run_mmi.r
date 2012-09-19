source("r/classes.r")
bugs <- read.csv("validation/bugs_mmi.csv")
pred <- read.csv("validation/predictors_mmi.csv")

validation_data <- new("mmi", bugdata = bugs, predictors = pred)

results <- mean(validation_data)
results.table <- summary(results)
View(results.table)
