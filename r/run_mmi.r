source("r/classes.r")
bugs <- read.csv("validation/bugs_mmi.csv")
pred <- read.csv("validation/predictors_mmi.csv")

mmi_data <- new("mmi", bugdata = bugs, predictors = pred)

mmiresults <- score(mmi_data)
mmi.results.table <- summary(mmiresults)
View(mmi.results.table)

oe_data <- new("oe", bugdata = bugs, predictors = pred)

oeresults <- score(oe_data)
oe.results.table <- summary(oeresults)
View(oe.results.table)
