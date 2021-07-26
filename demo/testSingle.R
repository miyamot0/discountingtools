rm(list = ls())

library(dplyr)
library(discountingtools)
library(tidyr)

set.seed(65535)

dataFrame = data.frame(
  ids = 1:30,
  ks  = NA
)

dataFrame$ks  = rnorm(length(dataFrame$ids), 0.125, 0.05)
dataFrame$ks  = log(dataFrame$ks)

delays = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)

for (row in 1:nrow(dataFrame)) {
  ys = hyperbolicDiscountFunc(delays, dataFrame[row, "ks"]) + rnorm(length(delays),
                                                                    0,
                                                                    0.05)

  dataFrame[row, as.character(delays)] = ys
}

dataFrame.long = dataFrame %>%
  gather(Delay, Value, -ids, -ks) %>%
  mutate(Delay = as.numeric(Delay))

dataFrame.long[,"Value"] = ifelse(dataFrame.long[,"Value"] > 1, 1, dataFrame.long[,"Value"])
dataFrame.long[,"Value"] = ifelse(dataFrame.long[,"Value"] < 0, 0, dataFrame.long[,"Value"])

results = fitDDCurves(data = dataFrame.long,
            settings = list(Delays     = Delay,
                            Values     = Value,
                            Individual = ids),
            maxValue = 1) %>%
  dd_modelOptions(plan = c("mazur")) %>%
  dd_analyze(modelSelection = FALSE)

summary(results)

png(filename = "SingleModelEvaluation.png", width = 6, height = 6, res = 300, units = "in")

plot(results,
     logAxis = "x",
     position = "topright")

dev.off()
