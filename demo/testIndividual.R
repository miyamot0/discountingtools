rm(list = ls())

library(dplyr)
library(discountingtools)
library(tidyr)

set.seed(65535)

dataFrame = data.frame(
  ids = 1:10,
  ks  = NA
)

dataFrame$ks  = rnorm(length(dataFrame$ids), 0.09, 0.05)
dataFrame$ks  = log(dataFrame$ks)

delays = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)

for (row in 1:nrow(dataFrame)) {
  ys = hyperbolicDiscountFunc(delays, dataFrame[row, "ks"]) + rnorm(length(delays),
                                                                    0,
                                                                    0.05)

  dataFrame[row, as.character(delays)] = ys
}

dataFrame[,"1"]    = ifelse(dataFrame[,"1"] > 1,    1, dataFrame[,"1"])
dataFrame[,"8640"] = ifelse(dataFrame[,"8640"] < 0, 0, dataFrame[,"8640"])

dataFrame.long = dataFrame %>%
  gather(Delay, Value, -ids, -ks) %>%
  mutate(Delay = as.numeric(Delay))

str(dataFrame.long)

results = fitDDCurves(data = dataFrame.long,
            settings = list(Delays     = Delay,
                            Values     = Value,
                            Individual = ids),
            maxValue = 1) %>%
  dd_modelOptions(plan = c("mazur",
                           "bleichrodt",
                           "ebertprelec",
                           "exponential",
                           "greenmyerson",
                           "laibson",
                           "noise",
                           "rachlin",
                           "rodriguezlogue")) %>%
  dd_metricOptions(metrics = c("lned50",
                               "mbauc",
                               "logmbauc")) %>%
  dd_analyze()

png(filename = "MultiModelEvaluation.png", width = 6, height = 6, res = 300, units = "in")

plot(results,
     logAxis = "x",
     position = "topright")

dev.off()
