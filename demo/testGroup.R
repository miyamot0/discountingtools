rm(list = ls())

library(dplyr)
library(discountingtools)
library(tidyr)

set.seed(65535)

dataFrame = data.frame(
  ids = 1:10,
  ks  = NA,
  grp = "A"
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

dataFrame2 = data.frame(
  ids = 11:20,
  ks  = NA,
  grp = "B"
)

dataFrame2$ks  = rnorm(length(dataFrame2$ids), 0.018, 0.005)
dataFrame2$ks  = log(dataFrame2$ks)

for (row in 1:nrow(dataFrame2)) {
  ys = hyperbolicDiscountFunc(delays, dataFrame2[row, "ks"]) + rnorm(length(delays),
                                                                    0,
                                                                    0.05)

  dataFrame2[row, as.character(delays)] = ys
}

dataFrame = rbind(dataFrame,
                  dataFrame2)

dataFrame.long = dataFrame %>%
  gather(Delay, Value, -ids, -ks, -grp) %>%
  mutate(Delay = as.numeric(Delay))

dataFrame.long[,"Value"] = ifelse(dataFrame.long[,"Value"] > 1, 1, dataFrame.long[,"Value"])
dataFrame.long[,"Value"] = ifelse(dataFrame.long[,"Value"] < 0, 0, dataFrame.long[,"Value"])

str(dataFrame.long)

results = fitDDCurves(data = dataFrame.long,
            settings = list(Delays     = Delay,
                            Values     = Value,
                            Individual = ids,
                            Group      = grp),
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

vecGroups = unique(results$data[,as.character(results$settings['Group'])])

plot(results, logAxis = "x", position = "topright", which = "group")
