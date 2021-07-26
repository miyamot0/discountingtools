rm(list = ls())

library(dplyr)
library(discountingtools)
library(tidyr)

set.seed(65535)

dataFrame = data.frame(
  ids = 1:30,
  ks  = NA,
  grp = "Group A"
)

dataFrame$ks  = rnorm(length(dataFrame$ids), 0.125, 0.08)
dataFrame$ks  = log(dataFrame$ks)

delays = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)

for (row in 1:nrow(dataFrame)) {
  ys = hyperbolicDiscountFunc(delays, dataFrame[row, "ks"]) + rnorm(length(delays),
                                                                    0,
                                                                    0.05)

  dataFrame[row, as.character(delays)] = ys
}

dataFrame2 = data.frame(
  ids = 31:60,
  ks  = NA,
  grp = "Group B"
)

dataFrame2$ks  = rnorm(length(dataFrame2$ids), 0.225, 0.035)
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

results = fitDDCurves(data = dataFrame.long,
            settings = list(Delays     = Delay,
                            Values     = Value,
                            Individual = ids,
                            Group      = grp),
            maxValue = 1,
            verbose  = TRUE) %>%
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

summary(results)

png(filename = "ED50.png", width = 6, height = 6, res = 300, units = "in")

plot(results, which = "ED50")

dev.off()

png(filename = "MBAUC.png", width = 6, height = 6, res = 300, units = "in")

plot(results, which = "MBAUC")

dev.off()

png(filename = "Log10MBAUC.png", width = 6, height = 6, res = 300, units = "in")

plot(results, which = "Log10MBAUC")

dev.off()

