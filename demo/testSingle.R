rm(list = ls())

library(dplyr)
library(discountingtools)
library(tidyr)

set.seed(65535)

dataFrame = data.frame(
  ids = 1:100,
  ks  = NA
)

dataFrame$ks  = rnorm(length(dataFrame$ids), 0.35, 0.125)
dataFrame$ks  = log(dataFrame$ks)

delays = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)

for (row in 1:nrow(dataFrame)) {
  ys = hyperbolicDiscountFunc(delays, dataFrame[row, "ks"]) + rnorm(length(delays),
                                                                    0,
                                                                    0.025)

  dataFrame[row, as.character(delays)] = ys
}

dataFrame.long = dataFrame %>%
  gather(Delay, Value, -ids, -ks) %>%
  mutate(Delay = as.numeric(Delay)) %>%
  mutate(Value = ifelse(Value < 0, 0, Value)) %>%
  mutate(Value = ifelse(Value > 1, 0, Value))

results = fitDDCurves(data = dataFrame.long,
            settings = list(Delays     = Delay,
                            Values     = Value,
                            Individual = ids),
            maxValue = 1,
            verbose  = TRUE) %>%
  dd_modelOptions(plan   = c("mazur")) %>%
  dd_screenOption(screen = FALSE) %>%
  dd_analyze(modelSelection = FALSE)

summary(results)

png(filename = "SingleModelEvaluation.png", width = 6, height = 6, res = 300, units = "in")

plot(results,
     logAxis = "x",
     position = "topright")

dev.off()
