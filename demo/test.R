rm(list = ls())

library(dplyr)
library(discountingtools)

dat <- data.frame(X = c(1, 30, 180, 540, 1080, 2160, 4320, 8640,
                        1, 30, 180, 540, 1080, 2160, 4320, 8640,
                        1, 30, 180, 540, 1080, 2160, 4320, 8640),
                  Y = c(1000, 880, 680, 450, 400, 350, 200, 100,
                        1000, 980, 780, 550, 500, 450, 100, 50,
                        1000, 780, 480, 350, 200, 150, 100, 50),
                  ids = c(2, 2, 2, 2, 2, 2, 2, 2,
                          4, 4, 4, 4, 4, 4, 4, 4,
                          6, 6, 6, 6, 6, 6, 6, 6))

results = fitDDCurves(data = dat,
            settings = list(Delays     = X,
                            Values     = Y,
                            Individual = ids),
            maxValue = 1000) %>%
  dd_modelOptions(plan = c("bleichrodt",
                           "ebertprelec",
                           "exponential",
                           "greenmyerson",
                           "laibson",
                           "mazur",
                           "noise",
                           "rachlin",
                           "rodriguezlogue")) %>%
  dd_metricOptions(metrics = c("lned50",
                               "mbauc",
                               "logmbauc")) %>%
  dd_analyze()

summary(results)

