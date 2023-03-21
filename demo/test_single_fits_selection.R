# Example: Simulated individual fits and parameter recovery

rm(list = ls())

set.seed(65535)

library(tidyverse)
library(discountingtools)

dataFrame = data.frame(
  ids = 1:100,
  ks  = NA
)

dataFrame$ks  = rnorm(length(dataFrame$ids), 0.07, 0.03)
dataFrame$ks  = log(dataFrame$ks)

delays = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)

for (row in seq_len(nrow(dataFrame))) {
  ys = dd_discount_func_mazur(delays, dataFrame[row, "ks"]) + rnorm(length(delays),
                                                                    0,
                                                                    0.025)

  dataFrame[row, as.character(delays)] = ys
}

dataFrame.long = dataFrame %>%
  gather(Delay, Value, -ids, -ks) %>%
  mutate(Delay = as.numeric(Delay)) %>%
  mutate(Value = ifelse(Value < 0, 0, Value)) %>%
  mutate(Value = ifelse(Value > 1, 0, Value))

results = fit_dd_curves(
  data = dataFrame.long,
  settings = list(Delays     = Delay,
                  Values     = Value,
                  Individual = ids),
  maxValue = 1,
  plan = c('mazur', 'exponential', 'rachlin'),
  verbose  = TRUE) |>
dd_analyze(modelSelection = TRUE)

data_frame_results <- summary(results)

png(filename = "../man/figures/single_fits_selection.png",
    width = 6,
    height = 4,
    res = 600,
    units = "in")

par(mfrow = c(1, 1))

plot(results,
     logAxis = "x",
     position = "topright")

dev.off()
