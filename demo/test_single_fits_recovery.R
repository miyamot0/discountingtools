# Example: Simulated individual fits and parameter recovery

rm(list = ls())

set.seed(65535)

library(dplyr)
library(discountingtools)

dataFrame = data.frame(
  ids = 1:100,
  ks  = NA
)

dataFrame$ks  = rnorm(length(dataFrame$ids), 0.35, 0.125)
dataFrame$ks  = log(dataFrame$ks)

delays = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)

for (row in seq_len(nrow(dataFrame))) {
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
            verbose  = FALSE) |>
  dd_modelOptions(plan   = c("mazur")) |>
  dd_screenOption(screen = FALSE) |>
  dd_analyze(modelSelection = FALSE)

summary(results)

data_frame_results <- summary(results)

par(mfrow = c(1, 2))

plot(results,
     logAxis = "x",
     position = "topright")

plot(data_frame_results$Mazur.Lnk,
     dataFrame$ks,
     main = "Fitted vs. Simulated",
     ylab = "Fitted",
     xlab = "Simulated",
     ylim = c(-2.5, 0),
     xlim = c(-2.5, 0))

lines(x = c(-2.5, 0),
      y = c(-2.5, 0))
