# Example: Simulated group fits and parameter recovery

# Note: Will take a minute to run

rm(list = ls())

library(tidyverse)
library(discountingtools)

set.seed(65535)

n_per_group <- 50

data_frame = data.frame(
  ids = seq_len(n_per_group),
  ks  = NA,
  grp = "Group A"
)

data_frame$ks  = rnorm(length(data_frame$ids), 0.35, 0.125)
data_frame$ks  = log(data_frame$ks)

delays = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)

data_frame$auc = dd_mbauc_mazur(1, data_frame$ks, min(delays), max(delays))

for (row in 1:nrow(data_frame)) {
  ys = dd_discount_func_mazur(delays, data_frame[row, "ks"]) + rnorm(length(delays),
                                                                    0,
                                                                    0.05)

  data_frame[row, as.character(delays)] = ys
}

data_frame2 = data.frame(
  ids = 50 + seq_len(n_per_group),
  ks  = NA,
  grp = "Group B"
)

data_frame2$ks  = rnorm(length(data_frame2$ids), 0.075, 0.035)
data_frame2$ks  = log(data_frame2$ks)

data_frame2$auc = dd_mbauc_mazur(1, data_frame2$ks, min(delays), max(delays))

for (row in 1:nrow(data_frame2)) {
  ys = dd_discount_func_mazur(delays, data_frame2[row, "ks"]) + rnorm(length(delays),
                                                                    0,
                                                                    0.025)

  data_frame2[row, as.character(delays)] = ys
}

data_frame = rbind(data_frame,
                  data_frame2)

data_frame_long = data_frame %>%
  gather(Delay, Value, -ids, -ks, -grp, -auc) %>%
  mutate(Delay = as.numeric(Delay)) %>%
  mutate(Value = ifelse(Value < 0, 0, Value)) %>%
  mutate(Value = ifelse(Value > 1, 0, Value))

results = fit_dd_curves(data = data_frame_long,
            settings = list(Delays     = Delay,
                            Values     = Value,
                            Individual = ids,
                            Group      = grp),
            plan = c("mazur", "exponential"),
            maxValue = 1,
            verbose  = TRUE) |>
  dd_analyze(modelSelection = TRUE)

data_frame_results = summary(results)

png(filename = "../man/figures/single_fits_grouped_mbauc_log10.png",
    width = 8,
    height = 6,
    res = 600,
    units = "in")

plot(results, logAxis = "x", position = "topright", which = "Log10MBAUC")

dev.off()
