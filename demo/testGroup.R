# # Example: Simulated group fits and parameter recovery
#
# # Note: Will take a minute to run
#
# rm(list = ls())
#
# library(tidyverse)
# library(discountingtools)
#
# set.seed(65535)
#
# n_per_group <- 10
#
# data_frame = data.frame(
#   ids = seq_len(n_per_group),
#   ks  = NA,
#   grp = "Group A"
# )
#
# data_frame$ks  = rnorm(length(data_frame$ids), 0.35, 0.125)
# data_frame$ks  = log(data_frame$ks)
#
# delays = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)
#
# data_frame$auc = dd_mbauc_mazur_exact(1, data_frame$ks, min(delays), max(delays))
#
# for (row in 1:nrow(data_frame)) {
#   ys = hyperbolicDiscountFunc(delays, data_frame[row, "ks"]) + rnorm(length(delays),
#                                                                     0,
#                                                                     0.05)
#
#   data_frame[row, as.character(delays)] = ys
# }
#
# data_frame2 = data.frame(
#   ids = 50 + seq_len(n_per_group),
#   ks  = NA,
#   grp = "Group B"
# )
#
# data_frame2$ks  = rnorm(length(data_frame2$ids), 0.225, 0.035)
# data_frame2$ks  = log(data_frame2$ks)
#
# data_frame2$auc = dd_mbauc_mazur_exact(1, data_frame2$ks, min(delays), max(delays))
#
# for (row in 1:nrow(data_frame2)) {
#   ys = hyperbolicDiscountFunc(delays, data_frame2[row, "ks"]) + rnorm(length(delays),
#                                                                     0,
#                                                                     0.025)
#
#   data_frame2[row, as.character(delays)] = ys
# }
#
# data_frame = rbind(data_frame,
#                   data_frame2)
#
# data_frame_long = data_frame %>%
#   gather(Delay, Value, -ids, -ks, -grp, -auc) %>%
#   mutate(Delay = as.numeric(Delay)) %>%
#   mutate(Value = ifelse(Value < 0, 0, Value)) %>%
#   mutate(Value = ifelse(Value > 1, 0, Value))
#
# results = fit_dd_curves(data = data_frame_long,
#             settings = list(Delays     = Delay,
#                             Values     = Value,
#                             Individual = ids,
#                             Group      = grp),
#             maxValue = 1,
#             verbose  = TRUE) |>
#   dd_model_options(plan = c("mazur", "exponential"))  |>
#   dd_metric_options(metrics = c("mbauc"))  |>
#   dd_screen_options(screen = FALSE)  |>
#   dd_analyze(modelSelection = FALSE)
#
# #summary(results)
#
# # data_frame_results <- summary(results)
# #
# # par(mfrow = c(1, 2))
# #
# # plot(results, logAxis = "x", position = "topright", which = "group")
# #
# # vecColors = rainbow(2, alpha = 1)
# #
# # plot(data_frame_results[data_frame_results$Group == 'Group A','MBAUC'],
# #      data_frame[data_frame$grp == "Group A", 'auc'],
# #      main = "Fitted vs. Simulated",
# #      ylab = "Fitted",
# #      xlab = "Simulated",
# #      col = vecColors[1],
# #      ylim = c(0, 0.01),
# #      xlim = c(0, 0.01))
# #
# # points(data_frame_results[data_frame_results$Group == 'Group B','MBAUC'],
# #       data_frame[data_frame$grp == "Group B", 'auc'],
# #       main = "Fitted vs. Simulated",
# #       ylab = "Fitted AUC",
# #       xlab = "Simulated AUC",
# #       col = vecColors[2])
