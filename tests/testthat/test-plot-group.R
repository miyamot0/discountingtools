
describe("dd_plot: Various Individuals", {
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(dplyr))
  library(discountingtools)

  n_per_group <- 5

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
                          verbose  = FALSE) |>
    dd_analyze(modelSelection = TRUE)

  it("Plots individually in groups: Predictions", {
    expect_no_error(
      plot(results, logAxis = "x", position = "topright", which = "group")
    )
  })

  it("Plots individually in groups: Predictions [plotit = false]", {
    expect_no_error(
      plot(results,
           which = "group",
           plotit = FALSE)
    )
  })

  it("Plots individually in groups: Single Predictions", {
    expect_no_error(
      plot(results, logAxis = "x", position = "topright",
           which = "group", id = "1")
    )
  })

  it("Plots individually in groups: Single Predictions [plotit = false]", {
    expect_no_error(
      plot(results,
           which = "group",
           id = "1",
           plotit = FALSE)
    )
  })

  it("Plots individually in groups: models", {
    expect_no_error(
      plot(results, logAxis = "x", position = "topright", which = "model")
    )
  })

  it("Plots individually in groups: models [plotit = false]", {
    expect_no_error(
      plot(results,
           which = "model",
           plotit = FALSE)
    )
  })

  it("Plots individually: ED50", {
    expect_no_error(
      plot(results, logAxis = "x", position = "topright", which = "ED50")
    )
  })

  it("Plots individually: ED50 [plotit = false]", {
    expect_no_error(
      plot(results,
           which = "ED50",
           plotit = FALSE)
    )
  })

  it("Plots individually: MBAUC", {
    expect_no_error(
      plot(results, logAxis = "x", position = "topright", which = "MBAUC")
    )
  })

  it("Plots individually: MBAUC [plotit = false]", {
    expect_no_error(
      plot(results,
           which = "MBAUC",
           plotit = FALSE)
    )
  })

  it("Plots individually: MBAUC Log10 Scaled", {
    expect_no_error(
      plot(results, logAxis = "x", position = "topright", which = "Log10MBAUC")
    )
  })

  it("Plots individually: MBAUC Log10 Scaled [plotit = false]", {
    expect_no_error(
      plot(results,
           which = "Log10MBAUC",
           plotit = FALSE)
    )
  })
})

