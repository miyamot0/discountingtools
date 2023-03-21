rm(list = ls())

library(discountingtools)

og_k <- -0.7549554
og_ln_ed50 <- 0.8379741
og_mb_auc <- 0.002105094
og_mb_auc_log <- 0.1320793

data_frame = data.frame(
  ids = 1,
  ks  = og_k,
  delay = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)
)

data_frame[, 'value'] <- dd_discount_func_mazur(
  data_frame$delay,
  data_frame$ks) +
  c(0.0109249702,
    0.0305477299,
    -0.0118369887,
    0.0052715790,
    -0.0454706108,
    0.0068767728,
    0.0008717358,
    -0.0044386829)

## TODO: need good tests here

describe("dd_fit_curves", {
  it("Should pass", {
    testthat::expect_no_error(
      fit_dd_curves(
        data = data_frame,
        settings = list(Delays     = delay,
                        Values     = value,
                        Individual = ids),
        maxValue = 1,
        plan = c('mazur'),
        verbose = TRUE) |>
        dd_screen_options() |>
        dd_analyze(modelSelection = FALSE)
    )
  })

  it("Should pass - with model selection", {
    testthat::expect_no_error(
      fit_dd_curves(
        data = data_frame,
        settings = list(Delays     = delay,
                        Values     = value,
                        Individual = ids),
        maxValue = 1,
        plan = c('mazur'),
        verbose = TRUE) |>
        dd_screen_options() |>
        dd_analyze(modelSelection = TRUE)
    )
  })

  it("Should fail: Missing Delays", {
    testthat::expect_error(
      fit_dd_curves(
        data = data_frame,
        settings = list(Values     = value,
                        Individual = ids),
        maxValue = 1,
        plan = c('mazur')) |>
        dd_analyze(modelSelection = FALSE),
      "No Delays aesthetic specified"
    )
  })

  it("Should fail: Missing Values", {
    testthat::expect_error(
      fit_dd_curves(
        data = data_frame,
        settings = list(Delays     = delay,
                        Individual = ids),
        maxValue = 1,
        plan = c('mazur')) |>
        dd_analyze(modelSelection = FALSE),
      "No Values aesthetic specified"
    )
  })

  it("Should fail: Missing Individuals", {
    testthat::expect_error(
      fit_dd_curves(
        data = data_frame,
        settings = list(Delays     = delay,
                        Values     = value),
        maxValue = 1,
        plan = c('mazur')) |>
        dd_analyze(modelSelection = FALSE),
      "No Individual aesthetic specified"
    )
  })

  it("Should fail: No models", {
    testthat::expect_error(
      fit_dd_curves(
        data = data_frame,
        settings = list(Delays     = delay,
                        Values     = value,
                        Individual = ids),
        maxValue = 1) |>
        dd_analyze(modelSelection = FALSE),
      "No models specified"
    )
  })

  it("Should fail: No max value specified", {
    testthat::expect_error(
      fit_dd_curves(
        data = data_frame,
        settings = list(Delays     = delay,
                        Values     = value,
                        Individual = ids),
        plan = c('mazur')) |>
        dd_analyze(modelSelection = FALSE),
      "No maximum value specified"
    )
  })

  it("Should fail: bad strategy specified", {
    testthat::expect_error(
      fit_dd_curves(
        data = data_frame,
        settings = list(Delays     = delay,
                        Values     = value,
                        Individual = ids),
        plan = c('mazur'),
        maxValue = 1,
        strategy = "individualized") |>
        dd_analyze(modelSelection = FALSE),
      "Only `ind` or `group` strategies supported"
    )
  })

})

