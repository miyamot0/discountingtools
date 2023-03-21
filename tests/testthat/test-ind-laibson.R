rm(list = ls())

library(discountingtools)

og_d <- 0.7373318
og_b <- 1
og_ln_ed50 <- 0.8711429
og_mb_auc <- 0.0003574047
og_mb_auc_log <- 0.1026852

data_frame = data.frame(
  ids = 1,
  ds  = og_d,
  bs  = og_b,
  delay = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)
)

data_frame[, 'value'] <- dd_discount_func_laibson(
  data_frame$delay,
  data_frame$bs,
  data_frame$ds) +
  c(0.0109249702,
    0.0305477299,
    -0.0118369887,
    0.0052715790,
    -0.0454706108,
    0.0068767728,
    0.0008717358,
    -0.0044386829)

describe("dd_fit: Laibson Model", {

  cached_results = fit_dd_curves(
    data = data_frame,
    settings = list(Delays     = delay,
                    Values     = value,
                    Individual = ids),
    maxValue = 1,
    plan = c('laibson')) |>
    dd_analyze(modelSelection = FALSE) |>
    summary()

  it("Should not fail with simple data", {
    expect_no_error(
      fit_dd_curves(
        data = data_frame,
        settings = list(Delays     = delay,
                        Values     = value,
                        Individual = ids),
        maxValue = 1,
        plan = c('laibson')) |>
        dd_analyze(modelSelection = FALSE)
    )
  })

  it("Should be close to simulated parameter (15%)", {
    testthat::expect_equal(
      cached_results[1, 'Laibson.Delta'],
      og_d,
      tolerance = 0.15
    )
  })

  it("Should be close to simulated parameter (15%)", {
    testthat::expect_equal(
      cached_results[1, 'Laibson.Beta'],
      og_b,
      tolerance = 0.15
    )
  })

  it("Should be close to expected LnED50", {
    testthat::expect_equal(
      cached_results[1, 'Laibson.LnED50'],
      og_ln_ed50,
      tolerance = 0.1
    )
  })

  it("Should be close to expected MBAUC", {
    testthat::expect_equal(
      cached_results[1, 'Laibson.MBAUC'],
      og_mb_auc,
      tolerance = 0.05
    )
  })

  it("Should be close to expected Log10 MBAUC", {
    testthat::expect_equal(
      cached_results[1, 'Laibson.Log10MBAUC'],
      og_mb_auc_log,
      tolerance = 0.05
    )
  })
})

