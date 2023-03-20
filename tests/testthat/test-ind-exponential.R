rm(list = ls())

library(discountingtools)

og_k <- -0.7925272
og_ln_ed50 <- 0.4260143
og_mb_auc <- 0.0001626005
og_mb_auc_log <- 0.06856932

data_frame = data.frame(
  ids = 1,
  ks  = og_k,
  delay = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)
)

data_frame[, 'value'] <- dd_discount_func_exponential(
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

describe("dd_fit: Exponential Model", {

  cached_results = fit_dd_curves(
    data = data_frame,
    settings = list(Delays     = delay,
                    Values     = value,
                    Individual = ids),
    maxValue = 1,
    plan = c('exponential')) |>
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
        plan = c('exponential')) |>
        dd_analyze(modelSelection = FALSE)
    )
  })

  it("Should be close to simulated parameter (15%)", {
    testthat::expect_equal(
      cached_results[1, 'Exponential.Lnk'],
      og_k,
      tolerance = 0.15
    )
  })

  it("Should be close to expected LnED50", {
    testthat::expect_equal(
      cached_results[1, 'Exponential.LnED50'],
      og_ln_ed50,
      tolerance = 0.1
    )
  })

  it("Should be close to expected MBAUC", {
    testthat::expect_equal(
      cached_results[1, 'Exponential.MBAUC'],
      og_mb_auc,
      tolerance = 0.05
    )
  })

  it("Should be close to expected Log10 MBAUC", {
    testthat::expect_equal(
      cached_results[1, 'Exponential.Log10MBAUC'],
      og_mb_auc_log,
      tolerance = 0.05
    )
  })
})

