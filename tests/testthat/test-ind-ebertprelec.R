rm(list = ls())

library(discountingtools)

og_k <- -2
og_s <- 0.9
og_ln_ed50 <- 1.77213
og_mb_auc <- 0.001022831
og_mb_auc_log <- 0.1871914

data_frame = data.frame(
  ids = 1,
  ks  = og_k,
  ss  = og_s,
  delay = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)
)

data_frame[, 'value'] <- dd_discount_func_ebertprelec(
  data_frame$delay,
  data_frame$ks,
  data_frame$ss) +
  c(0.0109249702,
    0.0305477299,
    -0.0118369887,
    0.0052715790,
    -0.0454706108,
    0.0068767728,
    0.0008717358,
    -0.0044386829)

describe("dd_fit: Ebert and Prelec Model", {

  cached_results = fit_dd_curves(
    data = data_frame,
    settings = list(Delays     = delay,
                    Values     = value,
                    Individual = ids),
    maxValue = 1,
    plan = c('ebertprelec')) |>
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
        plan = c('ebertprelec')) |>
        dd_analyze(modelSelection = FALSE)
    )
  })

  it("Should be close to simulated parameter (15%)", {
    testthat::expect_equal(
      cached_results[1, 'EbertPrelec.Lnk'],
      og_k,
      tolerance = 0.15
    )
  })

  it("Should be close to simulated parameter (15%)", {
    testthat::expect_equal(
      cached_results[1, 'EbertPrelec.S'],
      og_s,
      tolerance = 0.15
    )
  })

  it("Should be close to expected LnED50", {
    testthat::expect_equal(
      cached_results[1, 'EbertPrelec.LnED50'],
      og_ln_ed50,
      tolerance = 0.1
    )
  })

  it("Should be close to expected MBAUC", {
    testthat::expect_equal(
      cached_results[1, 'EbertPrelec.MBAUC'],
      og_mb_auc,
      tolerance = 0.1
    )
  })

  it("Should be close to expected Log10 MBAUC", {
    testthat::expect_equal(
      cached_results[1, 'EbertPrelec.Log10MBAUC'],
      og_mb_auc_log,
      tolerance = 0.1
    )
  })
})

