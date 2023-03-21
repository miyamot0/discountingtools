
describe("dd_plot: Various Individuals", {
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(dplyr))
  library(discountingtools)

  dataFrame = data.frame(
    ids = 1:10,
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

  data_frame = dataFrame %>%
    tidyr::gather(Delay, Value, -ids, -ks) %>%
    dplyr::mutate(Delay = as.numeric(Delay)) %>%
    dplyr::mutate(Value = ifelse(Value < 0, 0, Value)) %>%
    dplyr::mutate(Value = ifelse(Value > 1, 0, Value))

  results = fit_dd_curves(
    data = data_frame,
    settings = list(Delays     = Delay,
                    Values     = Value,
                    Individual = ids),
    maxValue = 1,
    plan = c('mazur')) |>
    dd_analyze(modelSelection = TRUE)

  it("Plots individually: Predictions", {
    expect_no_error(
      plot(results,
           logAxis = "x",
           position = "topright")
    )
  })

  it("Plots individually: ED50", {
    expect_no_error(
      plot(results, which = "ED50")
    )
  })

  it("Plots individually: MBAUC", {
    expect_no_error(
      plot(results, which = "MBAUC")
    )
  })

  it("Plots individually: MBAUC Log10 Scaled", {
    expect_no_error(
      plot(results, which = "Log10MBAUC")
    )
  })
})

