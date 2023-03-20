## ----dataFormatsLongWide, echo=FALSE, message=FALSE, cache=TRUE, eval=FALSE----
#  
#  library(discountingtools)
#  library(dplyr)
#  library(knitr)
#  library(kableExtra)
#  library(tidyr)
#  
#  set.seed(65535)
#  
#  dataFrame = data.frame(
#    ids = 1:3,
#    ks  = NA
#  )
#  
#  dataFrame$ks  = rnorm(length(dataFrame$ids), 0.35, 0.125)
#  dataFrame$ks  = log(dataFrame$ks)
#  
#  delays = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)
#  
#  for (row in 1:nrow(dataFrame)) {
#    ys = hyperbolicDiscountFunc(delays, dataFrame[row, "ks"]) + rnorm(length(delays),
#                                                                      0,
#                                                                      0.025)
#  
#    dataFrame[row, as.character(delays)] = ys
#  }
#  
#  dataFrameShow = dataFrame
#  dataFrameShow[dataFrameShow < 0] = 0
#  
#  dataFrameShow %>%
#    select(-ks) %>%
#    kable(., caption = "Wide Data Frame") %>%
#    kable_styling(full_width = TRUE)
#  
#  dataFrame.long = dataFrame %>%
#    select(-ks) %>%
#    gather(Delay, Value, -ids) %>%
#    mutate(Delay = as.numeric(Delay)) %>%
#    mutate(Value = ifelse(Value < 0, 0, Value)) %>%
#    mutate(Value = ifelse(Value > 1, 0, Value)) %>%
#    arrange(ids)
#  
#  kable(dataFrame.long, caption = "Long Data Frame") %>%
#    kable_styling(full_width = TRUE)
#  

## ----loadData, echo=TRUE, cache=TRUE, eval=FALSE------------------------------
#  
#                        # Data argument takes long data frame
#  results = fitDDCurves(data = dataFrame.long,
#                        # Settings argument takes a named list
#                        # Delays is linked the Delay column in data frame
#                        # Values is linked the Value column in data frame
#                        # Individual is linked the ids column in data frame
#                        settings = list(Delays     = Delay,
#                                        Values     = Value,
#                                        Individual = ids),
#                        # The max value of the commodity (A) is listed here
#                        maxValue = 1,
#                        # Modeling output can be enabled/disable with Verbose
#                        verbose  = TRUE)
#  

## ----analyticalStrategySingle, message=FALSE, eval=FALSE, cache=TRUE, eval=FALSE----
#  
#  results = fitDDCurves(data = dataFrame.long,
#                        settings = list(Delays     = Delay,
#                                        Values     = Value,
#                                        Individual = ids),
#                        maxValue = 1,
#                        verbose  = TRUE) %>%
#    # A single model is listed below
#    dd_modelOptions(plan = c("mazur")) %>%
#    # model selection procedures are switched off
#    dd_analyze(modelSelection = FALSE)
#  
#  summary(results)
#  

## ----analyticalStrategyMultiple, message=FALSE, eval=FALSE, cache=TRUE, eval=FALSE----
#  
#  results = fitDDCurves(data = dataFrame.long,
#                        settings = list(Delays     = Delay,
#                                        Values     = Value,
#                                        Individual = ids),
#                        maxValue = 1,
#                        verbose  = TRUE) %>%
#    # a wide range of models are specified
#    dd_modelOptions(plan = c("mazur",
#                             "bleichrodt",
#                             "ebertprelec",
#                             "exponential",
#                             "greenmyerson",
#                             "laibson",
#                             "noise",
#                             "rachlin",
#                             "rodriguezlogue")) %>%
#    # dd_analyze() proceeds with defaults (model selection is enabled)
#    dd_analyze()
#  
#  summary(results)
#  

## ----ed50Example, message=FALSE, eval=FALSE, cache=TRUE, eval=FALSE-----------
#  results = fitDDCurves(data = dataFrame.long,
#                        settings = list(Delays     = Delay,
#                                        Values     = Value,
#                                        Individual = ids),
#                        maxValue = 1,
#                        verbose  = TRUE) %>%
#    dd_modelOptions(plan = c("mazur",
#                             "bleichrodt",
#                             "ebertprelec",
#                             "exponential",
#                             "greenmyerson",
#                             "laibson",
#                             "noise",
#                             "rachlin",
#                             "rodriguezlogue")) %>%
#    # dd_metricOptions is supplied and the desired metric is specified
#    dd_metricOptions(metrics = c("lned50")) %>%
#    dd_analyze()
#  
#  summary(results)

## ----mbaucExample, message=FALSE, eval=FALSE, cache=TRUE, eval=FALSE----------
#  results = fitDDCurves(data = dataFrame.long,
#                        settings = list(Delays     = Delay,
#                                        Values     = Value,
#                                        Individual = ids),
#                        maxValue = 1,
#                        verbose  = TRUE) %>%
#    dd_modelOptions(plan = c("mazur",
#                             "bleichrodt",
#                             "ebertprelec",
#                             "exponential",
#                             "greenmyerson",
#                             "laibson",
#                             "noise",
#                             "rachlin",
#                             "rodriguezlogue")) %>%
#    # dd_metricOptions is supplied and the desired metric is specified
#    dd_metricOptions(metrics = c("mbauc")) %>%
#    dd_analyze()
#  
#  summary(results)

## ----log10MbaucExample, message=FALSE, eval=FALSE, cache=TRUE, eval=FALSE-----
#  results = fitDDCurves(data = dataFrame.long,
#                        settings = list(Delays     = Delay,
#                                        Values     = Value,
#                                        Individual = ids),
#                        maxValue = 1,
#                        verbose  = TRUE) %>%
#    dd_modelOptions(plan = c("mazur",
#                             "bleichrodt",
#                             "ebertprelec",
#                             "exponential",
#                             "greenmyerson",
#                             "laibson",
#                             "noise",
#                             "rachlin",
#                             "rodriguezlogue")) %>%
#    # dd_metricOptions is supplied and the desired metric is specified
#    dd_metricOptions(metrics = c("logmbauc")) %>%
#    dd_analyze()
#  
#  summary(results)

## ----screeningMethods, message=FALSE, eval=FALSE, cache=TRUE, eval=FALSE------
#  
#  results = fitDDCurves(data = dataFrame.long,
#              settings = list(Delays     = Delay,
#                              Values     = Value,
#                              Individual = ids),
#              maxValue = 1,
#              verbose  = TRUE) %>%
#    dd_modelOptions(plan = c("mazur",
#                             "bleichrodt",
#                             "ebertprelec",
#                             "exponential",
#                             "greenmyerson",
#                             "laibson",
#                             "noise",
#                             "rachlin",
#                             "rodriguezlogue")) %>%
#    dd_metricOptions(metrics = c("lned50",
#                                 "mbauc",
#                                 "logmbauc")) %>%
#    # flag to enable screening is set to TRUE in dd_screenOption
#    dd_screenOption(screen   = TRUE) %>%
#    dd_analyze()
#  

## ----screeningAndFilteringMethods, message=FALSE, eval=FALSE, cache=TRUE, eval=FALSE----
#  
#  results = fitDDCurves(data = dataFrame.long,
#              settings = list(Delays     = Delay,
#                              Values     = Value,
#                              Individual = ids),
#              maxValue = 1,
#              verbose  = TRUE) %>%
#    dd_modelOptions(plan = c("mazur",
#                             "bleichrodt",
#                             "ebertprelec",
#                             "exponential",
#                             "greenmyerson",
#                             "laibson",
#                             "noise",
#                             "rachlin",
#                             "rodriguezlogue")) %>%
#    dd_metricOptions(metrics = c("lned50",
#                                 "mbauc",
#                                 "logmbauc")) %>%
#    # filterPassing is set to require JB1 and JB2 to pass to include in analysis
#    dd_screenOption(screen        = TRUE,
#                    filterPassing = c("JB1", "JB2")) %>%
#    dd_analyze()
#  

## ----visualizationData, echo=FALSE, message=FALSE, warning=FALSE, cache=TRUE, eval=FALSE----
#  
#  library(dplyr)
#  library(discountingtools)
#  library(tidyr)
#  
#  set.seed(65535)
#  
#  dataFrame = data.frame(
#    ids = 1:50,
#    ks  = NA,
#    grp = "Group A"
#  )
#  
#  dataFrame$ks  = rnorm(length(dataFrame$ids), 0.35, 0.125)
#  dataFrame$ks  = log(dataFrame$ks)
#  
#  delays = c(1, 30, 180, 540, 1080, 2160, 4320, 8640)
#  
#  for (row in 1:nrow(dataFrame)) {
#    ys = hyperbolicDiscountFunc(delays, dataFrame[row, "ks"]) + rnorm(length(delays),
#                                                                      0,
#                                                                      0.05)
#  
#    dataFrame[row, as.character(delays)] = ys
#  }
#  
#  dataFrame2 = data.frame(
#    ids = 51:100,
#    ks  = NA,
#    grp = "Group B"
#  )
#  
#  dataFrame2$ks  = rnorm(length(dataFrame2$ids), 0.225, 0.035)
#  dataFrame2$ks  = log(dataFrame2$ks)
#  
#  for (row in 1:nrow(dataFrame2)) {
#    ys = hyperbolicDiscountFunc(delays, dataFrame2[row, "ks"]) + rnorm(length(delays),
#                                                                      0,
#                                                                      0.05)
#  
#    dataFrame2[row, as.character(delays)] = ys
#  }
#  
#  dataFrame = rbind(dataFrame,
#                    dataFrame2)
#  
#  dataFrame.long = dataFrame %>%
#    gather(Delay, Value, -ids, -ks, -grp) %>%
#    mutate(Delay = as.numeric(Delay))
#  
#  dataFrame.long[,"Value"] = ifelse(dataFrame.long[,"Value"] > 1, 1, dataFrame.long[,"Value"])
#  dataFrame.long[,"Value"] = ifelse(dataFrame.long[,"Value"] < 0, 0, dataFrame.long[,"Value"])
#  
#  results = fitDDCurves(data = dataFrame.long,
#              settings = list(Delays     = Delay,
#                              Values     = Value,
#                              Individual = ids,
#                              Group      = grp),
#              maxValue = 1,
#              verbose  = TRUE) %>%
#    dd_modelOptions(plan = c("mazur",
#                             "bleichrodt",
#                             "ebertprelec",
#                             "exponential",
#                             "greenmyerson",
#                             "laibson",
#                             "noise",
#                             "rachlin",
#                             "rodriguezlogue")) %>%
#    dd_metricOptions(metrics = c("lned50",
#                                "mbauc",
#                                "logmbauc")) %>%
#    dd_analyze()
#  

## ----visInd, fig.width=6, fig.height=6, warning=FALSE, message=FALSE, eval=FALSE----
#  
#  plot(results, which    = "ind",
#                logAxis  = "x",
#                position = "topright")
#  

## ----visGrp, fig.width=6, fig.height=6, warning=FALSE, message=FALSE, eval=FALSE----
#  
#  plot(results, which    = "group",
#                logAxis  = "x",
#                position = "topright")
#  

## ----visGrpED50, fig.width=6, fig.height=6, warning=FALSE, message=FALSE, eval=FALSE----
#  
#  plot(results, which    = "ED50",
#                logAxis  = "x",
#                position = "topright")
#  

## ----visGrpMBAUC, fig.width=6, fig.height=6, warning=FALSE, message=FALSE, eval=FALSE----
#  
#  plot(results, which    = "MBAUC",
#                logAxis  = "x",
#                position = "topright")
#  

## ----visGrpLog10MBAUC, fig.width=6, fig.height=6, warning=FALSE, message=FALSE, eval=FALSE----
#  
#  plot(results, which    = "Log10MBAUC",
#                logAxis  = "x",
#                position = "topright")
#  

