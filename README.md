# discountingtools

The *discountingtools* R package was designed to support researchers in conducting delay discounting analyses. This approach and methodology features considerable variability and this package has been designed for use as a modular toolkit for approaching questions related to delay discounting. For example, this package can be used to evaluate individual, specific models (e.g., hyperbolic model), explore and compare competing models (e.g., cross-model questions), or simply as a tool to characterize data (e.g., systematic vs. non-systematic responding). Specific features and usage scenarios are illustrated in the sections below.

## Version

0.0.2.1 (beta)

## Installation (pre-CRAN state)

1\) Install and load the devtools package. 2) Install and load the digest package. 3) Run install_github to install the package.

``` r
install.packages("devtools")
install.packages("digest")

devtools::install_github("miyamot0/discountingtools")

library(discountingtools)
```

## Model Candidates

-   Noise Model: Intercept-only comparison model (included by default)

-   Exponential Model: Samuelson, P. A. (1937). A note on measurement of utility. The Review of Economic Studies, 4(2), 155–161. <https://doi.org/10.2307/2967612>

-   Hyperbolic Model: Mazur, J. E. (1987). An adjusting procedure for studying delayed reinforcement. In Quantitative analysis of behavior: Vol. 5. The effect of delay and intervening events on reinforcement value (pp. 55–73). Hillsdale, NJ: Erlbaum.

-   Rodriguez & Logue Model: Rodriguez, M. L., & Logue, A. W. (1988). Adjusting delay to reinforcement: Comparing choice in pigeons and humans. Journal of Experimental Psychology: Animal Behavior Processes, 14(1), 105–117. <https://doi.org/10.1037/0097-7403.14.1.105>

-   Beta Delta Model: Laibson, D. (1997). Golden eggs and hyperbolic discounting. The Quarterly Journal of Economics, 112(2), 443–478. <https://doi.org/10.1162/003355397555253>

-   Green & Myerson Model: Green, L., & Myerson, J. (2004). A discounting framework for choice with delayed and probabilistic rewards. Psychological Bulletin, 130(5), 769–792. <https://doi.org/10.1037/0033-2909.130.5.769>

-   Rachlin Model: Rachlin, H. (2006). Notes on discounting. Journal of the Experimental Analysis of Behavior, 85(3), 425–435. <https://doi.org/10.1901/jeab.2006.85-05>

-   Ebert & Prelec Model: Ebert, J. E. J., & Prelec, D. (2007). The Fragility of Time: Time-Insensitivity and Valuation of the Near and Far Future. Management Science, 53(9), 1423–1438. <https://doi.org/10.1287/mnsc.1060.0671>

-   Bleichrodt et al. Constant Relative Decreasing Impatience: Bleichrodt, H., Rohde, K. I. M., & Wakker, P. P. (2009). Non-hyperbolic time inconsistency. Games and Economic Behavior, 66(1), 27–38. <https://doi.org/10.1016/j.geb.2008.05.007>

## Usage Scenarios

### Single-model Evaluation

Various users may wish to explore delay discounting phenomena with an a priori assumption regarding the underlying data-generating process. The *discountingtools* package can be used to load a data frame (i.e., data = ...), map settings to the data frame (i.e., Delays, Values, and Individuals have references mapped via settings = ...). Modeling options are specified using *dd_modelOptions* and *dd_analyze* (with model selection disabled) can be used to evaluate a single model.

A short snippet is illustrated below and a complete example of this approach is illustrated in demo/testSingle.R.

``` r
results = fitDDCurves(data = dataFrame.long,
                      settings = list(Delays     = Delay,
                                      Values     = Value,
                                      Individual = ids),
                      maxValue = 1) %>%
  dd_modelOptions(plan = c("mazur")) %>%
  dd_analyze(modelSelection = FALSE)

summary(results)
```

![Figure of Single Model Evaluation](figures/SingleModelEvaluation.png "Single Model Evaluation")

### Multi-Model Evaluation

More recent discussions on delay discounting patterns and processes have questioned whether *any* a priori assumptions regarding discounting models are tenable. As such, it is now more common for investigators to explore competing models before conducting terminal analyses. The *discountingtools* package can be used to load delay discounting data and specify a range of modeling options (e.g., hyperbolic, exponential). Specifically, a range of modeling options can be specified using *dd_modelOptions* and *dd_analyze* (with model selection *enabled*) and approximate Bayesian model selection will be perform to identify the best-performing candidate (at the individual level).

A short snippet is illustrated below and a complete example of this approach is illustrated in demo/testIndividual.R.

``` r
results = fitDDCurves(data = dataFrame.long,
            settings = list(Delays     = Delay,
                            Values     = Value,
                            Individual = ids),
            maxValue = 1) %>%
  dd_modelOptions(plan = c("mazur",
                           "bleichrodt",
                           "ebertprelec",
                           "exponential",
                           "greenmyerson",
                           "laibson",
                           "noise",
                           "rachlin",
                           "rodriguezlogue")) %>%
  dd_metricOptions(metrics = c("lned50",
                               "mbauc",
                               "logmbauc")) %>%
  dd_analyze()
  
summary(results)
```

![Figure of Multi Model Evaluation](figures/MultiModelEvaluation.png "Multi Model Evaluation")

#### Effective Delay 50 (ED50)

The multi-model approach is frustrated by the presence of distinct parameters. As an alternative, researchers have suggested a metric based on the rate of decay (i.e., time until decay to 50%). The multi-model evaluation provided above is re-evaluated in terms of ED50 below.

The full code necessary to re-create this result is provided in demo/testIndividualED50.R.

``` R
results = fitDDCurves(data = dataFrame.long,
            settings = list(Delays     = Delay,
                            Values     = Value,
                            Individual = ids),
            maxValue = 1) %>%
  dd_modelOptions(plan = c("mazur",
                           "bleichrodt",
                           "ebertprelec",
                           "exponential",
                           "greenmyerson",
                           "laibson",
                           "noise",
                           "rachlin",
                           "rodriguezlogue")) %>%
  dd_metricOptions(metrics = c("lned50")) %>%
  dd_analyze()

plot(results, which = "ED50")
```

![Figure of Multi Model ED50](figures/MultiModelEvaluationED50.png "Multi Model ED50")

#### Numerical Integration Area (MB-AUC)

As an alternative to ED50, researchers have suggested a metric based on the model area. The multi-model evaluation provided above is re-evaluated in terms of MB-AUC below.

The full code necessary to re-create this result is provided in demo/testIndividualMBAUC.R.

``` R
results = fitDDCurves(data = dataFrame.long,
            settings = list(Delays     = Delay,
                            Values     = Value,
                            Individual = ids),
            maxValue = 1) %>%
  dd_modelOptions(plan = c("mazur",
                           "bleichrodt",
                           "ebertprelec",
                           "exponential",
                           "greenmyerson",
                           "laibson",
                           "noise",
                           "rachlin",
                           "rodriguezlogue")) %>%
  dd_metricOptions(metrics = c("mbauc")) %>%
  dd_analyze()

plot(results, which = "MBAUC")
```

![Figure of Model-based AUC](figures/MultiModelEvaluationMBAUC.png "Model-based AUC")

#### Log10-Scaled Numerical Integration Area (Log10 MB-AUC)

Researchers using area is a metric of discounting have suggested re-scaling delays. The multi-model evaluation provided above is re-evaluated in terms of log10-scaled MB-AUC below.

The full code necessary to re-create this result is provided in demo/testIndividualLog10MBAUC.R.

``` R
results = fitDDCurves(data = dataFrame.long,
            settings = list(Delays     = Delay,
                            Values     = Value,
                            Individual = ids),
            maxValue = 1) %>%
  dd_modelOptions(plan = c("mazur",
                           "bleichrodt",
                           "ebertprelec",
                           "exponential",
                           "greenmyerson",
                           "laibson",
                           "noise",
                           "rachlin",
                           "rodriguezlogue")) %>%
  dd_metricOptions(metrics = c("logmbauc")) %>%
  dd_analyze()

plot(results, which = "Log10MBAUC")
```

![Figure of Log10 Scaled MBAUC](figures/MultiModelEvaluationLog10MBAUC.png "Log10 Scaled MBAUC")

### Multi-Model Evaluation (Grouped)

...

A short snippet is illustrated below and a complete example of this approach is illustrated in demo/testGroup.R.

![Figure of Multi Model (Group) Method](figures/MultiModelEvaluationGroup.png "Multi Model (Group) Method")

#### ED50

...

The full code necessary to re-create this result is provided in demo/testGroupED50.R.

![Figure of Multi Model ED50 (Group)](figures/MultiModelEvaluationGroupED50.png "Multi Model ED50 (Group)")

#### MB-AUC

...

The full code necessary to re-create this result is provided in demo/testGroupMBAUC.R.

![Figure of Multi Model MBAUC (Group)](figures/MultiModelEvaluationGroupMBAUC.png "Multi Model MBAUC (Group)")

#### Log10 MB-AUC

...

The full code necessary to re-create this result is provided in demo/testGroupLog10MBAUC.R.

![Figure of Multi-Model Log10 MBAUC (Group)](figures/MultiModelEvaluationGroupLog10MBAUC.png "Multi-Model Log10 MBAUC (Group)")

## Referenced Works (academic works)

The Small N Stats Discounting Model Selector is based on the following academic works:

-   Borges, A. M., Kuang, J., Milhorn, H., & Yi, R. (2016). An alternative approach to calculating Area-Under-the-Curve (AUC) in delay discounting research. Journal of the Experimental Analysis of Behavior, 106(2), 145–155. <https://doi.org/10.1002/jeab.219>

-   Franck, C. T., Koffarnus, M. N., House, L. L., & Bickel, W. K. (2015). Accurate characterization of delay discounting: a multiple model approach using approximate Bayesian model selection and a unified discounting measure. Journal of the Experimental Analysis of Behavior, 103(1), 218–233. <https://doi.org/10.1002/jeab.128>

-   Gilroy, S. P., Franck, C. T., & Hantula, D. A. (2017). The discounting model selector: Statistical software for delay discounting applications. Journal of the Experimental Analysis of Behavior, 107(3), 388–401. <https://doi.org/10.1002/jeab.257>

-   Myerson, J., Green, L., & Warusawitharana, M. (2001). Area under the curve as a measure of discounting. Journal of the Experimental Analysis of Behavior, 76(2), 235–243. <https://doi.org/10.1901/jeab.2001.76-235>

-   Yoon, J. H., & Higgins, S. T. (2008). Turning k on its head: Comments on use of an ED50 in delay discounting research. Drug and Alcohol Dependence, 95(1–2), 169–172. <https://doi.org/10.1016/j.drugalcdep.2007.12.011>

## Acknowledgments

-   Donald A. Hantula, Decision Making Laboratory, Temple University [Site](http://astro.temple.edu/~hantula/)

-   Chris Franck, Laboratory for Interdisciplinary Statistical Analysis - Virginia Tech

## Questions, Suggestions, and Contributions

Questions? Suggestions for features? [sgilroy1\@lsu.edu](mailto:sgilroy1@lsu.edu).

## To Do

-   ~~Modeling options: individual (default)~~
-   Modeling options: pooled (fit to group)
-   ~~Plotting tools (corresponding with modeling options)~~

## License

GPL Version \>= 2
