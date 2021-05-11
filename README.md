# discountingtools
-------------------
R package for performing delay discounting analyses.  Features include approximate Bayesian model selection, derivation of the Effective Delay 50 (ED50) and both normal and log10 transformation model area using Numerical Integration

### Version
-------------------
0.0.2.0 (beta)

### Installation
----------------
1) Install and load the devtools package. 
2) Install and load the digest package. 
3) Run install_github to install the package.

```r
install.packages("devtools")
install.packages("digest")

devtools::install_github("miyamot0/discountingtools")

library(discountingtools)
```

### Model Candidates (Possible models within the model selection approach)
-------------------
* Noise Model: Intecept-only comparison model (included by default)

* Exponential Model: Samuelson, P. A. (1937). A note on measurement of utility. The Review of Economic Studies, 4(2), 155–161. https://doi.org/10.2307/2967612

* Hyperbolic Model: Mazur, J. E. (1987). An adjusting procedure for studying delayed reinforcement. In Quantitative analysis of behavior: Vol. 5. The effect of delay and intervening events on reinforcement value (pp. 55–73). Hillsdale, NJ: Erlbaum.

* Rodriguez & Logue Model: Rodriguez, M. L., & Logue, A. W. (1988). Adjusting delay to reinforcement: Comparing choice in pigeons and humans. Journal of Experimental Psychology: Animal Behavior Processes, 14(1), 105–117. https://doi.org/10.1037/0097-7403.14.1.105

* Beta Delta Model: Laibson, D. (1997). Golden eggs and hyperbolic discounting. The Quarterly Journal of Economics, 112(2), 443–478. https://doi.org/10.1162/003355397555253

* Green & Myerson Model: Green, L., & Myerson, J. (2004). A discounting framework for choice with delayed and probabilistic rewards. Psychological Bulletin, 130(5), 769–792. https://doi.org/10.1037/0033-2909.130.5.769

* Rachlin Model: Rachlin, H. (2006). Notes on discounting. Journal of the Experimental Analysis of Behavior, 85(3), 425–435. https://doi.org/10.1901/jeab.2006.85-05

* Ebert & Prelec Model: Ebert, J. E. J., & Prelec, D. (2007). The Fragility of Time: Time-Insensitivity and Valuation of the Near and Far Future. Management Science, 53(9), 1423–1438. https://doi.org/10.1287/mnsc.1060.0671

* Bleichrodt et al. Constant Relative Decreasing Impatience: Bleichrodt, H., Rohde, K. I. M., & Wakker, P. P. (2009). Non-hyperbolic time inconsistency. Games and Economic Behavior, 66(1), 27–38. https://doi.org/10.1016/j.geb.2008.05.007

### Examples
-------------------------

#### Effective Delay Plotting

```r

# TODO

```

![Alt text](...)

#### Model-based Area

```r

# TODO

```

![Alt text](...)

#### Model-based Area (logarithmic scaling)

```r

# TODO

```

![Alt text](...)

### Referenced Works (academic works)
-------------------
The Small N Stats Discounting Model Selector is based on the following academic works:
* Borges, A. M., Kuang, J., Milhorn, H., & Yi, R. (2016). An alternative approach to calculating Area-Under-the-Curve (AUC) in delay discounting research. Journal of the Experimental Analysis of Behavior, 106(2), 145–155. https://doi.org/10.1002/jeab.219

* Franck, C. T., Koffarnus, M. N., House, L. L., & Bickel, W. K. (2015). Accurate characterization of delay discounting: a multiple model approach using approximate Bayesian model selection and a unified discounting measure. Journal of the Experimental Analysis of Behavior, 103(1), 218–233. https://doi.org/10.1002/jeab.128

* Gilroy, S. P., Franck, C. T., & Hantula, D. A. (2017). The discounting model selector: Statistical software for delay discounting applications. Journal of the Experimental Analysis of Behavior, 107(3), 388–401. https://doi.org/10.1002/jeab.257

* Myerson, J., Green, L., & Warusawitharana, M. (2001). Area under the curve as a measure of discounting. Journal of the Experimental Analysis of Behavior, 76(2), 235–243. https://doi.org/10.1901/jeab.2001.76-235

* Yoon, J. H., & Higgins, S. T. (2008). Turning k on its head: Comments on use of an ED50 in delay discounting research. Drug and Alcohol Dependence, 95(1–2), 169–172. https://doi.org/10.1016/j.drugalcdep.2007.12.011

### Acknowledgments
-------------------
* Donald A. Hantula, Decision Making Laboratory, Temple University [Site](http://astro.temple.edu/~hantula/)

* Chris Franck, Laboratory for Interdisciplinary Statistical Analysis - Virginia Tech

### Questions, Suggestions, and Contributions
---------------------------------------------

Questions? Suggestions for features? <sgilroy1@lsu.edu>.

### License
-----------

GPL Version >= 2
