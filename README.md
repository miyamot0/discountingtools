# discountingtools
-------------------
R package for performing delay discounting analyses.  Features include approximate Bayesian model selection, derivation of the Effective Delay 50 (ED50) and both normal and log10 transformation model area using Numerical Integration

### Version
-------------------
0.0.0.3 (alpha)

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

### Examples
-------------------------

#### Effective Delay Plotting

```r
dat <- data.frame(X=  c(1, 30,  180, 540, 1080, 2160),
                  Y=  c(1, 0.9, 0.8, 0.7, 0.6,  0.4),
                  ids= c(1, 1,   1,   1,   1,   1))

dat2<- data.frame(X=  c(1, 30,  180, 540, 1080, 2160),
                  Y=  c(1, 0.9, 0.7, 0.6, 0.5,  0.4),
                  ids= c(2, 2,   2,   2,   2,   2))

dat <- rbind(dat, dat2)

dat$Y <- dat$Y * 100

library(discountingtools)

results <- discountingModelSelection(dat,
                                     summarize = TRUE,
                                     figures = "ed50",
                                     idCol = "ids",
                                     A = 100)

Calculating series id: 1
Model Comparison Summary:
Rank #:  1  =  ep
 Probability  =  0.785999
 BIC  =  -22.2279
 Most Probable ln(ED50)  =  7.334051
 Most Probable Model AUC  =  0.597619
 Most Probable Model AUC (log)  =  0.843852

Rank #:  2  =  Rachlin
 Probability  =  0.186338
 BIC  =  -19.349116

Rank #:  3  =  Mazur
 Probability  =  0.023587
 BIC  =  -15.215414

Rank #:  4  =  exp
 Probability  =  0.003802
 BIC  =  -11.565193

Rank #:  5  =  BD
 Probability  =  0.000267
 BIC  =  -6.252941

Rank #:  6  =  noise
 Probability  =  7e-06
 BIC  =  1.128501
 
print(results)
                                     
  id JB.C1 JB.C2 noise.mean noise.RMSE noise.BIC noise.AIC   exp.lnk   exp.RMSE    exp.BIC    exp.AIC Mazur.lnk Mazur.RMSE
1  1  TRUE  TRUE  0.7333333  0.2160247  1.128501  1.544982 -7.615692 0.07500726 -11.565193 -11.148712 -7.217119 0.05533467
2  2  TRUE  TRUE  0.6833333  0.2316607  1.967073  2.383554 -7.369143 0.12400210  -5.532629  -5.116148 -6.846340 0.08482117
  Mazur.BIC  Mazur.AIC   BD.beta  BD.delta   BD.RMSE    BD.BIC    BD.AIC Rachlin.lnk Rachlin.s Rachlin.RMSE Rachlin.BIC
1 -15.21541 -14.798933 0.9999994 0.9993449 0.1124525 -6.252941 -5.628220   -5.000000 0.6800000   0.03775746   -19.34912
2 -10.08967  -9.673188 0.9999967 0.9991767 0.1542008 -2.464244 -1.839522   -3.941156 0.5676217   0.02057374   -26.63513
  Rachlin.AIC    ep.lnk      ep.s    ep.RMSE    ep.BIC    ep.AIC noise.BF    exp.BF  Mazur.BF     BD.BF Rachlin.BF    ep.BF
1   -18.72439 -7.976814 0.5702152 0.02970408 -22.22790 -21.60318        1 570.69054 3540.2682 40.073737   27967.79 117971.8
2   -26.01040 -7.796464 0.4521181 0.02740064 -23.19652 -22.57180        1  42.51475  415.0383  9.167443 1625131.21 291208.5
    noise.prob     exp.prob   Mazur.prob      BD.prob Rachlin.prob   ep.prob probable.model probable.ED50 probable.AUC
1 6.662599e-06 3.802282e-03 0.0235873880 2.669952e-04    0.1863382 0.7859985             ep      7.334051    0.5976186
2 5.217008e-07 2.217998e-05 0.0002165258 4.782663e-06    0.8478323 0.1519237        Rachlin      6.943280    0.5302546
  probable.Log10AUC
1         0.8438519
2         0.7925787
```

![Alt text](Figure_ED50_1.png?raw=true "ED50 Visuals")

![Alt text](Figure_ED50_2.png?raw=true "ED50 Visuals")

#### Model-based Area

```r
dat <- data.frame(X=  c(1, 30,  180, 540, 1080, 2160),
                  Y=  c(1, 0.9, 0.8, 0.7, 0.6,  0.4),
                  ids= c(1, 1,   1,   1,   1,   1))

dat2<- data.frame(X=  c(1, 30,  180, 540, 1080, 2160),
                  Y=  c(1, 0.9, 0.7, 0.6, 0.5,  0.4),
                  ids= c(2, 2,   2,   2,   2,   2))

dat <- rbind(dat, dat2)

dat$Y <- dat$Y * 100

library(discountingtools)

results <- discountingModelSelection(dat,
                                     summarize = TRUE,
                                     figures = "ed50",
                                     idCol = "ids",
                                     A = 100)

Calculating series id: 1
Model Comparison Summary:
Rank #:  1  =  ep
 Probability  =  0.785999
 BIC  =  -22.2279
 Most Probable ln(ED50)  =  7.334051
 Most Probable Model AUC  =  0.597619
 Most Probable Model AUC (log)  =  0.843852

Rank #:  2  =  Rachlin
 Probability  =  0.186338
 BIC  =  -19.349116

Rank #:  3  =  Mazur
 Probability  =  0.023587
 BIC  =  -15.215414

Rank #:  4  =  exp
 Probability  =  0.003802
 BIC  =  -11.565193

Rank #:  5  =  BD
 Probability  =  0.000267
 BIC  =  -6.252941

Rank #:  6  =  noise
 Probability  =  7e-06
 BIC  =  1.128501
 
print(results)
                                     
  id JB.C1 JB.C2 noise.mean noise.RMSE noise.BIC noise.AIC   exp.lnk   exp.RMSE    exp.BIC    exp.AIC Mazur.lnk Mazur.RMSE
1  1  TRUE  TRUE  0.7333333  0.2160247  1.128501  1.544982 -7.615692 0.07500726 -11.565193 -11.148712 -7.217119 0.05533467
2  2  TRUE  TRUE  0.6833333  0.2316607  1.967073  2.383554 -7.369143 0.12400210  -5.532629  -5.116148 -6.846340 0.08482117
  Mazur.BIC  Mazur.AIC   BD.beta  BD.delta   BD.RMSE    BD.BIC    BD.AIC Rachlin.lnk Rachlin.s Rachlin.RMSE Rachlin.BIC
1 -15.21541 -14.798933 0.9999994 0.9993449 0.1124525 -6.252941 -5.628220   -5.000000 0.6800000   0.03775746   -19.34912
2 -10.08967  -9.673188 0.9999967 0.9991767 0.1542008 -2.464244 -1.839522   -3.941156 0.5676217   0.02057374   -26.63513
  Rachlin.AIC    ep.lnk      ep.s    ep.RMSE    ep.BIC    ep.AIC noise.BF    exp.BF  Mazur.BF     BD.BF Rachlin.BF    ep.BF
1   -18.72439 -7.976814 0.5702152 0.02970408 -22.22790 -21.60318        1 570.69054 3540.2682 40.073737   27967.79 117971.8
2   -26.01040 -7.796464 0.4521181 0.02740064 -23.19652 -22.57180        1  42.51475  415.0383  9.167443 1625131.21 291208.5
    noise.prob     exp.prob   Mazur.prob      BD.prob Rachlin.prob   ep.prob probable.model probable.ED50 probable.AUC
1 6.662599e-06 3.802282e-03 0.0235873880 2.669952e-04    0.1863382 0.7859985             ep      7.334051    0.5976186
2 5.217008e-07 2.217998e-05 0.0002165258 4.782663e-06    0.8478323 0.1519237        Rachlin      6.943280    0.5302546
  probable.Log10AUC
1         0.8438519
2         0.7925787
```

![Alt text](Figure_Model_AUC_1.png?raw=true "Model AUC Visuals")

![Alt text](Figure_Model_AUC_2.png?raw=true "Model AUC Visuals")

#### Model-based Area (logarithmic scaling)

```r
dat <- data.frame(X=  c(1, 30,  180, 540, 1080, 2160),
                  Y=  c(1, 0.9, 0.8, 0.7, 0.6,  0.4),
                  ids= c(1, 1,   1,   1,   1,   1))

dat2<- data.frame(X=  c(1, 30,  180, 540, 1080, 2160),
                  Y=  c(1, 0.9, 0.7, 0.6, 0.5,  0.4),
                  ids= c(2, 2,   2,   2,   2,   2))

dat <- rbind(dat, dat2)

dat$Y <- dat$Y * 100

library(discountingtools)

results <- discountingModelSelection(dat,
                                     summarize = TRUE,
                                     figures = "ed50",
                                     idCol = "ids",
                                     A = 100)

Calculating series id: 1
Model Comparison Summary:
Rank #:  1  =  ep
 Probability  =  0.785999
 BIC  =  -22.2279
 Most Probable ln(ED50)  =  7.334051
 Most Probable Model AUC  =  0.597619
 Most Probable Model AUC (log)  =  0.843852

Rank #:  2  =  Rachlin
 Probability  =  0.186338
 BIC  =  -19.349116

Rank #:  3  =  Mazur
 Probability  =  0.023587
 BIC  =  -15.215414

Rank #:  4  =  exp
 Probability  =  0.003802
 BIC  =  -11.565193

Rank #:  5  =  BD
 Probability  =  0.000267
 BIC  =  -6.252941

Rank #:  6  =  noise
 Probability  =  7e-06
 BIC  =  1.128501
 
print(results)

  id JB.C1 JB.C2 noise.mean noise.RMSE noise.BIC noise.AIC   exp.lnk   exp.RMSE    exp.BIC    exp.AIC Mazur.lnk Mazur.RMSE
1  1  TRUE  TRUE  0.7333333  0.2160247  1.128501  1.544982 -7.615692 0.07500726 -11.565193 -11.148712 -7.217119 0.05533467
2  2  TRUE  TRUE  0.6833333  0.2316607  1.967073  2.383554 -7.369143 0.12400210  -5.532629  -5.116148 -6.846340 0.08482117
  Mazur.BIC  Mazur.AIC   BD.beta  BD.delta   BD.RMSE    BD.BIC    BD.AIC Rachlin.lnk Rachlin.s Rachlin.RMSE Rachlin.BIC
1 -15.21541 -14.798933 0.9999994 0.9993449 0.1124525 -6.252941 -5.628220   -5.000000 0.6800000   0.03775746   -19.34912
2 -10.08967  -9.673188 0.9999967 0.9991767 0.1542008 -2.464244 -1.839522   -3.941156 0.5676217   0.02057374   -26.63513
  Rachlin.AIC    ep.lnk      ep.s    ep.RMSE    ep.BIC    ep.AIC noise.BF    exp.BF  Mazur.BF     BD.BF Rachlin.BF    ep.BF
1   -18.72439 -7.976814 0.5702152 0.02970408 -22.22790 -21.60318        1 570.69054 3540.2682 40.073737   27967.79 117971.8
2   -26.01040 -7.796464 0.4521181 0.02740064 -23.19652 -22.57180        1  42.51475  415.0383  9.167443 1625131.21 291208.5
    noise.prob     exp.prob   Mazur.prob      BD.prob Rachlin.prob   ep.prob probable.model probable.ED50 probable.AUC
1 6.662599e-06 3.802282e-03 0.0235873880 2.669952e-04    0.1863382 0.7859985             ep      7.334051    0.5976186
2 5.217008e-07 2.217998e-05 0.0002165258 4.782663e-06    0.8478323 0.1519237        Rachlin      6.943280    0.5302546
  probable.Log10AUC
1         0.8438519
2         0.7925787

```

![Alt text](Figure_Model_AUC_Log10_1.png?raw=true "Log10 Model AUC Visuals")

![Alt text](Figure_Model_AUC_Log10_2.png?raw=true "Log10 Model AUC Visuals")

### Model Candidates (Possible models within the model selection approach)
-------------------
* Noise Model: Intecept-only comparison model (included by default)
* Exponential Model: Samuelson, P. A. (1937). A note on measurement of utility. The Review of Economic Studies, 4(2), 155–161. https://doi.org/10.2307/2967612
* Hyperbolic Model: Mazur, J. E. (1987). An adjusting procedure for studying delayed reinforcement. In Quantitative analysis of behavior: Vol. 5. The effect of delay and intervening events on reinforcement value (pp. 55–73). Hillsdale, NJ: Erlbaum.
* Beta Delta Model: Laibson, D. (1997). Golden eggs and hyperbolic discounting. The Quarterly Journal of Economics, 112(2), 443–478. https://doi.org/10.1162/003355397555253
* Green & Myerson Model: Green, L., & Myerson, J. (2004). A discounting framework for choice with delayed and probabilistic rewards. Psychological Bulletin, 130(5), 769–792. https://doi.org/10.1037/0033-2909.130.5.769
* Rachlin Model: Rachlin, H. (2006). Notes on discounting. Journal of the Experimental Analysis of Behavior, 85(3), 425–435. https://doi.org/10.1901/jeab.2006.85-05
* Ebert & Prelec Model: Ebert, J. E. J., & Prelec, D. (2007). The Fragility of Time: Time-Insensitivity and Valuation of the Near and Far Future. Management Science, 53(9), 1423–1438. https://doi.org/10.1287/mnsc.1060.0671

### Referenced Works (academic works)
-------------------
The Small N Stats Discounting Model Selector is based on the following academic works:
* Franck, C. T., Koffarnus, M. N., House, L. L., & Bickel, W. K. (2015). Accurate characterization of delay discounting: a multiple model approach using approximate Bayesian model selection and a unified discounting measure. Journal of the Experimental Analysis of Behavior, 103(1), 218–233. https://doi.org/10.1002/jeab.128
* Gilroy, S. P., Franck, C. T., & Hantula, D. A. (2017). The discounting model selector: Statistical software for delay discounting applications. Journal of the Experimental Analysis of Behavior, 107(3), 388–401. https://doi.org/10.1002/jeab.257
* Borges, A. M., Kuang, J., Milhorn, H., & Yi, R. (2016). An alternative approach to calculating Area-Under-the-Curve (AUC) in delay discounting research. Journal of the Experimental Analysis of Behavior, 106(2), 145–155. https://doi.org/10.1002/jeab.219

### Acknowledgments
-------------------
* Donald A. Hantula, Decision Making Laboratory, Temple University [Site](http://astro.temple.edu/~hantula/)
* Chris Franck, Laboratory for Interdisciplinary Statistical Analysis - Virginia Tech

### Questions, Suggestions, and Contributions
---------------------------------------------

Questions? Suggestions for features? <shawn.gilroy@temple.edu>.

### License
-----------

GPL Version 2+
