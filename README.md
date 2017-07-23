# discountingtools
R package for performing delay discounting analyses.  Features include approximate Bayesian model selection, derivation of the Effective Delay 50 (ED50) and both normal and log10 transformation model area using Numerical Integration

### Version
0.0.0.0001 (alpha)

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
dat <- data.frame(X=c(1,30,180,540,1080, 2160), Y=c(1,0.9,0.8,0.7,0.6, 0.4))

results <- discountingModelSelection(dat, 
                                     models = c("exponential", "hyperbolic", "bd", "mg", "rachlin"), 
                                     figures = "ed50",
                                     lineSize = 1.5)
                                     
  noise.mean noise.RMSE noise.BIC noise.AIC   exp.lnk   exp.RMSE   exp.BIC   exp.AIC Mazur.lnk Mazur.RMSE Mazur.BIC Mazur.AIC
1  0.7333333  0.2160247  1.128501  1.544982 -7.615692 0.07500726 -11.56519 -11.14871 -7.217119 0.05533467 -15.21541 -14.79893
    BD.beta  BD.delta   BD.RMSE    BD.BIC   BD.AIC    MG.lnk     MG.s    MG.RMSE    MG.BIC    MG.AIC Rachlin.lnk Rachlin.s
1 0.9999994 0.9993449 0.1124525 -6.252941 -5.62822 -5.657531 0.378083 0.04995242 -15.99046 -15.36574          -5      0.68
  Rachlin.RMSE Rachlin.BIC Rachlin.AIC noise.BF   exp.BF Mazur.BF    BD.BF    MG.BF Rachlin.BF   noise.prob   exp.prob
1   0.03775746   -19.34912   -18.72439        1 570.6905 3540.268 40.07374 5215.976   27967.79 2.678395e-05 0.01528535
  Mazur.prob     BD.prob   MG.prob Rachlin.prob probable.model probable.ED50 probable.AUC probable.Log10AUC
1 0.09482236 0.001073333 0.1397044    0.7490878        Rachlin      7.352941    0.5972851         0.8461974
```

![Alt text](Figure_ED50.png?raw=true "ED50 Visuals")

#### Model-based Area

```r
dat <- data.frame(X=c(1,30,180,540,1080, 2160), Y=c(1,0.9,0.8,0.7,0.6, 0.4))

results <- discountingModelSelection(dat, 
                                     models = c("exponential", "hyperbolic", "bd", "mg", "rachlin"), 
                                     figures = "auc",
                                     lineSize = 1.5)
                                     
  noise.mean noise.RMSE noise.BIC noise.AIC   exp.lnk   exp.RMSE   exp.BIC   exp.AIC Mazur.lnk Mazur.RMSE Mazur.BIC Mazur.AIC
1  0.7333333  0.2160247  1.128501  1.544982 -7.615692 0.07500726 -11.56519 -11.14871 -7.217119 0.05533467 -15.21541 -14.79893
    BD.beta  BD.delta   BD.RMSE    BD.BIC   BD.AIC    MG.lnk     MG.s    MG.RMSE    MG.BIC    MG.AIC Rachlin.lnk Rachlin.s
1 0.9999994 0.9993449 0.1124525 -6.252941 -5.62822 -5.657531 0.378083 0.04995242 -15.99046 -15.36574          -5      0.68
  Rachlin.RMSE Rachlin.BIC Rachlin.AIC noise.BF   exp.BF Mazur.BF    BD.BF    MG.BF Rachlin.BF   noise.prob   exp.prob
1   0.03775746   -19.34912   -18.72439        1 570.6905 3540.268 40.07374 5215.976   27967.79 2.678395e-05 0.01528535
  Mazur.prob     BD.prob   MG.prob Rachlin.prob probable.model probable.ED50 probable.AUC probable.Log10AUC
1 0.09482236 0.001073333 0.1397044    0.7490878        Rachlin      7.352941    0.5972851         0.8461974
```

![Alt text](Figure_Model_AUC.png?raw=true "Model AUC Visuals")

#### Model-based Area (logarithmic scaling)

```r
dat <- data.frame(X=c(1,30,180,540,1080, 2160), Y=c(1,0.9,0.8,0.7,0.6, 0.4))

results <- discountingModelSelection(dat, 
                                     models = c("exponential", "hyperbolic", "bd", "mg", "rachlin"), 
                                     figures = "logauc",
                                     lineSize = 1.5)

  noise.mean noise.RMSE noise.BIC noise.AIC   exp.lnk   exp.RMSE   exp.BIC   exp.AIC Mazur.lnk Mazur.RMSE Mazur.BIC Mazur.AIC
1  0.7333333  0.2160247  1.128501  1.544982 -7.615692 0.07500726 -11.56519 -11.14871 -7.217119 0.05533467 -15.21541 -14.79893
    BD.beta  BD.delta   BD.RMSE    BD.BIC   BD.AIC    MG.lnk     MG.s    MG.RMSE    MG.BIC    MG.AIC Rachlin.lnk Rachlin.s
1 0.9999994 0.9993449 0.1124525 -6.252941 -5.62822 -5.657531 0.378083 0.04995242 -15.99046 -15.36574          -5      0.68
  Rachlin.RMSE Rachlin.BIC Rachlin.AIC noise.BF   exp.BF Mazur.BF    BD.BF    MG.BF Rachlin.BF   noise.prob   exp.prob
1   0.03775746   -19.34912   -18.72439        1 570.6905 3540.268 40.07374 5215.976   27967.79 2.678395e-05 0.01528535
  Mazur.prob     BD.prob   MG.prob Rachlin.prob probable.model probable.ED50 probable.AUC probable.Log10AUC
1 0.09482236 0.001073333 0.1397044    0.7490878        Rachlin      7.352941    0.5972851         0.8461974

```

![Alt text](Figure_Model_AUC_Log10.png?raw=true "Log10 Model AUC Visuals")

### Referenced Works (academic works)
The Small N Stats Discounting Model Selector is based on the following academic works:
* Franck, C. T., Koffarnus, M. N., House, L. L. & Bickel W. K. (2015). Accurate characterization of delay discounting: a multiple model approach using approximate Bayesian model selection and a unified discounting measure. Journal of the Experimental Analysis of Behavior, 103, 218-33.
* Gilroy, S. P., Franck, C. T. & Hantula, D. A. (2017). The Discounting Model Selector: Statistical software for delay discounting applications. Journal of the Experimental Analysis of Behavior.

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
