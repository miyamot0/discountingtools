# discountingtools
R package for performing delay discounting analyses.  Features include approximate Bayesian model selection, derivation of the Effective Delay 50 (ED50) and both normal and log10 transformation model area using Numerical Integration

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
```

#### Model-based Area

```r
dat <- data.frame(X=c(1,30,180,540,1080, 2160), Y=c(1,0.9,0.8,0.7,0.6, 0.4))

results <- discountingModelSelection(dat, 
                                     models = c("exponential", "hyperbolic", "bd", "mg", "rachlin"), 
                                     figures = "auc",
                                     lineSize = 1.5)
```

#### Model-based Area (logarithmic scaling)

```r
dat <- data.frame(X=c(1,30,180,540,1080, 2160), Y=c(1,0.9,0.8,0.7,0.6, 0.4))

results <- discountingModelSelection(dat, 
                                     models = c("exponential", "hyperbolic", "bd", "mg", "rachlin"), 
                                     figures = "logauc",
                                     lineSize = 1.5)
```
