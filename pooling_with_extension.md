Pooling with extension
================
Erik Ø. Sørensen
August 7, 2015

We have a table in the paper in which we pool the preferred data with the extension. This file generates the estimates and tests for equality of parameters across these experiments. The document was changed in a minor way July 20th 2017 to 1) remove internal discussions about response to the editor and referees in the revision process and 2) addition of version information at the end. Minor change October 13, 2018, in order to render markdown output in GitHub format.

First I need to extract and run the library, load the data

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: parallel

Estimation on the extension data
================================

``` r
ll.2 <- function(b1.c=0, b2.c=-0.06, b3.c=-3, b4.c=0.1, b5.c=-1000, mu=-0.85, lsigma=0.8) {
    pvector <- c(b1.c, b2.c, b3.c, b4.c, b5.c, mu, lsigma)
    ll.luck(ds2, pvector)
}
```

Without constraints:

``` r
fit.later1 <- mle(minuslog = ll.2, fixed=list(b1.c=0, b5.c=-1000), 
              method="BFGS")
summary(fit.later1)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll.2, method = "BFGS", fixed = list(b1.c = 0, 
    ##     b5.c = -1000))
    ## 
    ## Coefficients:
    ##          Estimate   Std. Error
    ## b2.c    0.7508749    0.3562485
    ## b3.c   -8.4538147   22.9255628
    ## b4.c    0.1495024    0.4093543
    ## mu     -0.5802732    0.1871551
    ## lsigma  2.9969365 1977.2139337
    ## 
    ## -2 log L: 870.7355

``` r
knitr::kable(delta.method(grtrans, vcov(fit.later1), coef(fit.later1)))
```

|           |    Estimate|    Std. error|
|-----------|-----------:|-------------:|
| share.se  |   0.2336273|  6.016000e-02|
| share.l   |   0.4950218|  7.236300e-02|
| share.le  |   0.0000498|  1.141200e-03|
| share.cc  |   0.2713011|  6.570920e-02|
| share.scc |   0.0000000|  0.000000e+00|
| mu        |  -0.5802732|  1.871551e-01|
| sigma     |  20.0240985|  3.959193e+04|

``` r
logLik(fit.later1)
```

    ## 'log Lik.' -435.3677 (df=5)

Now constraining the LE share to be zero:

``` r
fit.later2 <- mle(minuslog = ll.2, fixed=list(b1.c=0,b3.c=-1000, b5.c=-1000), 
              method="BFGS")
summary(fit.later2)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll.2, method = "BFGS", fixed = list(b1.c = 0, 
    ##     b3.c = -1000, b5.c = -1000))
    ## 
    ## Coefficients:
    ##          Estimate Std. Error
    ## b2.c    0.7504445  0.3583774
    ## b4.c    0.1491543  0.4094136
    ## mu     -0.5803687  0.1871673
    ## lsigma  2.0763186 10.6463332
    ## 
    ## -2 log L: 870.732

``` r
knitr::kable(delta.method(grtrans, vcov(fit.later2), coef(fit.later2)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.2337108|   0.0603795|
| share.l   |   0.4949856|   0.0728975|
| share.le  |   0.0000000|   0.0000000|
| share.cc  |   0.2713036|   0.0658208|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -0.5803687|   0.1871673|
| sigma     |   7.9750553|  84.9050964|

``` r
logLik(fit.later2)
```

    ## 'log Lik.' -435.366 (df=4)

Estimation on the pooled data.
==============================

``` r
ll.all <- function(b1.c=0, b2.c=-0.06, b3.c=-3, b4.c=0.1, b5.c=-1000, mu=-0.85, lsigma=0.8) {
    pvector <- c(b1.c, b2.c, b3.c, b4.c, b5.c, mu, lsigma)
    ll.luck(ds.all, pvector)
}
```

Without constraints:

``` r
fit.all1 <- mle(minuslog = ll.all, fixed=list(b1.c=0, b5.c=-1000), 
              method="BFGS")
summary(fit.all1)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll.all, method = "BFGS", fixed = list(b1.c = 0, 
    ##     b5.c = -1000))
    ## 
    ## Coefficients:
    ##           Estimate Std. Error
    ## b2.c    0.18517959  0.2001821
    ## b3.c   -3.82245838  1.3100968
    ## b4.c    0.03777019  0.2190664
    ## mu     -0.72003172  0.1432765
    ## lsigma  0.84354981  0.1858683
    ## 
    ## -2 log L: 2510.506

``` r
knitr::kable(delta.method(grtrans, vcov(fit.all1), coef(fit.all1)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.3063912|   0.0386808|
| share.l   |   0.3687218|   0.0398266|
| share.le  |   0.0067020|   0.0086722|
| share.cc  |   0.3181850|   0.0406907|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -0.7200317|   0.1432765|
| sigma     |   2.3246043|   0.4320704|

``` r
logLik(fit.all1)
```

    ## 'log Lik.' -1255.253 (df=5)

Now constraining the LE share to be zero:

``` r
fit.all2 <- mle(minuslog = ll.all, fixed=list(b1.c=0,b3.c=-1000, b5.c=-1000), 
              method="BFGS")
summary(fit.all2)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll.all, method = "BFGS", fixed = list(b1.c = 0, 
    ##     b3.c = -1000, b5.c = -1000))
    ## 
    ## Coefficients:
    ##           Estimate Std. Error
    ## b2.c    0.17897367  0.1997244
    ## b4.c    0.03708869  0.2187981
    ## mu     -0.74177234  0.1378824
    ## lsigma  0.84720084  0.1840516
    ## 
    ## -2 log L: 2511.85

``` r
knitr::kable(delta.method(grtrans, vcov(fit.all2), coef(fit.all2)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.3092362|   0.0388077|
| share.l   |   0.3698431|   0.0398915|
| share.le  |   0.0000000|   0.0000000|
| share.cc  |   0.3209207|   0.0408285|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -0.7417723|   0.1378824|
| sigma     |   2.3331070|   0.4294121|

``` r
logLik(fit.all2)
```

    ## 'log Lik.' -1255.925 (df=4)

Testing pooling
===============

First testing for the unconstrained model:

``` r
diff1 <- -2*( logLik(fit.all1) - (-817.0486 + logLik(fit.later1)))
diff1
```

    ## 'log Lik.' 5.673726 (df=5)

``` r
pchisq(diff1, df=6, lower.tail=FALSE)
```

    ## 'log Lik.' 0.4607131 (df=5)

Now testing with the share of LE constrained to zero:

``` r
diff2 <- -2*(logLik(fit.all2)  - (-817.836 + logLik(fit.later2)))
diff2
```

    ## 'log Lik.' 5.446122 (df=4)

``` r
pchisq(diff2, df=5, lower.tail=FALSE)
```

    ## 'log Lik.' 0.36389 (df=4)

For footnote to paper: We cannot reject the restricted models in columns (3) and (4) against the unrestricted one in column (1) and (2) of table 9 and (1) and (2) of table 6 (*χ*<sub>6</sub><sup>2</sup> = 5.67, *p* = 0.46 for the model with LE included, *χ*<sub>5</sub><sup>2</sup> = 5.44, *p* = 0.36 with the lambda^LE=0).

Session Info
============

This documents the versions used of R and the involved packages to calculate these numbers:

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.1 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=nb_NO.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=nb_NO.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=nb_NO.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=nb_NO.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ## [1] doMC_1.3.5        iterators_1.0.10  foreach_1.4.4     numDeriv_2016.8-1
    ## [5] glmmML_1.0.3      foreign_0.8-70   
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_0.12.19     codetools_0.2-15 digest_0.6.17    rprojroot_1.3-2 
    ##  [5] backports_1.1.2  magrittr_1.5     evaluate_0.12    highr_0.7       
    ##  [9] stringi_1.2.4    rmarkdown_1.10   tools_3.5.1      stringr_1.3.1   
    ## [13] yaml_2.2.0       compiler_3.5.1   htmltools_0.3.6  knitr_1.20
