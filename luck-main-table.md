Main table
================
Erik Ø. Sørensen
1.  feb. 2015

This document is to estimate the specifications in the main table. Make sure to first extract the R-code from `luck-model.Rmd` (using `knitr::purl('luck-model.Rmd')`, unfortunately, knitr cannot call itself to do it automatically). Minor change July 20th, 2017 to document the version of packages used at the end of the document. Minor change October 13, 2018, in order to render markdown output in GitHub format.

Setup - loading library and data
================================

First I need to extract and run the library

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: parallel

For estimation, first look at the first round of data.

``` r
ds.luck1 <- datastruct(players[players$round==1,], situations)
```

Next, define the likelihood function.

``` r
ll <- function(b1.c=0, b2.c=-0.06, b3.c=-3, b4.c=0.1, b5.c=-1000, mu=-0.85, lsigma=0.8) {
     pvector <- c(b1.c, b2.c, b3.c, b4.c, b5.c, mu, lsigma)
    ll.luck(ds.luck1, pvector)
}
```

Specifications
==============

First, estimate the main specification.

``` r
fit.all <- mle(minuslog = ll, fixed=list(b1.c=0, b5.c=-1000), 
               method="BFGS", control=list(trace=TRUE, REPORT=2))
```

    ## initial  value 817.245336 
    ## iter   2 value 817.162989
    ## iter   4 value 817.130089
    ## iter   6 value 817.071719
    ## iter   8 value 817.048638
    ## final  value 817.048588 
    ## converged

``` r
summary(fit.all)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll, method = "BFGS", fixed = list(b1.c = 0, b5.c = -1000), 
    ##     control = list(trace = TRUE, REPORT = 2))
    ## 
    ## Coefficients:
    ##            Estimate Std. Error
    ## b2.c   -0.075838630  0.2456603
    ## b3.c   -3.501276788  1.2602163
    ## b4.c    0.001694048  0.2598799
    ## mu     -0.872913927  0.1789912
    ## lsigma  0.765249182  0.1661261
    ## 
    ## -2 log L: 1634.097

``` r
logLik(fit.all)
```

    ## 'log Lik.' -817.0486 (df=5)

``` r
knitr::kable(delta.method(grtrans, vcov(fit.all), coef(fit.all)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.3379726|   0.0487091|
| share.l   |   0.3132890|   0.0463072|
| share.le  |   0.0101929|   0.0126138|
| share.cc  |   0.3385456|   0.0511019|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -0.8729139|   0.1789912|
| sigma     |   2.1495299|   0.3570929|

The ideals are, in order: SE, L, LE, CC (and sCC, not estimated in this document).

Remove LE:

``` r
fit.2 <- mle(minuslog = ll, fixed=list(b1.c=0, b3.c=-1000, b5.c=-1000), method="BFGS",
             control=list(trace=TRUE, REPORT=2))
```

    ## initial  value 818.010590 
    ## iter   2 value 817.888177
    ## iter   4 value 817.858940
    ## iter   6 value 817.838863
    ## iter   8 value 817.836046
    ## iter   8 value 817.836046
    ## iter   8 value 817.836046
    ## final  value 817.836046 
    ## converged

``` r
summary(fit.2)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll, method = "BFGS", fixed = list(b1.c = 0, b3.c = -1000, 
    ##     b5.c = -1000), control = list(trace = TRUE, REPORT = 2))
    ## 
    ## Coefficients:
    ##            Estimate Std. Error
    ## b2.c   -0.087386725  0.2447495
    ## b4.c    0.001767553  0.2595727
    ## mu     -0.911380255  0.1726082
    ## lsigma  0.760964914  0.1559004
    ## 
    ## -2 log L: 1635.672

``` r
knitr::kable(delta.method(grtrans, vcov(fit.2), coef(fit.2)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.3426897|   0.0489363|
| share.l   |   0.3140143|   0.0463975|
| share.le  |   0.0000000|   0.0000000|
| share.cc  |   0.3432960|   0.0514308|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -0.9113803|   0.1726082|
| sigma     |   2.1403405|   0.3336800|

``` r
logLik(fit.2)
```

    ## 'log Lik.' -817.836 (df=4)

Now, remove LE and SE. When I remove one of the substantial ideals, the likelihood becomes quite ill-behaved and sensitive to starting values. These starting values converge at the best (local) point of convergence that I have found.

``` r
fit.3 <- mle(minuslog = ll, start=list(b4.c=-0.292, mu=-10.8, lsigma=2.398), 
             fixed=list(b1.c=-1000, b2.c=0, b3.c=-1000, b5.c=-1000), method="BFGS",
             control=list(trace=TRUE, REPORT=2))
```

    ## initial  value 937.853972 
    ## iter   2 value 937.845896
    ## iter   4 value 937.835243
    ## iter   6 value 937.830871
    ## iter   6 value 937.830865
    ## iter   6 value 937.830865
    ## final  value 937.830865 
    ## converged

``` r
summary(fit.3)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll, start = list(b4.c = -0.292, mu = -10.8, lsigma = 2.398), 
    ##     method = "BFGS", fixed = list(b1.c = -1000, b2.c = 0, b3.c = -1000, 
    ##         b5.c = -1000), control = list(trace = TRUE, REPORT = 2))
    ## 
    ## Coefficients:
    ##           Estimate Std. Error
    ## b4.c    -0.2883506  0.2817704
    ## mu     -10.8043151 17.7347788
    ## lsigma   2.3955079  1.5797253
    ## 
    ## -2 log L: 1875.662

``` r
knitr::kable(delta.method(grtrans, vcov(fit.3), coef(fit.3)))
```

|           |     Estimate|  Std. error|
|-----------|------------:|-----------:|
| share.se  |    0.0000000|   0.0000000|
| share.l   |    0.5715923|   0.0689984|
| share.le  |    0.0000000|   0.0000000|
| share.cc  |    0.4284077|   0.0689984|
| share.scc |    0.0000000|   0.0000000|
| mu        |  -10.8043151|  17.7347788|
| sigma     |   10.9737699|  17.3355417|

``` r
logLik(fit.3)
```

    ## 'log Lik.' -937.8309 (df=3)

Remove LE and L:

``` r
fit.4 <- mle(minuslog = ll, start=list(b4.c=0.48, mu=-1.3, lsigma=1.9), 
             fixed=list(b1.c=0, b2.c=-1000, b3.c=-1000, b5.c=-1000),
             method="BFGS", control=list(trace=TRUE, REPORT=2))
```

    ## initial  value 976.044768 
    ## iter   2 value 976.036847
    ## iter   4 value 976.036297
    ## iter   6 value 976.035778
    ## iter   8 value 969.353285
    ## iter  10 value 965.461492
    ## iter  12 value 965.218008
    ## iter  14 value 965.162087
    ## iter  14 value 965.162086
    ## iter  14 value 965.162086
    ## final  value 965.162086 
    ## converged

``` r
summary(fit.4)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll, start = list(b4.c = 0.48, mu = -1.3, lsigma = 1.9), 
    ##     method = "BFGS", fixed = list(b1.c = 0, b2.c = -1000, b3.c = -1000, 
    ##         b5.c = -1000), control = list(trace = TRUE, REPORT = 2))
    ## 
    ## Coefficients:
    ##          Estimate Std. Error
    ## b4.c    0.6525592  0.2371877
    ## mu     -1.5711072  0.1418440
    ## lsigma  0.6536483  0.1036459
    ## 
    ## -2 log L: 1930.324

``` r
knitr::kable(delta.method(grtrans, vcov(fit.4), coef(fit.4)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.3424131|   0.0534067|
| share.l   |   0.0000000|   0.0000000|
| share.le  |   0.0000000|   0.0000000|
| share.cc  |   0.6575869|   0.0534067|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -1.5711072|   0.1418440|
| sigma     |   1.9225420|   0.1992637|

``` r
logLik(fit.4)
```

    ## 'log Lik.' -965.1621 (df=3)

Remove LE and CC:

``` r
fit.5 <- mle(minuslog = ll, start=list(b2.c=0.034, mu=-1.6, lsigma=1.0),
             fixed=list(b1.c=0, b3.c=-1000, b4.c=-1000, b5.c=-1000), method="BFGS",
             control=list(trace=TRUE, REPORT=2))
```

    ## initial  value 905.730338 
    ## iter   2 value 904.810000
    ## iter   4 value 904.050351
    ## iter   6 value 904.039198
    ## final  value 904.039096 
    ## converged

``` r
summary(fit.5)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll, start = list(b2.c = 0.034, mu = -1.6, lsigma = 1), 
    ##     method = "BFGS", fixed = list(b1.c = 0, b3.c = -1000, b4.c = -1000, 
    ##         b5.c = -1000), control = list(trace = TRUE, REPORT = 2))
    ## 
    ## Coefficients:
    ##          Estimate Std. Error
    ## b2.c    0.1579661  0.2277637
    ## mu     -1.5998692  0.1566917
    ## lsigma  1.8130024  8.2446208
    ## 
    ## -2 log L: 1808.078

``` r
bhhh.fit5 <- vcov.bhhh(ds.luck1, fit.5)
knitr::kable(delta.method(grtrans, vcov(fit.5), coef(fit.5)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.4605904|   0.0565872|
| share.l   |   0.5394096|   0.0565872|
| share.le  |   0.0000000|   0.0000000|
| share.cc  |   0.0000000|   0.0000000|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -1.5998692|   0.1566917|
| sigma     |   6.1288211|  50.5298060|

``` r
bhhh.fit5 <- vcov.bhhh(ds.luck1, fit.5)
logLik(fit.5)
```

    ## 'log Lik.' -904.0391 (df=3)

Different way of calculating standard errors
============================================

The above estimates of standard errors are based on the defaults of the stats4 library, which we cite in the paper. Some standard errors (for *σ* are not well estimated). One could think that this results from how the stats4 library calculates standard errors: Based, as far as I understand, on taking the direct numerical Hessian of the log likelihood function. An alternative is the BHHH outer product of the gradient. I do that here (and conclude from this that it doesn't really matter which way we do it).

``` r
bhhh.all <- vcov.bhhh(ds.luck1, fit.all)
knitr::kable(delta.method(grtrans, bhhh.all, coef(fit.all)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.3379726|   0.0494650|
| share.l   |   0.3132890|   0.0489238|
| share.le  |   0.0101929|   0.0125590|
| share.cc  |   0.3385456|   0.0560461|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -0.8729139|   0.1257640|
| sigma     |   2.1495299|   0.4352818|

``` r
bhhh.fit2 <- vcov.bhhh(ds.luck1, fit.2)
knitr::kable(delta.method(grtrans, bhhh.fit2, coef(fit.2)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.3426897|   0.0498517|
| share.l   |   0.3140143|   0.0491362|
| share.le  |   0.0000000|   0.0000000|
| share.cc  |   0.3432960|   0.0563952|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -0.9113803|   0.1264824|
| sigma     |   2.1403405|   0.3831409|

``` r
bhhh.fit3 <- vcov.bhhh(ds.luck1, fit.3)
knitr::kable(delta.method(grtrans, bhhh.fit3, coef(fit.3)))
```

|           |     Estimate|    Std. error|
|-----------|------------:|-------------:|
| share.se  |    0.0000000|     0.0000000|
| share.l   |    0.5715923|     0.0740776|
| share.le  |    0.0000000|     0.0000000|
| share.cc  |    0.4284077|     0.0740776|
| share.scc |    0.0000000|     0.0000000|
| mu        |  -10.8043151|  1083.8924759|
| sigma     |   10.9737699|  1059.2424666|

``` r
bhhh.fit4 <- vcov.bhhh(ds.luck1, fit.4)
knitr::kable(delta.method(grtrans, bhhh.fit4, coef(fit.4)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.3424131|   0.0563150|
| share.l   |   0.0000000|   0.0000000|
| share.le  |   0.0000000|   0.0000000|
| share.cc  |   0.6575869|   0.0563150|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -1.5711072|   0.2125927|
| sigma     |   1.9225420|   0.2188782|

``` r
bhhh.fit5 <- vcov.bhhh(ds.luck1, fit.5)
knitr::kable(delta.method(grtrans, bhhh.fit5, coef(fit.5)))
```

|           |    Estimate|    Std. error|
|-----------|-----------:|-------------:|
| share.se  |   0.4605904|     0.0694511|
| share.l   |   0.5394096|     0.0694511|
| share.le  |   0.0000000|     0.0000000|
| share.cc  |   0.0000000|     0.0000000|
| share.scc |   0.0000000|     0.0000000|
| mu        |  -1.5998692|     0.1694018|
| sigma     |   6.1288211|  1889.0846302|

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
