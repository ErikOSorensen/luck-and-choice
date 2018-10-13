Insurance and no insurance
================
Erik Ø. Sørensen
August 4, 2015

This document is to estimate on the insurance-no-insurance data. The editor wants us to comment on whether sigma is the same or different between those who bought and those who didn't buy insurance. The document was changed in a minor way July 20th 2017 to 1) remove internal discussions about response to the editor in the revision process and 2) addition of version information at the end. Minor change October 13, 2018, in order to render markdown output in GitHub format.

Setup - loading library and subsetting data
===========================================

First I need to extract and run the library

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: parallel

Data subsetting (the global datasets are defined in `luck-model`):

``` r
ds.noin <- datastruct(players[players$round==1 & players$insurance==0,], situations)
ds.in   <- datastruct(players[players$round==1 & players$insurance==1,], situations)
```

For the likelihood ratio test we also need the full dataset of round 1.

``` r
ds <- datastruct(players[players$round==1,], situations)
```

On the insured
==============

Defining the likelihood function on this subset of data:

``` r
ll.in <- function(b1.c=0, b2.c=-0.06, b3.c=-3, b4.c=0.1, b5.c=-1000, mu=-0.85, lsigma=0.8) {
    pvector <- c(b1.c, b2.c, b3.c, b4.c, b5.c, mu, lsigma)
    ll.luck(ds.in, pvector)
}
```

Without constraints:

``` r
fit.in1 <- mle(minuslog = ll.in, fixed=list(b1.c=0, b5.c=-1000), 
              method="BFGS")
summary(fit.in1)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll.in, method = "BFGS", fixed = list(b1.c = 0, 
    ##     b5.c = -1000))
    ## 
    ## Coefficients:
    ##           Estimate Std. Error
    ## b2.c   -0.55602277  0.2899368
    ## b3.c   -3.75836684  1.3353361
    ## b4.c   -0.04783383  0.2718853
    ## mu     -0.90150077  0.1996237
    ## lsigma  0.68738103  0.1698333
    ## 
    ## -2 log L: 1293.559

``` r
knitr::kable(delta.method(grtrans, vcov(fit.in1), coef(fit.in1)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.3921416|   0.0559882|
| share.l   |   0.2248875|   0.0463762|
| share.le  |   0.0091454|   0.0120239|
| share.cc  |   0.3738255|   0.0583107|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -0.9015008|   0.1996237|
| sigma     |   1.9885009|   0.3377136|

``` r
logLik(fit.in1)
```

    ## 'log Lik.' -646.7797 (df=5)

Now constraining the LE share to be zero:

``` r
fit.in2 <- mle(minuslog = ll.in, fixed=list(b1.c=0,b3.c=-1000, b5.c=-1000), 
              method="BFGS")
summary(fit.in2)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll.in, method = "BFGS", fixed = list(b1.c = 0, 
    ##     b3.c = -1000, b5.c = -1000))
    ## 
    ## Coefficients:
    ##           Estimate Std. Error
    ## b2.c   -0.56663517  0.2895356
    ## b4.c   -0.04693526  0.2716383
    ## mu     -0.93486489  0.1927251
    ## lsigma  0.68340761  0.1647060
    ## 
    ## -2 log L: 1294.851

``` r
knitr::kable(delta.method(grtrans, vcov(fit.in2), coef(fit.in2)))
```

|           |    Estimate|  Std. error|
|-----------|-----------:|-----------:|
| share.se  |   0.3965766|   0.0562558|
| share.l   |   0.2250301|   0.0464549|
| share.le  |   0.0000000|   0.0000000|
| share.cc  |   0.3783933|   0.0586561|
| share.scc |   0.0000000|   0.0000000|
| mu        |  -0.9348649|   0.1927251|
| sigma     |   1.9806154|   0.3262191|

``` r
logLik(fit.in2)
```

    ## 'log Lik.' -647.4254 (df=4)

On the non-insured
==================

Defining the likelihood function on this subset of data:

``` r
ll.noin <- function(b1.c=0, b2.c=-0.06, b3.c=-3, b4.c=0.1, b5.c=-1000, mu=-0.85, lsigma=0.8) {
    pvector <- c(b1.c, b2.c, b3.c, b4.c, b5.c, mu, lsigma)
    ll.luck(ds.noin, pvector)
}
```

First the model in without constraints:

``` r
fit.noin1 <- mle(minuslog = ll.noin, 
               fixed=list(b1.c=0, b5.c=-1000), method="BFGS")
summary(fit.noin1)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll.noin, method = "BFGS", fixed = list(b1.c = 0, 
    ##     b5.c = -1000))
    ## 
    ## Coefficients:
    ##          Estimate Std. Error
    ## b2.c    1.8200531  0.7575301
    ## b3.c   -0.5142790  1.6380698
    ## b4.c    0.5478258  0.9107973
    ## mu     -0.6818047  0.3879309
    ## lsigma  2.3274832 25.7853882
    ## 
    ## -2 log L: 325.6523

``` r
knitr::kable(delta.method(grtrans, vcov(fit.noin1), coef(fit.noin1)))
```

|           |    Estimate|   Std. error|
|-----------|-----------:|------------:|
| share.se  |   0.1052675|    0.0708836|
| share.l   |   0.6497307|    0.1160929|
| share.le  |   0.0629428|    0.0871360|
| share.cc  |   0.1820590|    0.0956431|
| share.scc |   0.0000000|    0.0000000|
| mu        |  -0.6818047|    0.3879309|
| sigma     |  10.2521067|  264.3545505|

``` r
logLik(fit.noin1)
```

    ## 'log Lik.' -162.8261 (df=5)

Now constraining the LE share to be zero:

``` r
fit.noin2 <- mle(minuslog = ll.noin, 
               fixed=list(b1.c=0, b3.c=-1000, b5.c=-1000), method="BFGS")
summary(fit.noin2)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll.noin, method = "BFGS", fixed = list(b1.c = 0, 
    ##     b3.c = -1000, b5.c = -1000))
    ## 
    ## Coefficients:
    ##          Estimate Std. Error
    ## b2.c    1.7462899  0.7471303
    ## b4.c    0.5343893  0.9048596
    ## mu     -0.8706983  0.3877741
    ## lsigma  2.4843519 32.7602327
    ## 
    ## -2 log L: 326.2609

``` r
knitr::kable(delta.method(grtrans, vcov(fit.noin2), coef(fit.noin2)))
```

|           |    Estimate|   Std. error|
|-----------|-----------:|------------:|
| share.se  |   0.1184877|    0.0772937|
| share.l   |   0.6793243|    0.1124415|
| share.le  |   0.0000000|    0.0000000|
| share.cc  |   0.2021880|    0.1007194|
| share.scc |   0.0000000|    0.0000000|
| mu        |  -0.8706983|    0.3877741|
| sigma     |  11.9933454|  392.9047854|

``` r
logLik(fit.noin2)
```

    ## 'log Lik.' -163.1305 (df=4)

Constraining *σ* to be the same
===============================

In order to use the variant likelihood functions, we need a function that indicates which group we are are in:

``` r
setGeneric("insured", function(s) standardGeneric("insured"))
```

    ## [1] "insured"

``` r
setMethod("insured", "individual",
          function(s) {
            if (s@insurance==1)
              1
            else
              2
          }
)
```

Defining the likelihood function on all data and with the group indicator function:

``` r
ll.in.2K <- function(b1.c=0, b2.c=-0.06, b3.c=-3, b4.c=0.1, b5.c=-1000, mu=-0.85, lsigma=0.8,
                     b1.d=0.01, b2.d=0.01, b3.d=0.01, b4.d=0.01, b5.d=0.01, mu.d=0.01, lsigma.d=0.01) {
    pvector <- c(b1.c, b2.c, b3.c, b4.c, b5.c, mu, lsigma, b1.d, b2.d, b3.d, b4.d, b5.d, mu.d, lsigma.d)
    ll.luck.2K(ds, pvector, insured)
}
```

Now estimating on the unconstrained model to confirm that the variant estimation command works:

``` r
fit.2K <- mle(minuslog = ll.in.2K, 
              start=list(b2.c=-0.57, b4.c=-0.05, mu=-0.93, lsigma=0.68,
                         b2.d=2.3, b4.d=0.58, mu.d=0.06, lsigma.d=1.8),
              fixed=list(b1.c=0, b1.d=0, b3.c=-1000, b3.d=0, b5.c=-1000, b5.d=0),
              method="BFGS")
summary(fit.2K)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll.in.2K, start = list(b2.c = -0.57, b4.c = -0.05, 
    ##     mu = -0.93, lsigma = 0.68, b2.d = 2.3, b4.d = 0.58, mu.d = 0.06, 
    ##     lsigma.d = 1.8), method = "BFGS", fixed = list(b1.c = 0, 
    ##     b1.d = 0, b3.c = -1000, b3.d = 0, b5.c = -1000, b5.d = 0))
    ## 
    ## Coefficients:
    ##             Estimate Std. Error
    ## b2.c     -0.56663491  0.2895351
    ## b4.c     -0.04694588  0.2716387
    ## mu       -0.93487556  0.1927229
    ## lsigma    0.68341408  0.1647041
    ## b2.d      2.31288732  0.8012553
    ## b4.d      0.58122766  0.9447593
    ## mu.d     -0.87073766  0.3877679
    ## lsigma.d  1.80023781 32.6438629
    ## 
    ## -2 log L: 1621.112

``` r
logLik(fit.2K)
```

    ## 'log Lik.' -810.5559 (df=8)

We see that these numbers conform to those estimated seperately above, so this seems like a success.

Estimating with constrained *σ*
-------------------------------

First with non-constrained share of luck egalitarians:

``` r
fit.2K.1 <- mle(minuslog = ll.in.2K, 
              start=list(b2.c=-0.57, b3.c=-0.9, b4.c=-0.05, mu=-0.93, lsigma=1.2,
                         b2.d=2.3, b3.d=0.1, b4.d=0.58, mu.d=0.06),
              fixed=list(b1.c=0, b1.d=0, b5.c=-1000, b5.d=0, lsigma.d=0),
              method="BFGS")
summary(fit.2K.1)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll.in.2K, start = list(b2.c = -0.57, b3.c = -0.9, 
    ##     b4.c = -0.05, mu = -0.93, lsigma = 1.2, b2.d = 2.3, b3.d = 0.1, 
    ##     b4.d = 0.58, mu.d = 0.06), method = "BFGS", fixed = list(b1.c = 0, 
    ##     b1.d = 0, b5.c = -1000, b5.d = 0, lsigma.d = 0))
    ## 
    ## Coefficients:
    ##           Estimate Std. Error
    ## b2.c   -0.55094280  0.2904843
    ## b3.c   -3.74982848  1.3303203
    ## b4.c   -0.04848017  0.2723978
    ## mu     -0.91312685  0.1927261
    ## lsigma  0.73093619  0.1578585
    ## b2.d    2.23641411  0.7709021
    ## b3.d    3.33861650  1.9246491
    ## b4.d    0.53474284  0.9165443
    ## mu.d   -0.59880325  0.3764450
    ## 
    ## -2 log L: 1619.776

``` r
logLik(fit.2K.1)
```

    ## 'log Lik.' -809.8882 (df=9)

Now with the share of luck egalitarians constrained to zero:

``` r
fit.2K.2 <- mle(minuslog = ll.in.2K, 
              start=list(b2.c=-0.57, b4.c=-0.05, mu=-0.93, lsigma=1.2,
                         b2.d=2.3, b4.d=0.58, mu.d=0.06),
              fixed=list(b1.c=0, b1.d=0, b3.c=-1000, b3.d=0, 
                         b5.c=-1000, b5.d=0, lsigma.d=0),
              method="BFGS")
summary(fit.2K.2)
```

    ## Maximum likelihood estimation
    ## 
    ## Call:
    ## mle(minuslogl = ll.in.2K, start = list(b2.c = -0.57, b4.c = -0.05, 
    ##     mu = -0.93, lsigma = 1.2, b2.d = 2.3, b4.d = 0.58, mu.d = 0.06), 
    ##     method = "BFGS", fixed = list(b1.c = 0, b1.d = 0, b3.c = -1000, 
    ##         b3.d = 0, b5.c = -1000, b5.d = 0, lsigma.d = 0))
    ## 
    ## Coefficients:
    ##           Estimate Std. Error
    ## b2.c   -0.55882154  0.2903652
    ## b4.c   -0.04792119  0.2724303
    ## mu     -0.94991306  0.1864074
    ## lsigma  0.75082019  0.1496357
    ## b2.d    2.16321370  0.7559477
    ## b4.d    0.53617386  0.9066019
    ## mu.d   -0.83005681  0.4091584
    ## 
    ## -2 log L: 1622.13

``` r
logLik(fit.2K.2)
```

    ## 'log Lik.' -811.0652 (df=7)

Tests and pvalues of common *σ*
-------------------------------

First testing column 1 and 3 are the same.

``` r
diff1 <- -2*(logLik(fit.2K.1) - (logLik(fit.in1)  + logLik(fit.noin1)))
diff1 
```

    ## 'log Lik.' 0.5647638 (df=9)

``` r
( p_non_constrained <- pchisq(diff1, df=1, lower.tail=FALSE))
```

    ## 'log Lik.' 0.4523472 (df=9)

Now testing with the share of LE constrained to zero (2 and 4):

``` r
diff2 <- -2*(logLik(fit.2K.2) - (logLik(fit.in2)  + logLik(fit.noin2)))
diff2
```

    ## 'log Lik.' 1.018713 (df=7)

``` r
( p_constrained <- pchisq(diff2, df=1, lower.tail=FALSE))
```

    ## 'log Lik.' 0.3128246 (df=7)

Standard analytic inference on maximum likelihood estimates is based on asymptotic theory, and asymptotically, confidence intervals in the interior of the parameter space are symmetric. But using only the local information at the unconstrained estimates might be restrictive, and we now include a footnote with likelihood ratio test for the equality of *σ* between these specifications (which cannot be rejected). Since the p-values are far from being close to significant, we believe that developing small-sample properties of *σ* will not contribute much to the paper.

Footnote for paper: The constraint that the *σ*s are the same in column 1 and 3 and 2 and 4 cannot be rejected (likelihood ratio tests based on a pooled and restricted model: *χ*<sub>1</sub><sup>2</sup> = 0.56, *p* = 0.45 and *χ*<sub>1</sub><sup>2</sup> = 1.02, *p* = 0.31, respectively).

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
