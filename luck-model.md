Calculations with choice model: model
================
Erik Ø. Sørensen
1.  aug. 2015

This file contains library code for the choice model estimation, but not any actual estimation itself. Instead, all the sample and parameter restrictions are kept separate, one for each table. Minor change October 13, 2018, in order to render markdown output in GitHub format.

Preliminaries
=============

There are some dependencies: To load the S4 object system, to read Stata files, and to get a simple parallel computation on multicore machines.

``` r
require(methods)
library(foreign)
library(stats4)
library(glmmML)
library(numDeriv)
library(doMC)
```

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: parallel

``` r
registerDoMC(8)
getDoParWorkers()
```

    ## [1] 8

``` r
situations <- read.dta("data/DataErik_long.dta", 
                       convert.factors=FALSE)
players <- read.dta("data/DataErik_players.dta", 
                    convert.factors=FALSE)
```

The ideals and the choice model
===============================

The situations and the fairness ideals are enumerated:

|  situation| outcome1 | outcome2 |   e1|   e2|    X|  Fse|   Fl|  Fle|  Fcc|  Fscc|
|----------:|:---------|:---------|----:|----:|----:|----:|----:|----:|----:|-----:|
|          1| A        | B        |   24|    0|   24|   12|   24|   24|   24|    24|
|          2| A        | C        |   24|    0|   24|   12|   24|   18|   24|    24|
|          3| A        | B\*      |   24|   12|   36|   18|   24|   24|   18|    26|
|          4| B\*      | C\*      |   12|    0|   12|    6|   12|    6|    6|     6|
|          5| A        | C\*      |   24|    0|   24|   12|   24|   18|   12|    20|
|          6| B\*      | C        |   12|    0|   12|    6|   12|    6|   12|    10|
|          7| B\*      | B        |   12|    0|   12|    6|   12|   12|   12|    10|
|          8| A\*      | C        |   12|    0|   12|    6|   12|    6|   12|    10|
|          9| A        | A\*      |   24|   12|   36|   18|   24|   24|   18|    26|
|         10| A\*      | C\*      |   12|    0|   12|    6|   12|    6|    6|     6|
|         11| A\*      | B        |   12|    0|   12|    6|   12|   12|   12|    10|

The model we are thinking of is one in which an individual *i* is characterized by a fairness ideal *k*(*i*) and a *strictness factor* *γ*<sub>*i*</sub>, and compares the utility of the two alternatives available. Since all decisions are spectator decisions, there is no monetary self interest, only an internal cost of deviation from the fair distribution. There is also a random element to the utility of each alternative, representing behavioral randomness of the individual, so utility of an allocation *y* to the first individual, when the first individual under fairness ideal *k*(*i*) is entitled to *F*<sup>*k*(*i*)</sup> is
*U*(*y*; ⋅)= − *γ*<sub>*i*</sub>*V*(*y*, *F*<sup>*k*(*i*)</sup>; ⋅)+*ϵ*<sub>*y*</sub>.
 Here the *strictness factor* *γ*<sub>*i*</sub> express the weight on fairness relative to behavioral noise, and *V* is a loss-function corresponding to unfair choices. The limit case with *γ*<sub>*i*</sub> = 0 means that choice probabilities are always uniform (50-50), while as *γ* → ∞, choices converge to always being in line with the deterministic prediction of the fairness ideal. Utility is normalized by the assumption that the *ϵ*<sub>*y*</sub> are IID Extreme Value distributed (of type 1, such that *v**a**r*(*ϵ*<sub>*y*</sub>)=*π*<sup>2</sup>/6). Given that it is only the rate of substitution between the loss *V* and the noise *ϵ* that matters, normalizing the variance of *ϵ* or the scale of *γ* is necessary.

Different choises of the loss functions imply different behavioral models. Of special interest is a quadratic loss (such as previously estimated by Cappelen et al). (But note that a \`\`tremble'' model can be implemented by letting *V* be an indicator function, (*y* ≠ *F*<sup>*k*(*i*)</sup>). With the limited variation in possible deviations from fairness in the data from the luck experiment, there won't be much of a difference between quadratic- and indicator loss in practice.)

The quadratic loss formulation corresponding to that in Cappelen et al (AER, 2013) is used in the current estimates:

``` r
loss <- function(y, F, X) {
  (y-F)^2 / X
}
```

Datastructures
==============

In order to calculate the likelihood, data should be stored in a way that is easy to work with, and the default flat table output from Stata is not the easiest form to work with. Instead I construct new datatypes to organize the data: One to hold data from a single choice situation, and one to hold all information about an individual (including a list of the choice situations).

Situations
----------

To represent a situation, we keep track of the index of the situation type (for debugging purposes), the earnings of each individual (and the total), the income implemented for the first individual in the pair, a vector of what the ideals imply, and a vector of the alternatives available (for income to the first individual in the pair). Since all ideals are efficient, it suffices to keep track of the income allocated to one of the individuals, since the income of the other follows residually.

``` r
setClass("situation", representation(t="numeric",
                                     e1="numeric",
                                     e2="numeric",
                                     X="numeric",
                                     y1="numeric",
                                     ideals="numeric",
                                     Yset="numeric")
         )
```

There is a method to initialize an individual situation:

``` r
setMethod("initialize",
          signature(.Object = "situation"),
          function (.Object, type, e1, e2, y1, f1,f2,f3,f4,f5) {
            .Object@t  <- type
            .Object@e1 <- e1
            .Object@e2 <- e2
            .Object@X <- e1+e2
            .Object@y1 <- y1
            .Object@ideals <- c(f1,f2,f3,f4,f5)
            .Object@Yset <- c(e1, (e1+e2)/2)
            .Object
          }
          )
```

The helper function creates a situation from the natural arguments: The technology, claims for each individal, the endowment available, and amount allocated correspondingly to the two individuals.

``` r
situation <- function(t, e1,e2,y1, f1,f2,f3,f4,f5) {
    new("situation", t, e1,e2,y1, f1,f2,f3,f4,f5)
}
```

Individuals
-----------

Individuals hold all information collected about subjects, and information about the choices are kept as a list of situations.

``` r
setClass("individual", representation(pid="numeric",
                                      session="numeric",
                                      round="numeric",
                                      sex="numeric",
                                      spectator="numeric",
                                      insurance="numeric",
                                      risk="numeric",
                                      situations="list")
)
```

The initialize method takes named arguments with answers to survey questions and initializes to an empty list of situations.

``` r
setMethod("initialize",
          signature(.Object = "individual"),
          function (.Object, pid, session, round, sex, spectator, insurance, risk ) {
              .Object@pid <- pid
              .Object@session <- session
              .Object@round <- round
              .Object@sex <- sex
              .Object@spectator <- spectator
              .Object@insurance <- insurance
              .Object@risk <- risk
              .Object@situations <- list()
              .Object
          }
          )
```

The addsit method (a replace method) adds situations to an individual.

``` r
setGeneric("addsit<-", function(this, value) standardGeneric("addsit<-"))
setReplaceMethod("addsit", signature=signature("individual","situation"),
                 function(this, value) {
                     nsit <- length(this@situations)
                     this@situations[[nsit+1]] <- value
                     this
                 }
                 )
```

The helper function constructs an individual based on its personal id and answers to they survey questions:

``` r
individual <- function(pid,session,round,sex,spectator,insurance,risk) {
     new("individual",pid,session,round,sex,spectator,insurance,risk)
}
```

The helper function `n` on an individual returns the number of situations

``` r
setGeneric("n", function(s) standardGeneric("n"))
setMethod("n", "individual",
          function(s) {
              length(s@situations)
          }
          )
```

Constructing data as list of objects
------------------------------------

The total amount of data is held as a list of individual-objects with all the situations added in. The `datastruct` function takes two dataframes as arguments: One of players, and one of situstions (linked by having the common identifier `pid`).

``` r
datastruct <- function(pls, sits) {
    finalstructure <- list()
    for(i in pls$pid) {
        ind <- individual(i,
                          pls[pls$pid==i,]$Session,
                          pls[pls$pid==i,]$round,
                          pls[pls$pid==i,]$sex,
                          pls[pls$pid==i,]$spectator,
                          pls[pls$pid==i,]$insurance,
                          pls[pls$pid==i,]$risk)
        df <- subset(sits,pid==i)
        for (j in 1:length(df$pid)) {
            sit <- situation(df[j,]$situation,
                             df[j,]$e1,
                             df[j,]$e2,
                             df[j,]$inc1,
                             df[j,]$Fse,
                             df[j,]$Fl,
                             df[j,]$Fle,
                             df[j,]$Fcc,
                             df[j,]$Fscc)
            addsit(ind) <- sit
        }
        nind <- length(finalstructure)
        finalstructure[[nind+1]] <- ind
    }
    finalstructure
}
```

It should now work to load all data into a meaningful datastructure.

The likelihood function
=======================

Given a gamma and an ideal, it is straightforward to calculate choice probabilities. The generic function `pk.given.gamma` does this:

``` r
setGeneric("pk.given.gamma", function(s, k, gamma) standardGeneric("pk.given.gamma"))
```

There are two methods that implements this. On an individual situation, the choices are simply logit kernels.

``` r
setMethod("pk.given.gamma", 
          signature=signature(s="situation",k="numeric", gamma="numeric"), 
          function(s, k, gamma) {
            ds <- - gamma*loss(s@Yset, s@ideals[k], s@X)
            mds = max(ds)
            expV <-  exp(ds - mds)
            if (length(s@Yset)>1) {
              exp( - gamma*loss(s@y1, s@ideals[k], s@X) - mds) / sum(expV)
            } else {
              1.0
            }
          }
)
```

Note that the maximum is subtracted -- algebraically, this makes no difference, but it enhances numerical robustness.

There is a method for an individual which just multiplies the choice probabilities for all the situations an individual is involved in:

``` r
setMethod("pk.given.gamma", 
          signature=signature(s="individual", k="numeric", gamma="numeric"),
          function(s, k, gamma) {
              p <- 1.0
              for (j in 1:length(s@situations)) {
                  p <- p*pk.given.gamma(s@situations[[j]], k, gamma)
              }
              p
          }
          )
```

Conditional on ideal *k*, the likelihood takes the form
$$
L\_{ik} = \\int\_0^\\infty \\prod\_{s}\\left\[ \\frac{e^{V(y\_s;\\gamma, \\cdot)}}{\\sum\_{y\_j\\in Y\_s}
    e^{V(y\_j;\\gamma,\\cdot)}}\\right\] \\frac{1}{\\gamma \\sigma \\sqrt{2\\pi}}
   e^{- \\frac{1}{2}\\left(\\frac{\\log \\gamma - \\mu}{\\sigma}\\right)^2 }d \\gamma,
$$
 where the product is the choice probabilities given a *γ*, and the second expression is the log normal distribution assumed for *γ*. The integration is over the domain of *γ*. We want to evaluate this integral with Gauss-Hermite integration. For notation, let
$$
f(\\gamma) = \\prod\_{s}\\left\[ \\frac{e^{V(y\_s;\\gamma, \\cdot)}}{\\sum\_{y\_j\\in Y\_s}
    e^{V(y\_j;\\gamma,\\cdot)}}\\right\]. $$

Gauss-Hermite integration provides a way to choose nodes (*s*<sub>*i*</sub>)<sub>*i*</sub> and weights (*w*<sub>*i*</sub>)<sub>*i*</sub> such that
∫<sub>−∞</sub><sup>∞</sup>*f*(*s*)*e*<sup>−*s*<sup>2</sup></sup>*d**s* ≈ ∑<sub>*i*</sub>*w*<sub>*i*</sub>*f*(*s*<sub>*i*</sub>),
 in a way that is optimal for polynomial *f*s. In order to do so, we need to transform variables such that the kernel is the same,
$$ s^2 = \\frac{1}{2} \\left(\\frac{\\log \\gamma - \\mu}{\\sigma}\\right)^2.$$
 This implies that $s=(\\log \\gamma - \\mu)/\\sigma\\sqrt{2}$, and that $d\\gamma = \\gamma \\sqrt{2}\\sigma d s$, and the integral can be written
$$ L\_{ik} = \\frac{1}{\\sqrt{\\pi}} \\int\_0^\\infty f( \\exp(\\mu +
\\sqrt{2}\\sigma s)) e^{-s^2}\\approx \\frac{1}{\\sqrt{\\pi}}\\sum\_t w\_t f(\\exp(\\mu
+ \\sqrt{2}\\sigma s\_t)).$$
 The package `glmmML` provides a function to generate these nodes and weights. We set a global named list of vectors with nine nodes, which gives us the classical bell shape:

``` r
gh <- glmmML::ghq(9, FALSE)
```

We can now define a function that performs this calculation using the Gauss-Hermite nodes:

``` r
setGeneric("pk.given.theta",
           function(ind, k, mu, sigma) standardGeneric("pk.given.theta"))
setMethod("pk.given.theta", signature=signature(ind="individual", k="numeric",
                            mu="numeric", sigma="numeric"),
          function(ind, k, mu, sigma) {
              gammas <- exp( mu + sigma*sqrt(2)*gh$zeros)
              values <- rep(0, times=length(gammas))
              for (i in 1:length(gammas)) {
                  values[i] <- gh$weights[i]*pk.given.gamma(ind, k, gammas[i])
              }
              sum(values)/sqrt(pi)
          }
          )
```

In order to generate the probabilities that an individual has a given type, we define the function `p.of.k` that takes an individual and two parameters as arguments. The two arguments are transformed into the five-dimensional unit simplex in the standard way. This is to ensure that the likelihood maximization can proceed as if maximizing over **R**<sup>*D*</sup>, while in fact we are using functions on the inside of the likelihood function interface to transform these values to the parameter space *Θ*.

``` r
setGeneric("p.of.k", function(b1.c, b2.c, b3.c, b4.c, b5.c) standardGeneric("p.of.k"))
setMethod("p.of.k", signature=signature(b1.c="numeric", 
                                        b2.c="numeric", 
                                        b3.c="numeric", 
                                        b4.c="numeric"),
          function(b1.c, b2.c, b3.c, b4.c, b5.c) {
              t1 <- exp(b1.c)
              t2 <- exp(b2.c)
              t3 <- exp(b3.c)
              t4 <- exp(b4.c)
              t5 <- exp(b5.c)
              c( t1, t2, t3, t4, t5)/(t1+t2+t3+t4+t5)
          }
          )
```

The log likelihood of an individual is the type-probability weighted sum of the integral over choice probabilities given *k*,

``` r
setGeneric("log.li", function(ind, b1.c, b2.c, b3.c, b4.c, b5.c,
                              mu, sigma) standardGeneric("log.li"))
setMethod("log.li", signature(ind="individual", 
                              b1.c="numeric", 
                              b2.c="numeric", 
                              b3.c="numeric",
                              b4.c="numeric", 
                              b5.c="numeric",
                              mu="numeric", 
                              sigma="numeric"),
          function(ind, b1.c, b2.c, b3.c, b4.c, b5.c, mu, sigma) {
              lambdas <- p.of.k(b1.c, b2.c, b3.c, b4.c, b5.c)
              l<- 0.0
              for (k in 1:5) {
                  l <- l + lambdas[k]*pk.given.theta(ind,k,mu,sigma)
              }
              log(l)
          }
          )
```

The function to calculate log likelihood in fact searches for values for log*σ*, transformed back in this chunk:

``` r
llvector.luck  <- function(ds, pvector) {
    b1.c <- pvector[1]
    b2.c <- pvector[2]
    b3.c <- pvector[3]
    b4.c <- pvector[4]
    b5.c <- pvector[5]
    mu <- pvector[6]
    sigma <- exp(pvector[7])
    lls <- foreach(d=ds, .combine="c") %dopar% -log.li(d,b1.c,b2.c, b3.c, b4.c, b5.c, mu, sigma)
    lls
}
ll.luck  <- function(ds, pvector) {
    sum(llvector.luck(ds, pvector))
}
```

Variant likelihood functions with 2K parameters.
------------------------------------------------

For testing if a subset of individuals have different parameters, it is useful to be able to do likelihood ratio tests of arbitrary sets of variables. Requested in the revision, is a test of *σ* equal across those who did and did not buy insurance.

Such variants are not automatic with the code above, since that doesn't distinguish people by demographic group. Instead we want to create a variant that distinguishes between two groups of individuals based on demographics, and each of these can get (potentially) different parameter vectors. I create a variant of `ll.luck` that takes an `fnind` argument which is a mapping from an individual to {1, 2}, and these can get variant parameter vectors. For practical reasons, it is probably useful to make the two parameter vectors to be differences from each other, since then it is easy to impose that the vectors are the same along some dimensions.

``` r
setGeneric("log.li.2K", function(ind, b1.c.1, b2.c.1, b3.c.1, b4.c.1, b5.c.1, mu.1, sigma.1,
                                 b1.c.2, b2.c.2, b3.c.2, b4.c.2, b5.c.2, mu.2, sigma.2,
                                 fnind) standardGeneric("log.li.2K"))
```

    ## [1] "log.li.2K"

``` r
setMethod("log.li.2K", signature(ind="individual", 
                              b1.c.1="numeric", 
                              b2.c.1="numeric", 
                              b3.c.1="numeric",
                              b4.c.1="numeric", 
                              b5.c.1="numeric",
                              mu.1="numeric", 
                              sigma.1="numeric",
                              b1.c.2="numeric", 
                              b2.c.2="numeric", 
                              b3.c.2="numeric",
                              b4.c.2="numeric", 
                              b5.c.2="numeric",
                              mu.2="numeric", 
                              sigma.2="numeric",
                              fnind="standardGeneric"
                              ),
          function(ind, b1.c.1, b2.c.1, b3.c.1, b4.c.1, b5.c.1, mu.1, sigma.1,
                                 b1.c.2, b2.c.2, b3.c.2, b4.c.2, b5.c.2, mu.2, sigma.2,
                                 fnind) {
            if (fnind(ind)==1) {
              print("Type 1")
              lambdas <- p.of.k(b1.c.1, b2.c.1, b3.c.1, b4.c.1, b5.c.1)
              mu <- mu.1
              sigma <- sigma.1
            } else {
              print("Type 2")
              lambdas <- p.of.k(b1.c.2, b2.c.2, b3.c.2, b4.c.2, b5.c.2)
              mu <- mu.2
              sigma <- sigma.2          
            }
            l<- 0.0
            for (k in 1:5) {
              l <- l + lambdas[k]*pk.given.theta(ind,k,mu,sigma)
            }
            log(l)
          }
)
```

``` r
llvector.luck.2K <- function(ds, pvector, fnind) {
  b1.c.1 <- pvector[1]
  b2.c.1 <- pvector[2]
  b3.c.1 <- pvector[3]
  b4.c.1 <- pvector[4]
  b5.c.1 <- pvector[5]
  b1.c.2 <- pvector[1] + pvector[8]
  b2.c.2 <- pvector[2] + pvector[9]
  b3.c.2 <- pvector[3] + pvector[10]
  b4.c.2 <- pvector[4] + pvector[11]
  b5.c.2 <- pvector[5] + pvector[12]
  mu.1 <- pvector[6]
  mu.2 <- pvector[13]
  sigma.1 <- exp(pvector[7])
  sigma.2 <- exp(pvector[7] + pvector[14])
  lls <- foreach(d=ds, .combine="c") %dopar% -log.li.2K(d, b1.c.1, b2.c.1, b3.c.1, b4.c.1, b5.c.1, mu.1, sigma.1,
                                                        b1.c.2, b2.c.2, b3.c.2, b4.c.2, b5.c.2, mu.2, sigma.2,
                                                        fnind)
  lls
}


ll.luck.2K <- function(ds, pvector, fnind) {
  sum(llvector.luck.2K(ds, pvector, fnind))
}
```

I don't have a nice and tidy delta-transformed presentation view of this: All we are interested in is basically estimating the constrained model for likelihood comparisons.

Standard errors
===============

The `stats4` library by default attempts to take a numerical Hessian of the likelihood function. This does not seem to work well. An alternative is to use the outer product of the gradient Bernd74. The `numDeriv` package handles calculating the gradient.

The first function takes two parameter vectors: One long (x) and one short (y) and generate an index of where the names in the short vector appear in the names of the long vector. This is to pick out those variables that were actually actively optimized over.

``` r
vcovindices <- function(x,y) {
  z <- c()
  for (s in y) {
    z <- c( z, which(x==s)[1])
  }
  z
}
```

Now, for calculation of the BHHH covariance matrix with correct names.

``` r
vcov.bhhh <- function(ds, fit) {
  long_pv <- fit@fullcoef
  short_pv <- fit@coef
  f <- function(x) {
    - llvector.luck(ds, x)
  }
  g_long <- numDeriv::jacobian(f, long_pv)
  a_ind <- vcovindices(names(long_pv), names(short_pv))
  g <- g_long[, a_ind]
  vcv <- solve( t(g) %*% g)
  colnames(vcv) <- names(short_pv)
  rownames(vcv) <- names(short_pv)
  vcv
}
#VVV<-vcov.bhhh(ds.luck1, fit.all)
```

But what we really want is the covariance matrix of the transformed parameters. Let the `gtrans` function transform the parameters as they are free for estimation to the natural parameters.

``` r
grtrans <- function(x) {
    t1 <- exp(x[1])
    t2 <- exp(x[2])
    t3 <- exp(x[3])
    t4 <- exp(x[4])
    t5 <- exp(x[5])
    
    tsum <- t1+t2+t3+t4+t5
    se <- t1/tsum
    l  <- t2/tsum
    le <- t3/tsum
    cc <- t4/tsum
    scc <- t5/tsum
    mu <- x[6]
    sigma <- exp(x[7])
    fv <- c(se, l, le, cc, scc,  mu, sigma)
    names(fv) <- c("share.se", "share.l", "share.le", "share.cc", "share.scc", "mu", "sigma")
    fv
}
```

We calculate standard errors of the transformed parameters using the Delta method (approximating the function by a first order Taylor approximation).

``` r
delta.method <- function(func, vcov, x) {
    g <- numDeriv::jacobian(func, x)
    snames <- colnames(vcov)
    colnames(g) <- names(x)
    gm <- g[,snames]
    new.vcov <- gm %*% vcov %*% t(gm)
    std.errors <- sqrt(diag(new.vcov))
    params <- grtrans(x)
    out <- cbind(params,std.errors)
    colnames(out) <- c("Estimate", "Std. error")
    rownames(out) <- names(func(x))
    out
}
```

Actually calculating the standard errors of the overall would now be:

    delta.method(grtrans, vcov(fit.all), coef(fit.all))
