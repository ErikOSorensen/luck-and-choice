## ------------------------------------------------------------------------
require(methods)
library(foreign)
library(stats4)
library(glmmML)
library(numDeriv)
library(doMC)

registerDoMC(8)
getDoParWorkers()

situations <- read.dta("data/DataErik_long.dta", 
                       convert.factors=FALSE)
players <- read.dta("data/DataErik_players.dta", 
                    convert.factors=FALSE)

## ----echo=FALSE----------------------------------------------------------
ideals <- read.dta("data/situations.dta")
knitr::kable(ideals, caption="Ideals (allocation in individual 1) in different situations")

## ------------------------------------------------------------------------
loss <- function(y, F, X) {
  (y-F)^2 / X
}

## ----results='hide', tidy=FALSE------------------------------------------
setClass("situation", representation(t="numeric",
                                     e1="numeric",
                                     e2="numeric",
                                     X="numeric",
                                     y1="numeric",
                                     ideals="numeric",
                                     Yset="numeric")
         )

## ----results='hide', tidy=FALSE------------------------------------------
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

## ----situation, results='hide', tidy=FALSE-------------------------------
situation <- function(t, e1,e2,y1, f1,f2,f3,f4,f5) {
    new("situation", t, e1,e2,y1, f1,f2,f3,f4,f5)
}

## ----results='hide', tidy=FALSE------------------------------------------
setClass("individual", representation(pid="numeric",
                                      session="numeric",
                                      round="numeric",
                                      sex="numeric",
                                      spectator="numeric",
                                      insurance="numeric",
                                      risk="numeric",
                                      situations="list")
)

## ----results='hide', tidy=FALSE------------------------------------------
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

## ----results='hide', tidy=FALSE------------------------------------------
setGeneric("addsit<-", function(this, value) standardGeneric("addsit<-"))
setReplaceMethod("addsit", signature=signature("individual","situation"),
                 function(this, value) {
                     nsit <- length(this@situations)
                     this@situations[[nsit+1]] <- value
                     this
                 }
                 )

## ----tidy=FALSE----------------------------------------------------------
individual <- function(pid,session,round,sex,spectator,insurance,risk) {
     new("individual",pid,session,round,sex,spectator,insurance,risk)
}

## ----results='hide', tidy=FALSE------------------------------------------
setGeneric("n", function(s) standardGeneric("n"))
setMethod("n", "individual",
          function(s) {
              length(s@situations)
          }
          )

## ----datastruct, tidy=FALSE----------------------------------------------
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

## ----results='hide'------------------------------------------------------
setGeneric("pk.given.gamma", function(s, k, gamma) standardGeneric("pk.given.gamma"))

## ----results='hide', tidy=FALSE------------------------------------------
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

## ----results='hide', tidy=FALSE------------------------------------------
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

## ----gh------------------------------------------------------------------
gh <- glmmML::ghq(9, FALSE)

## ----results='hide', tidy=FALSE------------------------------------------
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

## ----results='hide', tidy=FALSE------------------------------------------
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

## ----results='hide', tidy=FALSE------------------------------------------
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

## ------------------------------------------------------------------------
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

## ----tidy=FALSE----------------------------------------------------------
setGeneric("log.li.2K", function(ind, b1.c.1, b2.c.1, b3.c.1, b4.c.1, b5.c.1, mu.1, sigma.1,
                                 b1.c.2, b2.c.2, b3.c.2, b4.c.2, b5.c.2, mu.2, sigma.2,
                                 fnind) standardGeneric("log.li.2K"))
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

## ----rtidy=FALSE---------------------------------------------------------
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

## ----tidy=FALSE----------------------------------------------------------
vcovindices <- function(x,y) {
  z <- c()
  for (s in y) {
    z <- c( z, which(x==s)[1])
  }
  z
}

## ------------------------------------------------------------------------
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

## ----results='hide'------------------------------------------------------
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

## ----results='hide'------------------------------------------------------
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

