
## install packages:

require(ggplot2)  ## for plotting
require(mgcv) ## for GAM
require(np)  ## for kernel estimation
require(foreach)  # for parallel computing
require(doParallel) # for parallel computing
require(Lmoments) # for L-kurtosis
require(lmtest) # for wald test

require(sandwich)
require(pcse)

## ggplot2
## mgcv (for GAM)
## np, foreach, doParallel

## simulated data

rm(list=ls())
source("hmx.R")


###simulated example of linear cont. interaction

set.seed(1234)
n<-200
d1<-sample(c(0,1),n,replace=TRUE)
d2<-rnorm(n,3,1)
x<-rnorm(n,3,1)
z<-rnorm(n,3,1)
e<-rnorm(n,0,1)
y1<-5 - 4 * x - 9 * d1 + 3 * x * d1 + 1 * z + 2 * e
y2<-5 - 4 * x - 9 * d2 + 3 * x * d2 + 1 * z + 2 * e
s1<-cbind.data.frame(Y = y1, D = d1, X = x, Z1 = z)
s2<-cbind.data.frame(Y = y2, D = d2, X = x, Z1 = z)
head(s1)
head(s2)

## Raw Plots
inter.rawplot(Y = "Y", D = "D", X = "X", data = s1, weights = NULL,
              Ylabel = "Outcome", Dlabel = "Treatment", Xlabel="Moderator")

inter.rawplot(Y = "Y", D = "D", X = "X", data = s2,
              Ylabel = "Outcome", Dlabel = "Treatment", Xlabel="Moderator")

inter.gam(Y="Y", D="D", X="X", data=s2)

## Binning Estimates

inter.binning(Y = "Y", D = "D", X = "X", Z = "Z1", data = s1, nbins = 3,
              Ylabel = "Outcome", Dlabel = "Treatment", Xlabel="Moderator")

inter.binning(Y = "Y", D = "D", X = "X", Z = "Z1", data = s2, 
              Ylabel = "Outcome", Dlabel = "Treatment", Xlabel="Moderator")

inter.binning(Y = "Y", D = "D", X = "X", Z = "Z1", data = s2, 
              Ylabel = "Outcome", Dlabel = "Treatment", Xlabel="Moderator",
              cuts = c(0,2,4,7))

## Kernel Estimates

inter.kernel(Y = "Y", D = "D", X = "X", Z="Z1", data=s1, nboot=500,
             cl=NULL, parallel=TRUE, cores = 4)

inter.kernel(Y = "Y", D = "D", X = "X", Z="Z1", data=s2, nboot=500, h=1.3)


## Wald test

inter.wald(Y = "Y", D = "D", X = "X", Z = "Z1", data=s1,
           vartype = "homoscedastic", cl = NULL,
           weights = NULL, nbins=3, cuts = NULL)

inter.wald(Y = "Y", D = "D", X = "X", Z = "Z1", data=s2)










