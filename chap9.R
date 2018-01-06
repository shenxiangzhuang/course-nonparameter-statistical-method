# 非参数密度估计和非参数回归简介

# 9.1非参数密度估计
library(KernSmooth)
x <- faithful$waiting
plot(x=c(30, 110), y=c(0, 0.04), type="n", bty="l",
     xlab="waiting time (minute)", ylab="density")
lines(bkde(x, bandwidth = dpik(x)))
lines(locpoly(x, bandwidth = dpik(x)), lty=3)

# 9.2非参数回归
# pic9.8
library(MASS)
par(mfrow=c(2,2))
X <- mcycle[, 1]
Y <- mcycle[, 2]
bw <- list("(a) h=1", "(b) h=2", "(c) h=3", "(d) h=5")
plot(X,Y,main=bw[[1]]);lines(ksmooth(X,Y,"normal", bandwidth = 1))
plot(X,Y,main=bw[[2]]);lines(ksmooth(X,Y,"normal", bandwidth = 2))
plot(X,Y,main=bw[[3]]);lines(ksmooth(X,Y,"normal", bandwidth = 3))
plot(X,Y,main=bw[[4]]);lines(ksmooth(X,Y,"normal", bandwidth = 5))

# pic9.9
par(mfrow=c(2,2))
library(MASS)
attach(mcycle)
X = mcycle[, 1]
Y = mcycle[, 2]
par(mfrow=c(2,2))
plot(accel~times, mcycle, main="Lowess")
lines(lowess(mcycle, f=.1))

# 这里和书上不太一致
a = loess(accel~times, mcycle, span=0.15)
plot(accel~times, mcycle, main="Loess")
lines(times, b$fit)

plot(accel~times, mcycle, main="Friedman's SuperSmoother")
lines(supsmu(X, Y))
plot(accel~times, mcycle, main="Smoothing Spline")
lines(ksmooth(X, Y, "normal", bandwidth = 2))
