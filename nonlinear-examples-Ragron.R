## Illustrate the use of nonlinear models
## 
## R agronomists: Why should I use a nonlinear model?
## Author: Fernando E. Miguez
## Date: 2020-05-01

## Load libraries
library(ggplot2)
library(nlraa)
library(nlme)
library(quantreg)

## Orange example
data(Orange)

ggplot(data = Orange, aes(x = age, y = circumference)) + 
  geom_point() + 
  ggtitle("Growth of Orange Trees") + 
  ylab("Circumference (cm)") + 
  xlab("age (days)")

## Polynomial fit
fit.poly <- lm(circumference ~ poly(age, 6), data = Orange)
ndat <- data.frame(age = 0:1700)
ndat$prd <- predict(fit.poly, newdata = ndat)

ggplot() + 
  geom_point(data = Orange, aes(x = age, y = circumference)) + 
  geom_line(data = ndat, aes(x = age, y = prd)) + 
  ggtitle("Growth of Orange Trees \n plus polynomial sixth degree")

## Short simple examples
data(barley, package = "nlraa")

barley2 <- subset(barley, year < 1975)

fit.bar <- nls(yield ~ SSlinp(NF, a, b, xs), data = barley2)

ndat <- data.frame(NF = seq(0, 12, 0.05))
ndat$prd <- predict(fit.bar, newdata = ndat)

ggplot() + 
  geom_point(data = barley2, aes(x = NF, y = yield)) + 
  geom_line(data = ndat,  aes(x = NF, y = prd)) + 
  ylab("Yield (g/m2)") + xlab("Nitrogen fertilizer (g/m2)") + 
  geom_vline(xintercept = coef(fit.bar)[3], linetype = 2) + 
  ggtitle("Nonlinear model: Linear-plateau")

## Fitting a logistic growth model
fit.nl <- nls(circumference ~ SSlogis(age, Asym, xmid, scal), 
              data = Orange)

ggplot(data = Orange, aes(x = age, y = circumference)) + 
  geom_point() + 
  geom_line(aes(y = fitted(fit.nl))) + 
  ggtitle("Growth of Orange Trees \n plus logistic")

## Extending the model
fit.gnls <- gnls(circumference ~ SSlogis(age, Asym, xmid, scal),
                 data = Orange, weights = varFixed(~age))

plot(fit.gnls, abline = 0)

## This code is for the alternative modeling of the variance
## Quantifying the uncertainty in the mean response using bootstrap
# fit.gnls <- gnls(circumference ~ SSlogis(age, Asym, xmid, scal),
#                  data = Orange, weights = varPower())

prd_fun <- function(x) predict(x, 
                               newdata = data.frame(age = 100:1600))

system.time(fit.gnls.bt <- boot_nlme(fit.gnls, prd_fun, R = 1e3, 
                                     parallel = "multicore",
                                     ncpus = 4))

ndat2 <- data.frame(age = 100:1600)
ndat2$prd.q5 <- apply(t(fit.gnls.bt$t), 1, quantile, probs = 0.05, na.rm = TRUE)
ndat2$prd.q50 <- apply(t(fit.gnls.bt$t), 1, quantile, probs = 0.50, na.rm = TRUE)
ndat2$prd.q95 <- apply(t(fit.gnls.bt$t), 1, quantile, probs = 0.95, na.rm = TRUE)

ggplot() + 
  geom_point(data = Orange, aes(x = age, y = circumference)) + 
  geom_line(data = ndat2, aes(x = age, y = prd.q5), color = "blue", linetype = 2) + 
  geom_line(data = ndat2, aes(x = age, y = prd.q50), color = "blue", linetype = 1) +
  geom_line(data = ndat2, aes(x = age, y = prd.q95), color = "blue", linetype = 2) +
  geom_text(aes(x = 1650, y = 150), label = "5%", color = "blue") + 
  geom_text(aes(x = 1650, y = 180), label = "50%", color = "blue") + 
  geom_text(aes(x = 1650, y = 210), label = "95%", color = "blue") + 
  ggtitle("Growth of Orange Trees \n bootstrapped mean function (varFixed)")

## Nonlinear quantile regression
fit.nlrq1 <- nlrq(circumference ~ SSlogis(age, Asym, xmid, scal),
                  data = Orange, tau = 0.05)
fit.nlrq2 <- nlrq(circumference ~ SSlogis(age, Asym, xmid, scal),
                  data = Orange, tau = 0.25)
fit.nlrq3 <- nlrq(circumference ~ SSlogis(age, Asym, xmid, scal),
                  data = Orange, tau = 0.5)
fit.nlrq4 <- nlrq(circumference ~ SSlogis(age, Asym, xmid, scal),
                  data = Orange, tau = 0.75)
fit.nlrq5 <- nlrq(circumference ~ SSlogis(age, Asym, xmid, scal),
                  data = Orange, tau = 0.95)

ndat <- data.frame(age = 100:1600)
ndat$prd.q05 <- predict(fit.nlrq1, newdata = ndat)
ndat$prd.q25 <- predict(fit.nlrq2, newdata = ndat)
ndat$prd.q50 <- predict(fit.nlrq3, newdata = ndat)
ndat$prd.q75 <- predict(fit.nlrq4, newdata = ndat)
ndat$prd.q95 <- predict(fit.nlrq5, newdata = ndat)

ggplot() + 
  geom_point(data = Orange, aes(x = age, y = circumference)) + 
  geom_line(data = ndat, aes(x = age, y = prd.q05), color = "blue", linetype = 3) + 
  geom_line(data = ndat, aes(x = age, y = prd.q25), color = "blue", linetype = 2) +
  geom_line(data = ndat, aes(x = age, y = prd.q50), color = "blue", linetype = 1) +
  geom_line(data = ndat, aes(x = age, y = prd.q75), color = "blue", linetype = 2) +
  geom_line(data = ndat, aes(x = age, y = prd.q95), color = "blue", linetype = 3) +
  geom_text(aes(x = 1650, y = 140), label = "5%", color = "blue") + 
  geom_text(aes(x = 1650, y = 150), label = "25%", color = "blue") + 
  geom_text(aes(x = 1650, y = 180), label = "50%", color = "blue") + 
  geom_text(aes(x = 1650, y = 205), label = "75%", color = "blue") + 
  geom_text(aes(x = 1650, y = 225), label = "95%", color = "blue") + 
  ggtitle("Growth of Orange Trees \n plus quantile nonlinear logistic regression")



