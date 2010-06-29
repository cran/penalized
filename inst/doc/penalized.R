###################################################
### chunk number 1: options
###################################################
options(continue = "  ")


###################################################
### chunk number 2: load
###################################################
library(penalized)
library(survival)
data(nki70)


###################################################
### chunk number 3: setseed
###################################################
set.seed(1)


###################################################
### chunk number 4: first
###################################################
fit <- penalized(ER, ~DIAPH3+NUSAP1, data=nki70, lambda2=1)
fit <- penalized(ER, nki70[,10:11], data=nki70, lambda2=1)
fit <- penalized(ER~DIAPH3+NUSAP1, data=nki70, lambda2=1)


###################################################
### chunk number 5: survival
###################################################
fit <- penalized(Surv(time,event)~DIAPH3+NUSAP1, data=nki70, lambda2=1)


###################################################
### chunk number 6: attach
###################################################
attach(nki70)


###################################################
### chunk number 7: extract
###################################################
residuals(fit)[1:10]
fitted(fit)[1:10]
basesurv(fit)


###################################################
### chunk number 8: coefficients
###################################################
coefficients(fit, "all")


###################################################
### chunk number 9: loglik_penalty
###################################################
loglik(fit)
penalty(fit)


###################################################
### chunk number 10: predict
###################################################
predict(fit, ~DIAPH3+NUSAP1, data=nki70[1:3,])
predict(fit, nki70[1:3,c("DIAPH3","NUSAP1")])


###################################################
### chunk number 11: predict_survival
###################################################
pred <- predict(fit, nki70[1:3,c("DIAPH3","NUSAP1")])
survival(pred, time=5)


###################################################
### chunk number 12: weights
###################################################
coefficients(fit)
coefficients(fit, standardize = TRUE)
weights(fit)


###################################################
### chunk number 13: unpenalized
###################################################
fit <- penalized(Surv(time,event), nki70[,8:77], ~ER, lambda2=1)
fit <- penalized(Surv(time,event)~ER, nki70[,8:77], lambda2=1)


###################################################
### chunk number 14: steps
###################################################
fit <- penalized(Surv(time,event), nki70[,8:77], lambda1=1,
    steps=50, trace = FALSE)
plotpath(fit, log="x")


###################################################
### chunk number 15: stepsplot
###################################################
plotpath(fit, log="x")


###################################################
### chunk number 16: positive
###################################################
fit <- penalized(Surv(time,event), nki70[,8:77], positive=TRUE)


###################################################
### chunk number 17: positivecoefficients
###################################################
coefficients(fit)


###################################################
### chunk number 18: positivestepsplot
###################################################
fit <- penalized(Surv(time,event), nki70[,8:77], positive=TRUE,
    steps=50)
plotpath(fit)


###################################################
### chunk number 19: globaltest_install eval=FALSE
###################################################
## source("http://bioconductor.org/biocLite.R")
## biocLite("globaltest")


###################################################
### chunk number 20: globaltest
###################################################
library(globaltest)
gt(Surv(time,event), nki70[,8:77])


###################################################
### chunk number 21: cvl1
###################################################
fit <- cvl(Surv(time,event), nki70[,8:77], lambda1=1, fold=10)


###################################################
### chunk number 22: cvl2
###################################################
fit$cvl
fit$fullfit


###################################################
### chunk number 23: cvl3
###################################################
fit <- cvl(Surv(time,event), nki70[,8:77], lambda1=2, fold=fit$fold)


###################################################
### chunk number 24: breslow
###################################################
fit$predictions
time(fit$predictions)
as.data.frame(basesurv(fit$fullfit))[1:10,]
plot(fit$predictions)


###################################################
### chunk number 25: breslowplot
###################################################
plot(fit$predictions)


###################################################
### chunk number 26: cv-survival
###################################################
survival(fit$predictions, 5)[1:10]


###################################################
### chunk number 27: prof
###################################################
fit1 <- profL1(Surv(time,event), nki70[,50:70], fold=10)
plot(fit1$lambda, fit1$cvl, type="l")
fit2 <- profL2(Surv(time,event), nki70[,50:70], fold=fit1$fold,
    minl = 0.01, maxl = 1000)
plot(fit2$lambda, fit2$cvl, type="l", log="x")


###################################################
### chunk number 28: profplot1
###################################################
plot(fit1$lambda, fit1$cvl, type="l")


###################################################
### chunk number 29: profplot2
###################################################
plot(fit2$lambda, fit2$cvl, type="l", log="x")


###################################################
### chunk number 30: profpath
###################################################
plotpath(fit2$fullfit, log="x")


###################################################
### chunk number 31: profpathplot
###################################################
plotpath(fit2$fullfit, log="x")


###################################################
### chunk number 32: opt1
###################################################
opt1 <- optL1(Surv(time,event), nki70[,50:70], fold=fit1$fold)


###################################################
### chunk number 33: optres
###################################################
opt1$lambda
opt1$cvl


###################################################
### chunk number 34: opt2
###################################################
opt2 <- optL2(Surv(time,event), nki70[,50:70], fold=fit2$fold)


