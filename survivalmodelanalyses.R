library(ggplot2)
library(dplyr)
library(survival)
library(ciTools)
library(tibble)
library(pscl)
library(splines)
library(lattice)
library(JM)


#### BASIC SURVIVAL ANALYSIS

head(germ.day.data)

str(germ.day.data)

### Estimators of the Survival Function and survival rate plot


KM.fit <- survfit(Surv(Germination.day, germinated.or.not)
                  ~ 1, data = germ.day.data) # Kaplan-Meier Estimator

KM.fit

KM.fit.trays <- sapply(trays, function(x) {
  list(survfit(Surv(germ.day.data[germ.day.data$germ.trays == x,"Germination.day"],
               germ.day.data[germ.day.data$germ.trays == x,"germinated.or.not"]) ~ 1))
})

head(KM.fit.trays)

traynames <- paste(rep(c("AA", "AN", "FS"), c(8,2,6)), toglm$Provenance, 
                   toglm$Moisture.treatment, sep = ".")

plot(KM.fit.trays[[1]], conf.int = FALSE, lwd = 2, 
     col = rainbow(16)[1], xlab = "Time(Days)", 
     ylab = "Probability of not germinating")

sapply(2:15, function(i) {
  lines(KM.fit.trays[[i]], conf.int = FALSE, lwd =2, col = rainbow(16)[i])
})
legend("bottomleft", traynames, col = rainbow(16) , 
       lwd =2 , lty =1, bty = "n" ,cex = 0.45)


### Statistical Tests

## Log-rank test (Nonparametric)

lg.rnk.all <- survdiff(Surv(Germination.day, germinated.or.not) ~ 
                           Species+Provenance+Treatment, data = germ.day.data)
lg.rnk.all

lg.rnk.sp.pr <- survdiff(Surv(Germination.day, germinated.or.not) ~ 
                         Species+Provenance, data = germ.day.data)
lg.rnk.sp.pr

lg.rnk.sp <- survdiff(Surv(Germination.day, germinated.or.not) ~ 
                     Species, data = germ.day.data)
lg.rnk.sp

lg.rnk.aa <- survdiff(Surv(Germination.day, germinated.or.not) ~ 
                        Provenance + Treatment, data = germ.day.data.aa)
lg.rnk.aa

lg.rnk.an <- survdiff(Surv(Germination.day, germinated.or.not) ~ 
                        Treatment, data = germ.day.data.an)
lg.rnk.an

lg.rnk.fs <- survdiff(Surv(Germination.day, germinated.or.not) ~ 
                        Provenance + Treatment, data = germ.day.data.fs)
lg.rnk.fs

## Cox proportional hazards model (Semi-Parametric)

cox.mod.all <- coxph(Surv(Germination.day, germinated.or.not) ~ 
                   Species + Provenance + Treatment, data = germ.day.data)

cox.mod.all

cox.mod.sp.tr <- coxph(Surv(Germination.day, germinated.or.not) ~ 
                         Species + Treatment, data = germ.day.data)
cox.mod.sp.tr

cox.mod.aa <- coxph(Surv(Germination.day, germinated.or.not) ~
                      Provenance + Treatment, data = germ.day.data.aa)
cox.mod.aa

cox.mod.an <- coxph(Surv(Germination.day, germinated.or.not) ~ 
                      Treatment, data = germ.day.data.an)

cox.mod.an

cox.mod.fs <- coxph(Surv(Germination.day, germinated.or.not) ~ 
                      Provenance + Treatment, data = germ.day.data.fs)
cox.mod.fs

## Accelerated Failure Time Model (Parametric)

aft.wei <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                        Treatment, data = germ.day.data) #weibull

aft.wei

aft.ln <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                        Treatment, data = germ.day.data, dist = "lognormal")

summary(aft.ln)

aft.ll <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                        Treatment, data = germ.day.data, dist = "loglogistic")
summary(aft.ll)

aft.e <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                       Treatment, data = germ.day.data, dist = "exponential")

summary(aft.e) 


aft.g <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                       Treatment, data = germ.day.data, dist = "gaussian")

summary(aft.g)


AIC(aft.e)
AIC(aft.ll)
AIC(aft.ln) # log normal is the lowest one so go on with this.
AIC(aft.wei)
AIC(aft.g)

aft.ln.aa <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                           Treatment, data = germ.day.data.aa, dist = "lognormal")

summary(aft.ln.aa)

aft.ln.an <- survreg(Surv(Germination.day, germinated.or.not) ~
                           Treatment, data = germ.day.data.an, dist = "lognormal")

summary(aft.ln.an)


aft.ln.fs <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                           Treatment, data = germ.day.data.fs, dist = "lognormal")

summary(aft.ln.fs)

# Checking Residuals

fitted.values <- aft.ln$linear.predictors

resids <- (log(aft.ln$y[, 1]) - fitted.values) / aft.ln$scale

res.KM <- survfit(Surv(resids, germinated.or.not) ~ 1, data = germ.day.data)

plot(res.KM, mark.time = FALSE, col="black", xlab = "AFT Residuals", 
     ylab = "Survival Probability")
xx <- seq(min(resids), max(resids), length.out = 35)
yy <- pnorm(xx, lower.tail = FALSE)
lines(xx, yy, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                       "Survival function of Extreme Value distribution"), 
       lty = c(1,2,1), col = c("black","black","red"), bty = "n")


### flexsurv

library(flexsurv)

flexsurvreg(Surv(Germination.day, germinated.or.not) ~ Species +
          Treatment, data = germ.day.data, dist = "lognormal")

flexsurvreg(Surv(Germination.day, germinated.or.not) ~
          Treatment, data = germ.day.data.an, dist = "lognormal")

?flexsurvreg

