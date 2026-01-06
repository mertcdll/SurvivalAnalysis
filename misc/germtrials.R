### Get overview of the germination data

germ <- read.csv("Germination Data/germdata.csv", header = TRUE)

germ.wo.seedloc <- germ[,!grepl("Seed.", colnames(germ))]

germ.wo.seedloc <- germ.wo.seedloc[,-14]

overview <- germ.wo.seedloc[germ.wo.seedloc[,"Position.in.chamber"] == "finish", ]

overview[,"Provenance"] = as.factor(overview[,"Provenance"])

overview[,"Moisture.treatment"] = as.factor(overview[,"Moisture.treatment"])

summary(overview)

### Seed Increase per Interval

seedmat <- matrix(germ.wo.seedloc[germ.wo.seedloc[,2] != "absolut","Summe"],
           nrow=16, ncol=23, byrow=TRUE)

seedmat <- as.data.frame(seedmat)

colnames(seedmat) <- as.Date(germ.wo.seedloc[,"Date"][1:23], tryFormats = "%m/%d/%Y")

head(seedmat)

metadata <- aggregate(Summe ~ Tray.number + Species + Country +
                        Provenance , data=germ.wo.seedloc, mean)[-5]

metadata <- metadata[order(metadata[,1]),]

seed.count.days <- cbind(metadata, seedmat)

library(tibble)

seed.count.days  <- add_column(seed.count.days, 
                                      Treatment = rep(c("Dry", "Moist"), times = 8), .after = 4)

seed.count.days$Final.percentage <- seed.count.days[,27] / 100

seed.increase <- matrix(data=NA, nrow = 16, ncol=22)
seed.increase <- sapply(1:16, function(j) {
  sapply(1:22, function (i){
    seed.increase[j,i] <- as.numeric(c(seed.count.days[,6:28][j,][i+1] - 
                                         seed.count.days[,6:28][j,][i]))
  })
})


row.names(seed.increase) <- paste0("Int",c(1:22))

colnames(seed.increase) <- seed.count.days[,1]

library(lubridate)

dates <- as.Date(colnames(seed.count.days)[6:28])

intervals <- sapply(1:22, function(i) {
  intervals <- list(interval(dates[i], dates[i+1]))
})


days <- sapply(1:22, function(x) {
  c <- trunc(time_length(intervals[[x]], "day"))
})

days

incrate.perday <- seed.increase / days

seed.increase.intervals <- cbind(seed.count.days[,1:5], t(seed.increase))

seed.increase.rate.intervals <- cbind(seed.count.days[,1:5], t(incrate.perday))

### Plots of Seed Increase along intervals

matplot(incrate.perday, type ="l", 
     xlab = "intervals", ylab ="Germination rate per day")

matplot(seed.increase, type ="l", 
        xlab = "intervals", ylab ="Increase of germinated seeds")

matplot(as.Date(colnames(seed.count.days)[6:28]), 
        t(seed.count.days[,6:28]), type="l" , xlab = "Date",
        ylab = "Number of Seeds Germinated")


par(mfrow=c(4,4))
sapply(1:16, function (i) {
  barplot(seed.increase[,i], main = seed.count.days[i,1] , breaks = 5)
}) #Increase per interval


par(mfrow=c(4,4))
sapply(1:16, function (i) {
  hist(seed.increase[,i], main = seed.count.days[i,1] , breaks = 5)
}) #Frequency of increase of number of germinated seeds


sapply(1:16, function(i){
  shapiro.test(seed.increase[,i])
}) # checking for normality


par(mfrow=c(4,4))
sapply(1:16, function (i) {
  hist(incrate.perday[,i], main = seed.count.days[i,1] , breaks = 5)
})  #Frequency of increase rate

sapply(1:16, function(i){
  shapiro.test(incrate.perday[,i])
}) # checking for normality

### Mean Increase Rates

apply(incrate.perday, 2, mean) # mean increase rate for each tray 
                               #(calculated from increase rates per interval)

duration <- interval(colnames(seed.count.days)[6], 
                     colnames(seed.count.days)[28])

ndays <- trunc(time_length(duration, "day"))

seed.count.days[,28] / ndays # mean increase rate for each tray 
                             # calculated directly.

### Preparing the interval increase data for linear model

Counts = as.vector(seed.increase)

metadata$Treatment <- rep(c("Dry", "Moist"), times = 8)

md <- sapply(1:5, function(i) {
  md <- cbind(rep(metadata[,i], each = 22))
})

md <- as.data.frame(md)

colnames(md) <- colnames(metadata)

md[,"Species"] = as.factor(md[,"Species"])
md[,"Treatment"] = as.factor(md[,"Treatment"])
md[,"Provenance"] = as.factor(md[,"Provenance"])

seedcounts <- cbind(md, Counts)

head(seedcounts)

increase.rate <- as.vector(incrate.perday)


seedcounts$Increase.rate <- increase.rate

head(seedcounts)

### Linear models for seed increase per interval data

library(lme4)


linearmodel <- lme4::lmer(Counts ~ 1 + (1|Provenance) + Species + Treatment, 
                   data = seedcounts, REML = TRUE)

summary(linearmodel)


glinearmodelpoisson <- lme4::glmer(Counts ~ 1 + Provenance * Treatment, 
                            data = seedcounts, family = "poisson")

glinearmodelpoisson

linearpois.rate <- lme4::lmer(Increase.rate ~ 1 + (1|Provenance) + Species + Treatment, 
                                data = seedcounts,REML = TRUE)

glinearqpois.rate <- lme4::glmer(Increase.rate ~ 1 + (1|Provenance) + Species + Treatment, 
                              data = seedcounts, family = "quasipoisson")


glinearmodelbinom <- lme4::glmer(Counts ~ 1 + (1|Provenance) + Species + Treatment, 
                                   data = seedcounts, family = binomial(link = "logit"))



summary(glm(Counts ~ Provenance*Treatment, family = "poisson", data=seedcounts, subset = Species=="Abies alba"))
summary(glm(Counts ~ Treatment, data=seedcounts,family = "poisson", subset = Species=="Abies nordmanniana"))
summary(lm(Counts ~ Provenance*Treatment, family = "poisson",data=seedcounts, subset = Species=="Fagus sylvatica"))

### Transforming the data for survival analysis

germ.wo.absolut <- germ[germ[,2] != "absolut",]

tray.day.seed <- germ.wo.absolut[,c(1,11,13:114)]

head(tray.day.seed)
dim(tray.day.seed)

library(dplyr)


### for all trays

germination.day <- c()

names <- sapply(1:16, function(i) {
  as.name(trays[i])
})

names.in <- sapply(1:16, function(i) {
  as.name(paste0(trays[i], ".ind"))
})

names.day <- sapply(1:16, function(i) {
  as.name(paste0(trays[i], "day"))
})

for (i in 1:16) {
  
  
  names[[i]] <- tray.day.seed %>% filter(Tray.number==trays[i])
  
  names.in[[i]] <- sapply(3:102, function(j){
    if (sum(names[[i]][,j]) >=1) {
      print(min(which(names[[i]][,j]==1)))
    } else {
      print(NA)
    }
  })
  names.day[[i]] <- names[[i]][names.in[[i]],2]
  germination.day <- c(germination.day, names.day[[i]])
}


### for each tray 

A1 <- tray.day.seed %>% filter(Tray.number=="A1")

a1<-sapply(3:102, function(j){
    if (sum(A1[,j]) >=1) {
      print(min(which(A1[,j]==1)))
    } else {
      print(NA)
    }
})

days.a1 <- A1[a1,2]

days.a1

###

A2 <- tray.day.seed %>% filter(Tray.number=="A2")

a2<-sapply(3:102, function(j){
  if (sum(A2[,j]) >=1) {
    print(min(which(A2[,j]==1)))
  } else {
    print(NA)
  }
})

days.a2 <- A2[a2,2]

days.a2

####

B1 <- tray.day.seed %>% filter(Tray.number=="B1")

b1<-sapply(3:102, function(j){
  if (sum(B1[,j]) >=1) {
    print(min(which(B1[,j]==1)))
  } else {
    print(NA)
  }
})

days.b1 <- B1[b1,2]

days.b1


###


B2 <- tray.day.seed %>% filter(Tray.number=="B2")

b2<-sapply(3:102, function(j){
  if (sum(B2[,j]) >=1) {
    print(min(which(B2[,j]==1)))
  } else {
    print(NA)
  }
})

days.b2 <- B2[b2,2]

days.b2

###

C1 <- tray.day.seed %>% filter(Tray.number=="C1")

c1 <-sapply(3:102, function(j){
  if (sum(C1[,j]) >=1) {
    print(min(which(C1[,j]==1)))
  } else {
    print(NA)
  }
})

days.c1 <- C1[c1,2]

days.c1

####

C2 <- tray.day.seed %>% filter(Tray.number=="C2")

c2 <-sapply(3:102, function(j){
  if (sum(C2 [,j]) >=1) {
    print(min(which(C2[,j]==1)))
  } else {
    print(NA)
  }
})

days.c2 <- C2[c2,2]

days.c2


####

D1 <- tray.day.seed %>% filter(Tray.number=="D1")

d1 <-sapply(3:102, function(j){
  if (sum(D1 [,j]) >=1) {
    print(min(which(D1[,j]==1)))
  } else {
    print(NA)
  }
})

days.d1 <- D1[d1,2]

days.d1


####


D2 <- tray.day.seed %>% filter(Tray.number=="D2")

d2 <-sapply(3:102, function(j){
  if (sum(D2[,j]) >=1) {
    print(min(which(D2[,j]==1)))
  } else {
    print(NA)
  }
})

days.d2 <- D2[d2,2]

days.d2


####

E1 <- tray.day.seed %>% filter(Tray.number=="E1")

e1 <-sapply(3:102, function(j){
  if (sum(E1[,j]) >=1) {
    print(min(which(E1[,j]==1)))
  } else {
    print(NA)
  }
})

days.e1 <- E1[e1,2]

days.e1


####

E2 <- tray.day.seed %>% filter(Tray.number=="E2")

e2 <-sapply(3:102, function(j){
  if (sum(E2[,j]) >=1) {
    print(min(which(E2[,j]==1)))
  } else {
    print(NA)
  }
})

days.e2 <- E2[e2,2]

days.e2


####

F1 <- tray.day.seed %>% filter(Tray.number=="F1")

f1 <-sapply(3:102, function(j){
  if (sum(F1[,j]) >=1) {
    print(min(which(F1[,j]==1)))
  } else {
    print(NA)
  }
})

days.f1 <- F1[f1,2]

days.f1

####

F2 <- tray.day.seed %>% filter(Tray.number=="F2")

f2 <-sapply(3:102, function(j){
  if (sum(F2[,j]) >=1) {
    print(min(which(F2[,j]==1)))
  } else {
    print(NA)
  }
})

days.f2 <- F2[f2,2]

days.f2

###

G1 <- tray.day.seed %>% filter(Tray.number=="G1")

g1 <-sapply(3:102, function(j){
  if (sum(G1[,j]) >=1) {
    print(min(which(G1[,j]==1)))
  } else {
    print(NA)
  }
})

days.g1 <- G1[g1,2]

days.g1

####

G2 <- tray.day.seed %>% filter(Tray.number=="G2")

g2 <-sapply(3:102, function(j){
  if (sum(G2[,j]) >=1) {
    print(min(which(G2[,j]==1)))
  } else {
    print(NA)
  }
})

days.g2 <- G2[g2,2]

days.g2

####


H1 <- tray.day.seed %>% filter(Tray.number=="H1")

h1 <-sapply(3:102, function(j){
  if (sum(H1[,j]) >=1) {
    print(min(which(H1[,j]==1)))
  } else {
    print(NA)
  }
})

days.h1 <- H1[h1,2]

days.h1

#### 

H2 <- tray.day.seed %>% filter(Tray.number=="H2")

h2 <-sapply(3:102, function(j){
  if (sum(H2[,j]) >=1) {
    print(min(which(H2[,j]==1)))
  } else {
    print(NA)
  }
})

days.h2 <- H2[h2,2]

days.h2


### Compiling the dataset


Germday <- c(days.a1, days.a2, days.b1, days.b2, days.c1, days.c2, days.d1,
             days.d2, days.e1, days.e2, days.f1, days.f2, days.g1, days.g2,
             days.h1, days.h2)

trays <- unique(tray.day.seed[,1])
trays

seedid <- sapply(1:16, function(x) {
  paste0(trays[x],"_Seed", c(1:100))
})



sps <- as.factor(c(rep("Abies alba", 800), rep("Abies nordmanniana", 200), 
             rep("Fagus Sylvatica", 600)))

prv <- unique(seedcounts[,"Provenance"])


prv <- rep(provenance, each=200)



trt <- unique(seedcounts[,"Treatment"])


trt <- rep(rep(treatment, each=100), times=8)

treatments

germ.day.data <- data.frame(SeedID = seedid, Species = sps, 
                           Provenance = prv, Treatment = trt, 
                           Germination.day = Germday) #dataset with factors

save(germ.day.data, file = "germinationdays.Rdata")

summary(germ.day.data)

germ.day.data.ch <- germ.day.data # dataset with no factor

germ.day.data.ch$Species <- as.character(germ.day.data.ch$Species)
germ.day.data.ch$Provenance <- as.character(germ.day.data.ch$Provenance)
germ.day.data.ch$Treatment <- as.character(germ.day.data.ch$Treatment)

str(germ.day.data.ch)


### Plots on germination data of each seed

par(mfrow=c(1,3))
hist(AN$Germination.day, main = "Abies nordmanniana", col="coral", xlab="Day of Germination")
hist(AA$Germination.day, main = "Abies alba", col="green", xlab="Day of Germination")
hist(FS$Germination.day, main = "Fagus Sylvatica", col="light blue", xlab="Day of Germination")


AA <- germ.day.data.ch %>% filter(Species == "Abies alba")
AN <- germ.day.data.ch %>% filter(Species == "Abies nordmanniana")
FS <- germ.day.data.ch %>% filter(Species == "Fagus Sylvatica")

par(mar=c(12,4,1,1))
par(mfrow=c(1,3))
boxplot(AA$Germination.day ~ AA$Treatment+ AA$Provenance , 
        col=2:6, las = 2, cex.axis= 1, xlab="")
boxplot(AN$Germination.day ~ AN$Treatment + AN$Provenance , 
        col=2:6, las = 2, cex.axis= 1, xlab="")
boxplot(FS$Germination.day ~ FS$Treatment + FS$Provenance , 
        col=2:6, las = 2, cex.axis= 1, xlab="")


library(dplyr)

library(ggplot2)

par(mfrow=c(1,3))

## Time spent until germination with zeros

ggplot(AA, aes(x=interaction(Treatment,Provenance), 
               y= Germination.day, fill=Provenance)) +
         geom_boxplot() +
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title="Abies alba", y= "Day of Germination", x = "Treatment & Provenance") 

ggplot(AN, aes(x=interaction(Treatment,Provenance), 
               y= Germination.day, fill=Provenance)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title="Abies nordmanniana",  "Day of Germination", x = "Treatment & Provenance") 


ggplot(FS, aes(x=interaction(Treatment,Provenance), 
               y= Germination.day, fill=Provenance)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title="Fagus sylvatica", y= "Day of Germination", x = "Treatment & Provenance") 

## Boxplots with interactions

ggplot(germ.day.data.ch, aes(x=interaction(Treatment,Provenance,Species), 
                             y=Germination.day, fill=Species)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title="All Species", y= "Day of Germination", x = "Species & Province & Provenance") 

ggplot(germ.day.data, aes(x=interaction(Treatment,Provenance,Species), 
                             y=Germday.nas, fill=Species)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title= "All Species", y= "Day of Germination", x = "Species & Province & Provenance")

####

library(dplyr)

ggarrange(bxp, dp, bp + rremove("x.text"), 
          labels = c("A", "B", "C"),
          ncol = 2, nrow = 2)

## Time spent until germination for each species

AA.n <- germ.day.data %>% filter(Species == "Abies alba")
AN.n <- germ.day.data %>% filter(Species == "Abies nordmanniana")
FS.n <- germ.day.data %>% filter(Species == "Fagus Sylvatica")

ggplot(AA.n, aes(x=interaction(Treatment,Provenance), 
               y= Germday.nas, fill=Provenance)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title="Abies alba", y= "Day of Germination", x = "Treatment & Provenance") 

ggplot(AN.n, aes(x=interaction(Treatment,Provenance), 
               y= Germination.day, fill=Provenance)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title="Abies nordmanniana",y= "Day of Germination", x = "Treatment & Provenance") 


ggplot(FS.n, aes(x=interaction(Treatment,Provenance), 
               y= Germination.day, fill=Provenance)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(title="Fagus sylvatica",y= "Day of Germination", x = "Treatment & Provenance") 

###

germ.day.data$Germday.zeros <- germ.day.data.ch$Germination.day
germ.day.data$Remaining.days <- germ.day.data.ch$Remaining.days

###

germ.trays <- rep(trays, each=100)
germ.trays

rm(Treatment)

library(tibble)

germ.day.data <- add_column(germ.day.data, germ.trays, .after = 1)


germ.day.data$germ.trays <- as.factor(germ.day.data$germ.trays)

a1tray <- germ.day.data %>% filter(germ.trays=="A1")

par(mfrow=c(4,4))
for (i in trays) {
  hist(germ.day.data.ch[germ.day.data.ch$germ.trays== i,]$Germination.day, 
       main = i, xlab = "Day")
}
  
library(tibble)  
  
###### zero-inflated model

library(pscl)

zeroinflatedpois <- zeroinfl(Germday.zeros ~ Species | Provenance, 
                             data = germ.day.data)


zeroinflatedpois.prtr <- zeroinfl(Germday.zeros ~ Provenance | Treatment, 
                                  data = germ.day.data)

summary(zeroinflatedpois)
summary(zeroinflatedpois.prtr)

zeroinflatedpois.all <- zeroinfl(Germday.zeros ~ Species, 
                                  data = germ.day.data)

summary(zeroinflatedpois.all)

?zeroinfl

###

meanzeros <- c()
meandays <- c()
totalzeros <- c()

for (i in trays) {
  meanzeros <- c(meanzeros, mean(germ.day.data[germ.day.data$germ.trays==i,]$Germday.zeros))
  meandays <- c(meandays, mean(germ.day.data[germ.day.data$germ.trays==i,]$Germination.day, na.rm = TRUE))
  totalzeros <- c(totalzeros, sum(germ.day.data[germ.day.data$germ.trays==i,]$Germday.zeros == 0))
}

sum(germ.day.data[germ.day.data$germ.trays=="A1",]$Germday.zeros == 0)

meanzeros
meandays

meangermdays <- data.frame(NonGerminatedIncluded= meanzeros, 
                           NonGerminatedExcluded=round(meandays,2),
                           total.zeros = totalzeros)

rownames(meangermdays) <- trays


meangermdays


germ.day.data$germinated.or.not <- ifelse(is.na(germ.day.data$Germination.day) == TRUE, 0,1)


####### cox semiparametric model

germ.day.data$germinated.or.not <- ifelse(is.na(germ.day.data$Germination.day) == TRUE, 0,1)

library(survival)

germ.day.data$Germday.nas <- germ.day.data$Germination.day

germ.day.data$Germination.day[is.na(germ.day.data$Germination.day)] = 72


attach(germ.day.data)

cox.mod <- coxph(Surv(Germination.day, germinated.or.not) ~ Species + Provenance + Treatment)

cox.mod.sp <- coxph(Surv(Germination.day, germinated.or.not) ~ Species + Treatment)


summary(cox.mod)
summary(cox.mod.sp)

germ.day.data.aa <- germ.day.data %>% filter(Species=="Abies alba")

germ.day.data.an <- germ.day.data %>% filter(Species=="Abies nordmanniana")

germ.day.data.fs <- germ.day.data %>% filter(Species=="Fagus Sylvatica")

germ.day.data.aa$germ.trays <- as.factor(as.character(germ.day.data.aa$germ.trays))
germ.day.data.aa$Provenance <- as.factor(as.character(germ.day.data.aa$Provenance))
germ.day.data.aa$Species <- as.factor(as.character(germ.day.data.aa$Species))

str(germ.day.data.aa)

attach(germ.day.data.aa)

cox.mod.aa <- coxph(Surv(Germination.day, germinated.or.not) ~ Provenance + Treatment)

summary(cox.mod.aa)

germ.day.data.an$germ.trays <- as.factor(as.character(germ.day.data.an$germ.trays))
germ.day.data.an$Provenance <- as.factor(as.character(germ.day.data.an$Provenance))
germ.day.data.an$Species <- as.factor(as.character(germ.day.data.an$Species))
germ.day.data.an$Treatment <- as.factor(as.character(germ.day.data.an$Treatment))

summary(germ.day.data.an)

attach(germ.day.data.an)
cox.mod.an <- coxph(Surv(Germination.day, germinated.or.not) ~ Treatment)

summary(cox.mod.an)

germ.day.data.fs$germ.trays <- as.factor(as.character(germ.day.data.fs$germ.trays))
germ.day.data.fs$Provenance <- as.factor(as.character(germ.day.data.fs$Provenance))
germ.day.data.fs$Species <- as.factor(as.character(germ.day.data.fs$Species))
germ.day.data.fs$Treatment <- as.factor(as.character(germ.day.data.fs$Treatment))


summary(germ.day.data.fs)

attach(germ.day.data.fs)

cox.mod.fs <- coxph(Surv(Germination.day, germinated.or.not) ~ Provenance + 
                      Treatment)

summary(cox.mod.fs)


##### accelerated failure time model


library(ciTools)

aft.wei.sp <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                 Treatment, data = germ.day.data)

aft.wei.ln <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                        Treatment, data = germ.day.data, dist = "lognormal")

summary(aft.wei.ln)

aft.wei.ll <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                        Treatment, data = germ.day.data, dist = "loglogistic")

summary(aft.wei.ll)

aft.wei.e <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                       Treatment, data = germ.day.data, dist = "exponential")

summary(aft.wei.e) 


aft.wei.g <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                       Treatment, data = germ.day.data, dist = "gaussian")
  

summary(aft.wei.g)


AIC(aft.wei.e)
AIC(aft.wei.ll)
AIC(aft.wei.ln)
AIC(aft.wei.sp)
AIC(aft.wei.g)


fitted_values <- aft.wei.ln$linear.predictors
resids <- (log(aft.wei.ln$y[, 1]) - fitted_values) / aft.wei.ln$scale


resKM <- survfit(Surv(resids, germinated.or.not) ~ 1, data = germ.day.data)


plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Probability of not germinating")
xx <- seq(min(resids), max(resids), length.out = 35)
yy <- pnorm(xx, lower.tail = FALSE)
lines(xx, yy, col = "red", lwd = 2)
legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                       "Survival function of Extreme Value distribution"), 
       lty = c(1,2,1), col = c(1,1,2), bty = "n")


aft.wei.ln.aa <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
          Treatment, data = germ.day.data.aa, dist = "lognormal")

summary(aft.wei.ln.aa)

aft.wei.ln.an <- survreg(Surv(Germination.day, germinated.or.not) ~
                           Treatment, data = germ.day.data.an, dist = "lognormal")

summary(aft.wei.ln.an)

aft.wei.ln.fs <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                         Treatment, data = germ.day.data.fs, dist = "lognormal")

summary(aft.wei.ln.fs)

save(germ.day.data, germ.day.data.aa, 
     germ.day.data.an, germ.day.data.fs, toglm, file = "tosurvanalysis.RData" )

write.csv(germ.day.data, file="survivalanalysis_dataformat.csv")

