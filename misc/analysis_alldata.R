library(tibble)
library(dplyr)
library(ggplot2)
library(survival)
library(flexsurv)
library(msm)
library(ggh4x)
library(drc)

### import data

germination.data <- read.csv("Germination Data/GermDataAll.csv", header = TRUE)

### give different names to trays for old and new experiment

germination.data$Tray.number = c(paste0(germination.data[1:368,]$Tray.number, "_Old"), 
                                 paste0(germination.data[369:720,]$Tray.number, "_New"))

### Total number of seeds germinated & Seed Increase per interval 

germ.data.wo.seeds <- germination.data[,!grepl("Seed.", colnames(germination.data))]

germ.data.overview <- germ.data.wo.seeds[germ.data.wo.seeds[,"Position.in.chamber"] == "finish", ]

# Old experiment seed increase

germination.data.old <- germination.data[1:368,]
germ.data.wo.seeds.old <- germination.data.old[,!grepl("Seed.", colnames(germination.data.old))]

seed.matrix.old <- matrix(germ.data.wo.seeds.old [,"Summe"],
                          nrow=16, ncol=23, byrow=TRUE)

seed.df.old <- as.data.frame(seed.matrix.old)

colnames(seed.df.old) <- germ.data.wo.seeds.old[,11][1:23]

metadata.old <- aggregate(Summe ~ Tray.number + Species + Country +
                        Provenance + Moisture.treatment , data=germ.data.wo.seeds.old, mean)[-6]

metadata.old <- metadata.old[order(metadata.old[,1]),]

seed.count.days.old <- cbind(metadata.old, seed.df.old)


# New experiment seed increase

germination.data.new <- germination.data[369:720,]

germ.data.wo.seeds.new <- germination.data.new[,!grepl("Seed.", colnames(germination.data.new))]

seed.matrix.new <- matrix(germ.data.wo.seeds.new [,"Summe"],
                          nrow=16, ncol=22, byrow=TRUE)
  
seed.df.new <- as.data.frame(seed.matrix.new)

colnames(seed.df.new) <- germ.data.wo.seeds.new[,11][1:22]

metadata.new <- aggregate(Summe ~ Tray.number + Species + Country +
                            Provenance + Moisture.treatment , data=germ.data.wo.seeds.new, mean)[-6]
metadata.new <- metadata.new[order(metadata.new[,1]),]


seed.count.days.new <- cbind(metadata.new, seed.df.new)

# New experiment seed increase (Long one)

germdata.new.long <- read.csv("Germination Data/germ_newexplong.csv", header = TRUE)

germdata.new.long.noabs <- germdata.new.long[germdata.new.long$Position.in.chamber != "absolut" ,]

germdata.new.long.noabs.woseed <- germdata.new.long.noabs[,!grepl("Seed.", colnames(germdata.new.long.noabs))]


### Plots of Seed Increase along intervals

library(ggplot2)
library(dplyr)

### Old experiment

#Abies alba & Abies nordmanniana

aa.germ.data.old <- germ.data.wo.seeds.old[germ.data.wo.seeds.old$Species 
                                           %in% c("Abies alba", "Abies nordmanniana"),]

ggplot(aa.germ.data.old
       , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment)) + 
  scale_color_manual(values = 1:5) +
  scale_x_continuous(breaks= seq(0,71, by=5)) +
  scale_y_continuous(breaks = seq (0,60,by = 5)) +
  geom_vline(xintercept = 57, linetype="dotted") +
  ggtitle("Abies Alba & Abies Nordmanniana (Ambrolauri Tlugi)") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds")


#Fagus sylvatica

fs.germ.data.old <- germ.data.wo.seeds.old[germ.data.wo.seeds.old$Species 
                                           %in% c("Fagus sylvatica"),]


ggplot(fs.germ.data.old
       , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment)) + 
  scale_color_manual(values = 1:5) +
  scale_x_continuous(breaks= seq(0,71, by=5)) +
  scale_y_continuous(breaks = seq (0,50,by = 5)) +
  geom_vline(xintercept = 57, linetype="dotted") +
  ggtitle("Fagus Sylvatica") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds")

### New Experiment

# Abies alba

aa.germ.data.new <- germ.data.wo.seeds.new[germ.data.wo.seeds.new$Species 
                                           %in% c("Abies alba"),]

ggplot(aa.germ.data.new
       , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment)) + 
  scale_color_manual(values = 1:5) +
  scale_x_continuous(breaks= seq(0,71, by=5)) +
  scale_y_continuous(breaks = seq (0,60,by = 5)) +
  geom_vline(xintercept = 61, linetype="dotted") +
  ggtitle("Abies Alba") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds")

# Fagus sylvatica


fs.germ.data.new <- germ.data.wo.seeds.new[germ.data.wo.seeds.new$Species 
                                           %in% c("Fagus sylvatica", "Fagus orientalis"),]


ggplot(fs.germ.data.new
       , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment)) + 
  scale_color_manual(values = 1:5) +
  scale_x_continuous(breaks= seq(0,71, by=5)) +
  scale_y_continuous(breaks = seq (0,60,by = 5)) +
  geom_vline(xintercept = 61, linetype="dotted") +
  ggtitle("Fagus Sylvatica & Fagus Orientalis") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds")

### New Experiment (Long one)

#Abies alba

aa.germ.data.newlong <- germdata.new.long.noabs.woseed[germdata.new.long.noabs.woseed$Species 
                                           %in% c("Abies alba"),]

ggplot(aa.germ.data.newlong
       , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment), size = 1) + 
  scale_color_manual(values = 1:5) +
  scale_x_continuous(breaks= seq(0,105, by=5)) +
  scale_y_continuous(breaks = seq (0,60,by = 5)) +
  geom_vline(xintercept = 61, linetype="dotted") +
  ggtitle("Abies Alba") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds")

# Fagus sylvatica & Fagus orientalis

fs.germ.data.newlong <- germdata.new.long.noabs.woseed[germdata.new.long.noabs.woseed$Species 
                                           %in% c("Fagus sylvatica", "Fagus orientalis"),]


ggplot(fs.germ.data.newlong
       , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment), size = 1) + 
  scale_color_manual(values = 1:5) +
  scale_x_continuous(breaks= seq(0,105, by=5)) +
  scale_y_continuous(breaks = seq (0,60,by = 5)) +
  geom_vline(xintercept = 61, linetype="dotted") +
  ggtitle("Fagus Sylvatica & Fagus Orientalis") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds")



### Compiling all

# Abies alba & Abies nordmanniana


aa.germ.data <- germ.data.wo.seeds[germ.data.wo.seeds$Species 
                                       %in% c("Abies alba", "Abies nordmanniana"),]

library(ggplot2)

ggplot(aa.germ.data
       , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment), size = 1) + 
  scale_color_manual(values = 1:9) +
  scale_x_continuous(breaks= seq(0,71, by=5)) +
  scale_y_continuous(breaks = seq (0,60,by = 5)) +
  geom_vline(xintercept = c(57,61), linetype="dotted") +
  ggtitle("Abies alba & Abies nordmanniana") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds")

# Fagus sylvatica & fagus orientalis

fs.germ.data <- germ.data.wo.seeds[germ.data.wo.seeds$Species 
                                   %in% c("Fagus sylvatica", "Fagus orientalis"),]

ggplot(fs.germ.data
       , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment), size = 1) + 
  scale_color_manual(values = 1:7) +
  scale_x_continuous(breaks= seq(0,71, by=5)) +
  scale_y_continuous(breaks = seq (0,60,by = 5)) +
  geom_vline(xintercept = c(57,61), linetype="dotted") +
  ggtitle("Fagus Sylvatica & Fagus orientalis") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds")

### Compiling all for old and new experiment(long one)

## Abies

aa.germ.data.all.long <- rbind(aa.germ.data.old, aa.germ.data.newlong)

abies.plot <- ggplot(aa.germ.data.all.long
       , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment), size = 1) + 
  scale_color_manual(values = 1:9) +
  scale_x_continuous(breaks= seq(0,105, by=5)) +
  scale_y_continuous(breaks = seq (0,60,by = 5)) +
  geom_vline(xintercept = c(57,61), linetype="dotted") +
  ggtitle("Abies alba & Abies nordmanniana") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds")

abies.plot
  
nsrmoist <- nsr[nsr$Moisture.treatment == "moist",]
nsrdry <-   nsr[nsr$Moisture.treatment == "dry",]
  
nsr <- aa.germ.data.all.long[aa.germ.data.all.long$Provenance == "North Slovakia Region",]
pdhp <- aa.germ.data.all.long[aa.germ.data.all.long$Provenance == "Prealpes de Haute-Provence",]


geom_ribbon(data=nsr, 
            aes(ymin=nsr[nsr$Moisture.treatment == "dry",]$Summe,
                ymax=nsr[nsr$Moisture.treatment == "moist",]$Summe),
            fill="blue", alpha=0.5)


nsr[nsr$Moisture.treatment == "moist",]$Summe
nsr[nsr$Moisture.treatment == "dry",]$Summe


## Fagus

fs.germ.data.all.long <- rbind(fs.germ.data.old, fs.germ.data.newlong)

fagus.plot <- ggplot(fs.germ.data.all.long
       , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment), size = 1) + 
  scale_color_manual(values = 1:9) +
  scale_x_continuous(breaks= seq(0,105, by=5)) +
  scale_y_continuous(breaks = seq (0,60,by = 5)) +
  geom_vline(xintercept = c(57,61), linetype="dotted") +
  ggtitle("Fagus Sylvatica & Fagus orientalis") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds")

fagus.plot

### DRM


model.aa <- drm(Summe ~ Time..d., fct = LN.2(),
             type = "event", curveid = Tray.number, 
             data = aa.germ.data.all.long, 
             pmodels=list(~ 1, ~ Tray.number - 1))


### Transforming the data for survival analysis

## Old experiment

tray.day.seed.old <- germination.data.old[,c(1,11,13:112)]

trays.old <- unique(germ.data.wo.seeds.old$Tray.number)

germination.day.old <- c()

names.old <- sapply(1:16, function(i) {
  as.name(trays.old[i])
})

names.in.old <- sapply(1:16, function(i) {
  as.name(paste0(trays.old[i], ".ind"))
})

names.day.old <- sapply(1:16, function(i) {
  as.name(paste0(trays.old[i], "day"))
})

for (i in 1:16) {
  
  
  names.old[[i]] <- tray.day.seed.old %>% filter(Tray.number==trays.old[i])
  
  names.in.old[[i]] <- sapply(3:102, function(j){
    if (sum(names.old[[i]][,j]) >=1) {
      print(min(which(names.old[[i]][,j]==1)))
    } else {
      print(NA)
    }
  })
  names.day.old[[i]] <- names.old[[i]][names.in.old[[i]],2]
  germination.day.old <- c(germination.day.old, names.day.old[[i]])
}

seedid.old <- sapply(1:16, function(x) {
  paste0(trays.old[x],"_Seed", c(1:100))
})

as.vector(seedid.old)

sps.old <- as.factor(c(rep("Abies alba", 800), rep("Abies nordmanniana", 200), 
                   rep("Fagus sylvatica", 600)))

provenance.old <- unique(germination.data.old$Provenance)

prv.old <- rep(provenance.old, each=200)

treatment.old <- unique(germination.data.old$Moisture.treatment)

trt.old <- rep(rep(treatment.old , each=100), times=8)

germ.day.data.old <- data.frame(SeedID = as.vector(seedid.old), Species = sps.old, 
                            Provenance = prv.old, Treatment = trt.old, 
                            Germination.day = germination.day.old)


germday.old.na <- germ.day.data.old$Germination.day

germ.day.data.old$germinated.or.not <- 
  ifelse(is.na(germ.day.data.old$Germination.day) == TRUE, 0,1)

germ.day.data.old$Germination.day[is.na(germ.day.data.old$Germination.day)] = 72

germ.day.data.old$germday.nas <- germday.old.na

library(tibble)

germ.trays <- rep(trays.old, each=100)

germ.day.data.old <- add_column(germ.day.data.old, germ.trays, .after = 1)


## New experiment

tray.day.seed.new <- germination.data.new[,c(1,11,13:112)]

trays.new <- unique(germ.data.wo.seeds.new$Tray.number)

germination.day.new <- c()

names.new <- sapply(1:16, function(i) {
  as.name(trays.new[i])
})

names.in.new <- sapply(1:16, function(i) {
  as.name(paste0(trays.new[i], ".ind"))
})

names.day.new <- sapply(1:16, function(i) {
  as.name(paste0(trays.new[i], "day"))
})

for (i in 1:16) {
  
  
  names.new[[i]] <- tray.day.seed.new %>% filter(Tray.number==trays.new[i])
  
  names.in.new[[i]] <- sapply(3:102, function(j){
    if (sum(names.new[[i]][,j]) >=1) {
      print(min(which(names.new[[i]][,j]==1)))
    } else {
      print(NA)
    }
  })
  names.day.new[[i]] <- names.new[[i]][names.in.new[[i]],2]
  germination.day.new <- c(germination.day.new, names.day.new[[i]])
}

germination.day.new


seedid.new <- sapply(1:16, function(x) {
  paste0(trays.new[x],"_Seed", c(1:100))
})

as.vector(seedid.new)

sps.new <- as.factor(c(rep("Abies alba", 800), 
                       rep("Fagus sylvatica", 600), rep("Fagus orientalis", 200)))

provenance.new <- unique(germination.data.new$Provenance)

prv.new <- rep(provenance.new, each=200)

treatment.new <- unique(germination.data.new$Moisture.treatment)

trt.new <- rep(rep(treatment.new , each=100), times=8)

germ.day.data.new <- data.frame(SeedID = as.vector(seedid.new), Species = sps.new, 
                                Provenance = prv.new, Treatment = trt.new, 
                                Germination.day = germination.day.new)


germday.new.na <- germ.day.data.new$Germination.day

germ.day.data.new$germinated.or.not <- 
  ifelse(is.na(germ.day.data.new$Germination.day) == TRUE, 0,1)

germ.day.data.new$Germination.day[is.na(germ.day.data.new$Germination.day)] = 72

germ.day.data.new$germday.nas <- germday.new.na

germ.trays <- rep(trays.new, each=100)

germ.day.data.new <- add_column(germ.day.data.new, germ.trays, .after = 1)


## New experiment (Long one)

tray.day.seed.newlong <- germdata.new.long.noabs[,c(1,11,13:112)]

trays.newlong <- unique(germdata.new.long.noabs$Tray.number)

germination.day.newlong <- c()

names.newlong <- sapply(1:16, function(i) {
  as.name(trays.newlong[i])
})

names.in.newlong <- sapply(1:16, function(i) {
  as.name(paste0(trays.newlong[i], ".ind"))
})

names.day.newlong <- sapply(1:16, function(i) {
  as.name(paste0(trays.newlong[i], "day"))
})

for (i in 1:16) {
  
  
  names.newlong[[i]] <- tray.day.seed.newlong %>% filter(Tray.number==trays.newlong[i])
  
  names.in.newlong[[i]] <- sapply(3:102, function(j){
    if (sum(names.newlong[[i]][,j]) >=1) {
      print(min(which(names.newlong[[i]][,j]==1)))
    } else {
      print(NA)
    }
  })
  names.day.newlong[[i]] <- names.newlong[[i]][names.in.newlong[[i]],2]
  germination.day.newlong <- c(germination.day.newlong, names.day.newlong[[i]])
}

germination.day.newlong


seedid.newlong <- sapply(1:16, function(x) {
  paste0(trays.newlong[x],"_Seed", c(1:100))
})

as.vector(seedid.newlong)

sps.newlong <- as.factor(c(rep("Abies alba", 800), 
                       rep("Fagus sylvatica", 600), rep("Fagus orientalis", 200)))

provenance.newlong <- unique(germination.data.new$Provenance)

prv.newlong <- rep(provenance.new, each=200)

treatment.newlong <- unique(germination.data.new$Moisture.treatment)

trt.newlong <- rep(rep(treatment.new , each=100), times=8)

germ.day.data.newlong <- data.frame(SeedID = as.vector(seedid.newlong), Species = sps.newlong, 
                                Provenance = prv.newlong, Treatment = trt.newlong, 
                                Germination.day = germination.day.newlong)


germday.newlong.na <- germ.day.data.newlong$Germination.day


germ.day.data.newlong$germinated.or.not <- 
  ifelse(is.na(germ.day.data.newlong$Germination.day) == TRUE, 0,1)

germ.day.data.newlong$Germination.day[is.na(germ.day.data.newlong$Germination.day)] = 104

germ.day.data.newlong$germday.nas <- germday.newlong.na


germ.trays <- rep(trays.newlong, each=100)

germ.day.data.newlong <- add_column(germ.day.data.newlong, germ.trays, .after = 1)

### Basic survival analysis of new experiment

##Estimators of the Survival Function and survival rate plot

KM.fit.new <- survfit(Surv(Germination.day, germinated.or.not)
                  ~ 1, data = germ.day.data.new) # Kaplan-Meier Estimator

KM.fit.new

KM.fit.trays.new <- sapply(trays.new, function(x) {
  list(survfit(Surv(germ.day.data.new[germ.day.data.new$germ.trays == x,"Germination.day"],
                    germ.day.data.new[germ.day.data.new$germ.trays == x,"germinated.or.not"]) ~ 1))
})

head(KM.fit.trays.new)

traynames.new <- paste(rep(c("AA", "FS", "FO"), c(8,6,2)), 
                       germ.data.overview$Provenance[17:32], 
                       germ.data.overview$Moisture.treatment[17:32], sep = ".")

plot(KM.fit.trays.new[[1]], conf.int = FALSE, lwd = 2, 
     col = rainbow(16)[1], xlab = "Time(Days)", 
     ylab = "Probability of not germinating")

sapply(2:16, function(i) {
  lines(KM.fit.trays.new[[i]], conf.int = FALSE, lwd =2, col = rainbow(16)[i])
})

legend("bottomleft", traynames.new, col = rainbow(16) , 
       lwd =2 , lty =1, bty = "n" ,cex = 0.45)


## Estimators of the Survival Function and survival rate plot (long one)

KM.fit.newlong <- survfit(Surv(Germination.day, germinated.or.not)
                      ~ 1, data = germ.day.data.newlong) # Kaplan-Meier Estimator

KM.fit.newlong

KM.fit.trays.newlong <- sapply(trays.newlong, function(x) {
  list(survfit(Surv(germ.day.data.newlong[germ.day.data.newlong$germ.trayslong == x,"Germination.day"],
                    germ.day.data.newlong[germ.day.data.newlong$germ.trayslong == x,"germinated.or.not"]) ~ 1))
})

head(KM.fit.trays.newlong)

traynames.newlong <- paste(rep(c("AA", "FS", "FO"), c(8,6,2)), 
                       germ.data.overview$Provenance[17:32], 
                       germ.data.overview$Moisture.treatment[17:32], sep = ".")

plot(KM.fit.trays.newlong[[1]], conf.int = FALSE, lwd = 2, 
     col = rainbow(16)[1], xlab = "Time(Days)", 
     ylab = "Probability of not germinating")

sapply(2:16, function(i) {
  lines(KM.fit.trays.newlong[[i]], conf.int = FALSE, lwd =2, col = rainbow(16)[i])
})

legend("bottomleft", traynames.newlong, col = rainbow(16) , 
       lwd =2 , lty =1, bty = "n" ,cex = 0.45)


## Accelerated Failure Time Model for new experiment (Parametric)

aft.wei.new <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                     Treatment, data = germ.day.data.new) #weibull

aft.wei.new

aft.ln.new <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                    Treatment, data = germ.day.data.new, dist = "lognormal")

summary(aft.ln.new)

aft.ll.new <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                    Treatment, data = germ.day.data.new, dist = "loglogistic")
summary(aft.ll.new)

aft.e.new <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                   Treatment, data = germ.day.data.new, dist = "exponential")

summary(aft.e.new) 


aft.g.new <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                   Treatment, data = germ.day.data.new, dist = "gaussian")

summary(aft.g.new)


AIC(aft.e.new)
AIC(aft.ll.new)# log logistic is the lowest one so go on with this.
AIC(aft.ln.new) 
AIC(aft.wei.new)
AIC(aft.g.new)

library(flexsurv)

# for species difference

flexsurvreg(Surv(Germination.day, germinated.or.not) ~ Species +
              Treatment, data = germ.day.data.new, dist = "llogis")


library(dplyr)

germ.day.data.new.aa <- germ.day.data.new %>% filter(Species == "Abies alba")
germ.day.data.new.fs <- germ.day.data.new %>% filter(Species == "Fagus sylvatica")
germ.day.data.new.fo <- germ.day.data.new %>% filter(Species == "Fagus orientalis")

# for abies alba

aft.ll.aa.new <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                       Treatment, data = germ.day.data.new.aa, dist = "loglogistic")

summary(aft.ll.aa.new)


flexsurvreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
              Treatment, data = germ.day.data.new.aa, dist = "llogis")

# for fagus sylvatica

aft.ll.fs.new <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                       Treatment, data = germ.day.data.new.fs, dist = "loglogistic")

summary(aft.ln.fs.new) ## all zero - cannot perform analysis

## for fagus orientalis

aft.ll.fo.new <- survreg(Surv(Germination.day, germinated.or.not) ~ 
                       Treatment, data = germ.day.data.fs, dist = "loglogistic")

summary(aft.ll.fo.new)

### AFT model for new experiment (long one)

aft.wei.newlong <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                         Treatment, data = germ.day.data.newlong) #weibull

aft.wei.newlong

aft.ln.newlong <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                        Treatment, data = germ.day.data.newlong, dist = "lognormal")

summary(aft.ln.newlong)

aft.ll.newlong <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                        Treatment, data = germ.day.data.newlong, dist = "loglogistic")
summary(aft.ll.newlong)

aft.e.newlong <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                       Treatment, data = germ.day.data.newlong, dist = "exponential")

summary(aft.e.newlong) 


aft.g.newlong <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                       Treatment, data = germ.day.data.newlong, dist = "gaussian")

summary(aft.g.newlong)


AIC(aft.e.newlong)
AIC(aft.ll.newlong)
AIC(aft.ln.newlong) # log normal is the lowest one so go on with this.
AIC(aft.wei.newlong)
AIC(aft.g.newlong)

library(flexsurv)

# for species difference

flexsurvreg(Surv(Germination.day, germinated.or.not) ~ Species +
              Treatment, data = germ.day.data.newlong, dist = "lognormal")


library(dplyr)

germ.day.data.new.aalong <- germ.day.data.newlong %>% filter(Species == "Abies alba")
germ.day.data.new.fslong <- germ.day.data.newlong %>% filter(Species == "Fagus sylvatica")
germ.day.data.new.folong <- germ.day.data.newlong %>% filter(Species == "Fagus orientalis")

# for abies alba

aft.ln.aa.newlong <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                           Treatment, data = germ.day.data.new.aalong, dist = "lognormal")

summary(aft.ln.aa.newlong)


flexsurvreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
              Treatment, data = germ.day.data.new.aalong, dist = "lognormal")

# for fagus sylvatica

aft.ln.fs.newlong <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                           Treatment, data = germ.day.data.new.fslong, dist = "lognormal")

summary(aft.ln.fs.newlong)

## for fagus orientalis

aft.ln.fo.newlong <- survreg(Surv(Germination.day, germinated.or.not) ~ 
                           Treatment, data = germ.day.data.new.folong, dist = "lognormal")

summary(aft.ln.fo.newlong)

### Basic survival analysis of two experiments combined

##Estimators of the Survival Function and survival rate plot

germ.day.data.all <- rbind(germ.day.data.old, germ.day.data.new)

KM.fit.all <- survfit(Surv(Germination.day, germinated.or.not)
                      ~ 1, data = germ.day.data.all) # Kaplan-Meier Estimator

KM.fit.all

trays.all <- unique(germ.day.data.all$germ.trays)

KM.fit.trays.all <- sapply(trays.all, function(x) {
  list(survfit(Surv(germ.day.data.all[germ.day.data.all$germ.trays == x,"Germination.day"],
                    germ.day.data.all[germ.day.data.all$germ.trays == x,"germinated.or.not"]) ~ 1))
})

head(KM.fit.trays.new)

tray.names.old <-paste(rep(c("AA", "AN", "FS"), c(8,2,6)), 
                       germ.data.overview$Provenance[1:16], 
                       germ.data.overview$Moisture.treatment[1:16], sep = ".")

traynames.new <- paste(rep(c("AA", "FS", "FO"), c(8,6,2)), 
                       germ.data.overview$Provenance[17:32], 
                       germ.data.overview$Moisture.treatment[17:32], sep = ".")

tray.names.all <- c(tray.names.old, traynames.new)


plot(KM.fit.trays.all[[1]], conf.int = FALSE, lwd = 2, 
     col = rainbow(32)[1], xlab = "Time(Days)", 
     ylab = "Probability of not germinating")

sapply(2:32, function(i) {
  lines(KM.fit.trays.all[[i]], conf.int = FALSE, lwd =2, col = rainbow(32)[i])
})

legend("bottomleft", tray.names.all, col = rainbow(32) , 
       lwd =2 , lty =1, bty = "n" ,cex = 0.23)

## Accelerated Failure Time Model for all data (Parametric)

aft.wei.all <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                     Treatment, data = germ.day.data.all) #weibull

summary(aft.wei.all)

aft.ln.all <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                    Treatment, data = germ.day.data.all, dist = "lognormal")

summary(aft.ln.all)

aft.ll.all <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                    Treatment, data = germ.day.data.all, dist = "loglogistic")
summary(aft.ll.all)

aft.e.all <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                   Treatment, data = germ.day.data.all, dist = "exponential")

summary(aft.e.all) 


aft.g.all <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                   Treatment, data = germ.day.data.all, dist = "gaussian")

summary(aft.g.all)


AIC(aft.e.all)
AIC(aft.ll.all)
AIC(aft.ln.all) # log normal is the lowest one so go on with this.
AIC(aft.wei.all)
AIC(aft.g.all)

germ.day.data.all.aa <- germ.day.data.all %>% filter(Species == "Abies alba")
germ.day.data.all.an <- germ.day.data.all %>% filter(Species == "Abies nordmanniana")
germ.day.data.all.fs <- germ.day.data.all %>% filter(Species == "Fagus sylvatica")
germ.day.data.all.fo <- germ.day.data.all %>% filter(Species == "Fagus orientalis")


aft.ln.aa.all <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                       Treatment, data = germ.day.data.all.aa, dist = "lognormal")

summary(aft.ln.aa.all)

aft.ln.an.all <- survreg(Surv(Germination.day, germinated.or.not) ~
                       Treatment, data = germ.day.data.all.an, dist = "lognormal")

summary(aft.ln.an.all)


aft.ln.fs.all <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                       Treatment, data = germ.day.data.all.fs, dist = "lognormal")

summary(aft.ln.fs.all)


aft.ln.fo.all <- survreg(Surv(Germination.day, germinated.or.not) ~ 
                           Treatment, data = germ.day.data.all.fo, dist = "lognormal")

summary(aft.ln.fo.all)

### Basic survival analysis of two experiments combined (with longer new one)

germ.day.data.all.long <- rbind(germ.day.data.old, germ.day.data.newlong)

KM.fit.alllong <- survfit(Surv(Germination.day, germinated.or.not)
                      ~ 1, data = germ.day.data.all.long)
KM.fit.alllong

trays.alllong <- unique(germ.day.data.all.long$germ.trays)

KM.fit.trays.alllong <- sapply(trays.alllong, function(x) {
  list(survfit(Surv(germ.day.data.all.long[germ.day.data.all.long$germ.trays == x,"Germination.day"],
                    germ.day.data.all.long[germ.day.data.all.long$germ.trays == x,"germinated.or.not"]) ~ 1))
})


head(KM.fit.trays.alllong)


plot(KM.fit.trays.alllong[[1]], conf.int = FALSE, lwd = 2, 
     col = rainbow(32)[1], xlab = "Time(Days)", 
     ylab = "Probability of not germinating")

sapply(2:32, function(i) {
  lines(KM.fit.trays.alllong[[i]], conf.int = FALSE, lwd =2, col = rainbow(32)[i])
})

legend("bottomleft", tray.names.all, col = rainbow(32) , 
       lwd =2 , lty =1, bty = "n" ,cex = 0.23)

## Accelerated Failure Time Model for all data (With longer new experiment)

aft.wei.all.long <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                         Treatment, data = germ.day.data.all.long) #weibull

summary(aft.wei.all.long)

aft.ln.all.long <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                        Treatment, data = germ.day.data.all.long, dist = "lognormal")

summary(aft.ln.all.long)

aft.ll.all.long <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                        Treatment, data = germ.day.data.all.long, dist = "loglogistic")
summary(aft.ll.all.long)

aft.e.all.long <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                       Treatment, data = germ.day.data.all.long, dist = "exponential")

summary(aft.e.all.long) 


aft.g.all.long <- survreg(Surv(Germination.day, germinated.or.not) ~ Species +
                       Treatment, data = germ.day.data.all.long, dist = "gaussian")

summary(aft.g.all.long)


AIC(aft.e.all.long)
AIC(aft.ll.all.long)
AIC(aft.ln.all.long) # log normal is the lowest one so go on with this.
AIC(aft.wei.all.long)
AIC(aft.g.all.long)

germ.day.data.all.aa.long <- germ.day.data.all.long %>% filter(Species == "Abies alba")
germ.day.data.all.an.long <- germ.day.data.all.long %>% filter(Species == "Abies nordmanniana")
germ.day.data.all.fs.long <- germ.day.data.all.long %>% filter(Species == "Fagus sylvatica")
germ.day.data.all.fo.long <- germ.day.data.all.long %>% filter(Species == "Fagus orientalis")


aft.ln.aa.all.long <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                           Treatment, data = germ.day.data.all.aa.long, dist = "lognormal")

summary(aft.ln.aa.all.long)

aft.ln.an.all.long <- survreg(Surv(Germination.day, germinated.or.not) ~
                           Treatment, data = germ.day.data.all.an.long, dist = "lognormal")

summary(aft.ln.an.all.long)


aft.ln.fs.all.long <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                           Treatment, data = germ.day.data.all.fs.long, dist = "lognormal")

summary(aft.ln.fs.all.long)


aft.ln.fo.all.long <- survreg(Surv(Germination.day, germinated.or.not) ~ 
                           Treatment, data = germ.day.data.all.fo.long
                           , dist = "lognormal")

summary(aft.ln.fo.all.long)


summary(aft.ln.aa.all)
summary(aft.ln.aa.all.long)


### Analysis of the germination stages (MSM)

# import data

germ.stages.old <- read.csv("Germination Data/germstagesold.csv", header = TRUE)

germ.stages.new <- read.csv("Germination Data/germstagesnew.csv", header = TRUE)

germ.stages.new.long <- read.csv("Germination Data/germstage_newexplong.csv", header = TRUE)

### convert the old stage data to a format which is compatible with msm

germ.stages.old$Tray.number = paste0(germ.stages.old$Tray.number,"_old")
germ.stages.new$Tray.number = paste0(germ.stages.new$Tray.number,"_new")


germ.stage.days.old <- germ.stages.old[,c(1,11,13:112)]
germ.stage.days.new <- germ.stages.new[,c(1,11,13:112)]
germ.stage.days.new.long <- germ.stages.new.long[,c(1,11,13:112)]

# old stage data

tray.stage.old <- unique(germ.stage.days.old$Tray.number)

names.stage.old <- sapply(1:16, function(i) {
  as.name(tray.stage.old[i])
})

stage.df.old <- sapply(1:16, function(i) {
  as.name(paste0(tray.stage.old [i], ".df"))
})

for (i in 1:16) {
  names.stage.old[[i]] <- germ.stage.days.old %>% filter(Tray.number == tray.stage.old[i])
  
  stage.df.old[[i]] <- data.frame (SeedID = paste(rep(colnames(names.stage.old[[i]][,3:102]), each = 23),
                                                  rep(names.stage.old[[i]]$Tray.number, 100), 
                                                  sep = "_"), 
                                   Time = rep(names.stage.old[[i]]$Time..d. , times= 100), 
                                   Stage = as.vector(as.matrix(names.stage.old[[i]][,3:102])))
  
}

stage.merged.old <- as.data.frame(do.call(rbind, stage.df.old))

head(stage.merged.old)

species.stages.old <- rep(germ.stages.old$Species, each=100)
provenance.stages.old <- rep(germ.stages.old$Provenance, each = 100)
treatment.stages.old <- rep(germ.stages.old$Moisture.treatment, each = 100)

stage.merged.old <-add_column(stage.merged.old, species.stages.old, .after = 1)
stage.merged.old <-add_column(stage.merged.old, provenance.stages.old, .after = 2)
stage.merged.old <-add_column(stage.merged.old, treatment.stages.old, .after = 3)


# new stage data

tray.stage.new <- unique(germ.stage.days.new$Tray.number)

names.stage.new <- sapply(1:10, function(i) {
  as.name(tray.stage.new[i])
})

stage.df.new <- sapply(1:10, function(i) {
  as.name(paste0(tray.stage.new [i], ".df"))
})

for (i in 1:10) {
  names.stage.new[[i]] <- germ.stage.days.new %>% filter(Tray.number == tray.stage.new[i])
  
  stage.df.new[[i]] <- data.frame (SeedID = paste(rep(colnames(names.stage.new[[i]][,3:102]), each = 22),
                                                  rep(names.stage.new[[i]]$Tray.number, 100),
                                                  sep = "_"), 
                                   Time = rep(names.stage.new[[i]]$Time..d. , times= 100), 
                                   Stage = as.vector(as.matrix(names.stage.new[[i]][,3:102])))
  
}

stage.merged.new <- as.data.frame(do.call(rbind, stage.df.new))


head(stage.merged.new)

species.stages.new <- rep(germ.stages.new$Species, each=100)
provenance.stages.new <- rep(germ.stages.new$Provenance, each = 100)
treatment.stages.new <- rep(germ.stages.new$Moisture.treatment, each = 100)

stage.merged.new <-add_column(stage.merged.new, species.stages.new, .after = 1)
stage.merged.new <-add_column(stage.merged.new, provenance.stages.new, .after = 2)
stage.merged.new <-add_column(stage.merged.new, treatment.stages.new, .after = 3)

head(stage.merged.new)

#New stage data (long)

tray.stage.new.long <- unique(germ.stage.days.new.long$Tray.number)

names.stage.new.long <- sapply(1:10, function(i) {
  as.name(tray.stage.new.long[i])
})

stage.df.new.long <- sapply(1:10, function(i) {
  as.name(paste0(tray.stage.new.long [i], ".df"))
})

for (i in 1:10) {
  names.stage.new.long[[i]] <- germ.stage.days.new.long %>% filter(Tray.number == tray.stage.new.long[i])
  
  stage.df.new.long[[i]] <- data.frame (SeedID = paste(rep(colnames(names.stage.new.long[[i]][,3:102]), each = 30),
                                                  rep(names.stage.new.long[[i]]$Tray.number, 100),
                                                  sep = "_"), 
                                   Time = rep(names.stage.new.long[[i]]$Time..d. , times= 100), 
                                   Stage = as.vector(as.matrix(names.stage.new.long[[i]][,3:102])))
  
}

stage.merged.new.long <- as.data.frame(do.call(rbind, stage.df.new.long))


head(stage.merged.new.long)

species.stages.new.long <- rep(germ.stages.new.long$Species, each=100)
provenance.stages.new.long <- rep(germ.stages.new.long$Provenance, each = 100)
treatment.stages.new.long <- rep(germ.stages.new.long$Moisture.treatment, each = 100)

stage.merged.new.long <-add_column(stage.merged.new.long, species.stages.new.long, .after = 1)
stage.merged.new.long <-add_column(stage.merged.new.long, provenance.stages.new.long, .after = 2)
stage.merged.new.long <-add_column(stage.merged.new.long, treatment.stages.new.long, .after = 3)


### MSM

library(msm)

statetable.msm(Stage, subject = SeedID,  data=stage.merged.old)


statetable.msm(Stage, subject = SeedID,  data=stage.merged.new.long)

colnames(stage.merged.new.long) <- c("SeedID", "Species", "Provenance", "Treatment", "Time", "Stage")
colnames(stage.merged.old) <- c("SeedID", "Species", "Provenance", "Treatment", "Time", "Stage")

stage.merged.all <- rbind (stage.merged.old, stage.merged.new.long)

# For each species

aa.an.stage <- stage.merged.all[stage.merged.all$Species %in% c("Abies alba", "Abies nordmanniana"),]

fs.fo.stage <- stage.merged.all[stage.merged.all$Species %in% c("Fagus sylvatica", "Fagus orientalis"),]

statetable.msm(Stage, subject = SeedID,  data=aa.an.stage)

dim(aa.an.stage)

statetable.msm(Stage, subject = SeedID,  data=fs.fo.stage)

dim(fs.fo.stage)

?statetable.msm

aa.an.stage$status <- aa.an.stage$Stage + 1

fs.fo.stage$status <- fs.fo.stage$Stage +1


q.abies <- rbind(c(0,0.5,0,0,0),
              c(0,0,0.5,0,0),
              c(0,0,0,0.5,0),
              c(0,0,0,0,0.5),
              c(0,0,0,0,0))


q.fagus <- rbind(c(0,0.5,0,0,0,0),
                 c(0,0,0.5,0,0,0),
                 c(0,0,0,0.5,0,0),
                 c(0,0,0,0,0.5,0),
                 c(0,0,0,0,0,0.5),
                 c(0,0,0,0,0,0))





q.abies.crude <- crudeinits.msm(status ~ Time, SeedID, data = aa.an.stage, qmatrix = q.abies)
q.abies.crude

q.abies.ratio <- statetable.msm(Stage,subject = SeedID, data=aa.an.stage)/46999

q.fagus.crude <- crudeinits.msm(status ~ Time, SeedID, data = fs.fo.stage, qmatrix = q.fagus)

?crudeinits.msm

abies.msm <- msm(status~Time, subject = SeedID, data=aa.an.stage,
                 qmatrix = q.abies.crude, deathexact = NULL, 
                 control=list(fnscale=4000,maxit = 10000))

abies.msm

pmatrix.msm(abies.msm, t= 100)
sojourn.msm(abies.msm)
totlos.msm(abies.msm)
plot(abies.msm, legend.pos = c(60, 0.3))

abies.prov.msm <- msm(status~Time, subject = SeedID, data=aa.an.stage,
                      qmatrix = q.abies.crude, deathexact = NULL, covariates = ~ Provenance, 
                      control=list(fnscale=4000,maxit = 10000))
abies.prov.msm

qmatrix.msm(abies.prov.msm, covariates = list(Provenance = "ProvenanceAmbrolauri Tlugi"))
qmatrix.msm(abies.prov.msm, covariates = list(Provenance = "ProvenanceSudety Region"))


abies.treat.msm <- msm(status~Time, subject = SeedID, data=aa.an.stage,
                       qmatrix = q.abies.crude, deathexact = NULL, covariates = ~ Treatment)

abies.treat.msm



fagus.msm <- msm(status~Time, subject = SeedID, data=fs.fo.stage,
                 qmatrix = q.fagus.crude, deathexact = NULL)
fagus.msm

pmatrix.msm(fagus.msm, t= 100)
sojourn.msm(fagus.msm)
totlos.msm(fagus.msm)
plot(fagus.msm, legend.pos = c(60, 0.77))

fagus.prov.msm <- msm(status~Time, subject = SeedID, data=fs.fo.stage,
                      qmatrix = q.fagus.crude, deathexact = NULL, covariates = ~ Provenance, 
                      control=list(fnscale=4000,maxit = 10000))

fagus.prov.msm

fagus.treat.msm <- msm(status~Time, subject = SeedID, data=fs.fo.stage,
                      qmatrix = q.fagus.crude, deathexact = NULL, covariates = ~ Treatment)
fagus.plot
fagus.treat.msm

qmatrix.msm(abies.treat.msm, covariates = list(Treatment="moist"))

qmatrix.msm(abies.treat.msm, covariates = list(Treatment="dry"))




er <- c()

for (i in 1:nrow(stage.merged.new.long)) {
  if (stage.merged.new.long$Stage[i+1] == 0 & stage.merged.new.long$Stage[i] == 1){
    er <- c(er, i)
  }
}

look <- stage.merged.new.long[er,]
look


abies.plot
