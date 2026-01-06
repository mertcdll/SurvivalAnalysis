library(survival)
library(ggplot2)
library(flexsurv)
library(tibble)
library(dplyr)
library(PerformanceAnalytics)
library(survminer)

### Analysis of the Germination Data

## Import data

# Phase 1

germination.data.ph1 <- read.csv("Germination Data/germdata_phs1.csv", header = TRUE)
head(germination.data.ph1) #germination data of phase 1

# Phase 2

germination.data.ph2 <- read.csv("Germination Data/germdata_phs2.csv", header = TRUE)
head(germination.data.ph2) #germination data of phase 2

# Change Tray Names

germination.data.ph1$Tray.number <- 
  paste0(germination.data.ph1$Tray.number, "_ph1") #to give each tray a unique name
                                                   #because two phases share same tray names.

germination.data.ph2$Tray.number <- 
  paste0(germination.data.ph2$Tray.number, "_ph2")

### Overview of the Germination data

germ.data.wo.seeds.ph1 <- germination.data.ph1[,!grepl("Seed.", colnames(germination.data.ph1))]
head(germ.data.wo.seeds.ph1) #Excludes individual observations on seeds and which time they germinated
                             #and only takes the metadata and how many seeds germinated at each observation
                             #point

germ.data.wo.seeds.ph2 <- germination.data.ph2[,!grepl("Seed.", colnames(germination.data.ph2))]
head(germ.data.wo.seeds.ph2)

# For phase 1

germ.data.overview.ph1 <- 
  germ.data.wo.seeds.ph1[germ.data.wo.seeds.ph1[,"Position.in.chamber"] == "finish", ]

germ.data.overview.ph1$Absolut <- 
  germ.data.wo.seeds.ph1$Summe[germ.data.wo.seeds.ph1$Position.in.chamber == "absolut"]

germ.data.overview.ph1$Percentsum <- germ.data.overview.ph1$Summe / 100 
germ.data.overview.ph1$Percentabs <- germ.data.overview.ph1$Absolut / 100

head(germ.data.overview.ph1) # this one shows the total number of seeds 
                             #that are germinated for each provenance

# For phase 2

germ.data.overview.ph2 <- 
  germ.data.wo.seeds.ph2[germ.data.wo.seeds.ph2[,"Position.in.chamber"] == "finish", ]

germ.data.overview.ph2$Absolut <- 
  germ.data.wo.seeds.ph2$Summe[germ.data.wo.seeds.ph2$Position.in.chamber == "absolut"]

germ.data.overview.ph2$Percentsum <- germ.data.overview.ph2$Summe / 100 #adding germination rates (Percentage of seeds that are germinated)
germ.data.overview.ph2$Percentabs <- germ.data.overview.ph2$Absolut / 100 #adding germination rates (Percentage of seeds that are germinated)

head(germ.data.overview.ph2) # this one shows the total number of seeds 
                             #that are germinated for each provenance


# Combine them

germ.data.overview.all <- rbind(germ.data.overview.ph1, germ.data.overview.ph2)
head(germ.data.overview.all)
write.csv(germ.data.overview.all, file = "Dataframes/Germination Data Overview.csv", row.names = FALSE)

### Combining moist and dry treatments within the data

sums <- c()
absolute <- c()
perc.sum <- c()
perc.abs <- c()

germ.data.overview.all$Identification[c(25,26,27,28,31,32)] = c("e", "e" , "f", "f", "h", "h")
#giving the temporary identification codes to those that do not have any, 
#to run the loop below

for(i in unique(germ.data.overview.all$Identification)){
  
  sums <- c(sums,
            sum(germ.data.overview.all[germ.data.overview.all$Identification == i,]$Summe))
  absolute <- c(absolute,
                sum(germ.data.overview.all[germ.data.overview.all$Identification== i,]$Absolut))
  perc.sum <- c(perc.sum,
                sum(germ.data.overview.all[germ.data.overview.all$Identification == i,]$Summe)/ 200)
  perc.abs <- c(perc.abs,
                sum(germ.data.overview.all[germ.data.overview.all$Identification== i,]$Absolut/ 200))
  
}

germ.data.overview.comb <- data.frame(Summe = sums, Absolut = absolute, 
                                      Percentsum = perc.sum, Percentabs = perc.abs)


germ.data.overview.comb <- 
  cbind(germ.data.overview.all[germ.data.overview.all$Moisture.treatment=="dry",][,3:7], 
        germ.data.overview.comb)

head(germ.data.overview.comb)



#Add seed weights and moisture content to combined data

seed.weight.ph1 <- read.csv("Germination Data/seedweight_phs1.csv", header = TRUE) 
seed.weight.ph2 <- read.csv("Germination Data/seedweight_phs2.csv", header = TRUE)

sw.ph1 <- seed.weight.ph1$weight.100.seeds[seed.weight.ph1$State == "dried"]
mc.ph1 <- seed.weight.ph1$moisture.content....[seed.weight.ph1$State == "after stratification"]

sw.ph2 <- seed.weight.ph1$weight.100.seeds[seed.weight.ph2$State == "dried"]
mc.ph2 <- seed.weight.ph1$moisture.content....[seed.weight.ph2$State == "after stratification"]

seed.weights <- c(sw.ph1, sw.ph2)
moist.cont <- c(mc.ph1, mc.ph2)

germ.data.overview.comb$DryWeight = seed.weights
germ.data.overview.comb$MoistureCont = moist.cont

head(germ.data.overview.comb)

write.csv(germ.data.overview.comb,
          file = "Dataframes/Germination Data Overview (Treatments combined).csv", 
          row.names = FALSE)

### Plots of Seed Increase along intervals

germ.data.wo.seeds.all <- rbind(germ.data.wo.seeds.ph1, germ.data.wo.seeds.ph2)
germ.data.wo.seeds.all.noabs <- 
  germ.data.wo.seeds.all[germ.data.wo.seeds.all$Position.in.chamber != "absolut",]

## Abies

aa.germ.data <- germ.data.wo.seeds.all.noabs[germ.data.wo.seeds.all.noabs$Species 
                                   %in% c("Abies alba", "Abies nordmanniana"),]

library(RColorBrewer)

abies.plot <- ggplot(aa.germ.data
                     , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment), size = 1.5) + 
  scale_color_brewer(palette="Paired") +
  scale_x_continuous(breaks= seq(0,105, by=5)) +
  scale_y_continuous(breaks = seq (0,60,by = 5)) +
  geom_vline(xintercept = c(57,61), linetype="dotted", size = 1.5) +
  ggtitle("Abies alba & Abies nordmanniana") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds") +
  theme(legend.position="top", 
        plot.title = element_text(hjust = 0.5), 
        text = element_text(size = 18, face = "bold"))

abies.plot

## Fagus

fs.germ.data <- germ.data.wo.seeds.all.noabs[germ.data.wo.seeds.all.noabs$Species 
                                   %in% c("Fagus sylvatica", "Fagus orientalis"),]

fagus.plot <- ggplot(fs.germ.data
                     , aes(x = Time..d., y = Summe)) + 
  geom_line(aes(color = Provenance, linetype = Moisture.treatment), size = 1.5) + 
  scale_color_manual(values = 1:7) +
  scale_x_continuous(breaks= seq(0,105, by=5)) +
  scale_y_continuous(breaks = seq (0,60,by = 5)) +
  geom_vline(xintercept = c(57,61), linetype="dotted", size = 1.5) +
  ggtitle("Fagus Sylvatica & Fagus orientalis") +
  xlab("Time (Days)") + ylab("Number of Germinated Seeds") +
  theme(legend.position="top", 
        plot.title = element_text(hjust = 0.5), 
        text = element_text(size = 20, face = "bold"))

fagus.plot

### Correlation between the moisture content/seed weight and germination rates

# A general look

pairs(data=germ.data.overview.comb,
      ~ Summe + Absolut + DryWeight + MoistureCont)

germ.data.overview.comb.a <- #Overview for abies
  germ.data.overview.comb[germ.data.overview.comb$Species %in% c("Abies alba", "Abies nordmanniana"),]

germ.data.overview.comb.a

germ.data.overview.comb.f <- #Overview for fagus
  germ.data.overview.comb[germ.data.overview.comb$Species %in% c("Fagus sylvatica", "Fagus sylvatica"),]
  
  
germ.data.comp.a = germ.data.overview.comb.a[c("Summe", "Absolut", "DryWeight", "MoistureCont")]

germ.data.comp.a

chart.Correlation(germ.data.comp.a, #correlation chart for abies
                  method="pearson",
                  histogram=TRUE,
                  pch=30)


?chart.Correlation
germ.data.comp.f <- germ.data.overview.comb.f[c("Summe", "Absolut", "DryWeight", "MoistureCont")]

germ.data.comp.f

chart.Correlation(germ.data.comp.f[c(1:3)], #correlation chart for fagus (species that are not germinated are excluded)
                  method="pearson",
                  histogram=TRUE,
                  pch=20)

# It seems that there is no correlation between the final counts and seed weight/moisture content
# performing pairwise correlation or linear regression is unnecessary.
# Seeds having less dry weight tended to have more final moisture content **

### BASIC SURVIVAL ANALYSIS

## Transforming the data

#For phase 1

tray.day.seed.ph1 <- germination.data.ph1[germination.data.ph1$Position.in.chamber != "absolut",c(1,11,13:112)]
#this object has only the tray number, germination time and individual observation on seeds.


trays.ph1 <- unique(germ.data.wo.seeds.ph1$Tray.number) # tray names

germination.day.ph1 <- c() #an empty vector to record germination times for all seeds. 

names.ph1 <- sapply(1:16, function(i) {
  as.name(trays.ph1[i])
}) #an empty list for each tray

names.in.ph1 <- sapply(1:16, function(i) {
  as.name(paste0(trays.ph1[i], ".ind"))
}) #an empty list for each tray to record the seeds that are germinated

names.day.ph1 <- sapply(1:16, function(i) {
  as.name(paste0(trays.ph1[i], "day"))
}) #an empty list for each tray to record the day of germination

for (i in 1:16) { #this loop looks for 1's in the dataframe and records both the seed id and germination time of the seed in separate lists.
  
  
  names.ph1[[i]] <- tray.day.seed.ph1 %>% filter(Tray.number==trays.ph1[i])
  
  names.in.ph1[[i]] <- sapply(3:ncol(names.ph1[[i]]), function(j){
    if (sum(names.ph1[[i]][,j]) >=1) {
      print(min(which(names.ph1[[i]][,j]==1)))
    } else {
      print(NA)
    }
  })
  names.day.ph1[[i]] <- names.ph1[[i]][names.in.ph1[[i]],2]
  germination.day.ph1 <- c(germination.day.ph1, names.day.ph1[[i]])
}


seedid.ph1 <- sapply(1:16, function(x) { # giving each seed a unique name.
  paste0(trays.ph1[x],"_Seed", c(1:100))
})


sps.ph1 <- as.factor(c(rep("Abies alba", 800), rep("Abies nordmanniana", 200), 
                       rep("Fagus sylvatica", 600))) #creating a species factor for the final dataframe

provenance.ph1 <- unique(germination.data.ph1$Provenance) 

prv.ph1 <- rep(provenance.ph1, each=200) #creating a provenance column for the final dataframe

treatment.ph1 <- unique(germination.data.ph1$Moisture.treatment)

trt.ph1 <- rep(rep(treatment.ph1 , each=100), times=8) #creating a treatment column for the final dataframe

germ.day.data.ph1 <- data.frame(SeedID = as.vector(seedid.ph1), Species = sps.ph1, 
                                Provenance = prv.ph1, Treatment = trt.ph1, 
                                Germination.day = germination.day.ph1) #merge them into one dataframe
 

germday.ph1.na <- germ.day.data.ph1$Germination.day #record germination day column with nas to a different object

germ.day.data.ph1$germinated.or.not <- 
  ifelse(is.na(germ.day.data.ph1$Germination.day) == TRUE, 0,1)#add a germination status column (1=germinated, 0=non-germinated)

germ.day.data.ph1$Germination.day[is.na(germ.day.data.ph1$Germination.day)] = 72 #add 72(for non-germinated seeds-means that it is not germinated in 71 days but it can germinate in or after 72. day) in the place of nas for germination days column 

germ.day.data.ph1$germday.nas <- germday.ph1.na #assign germination days vector with nas to a different column.

germ.trays <- rep(trays.ph1, each=100) 

germ.day.data.ph1 <- add_column(germ.day.data.ph1, germ.trays, .after = 1) #add tray names

head(germ.day.data.ph1)#this is the data input format for basic survival analysis

# For phase 2 (long observation - 103 days)

head(germination.data.ph2)

tray.day.seed.ph2 <- germination.data.ph2[germination.data.ph2$Position.in.chamber != "absolut",c(1,11,13:112)]
#this object has only the tray number, germination time and individual observation on seeds.

trays.ph2 <- unique(germ.data.wo.seeds.ph2$Tray.number) #tray names

germination.day.ph2 <- c()  #an empty vector to record germination times for all seeds. 

names.ph2 <- sapply(1:16, function(i) {
  as.name(trays.ph2[i])
}) #an empty list for each tray

names.in.ph2 <- sapply(1:16, function(i) {
  as.name(paste0(trays.ph2[i], ".ind"))
}) #an empty list for each tray to record the seeds that are germinated

names.day.ph2 <- sapply(1:16, function(i) {
  as.name(paste0(trays.ph2[i], "day"))
}) #an empty list for each tray to record the day of germination

for (i in 1:16) { #this loop looks for 1's in the dataframe and records both the seed id and germination time of the seed in separate lists.
  
  names.ph2[[i]] <- tray.day.seed.ph2 %>% filter(Tray.number==trays.ph2[i])
  
  names.in.ph2[[i]] <- sapply(3:102, function(j){
    if (sum(names.ph2[[i]][,j]) >=1) {
      print(min(which(names.ph2[[i]][,j]==1)))
    } else {
      print(NA)
    }
  })
  names.day.ph2[[i]] <- names.ph2[[i]][names.in.ph2[[i]],2]
  germination.day.ph2 <- c(germination.day.ph2, names.day.ph2[[i]])
}

germination.day.ph2


seedid.ph2 <- sapply(1:16, function(x) {
  paste0(trays.ph2[x],"_Seed", c(1:100))
}) # giving each seed a unique name.

sps.ph2 <- as.factor(c(rep("Abies alba", 800),  #creating a species factor for the final dataframe

                       rep("Fagus sylvatica", 600), rep("Fagus orientalis", 200)))

provenance.ph2 <- unique(germination.data.ph2$Provenance)

prv.ph2 <- rep(provenance.ph2, each=200)  #creating a provenance column for the final dataframe

treatment.ph2 <- unique(germination.data.ph2$Moisture.treatment)

trt.ph2 <- rep(rep(treatment.ph2 , each=100), times=8) #creating a treatment column for the final dataframe

germ.day.data.ph2 <- data.frame(SeedID = as.vector(seedid.ph2), Species = sps.ph2, 
                                Provenance = prv.ph2, Treatment = trt.ph2, 
                                Germination.day = germination.day.ph2) #merge them into a dataframe


germday.ph2.na <- germ.day.data.ph2$Germination.day #record germination day column with nas to a different object

germ.day.data.ph2$germinated.or.not <- 
  ifelse(is.na(germ.day.data.ph2$Germination.day) == TRUE, 0,1) #add a germination status column (1=germinated, 0=non-germinated)

germ.day.data.ph2$Germination.day[is.na(germ.day.data.ph2$Germination.day)] = 104 #add 104(for non-germinated seeds-means that it is not germinated in 104 days but it can germinate in or after 104. day) in the place of nas for germination days column 

germ.day.data.ph2$germday.nas <- germday.ph2.na #assign germination days vector with nas to a different column.

germ.trays <- rep(trays.ph2, each=100)

germ.day.data.ph2 <- add_column(germ.day.data.ph2, germ.trays, .after = 1) #add trays

head(germ.day.data.ph2) #this is the data format for basic survival analysis

# For phase 2 (short observation - 71 days) 

tray.day.seed.ph2.s <- germination.data.ph2[germination.data.ph2$Time..d. <= 71,c(1,11,13:112)]

trays.ph2.s <- unique(germ.data.wo.seeds.ph2$Tray.number)

germination.day.ph2.s <- c()

names.ph2.s <- sapply(1:16, function(i) {
  as.name(trays.ph2.s[i])
})

names.in.ph2.s <- sapply(1:16, function(i) {
  as.name(paste0(trays.ph2.s[i], ".ind"))
})

names.day.ph2.s <- sapply(1:16, function(i) {
  as.name(paste0(trays.ph2.s[i], "day"))
})

for (i in 1:16) {
  
  names.ph2.s[[i]] <- tray.day.seed.ph2.s %>% filter(Tray.number==trays.ph2.s[i])
  
  names.in.ph2.s[[i]] <- sapply(3:102, function(j){
    if (sum(names.ph2.s[[i]][,j]) >=1) {
      print(min(which(names.ph2.s[[i]][,j]==1)))
    } else {
      print(NA)
    }
  })
  names.day.ph2.s[[i]] <- names.ph2.s[[i]][names.in.ph2.s[[i]],2]
  germination.day.ph2.s <- c(germination.day.ph2.s, names.day.ph2.s[[i]])
}

germination.day.ph2.s

length(germination.day.ph2.s)

seedid.ph2.s <- sapply(1:16, function(x) {
  paste0(trays.ph2.s[x],"_Seed", c(1:100))
})

sps.ph2.s <- as.factor(c(rep("Abies alba", 800), 
                       rep("Fagus sylvatica", 600), rep("Fagus orientalis", 200)))

provenance.ph2.s <- unique(germination.data.ph2$Provenance)

prv.ph2.s <- rep(provenance.ph2.s, each=200)

treatment.ph2.s <- unique(germination.data.ph2$Moisture.treatment)

trt.ph2.s <- rep(rep(treatment.ph2.s , each=100), times=8)

germ.day.data.ph2.s <- data.frame(SeedID = as.vector(seedid.ph2.s), Species = sps.ph2.s, 
                                Provenance = prv.ph2.s, Treatment = trt.ph2.s, 
                                Germination.day = germination.day.ph2.s)


germday.ph2.s.na <- germ.day.data.ph2.s$Germination.day

germ.day.data.ph2.s$germinated.or.not <- 
  ifelse(is.na(germ.day.data.ph2.s$Germination.day) == TRUE, 0,1)

germ.day.data.ph2.s$Germination.day[is.na(germ.day.data.ph2.s$Germination.day)] = 72

germ.day.data.ph2.s$germday.nas <- germday.ph2.s.na

germ.trays <- rep(trays.ph2.s, each=100)

germ.day.data.ph2.s <- add_column(germ.day.data.ph2.s, germ.trays, .after = 1)

head(germ.day.data.ph2.s)

### Estimators of the Survival Function and survival rate plot 

germ.day.data.all <- rbind(germ.day.data.ph1, germ.day.data.ph2) # combine the germination data 
write.csv(germ.day.data.all, "Dataframes/basic survival analysis input data format.csv", row.names = FALSE)

germ.day.data.all.abies <- germ.day.data.all[germ.day.data.all$Species %in% c("Abies alba", "Abies nordmanniana"),]
#subset for abies
germ.day.data.all.fagus <- germ.day.data.all[germ.day.data.all$Species %in% c("Fagus sylvatica", "Fagus orientalis"),]
#subset for fagus

KM.fit.abies.prv <- survfit(Surv(Germination.day, germinated.or.not)
                          ~ Provenance , data = germ.day.data.all.abies) #fits the data into model

abies.surv.plot <- ggsurvplot(KM.fit.abies.prv, data = germ.day.data.all.abies ,
                              break.time.by = 10, ggtheme = theme_light(), size= 2)
                             #plots the "survival probability versus time" of each tray

abies.surv.plot$plot + 
  theme(text = element_text(size = 18, color = "black", face = "bold"))+
  ylab("Probability of not germinating") +
  xlab("Time (Days)")#survival probability(Probability of not germinating) versus time plot for abies

KM.fit.fagus.prv <- survfit(Surv(Germination.day, germinated.or.not)
                            ~ Provenance , data = germ.day.data.all.fagus) #fits the data into model

fagus.surv.plot <- ggsurvplot(KM.fit.fagus.prv, data = germ.day.data.all.fagus,
                              break.time.by = 10, ggtheme = theme_light(), size =2) #plots the "survival probability versus time" of each tray


fagus.surv.plot$plot + 
  theme(text = element_text(size = 18, color = "black", face = "bold"))+
  ylab("Probability of not germinating")+
  xlab("Time (Days)")  #survival probability(Probability of not germinating) versus time plot for fagus

### Accelerated Failure Time Model

## For all data (Experiment 1 and Experiment 2 combined)

'%!in%' <- function(x,y)!('%in%'(x,y)) #simple function for subsetting - does the opposite of %in%

#Exclude the provenances that have no or too low germination 

germ.day.data.wonogerm.fagus <- germ.day.data.all.fagus[germ.day.data.all.fagus$Provenance %!in% c("Switzerland Salenstein", 
                                                                                                   "Romania SE Carpathians", 
                                                                                                   "Italy Passo Fittanze"),]
## For abies - Model comparison

#Weibull
aft.wei.abies <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                           Treatment, data = germ.day.data.all.abies)

#log-normal
aft.ln.abies <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                             Treatment, data = germ.day.data.all.abies, dist = "lognormal")

#log logistic
aft.ll.abies <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                             Treatment, data = germ.day.data.all.abies, dist = "loglogistic")

#exponential
aft.e.abies <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                            Treatment, data = germ.day.data.all.abies, dist = "exponential")

#gaussian
aft.g.abies <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                            Treatment, data = germ.day.data.all.abies, dist = "gaussian")


AIC(aft.e.abies)
AIC(aft.ll.abies)
AIC(aft.ln.abies) # log normal is the lowest one so go on with this.
AIC(aft.wei.abies)
AIC(aft.g.abies)


summary(aft.ln.abies) #summary table of the aft model. includes parameter estimates and their standard errors along with Z and p values.

flex.surv.abies <- flexsurvreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
              Treatment, data = germ.day.data.all.abies, dist = "lognormal")

flex.surv.abies

#output of this function also includes the exp(parameter estimate) and its confidence interval.

## For fagus - Model comparison

#Weibull
aft.wei.fagus <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                           Treatment, data = germ.day.data.wonogerm.fagus)

#log-normal
aft.ln.fagus <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                          Treatment, data = germ.day.data.wonogerm.fagus, dist = "lognormal")

#log logistic
aft.ll.fagus <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                          Treatment, data = germ.day.data.wonogerm.fagus, dist = "loglogistic")

#exponential
aft.e.fagus <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                         Treatment, data = germ.day.data.wonogerm.fagus, dist = "exponential")

#gaussian
aft.g.fagus <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                         Treatment, data = germ.day.data.wonogerm.fagus, dist = "gaussian")


AIC(aft.e.fagus)
AIC(aft.ll.fagus)
AIC(aft.ln.fagus) # log normal is the lowest one so go on with this.
AIC(aft.wei.fagus)
AIC(aft.g.fagus)

summary(aft.ln.fagus) #summary table of the aft model. includes parameter estimates and their standard errors along with Z and p values.

flexsurvreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
              Treatment, data = germ.day.data.wonogerm.fagus, dist = "lognormal")

#output of this function also includes the exp(parameter estimate) and its confidence interval.

### Comparison between short(71 day) and long (103 day) phase two experiment

#Prepare data

germ.day.abies.ph2.short <- germ.day.data.ph2.s[germ.day.data.ph2.s$Species %in% c("Abies alba"),]
germ.day.abies.ph2.long <- germ.day.data.ph2[germ.day.data.ph2$Species %in% c("Abies alba"),]

germ.day.fagus.ph2.short <- germ.day.data.ph2.s[germ.day.data.ph2.s$Species %in% c("Fagus sylvatica", "Fagus orientalis"),]
germ.day.fagus.ph2.long <- germ.day.data.ph2[germ.day.data.ph2$Species %in% c("Fagus sylvatica", "Fagus orientalis"),]

germ.day.fagus.ph2.short <-germ.day.fagus.ph2.short[germ.day.fagus.ph2.short$Provenance %!in% c("Switzerland Salenstein", 
                                                                   "Romania SE Carpathians", 
                                                                   "Italy Passo Fittanze"),]

germ.day.fagus.ph2.long <-germ.day.fagus.ph2.long[germ.day.fagus.ph2.long$Provenance %!in% c("Switzerland Salenstein", 
                                                                                                "Romania SE Carpathians", 
                                                                                                "Italy Passo Fittanze"),]

##Compare with AFT

#For abies

aft.ln.abies.ph2.short <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                          Treatment, data = germ.day.abies.ph2.short, dist = "lognormal")

aft.ln.abies.ph2.long <- survreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
                               Treatment, data = germ.day.abies.ph2.long, dist = "lognormal")

summary(aft.ln.abies.ph2.short)
summary(aft.ln.abies.ph2.long)

flexsurvreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
              Treatment, data = germ.day.abies.ph2.short, dist = "lognormal")

flexsurvreg(Surv(Germination.day, germinated.or.not) ~ Provenance +
              Treatment, data = germ.day.abies.ph2.long, dist = "lognormal")


# For Fagus

aft.ln.fagus.ph2.short <- survreg(Surv(Germination.day, germinated.or.not) ~ Treatment +
                                    Treatment, data = germ.day.fagus.ph2.short, dist = "lognormal")

aft.ln.fagus.ph2.long <- survreg(Surv(Germination.day, germinated.or.not) ~ Treatment +
                                   Treatment, data = germ.day.fagus.ph2.long, dist = "lognormal")

summary(aft.ln.fagus.ph2.short)
summary(aft.ln.fagus.ph2.long)

flexsurvreg(Surv(Germination.day, germinated.or.not) ~ 
              Treatment, data = germ.day.fagus.ph2.short, dist = "lognormal")

flexsurvreg(Surv(Germination.day, germinated.or.not) ~
              Treatment, data = germ.day.fagus.ph2.long, dist = "lognormal")


#################### FIN #######################

