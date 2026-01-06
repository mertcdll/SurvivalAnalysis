library(msm)
library(ggplot2)
library(plyr)
library(ggpubr)
library(dplyr)
library(tibble)
### Analysis of the germination stages (MSM)

# import data

germ.stages.ph1 <- read.csv("Germination Data/germstage_phs1.csv", header = TRUE) #germination stages of phase 1

germ.stages.ph2 <- read.csv("Germination Data/germstage_phs2.csv", header = TRUE) #germination stages of phase 2

### convert the old stage data to a format which is compatible with msm

germ.stages.ph1$Tray.number = paste0(germ.stages.ph1$Tray.number,"_ph1") #give unique names to trays
germ.stages.ph2$Tray.number = paste0(germ.stages.ph2$Tray.number,"_ph2")


germ.stage.days.ph1 <- germ.stages.ph1[,c(1,11,13:112)] #these objects include only the tray name, time of stage transitions and individual observations on seeds. 
germ.stage.days.ph2 <- germ.stages.ph2[,c(1,11,13:112)]
germ.stage.days.ph1

# Phase 1 stage data

tray.stage.ph1 <- unique(germ.stage.days.ph1$Tray.number) #tray names

names.stage.ph1 <- sapply(1:16, function(i) {
  as.name(tray.stage.ph1[i])
}) #empty list to record the data from each tray

stage.df.ph1 <- sapply(1:16, function(i) {
  as.name(paste0(tray.stage.ph1 [i], ".df"))
}) #empty list to record the converted dataframes of each tray 

for (i in 1:16) { #this loop converts the seed stage matrix into a format that is compatible with msm
  names.stage.ph1[[i]] <- germ.stage.days.ph1 %>% filter(Tray.number == tray.stage.ph1[i])
  
  stage.df.ph1[[i]] <- data.frame (SeedID = paste(rep(colnames(names.stage.ph1[[i]][,3:102]), each = 23),
                                                  rep(names.stage.ph1[[i]]$Tray.number, 100), 
                                                  sep = "_"), 
                                   Time = rep(names.stage.ph1[[i]]$Time..d. , times= 100), 
                                   Stage = as.vector(as.matrix(names.stage.ph1[[i]][,3:102])))
  
}

stage.merged.ph1 <- as.data.frame(do.call(rbind, stage.df.ph1))

head(stage.merged.ph1)

species.stages.ph1 <- rep(germ.stages.ph1$Species, each=100) #add metadata to the dataframe
provenance.stages.ph1 <- rep(germ.stages.ph1$Provenance, each = 100)
treatment.stages.ph1 <- rep(germ.stages.ph1$Moisture.treatment, each = 100)

stage.merged.ph1 <-add_column(stage.merged.ph1, species.stages.ph1, .after = 1)
stage.merged.ph1 <-add_column(stage.merged.ph1, provenance.stages.ph1, .after = 2)
stage.merged.ph1 <-add_column(stage.merged.ph1, treatment.stages.ph1, .after = 3)

head(stage.merged.ph1) # this is the input format of MSM

# Phase 2 stage data

tray.stage.ph2 <- unique(germ.stage.days.ph2$Tray.number) #tray names

names.stage.ph2<- sapply(1:10, function(i) {
  as.name(tray.stage.ph2[i])
}) #empty list to record the data from each tray

stage.df.ph2 <- sapply(1:10, function(i) {
  as.name(paste0(tray.stage.ph2 [i], ".df"))
}) #empty list to record the converted dataframes of each tray

for (i in 1:10) { #this loop converts the seed stage matrix into a format that is compatible with msm
  names.stage.ph2[[i]] <- germ.stage.days.ph2 %>% filter(Tray.number == tray.stage.ph2[i])
  
  stage.df.ph2[[i]] <- data.frame (SeedID = paste(rep(colnames(names.stage.ph2[[i]][,3:102]), each = 30),
                                                       rep(names.stage.ph2[[i]]$Tray.number, 100),
                                                       sep = "_"), 
                                        Time = rep(names.stage.ph2[[i]]$Time..d. , times= 100), 
                                        Stage = as.vector(as.matrix(names.stage.ph2[[i]][,3:102])))
  
}

stage.merged.ph2<- as.data.frame(do.call(rbind, stage.df.ph2))

dim(stage.merged.ph2)
tail(stage.merged.ph2)

head(stage.merged.ph2)

species.stages.ph2 <- rep(germ.stages.ph2$Species, each=100) #add metadata
provenance.stages.ph2 <- rep(germ.stages.ph2$Provenance, each = 100)
treatment.stages.ph2 <- rep(germ.stages.ph2$Moisture.treatment, each = 100)

stage.merged.ph2 <-add_column(stage.merged.ph2, species.stages.ph2, .after = 1)
stage.merged.ph2 <-add_column(stage.merged.ph2, provenance.stages.ph2, .after = 2)
stage.merged.ph2 <-add_column(stage.merged.ph2, treatment.stages.ph2, .after = 3)

head(stage.merged.ph2)  # this is the input format of MSM

### MSM

library(msm)

colnames(stage.merged.ph1) <- c("SeedID", "Species", "Provenance", "Treatment", "Time", "Stage")
colnames(stage.merged.ph2) <- c("SeedID", "Species", "Provenance", "Treatment", "Time", "Stage")

stage.merged.all <- rbind (stage.merged.ph1, stage.merged.ph2) #merge the ph1 and ph2 seed stage data

write.csv(stage.merged.all, "Dataframes/MSM input data format.csv", row.names = FALSE)

stage.merged.ph2
## For each species

##For abies

aa.an.stage <- stage.merged.all[stage.merged.all$Species %in% c("Abies alba", "Abies nordmanniana"),]
#subset for abies

head(aa.an.stage)

aa.an.stage$status <- aa.an.stage$Stage + 1 #add new column by adding stages +1 because msm function does not accept 0 as a stage.

statetable.msm(status, subject = SeedID,  data=aa.an.stage) #stage transition matrix

q.abies <- rbind(c(0,0.5,0,0,0), #transition intensity matrix initial values
                 c(0,0,0.5,0,0),
                 c(0,0,0,0.5,0),
                 c(0,0,0,0,0.5),
                 c(0,0,0,0,0))

q.abies.crude <- crudeinits.msm(status ~ Time, SeedID, data = aa.an.stage, qmatrix = q.abies)
q.abies.crude #transition intensity matrix initial values calculated by crudeinits.msm function

abies.msm <- msm(status~Time, subject = SeedID, data=aa.an.stage,
                 qmatrix = q.abies.crude, deathexact = NULL,  
                 control=list(fnscale=4000,maxit = 10000)) #this function calculates transition intensities

abies.msm

pmatrix.msm(abies.msm, t= 100) # transition probability matrix
sojourn.msm(abies.msm) # mean sojourn times
totlos.msm(abies.msm) # total length of stay
plot(abies.msm, legend.pos = c(60, 0.4)) #survival plot

#Provenance as a covariate

abies.prov.msm <- msm(status~Time, subject = SeedID, data=aa.an.stage,
                      qmatrix = q.abies.crude, deathexact = NULL, covariates = ~ Provenance, 
                      control=list(fnscale=4000,maxit = 10000)) #this function calculates transition intensities and hazard ratios for provenance covariate.


abies.prov.msm # results
hazard.msm(abies.prov.msm) #hazard ratios of each abies provenance

# Heatmap of hazard ratios

abies.hazards <- sapply(1:8, function(i) {
  hazard.msm(abies.prov.msm)[[i]][,"HR"]
})

abies.hazards <- as.data.frame(abies.hazards)

crni.lazi.hazards <- exp(abies.prov.msm$estimates.t[1:4])

abies.hazards$crni.lazi <- crni.lazi.hazards #adding the intercept provenance

colnames(abies.hazards) <- c(gsub("Provenance", "",names(hazard.msm(abies.prov.msm))), "Croatia Crni Lazi")

abies.hazards
row.names(abies.hazards) <- c("Non-Germinated - Stage 1", "Stage 1 - Stage 2", "Stage 2 - Stage 3", "Stage 3 - Stage 4")
#hazard ratios of each provenance, as a dataframe
write.csv(abies.hazards, "Dataframes/Abies stage transition hazard ratios.csv")


rownames(abies.hazards)
colnames(abies.hazards)
as.vector(as.matrix(abies.hazards))

abies.hz <- data.frame(stage.tr = rep(rownames(abies.hazards), 9), 
                       Provenance = rep(colnames(abies.hazards), each=4), 
                       h.ratios = as.vector(as.matrix(abies.hazards))) #convert the dataframe to a plottable format

abies.hz 

ggplot(abies.hz, aes(x = stage.tr, y = Provenance, fill = h.ratios)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "blue") +
  ggtitle("Heatmap of Stage Transition Hazard Ratios for Abies alba & Abies nordmanniana") +
  xlab("Stage Transition") + ylab("Provenance") +
  theme(text = element_text(size=20))# heatmap of stage transition hazard ratios for abies


# Treatment as a covariate

abies.treat.msm <- msm(status~Time, subject = SeedID, data=aa.an.stage,
                       qmatrix = q.abies.crude, deathexact = NULL, covariates = ~ Treatment)

abies.treat.msm


## For fagus

fs.fo.stage <- stage.merged.all[stage.merged.all$Species %in% c("Fagus sylvatica", "Fagus orientalis"),]
#subset for fagus

fs.fo.stage$status <- fs.fo.stage$Stage +1 #add new column by adding stages +1 because msm function does not accept 0 as a stage.

statetable.msm(Stage, subject = SeedID,  data=fs.fo.stage) #stage transition matrix 

q.fagus <- rbind(c(0,0.5,0,0,0,0),  #transition intensity matrix initial values
                 c(0,0,0.5,0,0,0),
                 c(0,0,0,0.5,0,0),
                 c(0,0,0,0,0.5,0),
                 c(0,0,0,0,0,0.5),
                 c(0,0,0,0,0,0))

q.fagus.crude <- crudeinits.msm(status ~ Time, SeedID, data = fs.fo.stage, qmatrix = q.fagus)
q.fagus.crude #transition intensity matrix initial values calculated by crudeinits.msm function

fagus.msm <- msm(status~Time, subject = SeedID, data=fs.fo.stage,
                 qmatrix = q.fagus.crude, deathexact = NULL) #this function calculates transition intensities
fagus.msm

pmatrix.msm(fagus.msm, t= 100)  # transition probability matrix
sojourn.msm(fagus.msm) # mean sojourn times
totlos.msm(fagus.msm) # total length of stay
plot(fagus.msm, legend.pos = c(60, 0.77)) #survival plot


# Provenance as a covariate

fagus.prov.msm <- msm(status~Time, subject = SeedID, data=fs.fo.stage,
                      qmatrix = q.fagus.crude, deathexact = NULL, covariates = ~ Provenance, 
                      control=list(fnscale=4000,maxit = 10000))  #this function calculates transition intensities and hazard ratios for provenance covariate.

fagus.prov.msm

hazard.msm(fagus.prov.msm)  #hazard ratios of each abies provenance

# Heatmap of hazard ratios

fagus.hazards <- sapply(1:3, function(i) {
  hazard.msm(fagus.prov.msm)[[i]][,"HR"]
})

fagus.hazards <- as.data.frame(fagus.hazards) #hazard ratios of each provenance, as a dataframe

massif.hazards <- exp(abies.prov.msm$estimates.t[1:5])

fagus.hazards$massif <- massif.hazards #adding the intercept provenance

colnames(fagus.hazards) <- c(gsub("Provenance", "" ,names(hazard.msm(fagus.prov.msm))), 
                             "France Massif Armoricain") 

fagus.hazards
row.names(fagus.hazards) <- c("Non-Germinated - Stage 1", "Stage 1 - Stage 2", "Stage 2 - Stage 3", "Stage 3 - Stage 4", "Stage 4 - Stage 5")

write.csv(fagus.hazards, "Dataframes/Fagus stage transition hazard ratios.csv")

rownames(fagus.hazards)
colnames(fagus.hazards)
as.vector(as.matrix(fagus.hazards))

fagus.hz <- data.frame(stage.tr = rep(rownames(fagus.hazards), 4), 
                       Provenance = rep(colnames(fagus.hazards), each=5), 
                       h.ratios = as.vector(as.matrix(fagus.hazards))) #convert the dataframe to a plottable format

fagus.hz

ggplot(fagus.hz, aes(x = stage.tr, y = Provenance, fill = h.ratios)) +
  geom_tile(color = "black") +
  scale_fill_gradient(low = "white", high = "blue")+
  ggtitle("Heatmap of Stage Transition Hazard Ratios for Fagus sylvatica & Fagus orientalis") +
  xlab("Stage Transition") + ylab("Provenance") +
  theme(text = element_text(size=20)) # heatmap of stage transition hazard ratios for fagus


## Treatment as a covariate

fagus.treat.msm <- msm(status~Time, subject = SeedID, data=fs.fo.stage,
                       qmatrix = q.fagus.crude, deathexact = NULL, covariates = ~ Treatment)

fagus.treat.msm

#######################################################################################################


### Plotting stage transitions

germ.stages.all <- rbind(germ.stages.ph1, germ.stages.ph2)

abies.germ.stages <- germ.stages.all[germ.stages.all$Species %in% c("Abies alba", "Abies nordmanniana"),]
fagus.germ.stages <- germ.stages.all[germ.stages.all$Species %in% c("Fagus sylvatica", "Fagus orientalis"),]


## For abies - Converting the data

#Convert the data by counting the number of seeds that are in various stages for each observation time.

abies.tray.st <- unique(abies.germ.stages$Tray.number)

abies.prv.st <- paste(rep(unique(abies.germ.stages$Provenance), each = 2), rep(c("Dry", "Moist"),9), sep= " ")

names.ab.st <- sapply(1:18, function(i) {
  as.name(abies.prv.st[i])
})


names.ab.tray.st <- sapply(1:18, function(i) {
  as.name(abies.tray.st[i])
})


names.ab.df.st <- sapply(1:18, function(i) {
  as.name(abies.prv.st[i])
})


time.ph1 = abies.germ.stages$Time..d.[abies.germ.stages$Tray.number=="A1_ph1"]
time.ph2 = abies.germ.stages$Time..d.[abies.germ.stages$Tray.number=="A1_ph2"]

f.ab <- function(x) { # a simple counter function for abies stage data
  sapply(0:4, function(i){
    c(sum(x==i))
  })
}

stages.ab <- c("Non-Germinated", "Stage 1", "Stage 2", "Stage 3", "Stage 4")


for (i in 1:18) {
  
  names.ab.tray.st[[i]] <- abies.germ.stages %>% filter(Tray.number==abies.tray.st[i])
  
  names.ab.st[[i]] <- sapply(1:nrow(names.ab.tray.st[[i]]), function(j){
    f.ab(names.ab.tray.st[[i]][j,13:ncol(names.ab.tray.st[[i]])])
  })
  
  names.ab.df.st[[i]] <-  t(names.ab.st[[i]])
  
  if (nrow(names.ab.df.st[[i]]) == 23) {
    names.ab.df.st[[i]] <- data.frame(Count = as.vector(as.matrix(names.ab.df.st[[i]])),
                                      Time = rep(time.ph1,5), 
                                      Stage = rep(stages.ab, each = 23))
    
  }else if (nrow(names.ab.df.st[[i]]) == 30) {
    names.ab.df.st[[i]] <- data.frame(Count = as.vector(as.matrix(names.ab.df.st[[i]])),
                                      Time = rep(time.ph2,5), 
                                      Stage = rep(stages.ab, each = 30))
    
  }
    
}

head(names.ab.df.st) #stage count data in plottable format

plot.list.abies <- sapply(1:18, function(i) {
  as.name(abies.tray.st[i])
})

names(names.ab.df.st) <- abies.prv.st

for (i in seq(names.ab.df.st)) {
  plot.list.abies[[i]] <- ggplot(names.ab.df.st[[i]], aes(x = Time, y = Count))+
    geom_line(aes(color = Stage), size = 1.5) +
    ggtitle(names(names.ab.df.st[i])) +
    theme(text = element_text(size=15), legend.text = element_text(size = 20, face = "bold")) +
    xlab("Time (Days)")
}



ggarrange(plotlist = plot.list.abies, ncol = 4, nrow = 5, 
          common.legend = TRUE) #plots for abies, shows the number of seeds that are in various stages throughout the experiment.

## For Fagus - Converting the data

fagus.tray.st <- unique(fagus.germ.stages$Tray.number)

fagus.prv.st <- paste(rep(unique(fagus.germ.stages$Provenance), each = 2), rep(c("Dry", "Moist"),4), sep= " ")

names.fg.st <- sapply(1:8, function(i) {
  as.name(fagus.prv.st[i])
})


names.fg.tray.st <- sapply(1:8, function(i) {
  as.name(fagus.tray.st[i])
})


names.fg.df.st <- sapply(1:8, function(i) {
  as.name(fagus.prv.st[i])
})


f.fg <- function(x) { # a simple counter function for fagus stage data
  sapply(0:5, function(i){
    c(sum(x==i))
  })
}

stages.fg <- c("Non-Germinated", "Stage 1", "Stage 2", "Stage 3", "Stage 4", "Stage 5")


for (i in 1:8) {
  
  names.fg.tray.st[[i]] <- fagus.germ.stages %>% filter(Tray.number==fagus.tray.st[i])
  
  names.fg.st[[i]] <- sapply(1:nrow(names.fg.tray.st[[i]]), function(j){
    f.fg(names.fg.tray.st[[i]][j,13:ncol(names.fg.tray.st[[i]])])
  })
  
  names.fg.df.st[[i]] <-  t(names.fg.st[[i]])
  
  if (nrow(names.fg.df.st[[i]]) == 23) {
    names.fg.df.st[[i]] <- data.frame(Count = as.vector(as.matrix(names.fg.df.st[[i]])),
                                      Time = rep(time.ph1,6), 
                                      Stage = rep(stages.fg, each = 23))
    
  }else if (nrow(names.fg.df.st[[i]]) == 30) {
    names.fg.df.st[[i]] <- data.frame(Count = as.vector(as.matrix(names.fg.df.st[[i]])),
                                      Time = rep(time.ph2,6), 
                                      Stage = rep(stages.fg, each = 30))
    
  }
  
}

head(names.fg.df.st) #stage count data in plottable format

plot.list.fagus <- sapply(1:8, function(i) {
  as.name(fagus.tray.st[i])
})

names(names.fg.df.st) <- fagus.prv.st

for (i in seq(names.fg.df.st)) {
  plot.list.fagus[[i]] <- ggplot(names.fg.df.st[[i]], aes(x = Time, y = Count))+
    geom_line(aes(color = Stage), size = 1.5) +
    ggtitle(names(names.fg.df.st[i])) +
    xlab("Time (Days)")+
    theme(text = element_text(size=15), legend.text = element_text(size = 20, face = "bold")) +
    xlab("Time (Days)")
}



ggarrange(plotlist = plot.list.fagus, ncol = 4, nrow = 2, 
          common.legend = TRUE)  #plots for fagus, shows the number of seeds that are in various stages throughout the experiment.
