### Zero Inflated Model

## Frequency Distribution Plots

#Abies alba

germ.day.data$germzeros <- germ.day.data$Germday.nas
germ.day.data$germzeros[is.na(germ.day.data$germzeros)] = 0

trays <- levels(germ.day.data$germ.trays)
traynames <- paste(rep(c("AA", "AN", "FS"), c(8,2,6)), toglm$Provenance, 
                   toglm$Moisture.treatment, sep = ".")

par(mfrow=c(2,4))
for (i in trays[1:8]) {
  hist(germ.day.data[germ.day.data$germ.trays== i,]$germzeros, 
       main = traynames[which(trays==i)], xlab = "Day",
       col=rainbow(8)[which(trays==i)])
}

#Abies nordmanniana

par(mfrow=c(1,2))
for (i in trays[9:10]) {
  hist(germ.day.data[germ.day.data$germ.trays== i,]$germzeros, 
       main = traynames[which(trays==i)], xlab = "Day",
       col=rainbow(16)[which(trays==i)])
}

#Fagus sylvatica

par(mfrow=c(2,3))
for (i in trays[11:16]) {
  hist(germ.day.data[germ.day.data$germ.trays== i,]$germzeros, 
       main = traynames[which(trays==i)], xlab = "Day",
       col=rainbow(16)[which(trays==i)])
}


##Statistical Tests

library(pscl)

germ.day.data.aa$germzeros <- germ.day.data.aa$Germday.nas
germ.day.data.aa$germzeros[is.na(germ.day.data.aa$germzeros)] = 0

germ.day.data.an$germzeros <- germ.day.data.an$Germday.nas
germ.day.data.an$germzeros[is.na(germ.day.data.an$germzeros)] = 0

germ.day.data.fs$germzeros <- germ.day.data.fs$Germday.nas
germ.day.data.fs$germzeros[is.na(germ.day.data.fs$germzeros)] = 0


zeroinflatedpois.all <- zeroinfl(germzeros ~ Species + Treatment, 
                                 data = germ.day.data)

summary(zeroinflatedpois.all)

zeroinflatedpois.aa <- zeroinfl(germzeros ~ Provenance + Treatment, 
                                 data = germ.day.data.aa)

summary(zeroinflatedpois.aa)

zeroinflatedpois.an <- zeroinfl(germzeros ~ Treatment, 
                                data = germ.day.data.an)

summary(zeroinflatedpois.an)

zeroinflatedpois.fs <- zeroinfl(germzeros ~ Provenance + Treatment, 
                                data = germ.day.data.fs)

summary(zeroinflatedpois.fs)

