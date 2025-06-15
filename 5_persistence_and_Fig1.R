#this script explores the database for descriptive analysis.

#read the data-----
iucn <- read.csv(file = "data/Pollinators_IUCN.csv")
head(iucn)

dat <- read.csv(file = "data/Summary_species_presence.csv")
head(dat)
colnames(dat)[3:7] <- c("period_1921_1940",
                        "period_1941_1960",
                        "period_1961_1980",
                        "period_1981_2000",
                        "period_2001_2020")
#join
dim(dat)
dim(iucn)
iucn$Species <- gsub("_new", "", iucn$Species)
iucn$Species[which(!iucn$Species %in% dat$species)]
dat$species[which(!dat$species %in% iucn$Species)]
#IGNORE FOR NOW
dat <- merge(dat, iucn, by.x = "species", by.y = "Species")
dim(dat)
head(dat)

#plot the number of species that are in 0,1,2... time periods.
library(waffle)
#BEES
dat_bee <- subset(dat, taxo_group.x == "bees")
periods <- as.data.frame(table(dat_bee$Periods))
colnames(periods) <- c("group", "value")
periods$value <- round((periods$value/sum(periods$value))*100, digits = 0)
periods$value[6] <- 40 #manually adjust
sum(periods$value)
levels(periods$group) <- c("not in database", "1 time period", "2 time periods",
                           "3 time periods", "4 time periods", 
                           "5 time periods")
waffle(periods, flip = TRUE, size = 0, reverse = TRUE) #size = 2 for a waffle effect.
#We can add numbers to the first tile if we think is cool to show raw numbers.
#num <- table(dat$Periods)
#text(1,1, num[1], cex = .8)
#Not working, might be easier to edit the graph outside R.

#FLIES
dat_fly <- subset(dat, taxo_group.x == "hoverflies")
periods <- as.data.frame(table(dat_fly$Periods))
colnames(periods) <- c("group", "value")
periods$value <- round((periods$value/sum(periods$value))*100, digits = 0)
#periods$value[6] <- 40 #manually adjust
sum(periods$value)
levels(periods$group) <- c("not in database", "1 time period", "2 time periods",
                           "3 time periods", "4 time periods", 
                           "5 time periods")
waffle(periods, flip = TRUE, size = 0, reverse = TRUE) #size = 2 for a waffle effect.
#We can add numbers to the first tile if we think is cool to show raw numbers.
#num <- table(dat$Periods)
#text(1,1, num[1], cex = .8)
#Not working, might be easier to edit the graph outside R.












