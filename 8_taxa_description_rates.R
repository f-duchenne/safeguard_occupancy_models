#Load IUCN file with taxonomic authority per species
master <- read.csv(file = "Side_projects/description_rates_EU.csv")
head(master)
#CARLOS: maybe this is already elsewhere in the repo.

master$taxonomic_authority

#Load Our Sup Mat with species per region
sm <- read.csv(file = "Side_projects/Summary_records_and_trends.csv", sep = ";")
head(sm)

#CARLOS: I think its better if we have a list of species present per biogeoregion, as here I am only using species with trends.


#select Med species (Trend_mediterranean)
med <- sm[which(!is.na(sm$Trend_mediterranean)),]
head(med)
#select Cont species (Trend_atlantic)
atl <- sm[which(!is.na(sm$Trend_atlantic)),]
#select Med species (Trend_boreal)
bor <- sm[which(!is.na(sm$Trend_boreal)),]
#select Med species (Trend_continental)
cont <- sm[which(!is.na(sm$Trend_continental)),]
#select Med species (Trend_alpine)
alp <- sm[which(!is.na(sm$Trend_alpine)),]

#subset from master, med species ----
d_med <- subset(master, friendly_name %in% med$TAXON) #Assuming taxonomy matches!!
head(d_med)

med_year <- c()
matches <- regmatches(d_med$taxonomic_authority, gregexpr("[[:digit:]]+", d_med$taxonomic_authority))
med_year <- as.numeric(unlist(matches))

#subset from master, continental species ----
d_cont <- subset(master, friendly_name %in% cont$TAXON) #Assuming taxonomy matches!!

cont_year <- c()
matches <- regmatches(d_cont$taxonomic_authority, gregexpr("[[:digit:]]+", d_cont$taxonomic_authority))
cont_year <- as.numeric(unlist(matches))

#subset from master, alp species----
d_alp <- subset(master, friendly_name %in% alp$TAXON) #Assuming taxonomy matches!!

alp_year <- c()
matches <- regmatches(d_alp$taxonomic_authority, gregexpr("[[:digit:]]+", d_alp$taxonomic_authority))
alp_year <- as.numeric(unlist(matches))

#subset from master, bor species----
d_bor <- subset(master, friendly_name %in% bor$TAXON) #Assuming taxonomy matches!!

bor_year <- c()
matches <- regmatches(d_bor$taxonomic_authority, gregexpr("[[:digit:]]+", d_bor$taxonomic_authority))
bor_year <- as.numeric(unlist(matches))

#subset from master, atl species----
d_atl <- subset(master, friendly_name %in% atl$TAXON) #Assuming taxonomy matches!!

atl_year <- c()
matches <- regmatches(d_atl$taxonomic_authority, gregexpr("[[:digit:]]+", d_atl$taxonomic_authority))
atl_year <- as.numeric(unlist(matches))

#Put all lines toguether

plot(as.numeric(cumsum(table(med_year))) ~ as.numeric(names(table(med_year))),
     las = 1, 
     type = "l", 
     ylab = "Species Described", 
     xlab = "Year", 
     col = 1, 
     main = "Cumulative Species Descriptions Over Time")
lines(as.numeric(cumsum(table(bor_year))) ~ as.numeric(names(table(bor_year))),
      col = 2)
lines(as.numeric(cumsum(table(alp_year))) ~ as.numeric(names(table(alp_year))),
      col = 3)
lines(as.numeric(cumsum(table(cont_year))) ~ as.numeric(names(table(cont_year))),
      col = 4)
lines(as.numeric(cumsum(table(atl_year))) ~ as.numeric(names(table(atl_year))),
      col = 5)

#CARLOS: Match colors with the rest of the ms.
