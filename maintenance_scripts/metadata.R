#How to set up metadata follwoing https://annakrystalli.me/dataspice-tutorial/

# Get the last version of dataspice and load it
#install.packages("readr")
install.packages("devtools")
devtools::install_github("ropenscilabs/dataspice")
library(readr)
library(dataspice)

#load data (see demo.qmd for details)
githubURL <- "https://github.com/Joemillard/GlitrsDynametaLocal/raw/refs/heads/master/inst/shiny_data/prior_data.rds"
d <- readRDS(url(githubURL))

#create basic .csv metadata files and folders.
create_spice()

#Add creators
edit_creators()

#Add how to access the data
prep_access()
edit_access()

#Add metadata

#useful to get ranges before adding those:
range(d$Year) 
range(d$Latitude, na.rm = T)
range(d$Longitude, na.rm = T)

edit_biblio()

#Describe variables
prep_attributes()
colnames(d)
edit_attributes()

#create a json file
write_spice()

#look at the json:
jsonlite::read_json(here::here("data", "metadata", "dataspice.json")) %>% listviewer::jsonedit()

#build a webpage!
build_site()


