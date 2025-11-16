#How to set up metadata follwoing https://annakrystalli.me/dataspice-tutorial/

#Get last version on Github
install.packages("devtools")
devtools::install_github("ropenscilabs/dataspice")
library(dataspice)

#create basic .csv metadata files and folders.
create_spice()

#Add creators
edit_creators()

#Add how to access the data
prep_access()
edit_access()

#Add metadata

#useful to get ranges before adding those:
#species <- read.csv("data/final_and_intermediate_outputs/database_clean_filtered.csv", stringsAsFactors = T) ## 4158540 records.
range(species$endYear) 
range(species$decimalLatitude, na.rm = T)
range(species$decimalLongitude, na.rm = T)

edit_biblio()

#Describe variables
prep_attributes()
edit_attributes()

#create a json file
write_spice()

#look at the json:
jsonlite::read_json(here::here("data", "metadata", "dataspice.json")) %>% listviewer::jsonedit()

#build a webpage!
build_site()


