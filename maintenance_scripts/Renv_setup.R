#This script ensures R packages used are documented, including its version

#load Renv library
library(renv)

#to set up the R dependency management
renv::init()

#See which packages we use where:
renv::dependencies()

#to update the dependency management 
renv::snapshot()


#In case you need to restore the environment:
#renv::restore() #commented to avoid problems