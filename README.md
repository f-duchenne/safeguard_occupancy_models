# Scripts to fit occupancy models for the safeguard project.

Occupancy models were fitted using an approach based on what was done has been done in the *sparta* package (https://github.com/BiologicalRecordsCentre/sparta). However, as our list of species wsa too important to run the data preparation and models using *sparta*, we used customized code, taking advantage of the facilities provided by the *data.table* R package to improve performance.

Script 1 is filtering the species with enough data to be analysed and put the data in the right format to fit the model.

Script 2 is fitting the model for each species, independently, on a HPC.
