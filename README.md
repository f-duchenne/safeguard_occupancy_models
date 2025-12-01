[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17772148.svg)](https://doi.org/10.5281/zenodo.17772148) [![](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# [Code] Historical trends show more losers than winners among European pollinators

This Repo contains the code to reproduce all analysis and figures from "Historical trends show more losers than winners among European pollinators. 2026. C Martinez-Nuñez, F Duchenne N de Manincor, G Ghisbain, A Vujic, A Sentil, M Miličić, D Michez, I Bartomeus. Submitted"

You can cite as "Duchenne, Martinez-Nuñez and Bartomeus 2025. Code for Historical trends show more losers than winners among European pollinators. <https://zenodo.org/records/17772148>"

Note that the raw data is only available through the Zenodo repository due to its large size: [Link Zenodo]. In particular, the following files are not uploaded to github:

`/data/raw_data/Bee_DB_2025-09-03.rds`

`/data/raw_data/Hoverfly_DB_2025-09-03.rds`

`/data/raw_data/grids_shapefiles`

`/data/final_and_intermediate_outputs/database_clean_filtered.csv`

`/data/final_and_intermediate_outputs/traits_models_brms.RData`

Some large intermediate files are also not uploaded to GitHub, but will be created upon running the code.

Metadata (.json and .csv) for the three key datasets (inputs: `database_clean_filtered.csv, traits_table.csv` and output:`Table_S1.csv`) lives here: `/data/metadata` and can be visualized in html here: <https://f-duchenne.github.io/safeguard_occupancy_models/>

Scripts are numbered following the logical order to create first all needed inputs for subsequent analysis. For code requiring large computing time, intermediate outputs are stored in ``` /``final_and_intermediate_outputs ``` folder. Supplementary material scripts can be found in `/supplementary_codes`

If you want to use the exact same environment we used to run the analysis (i.e. same package versions) you can restore it using: `renv::restore()`. The folder `/maintenance_scrips` contains the scripts to set up `Renv` and create metadata.
