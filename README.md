# PALM-Crowdsourcing
Complementary scripts to the following paper:
van der Linden L, Hogan P, Maronga B, Hagemann R, Bechtel B (2023) Crowdsourcing air temperature data for the evaluation of the urban microscale model PALM—A case study in central Europe. PLOS Clim 2(8): e0000197. https://doi.org/10.1371/journal.pclm.0000197

Python script to process model results from PALM to display in QGIS and for further statistical analysis in R
- Python version 3.8.10
- required packages: xarray, numpy, geopandas, rioxarray, shapely.geometry, pyproj, pandas

R scripts to evaluate model results from PALM with professional weather station data and crowdsourced weather data
- R version 4.1.2
- required packages: ggplot2, lubridate, reshape2, hydroGOF, data.table, ggpubr, dplyr
