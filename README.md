# VanishingJumps
Code for Paper of pandemic jump modelling

## Models 
Own, as well as a Bayesian implementation of the Liu,Li(2015) models can be found `01_NimbleModels.R`
The models are fitted using **nimble** 

## Data 
Within the `Data` folder the following can be found: Deaths and Population data by Age, and time for multiple countries. The data was downloaded from both the Human Mortality Database ([HMD](https://www.mortality.org/)) as well as Eurostat (https://ec.europa.eu/eurostat/data/database). 

## Results
Within the `Results` folder the following can be found: 
  * Converged Samples from the posterior of both the Liu-Li as well as own model for all Countries.
  * Estimates of the WAIC and LOO-CV of both the Liu-Li as well as own model for all Countries.

## R Code
There are multiple R files each requiring different packages to run.  

* `01_NimbleModels.R` Includes all self written models that were used for the analysis.
* `02_Functions.R` Helper script that loads multiple self written functions needed for either of the R scripts. 
* `03_LoadingData.R` Script whichs loads the data into R, transforms it into long format as well as calculate the Mortality improvement rates for each Country. Finished Data for analysis is then saved as *CovidData.RData* and *UKWARData.RData* which can also be found in the `Data` folder.
* `04_UKWarEstimation.R` Code for the analysis of the England and Wales data including calculation of WAIC and LOO-CV. Results of the analysis was saved within the `Results` folder.
* `05_EstimationCovidJumps.R` Code for the analysis of data during the COVID pandemic for the US, Spain and Italy including calculation of WAIC and LOO-CV. Results of the analysis was saved within the `Results` folder.
* `06_VisualisationOfParameters` Code for recreation of all plots within the paper. 
* `07_Forecasts` Exemplary code for the generation of forecasts for Spain. In addition, comparison of own with the `StMoMo` forecasts.  

The estimation and analysis code requires the following packages: 
`nimble`, `tidyverse`, `rstan`, `openxlsx`, `reshape2`, `StMoMo`

In addition, the visualizations part needs 
`cowplot`, `ggdist`, `gghighlight`


