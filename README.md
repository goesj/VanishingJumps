# Bayesian mortality modelling with pandemics: a vanishing jump approach
Hereafter, you can find the code to reproduce the results from the paper “Bayesian mortality modelling with pandemics: a vanishing jump approach” by Julius Goes, Karim Barigou and Anne Leucht (https://arxiv.org/pdf/2311.04920.pdf).

## Models 
Our model as well as the Bayesian implementation of the Liu & Li model can be found in the R file `01_NimbleModels.R`
The models are fitted using the **nimble** R package 

## Data 
Within the `Data` folder the following can be found: Data on deaths and exposure categorized by age groups and time periods for multiple countries. 
The data was downloaded from both the Human Mortality Database (HMD) (https://www.mortality.org/) as well as Eurostat (https://ec.europa.eu/eurostat/data/database). 

## Results
Within the `Results` folder the following can be found: 
  * Converged Samples from the posterior of both the Liu-Li as well as own model for all countries.
  * Estimates of the WAIC and LOO-CV of both the Liu-Li as well as own model for all countries.

## R Code
There are multiple R files each requiring different packages to run.  

* `01_NimbleModels.R`: Includes all self-written models that were used for the analysis.
* `02_Functions.R`: Helper script that loads multiple self-written functions needed for either of the R scripts. 
* `03_LoadingData.R`: Script whichs loads the data into R, transforms it into long format as well as calculates mortality improvement rates for each country. Cleaned data for analysis is then saved as **CovidData.RData** and **UKWARData.RData** respectively. Both can be found in the `Data` folder.
* `04_UKWarEstimation.R`: Code for the analysis of the England and Wales data including calculation of WAIC and LOO-CV (section 5.1 of the paper). Results of the analysis was saved within the `Results` folder. 
* `05_EstimationCovidJumps.R`: Code for the analysis of data during the COVID pandemic for the US, Spain and Italy including calculation of WAIC and LOO-CV (section 5.2 of the paper). Results of the analysis was saved within the `Results` folder.
* `06_VisualisationOfParameters`: Code for recreation of all plots within the paper. 
* `07_Forecasts`: Exemplary code for the generation of forecasts for Spain. In addition, comparison of own forecasts with those of `StMoMo` (section 6.1 of the paper).  

The estimation and analysis code requires the following packages: 
`nimble`, `tidyverse`, `rstan`, `openxlsx`, `reshape2`, `StMoMo`,`loo`,`bayestestR`

In addition, the visualizations part needs 
`cowplot`, `ggdist`, `gghighlight`


