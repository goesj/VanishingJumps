### Script to Forecast future Death Rates ######################################
pacman::p_load("MCMCvis","scoringRules","tidyverse")

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 


#load Samples
load(file = file.path(getwd(),"Results/SamplesSp.RData")) # own Model
load(file = file.path(getwd(),"Results/SamplesUS.RData")) # own Model
load(file = file.path(getwd(),"Results/SamplesIt.RData")) # own Model


#load Data
load(file = file.path(getwd(),"Data/CovidData.RData")) # Load Data

#  Forecast Comparison for UK WAR DATA
########### 1. OOS COMPARISON ##################################################
load(file = file.path(getwd(),"Data/UKWARData.RData")) # Load Data
load(file.path(getwd(),"Results/Samples_UKWar_43.RData"))

pacman::p_load("MCMCvis","scoringRules")
S <- 10000
H <- 5
LastYearObs <- 1943

FutureZWAR_AR <- FutureZ(H=H, Mod = "AR", Samples = SamplesAR_War_43, S=S)
FutureZWAR_Liu <- FutureZ(H=H, Mod = "Liu", Samples = SamplesLiuLi_War_43, S=S)
FutureZWAR_MA <- FutureZ(H=H, Mod = "MA", Samples = SamplesMA_War_43, S=S)

LastYearObsInd <- which(1901:2010 == LastYearObs)

FutureLogMatArrayWar_AR <- FutureLogMatArrayWar_Liu <-  FutureLogMatArrayWar_MA <- 
  FutureLogMatArrayWar_LC <- 
  array(data = 0,dim = c(10,H,S),
        dimnames = list("NewAgeInd"=1:10, "FCHorizon"=1:H, "It"=1:S))

for(s in 1:S){
  FutureLogMatArrayWar_AR[,,s] <- log(LambdaMatWar)[,LastYearObsInd]+ #Last obs rate
    t(apply(FutureZWAR_AR$Rates[,,s],1,cumsum)) #Cumsum
  FutureLogMatArrayWar_Liu[,,s] <- log(LambdaMatWar)[,LastYearObsInd]+ #Last obs rate
    t(apply(FutureZWAR_Liu$Rates[,,s],1,cumsum)) #Cumsum
  FutureLogMatArrayWar_MA[,,s] <- log(LambdaMatWar)[,LastYearObsInd]+ #Last obs rate
    t(apply(FutureZWAR_MA$Rates[,,s],1,cumsum)) #Cumsum
}

### 1.1 Get Observed Values
YWar <- LambdaVecWar %>% filter(Year > LastYearObs & Year <= LastYearObs+H) %>% select(Rate) %>% pull() %>% log()
ZWar <- LambdaVecWar %>% filter(Year > LastYearObs & Year <= LastYearObs+H) %>% select(ZVal) %>% pull()


##### Out of Sample
Log_War_Z <- rbind("AR"=logs_sample(y = ZWar, 
                                    dat = (apply(FutureZWAR_AR$Rates,3, c))) %>% sum(),
                   "Liu"= logs_sample(y = ZWar, 
                                      dat = (apply(FutureZWAR_Liu$Rates,3, c))) %>% sum(),
                   "MA" = logs_sample(y = ZWar, 
                                      dat = (apply(FutureZWAR_MA$Rates,3, c))) %>% sum())

CRPS_War_Z <- rbind("AR"=crps_sample(y = ZWar, 
                                     dat = (apply(FutureZWAR_AR$Rates,3, c))) %>% sum(),
                    "Liu"= crps_sample(y = ZWar, 
                                       dat = (apply(FutureZWAR_Liu$Rates,3, c))) %>% sum(),
                    "MA" = crps_sample(y = ZWar, 
                                       dat = (apply(FutureZWAR_MA$Rates,3, c))) %>% sum())

ScoringFrame <- data.frame(#"DSS"=DSS_War,
  "LogS"=Log_War_Z,
  "CRPS"=CRPS_War_Z) %>% t()


# MSE and MAE 
MSEFrame <- data.frame(
  "MSE"=rbind("AR"=(YWar - apply((apply(FutureLogMatArrayWar_AR,3, c)),1,mean))^2 %>% mean(),
              "Liu"=(YWar - apply((apply(FutureLogMatArrayWar_Liu,3, c)),1,mean))^2 %>% mean(),
              "MA"=(YWar - apply((apply(FutureLogMatArrayWar_MA,3, c)),1,mean))^2 %>% mean(),
              "LC"=(YWar - apply((apply(FutureLogMatArrayWar_LC,3, c)),1,mean))^2 %>% mean()
  )*100,
  "MAE"=rbind("AR"=abs(YWar - apply((apply(FutureLogMatArrayWar_AR,3, c)),1,mean)) %>% mean(),
              "Liu"=abs(YWar - apply((apply(FutureLogMatArrayWar_Liu,3, c)),1,mean)) %>% mean(),
              "MA"=abs(YWar - apply((apply(FutureLogMatArrayWar_MA,3, c)),1,mean)) %>% mean(),
              "LC"=abs(YWar - apply((apply(FutureLogMatArrayWar_LC,3, c)),1,mean)) %>% mean()
  )*100
) %>% t()

