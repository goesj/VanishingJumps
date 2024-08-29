## Load necessary files and libraries ###
library(nimble);library(tidyverse); library(rstan);library(loo)

source("01_NimbleModels.R")
source("02_Functions.R") #needs rstan

############## 1. LOADING DATA #################################################
load(file = file.path(getwd(),"Data/UKWARData.RData")) # Load Data


############## 2. Estimating Parameters ########################################

#Estimating Parameters of both the Liu-Li model, as well as our Model 
#Route II Estimation Method in use (Mortality Improvements)

LastYearObs <- 1943 #Adjust here for different scenario
LastYearObsInd <- which(1901:2010 == LastYearObs)

#After World WAR I, 
ZMatWar_NOWWII <- ZMatWar[,1:LastYearObsInd]
NimbleConst_UK <- list("N_AgeGroups"=nrow(ZMatWar_NOWWII),
                       "N_Year"=ncol(ZMatWar_NOWWII))

NimbleData_UK <- list("ZMat"=ZMatWar_NOWWII) 

xlength <- nrow(ZMatWar_NOWWII)
klength <- ncol(ZMatWar_NOWWII)+1

#### ATTENTION ! CHANGE PRIORS IN NIMBLE CODE TO THE ONES SHOWN IN THE PAPER ###
### 2.1.Own Model AR  ##########################################################

### !!! Change constraint, that N_T = 0 !!!!!!!!!! 

OwnModel_AR_War <- nimbleModel(code=OwnModel_AR, 
                            constants = NimbleConst_UK, 
                            data=NimbleData_UK, 
                            buildDerivs = TRUE)

cOwnModel_AR <- configureMCMC(OwnModel_AR_War,monitors = c("sigma_eps",
                                                          "beta","betaJump",
                                                          "N_t","p",
                                                          "sigma_time",
                                                          "sdY","muY",
                                                          "drift",
                                                          "a",
                                                          "Y_t","J",
                                                          "b1","b2",
                                                          "mu","k"),
                                enableWAIC = TRUE, useConjugacy = FALSE)

params_to_remove <-  c(
  "muY","sdY","sigma_time",
  "a","p",
  paste0("b1[",1:xlength,"]"),
  paste0("b2[",1:xlength,"]"),
  paste0("Y_t[",1:(klength),"]")
)

cOwnModel_AR$removeSamplers(params_to_remove)
cOwnModel_AR$addSampler(target = "muY",type="slice")
cOwnModel_AR$addSampler(target = "a",type="slice")
cOwnModel_AR$addSampler(target = "sdY",type="slice")
cOwnModel_AR$addDefaultSampler(nodes="p", useConjugacy = TRUE)
cOwnModel_AR$addSampler(target = "sigma_time",type="slice")

for(j in 3:(klength)){
  cOwnModel_AR$addSampler(target=paste0("Y_t[",j,"]"),
                            type="slice")
}

nimbleHMC::addHMC(cOwnModel_AR, 
                  target = c("k[2:43]",
                             "b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.82, maxTreeDepth = 11))


bOwnModel_AR <- buildMCMC(cOwnModel_AR)
comOwnModel_AR <- compileNimble(OwnModel_AR_War, bOwnModel_AR)


SamplesAR_War_43 <- runMCMC(comOwnModel_AR$bOwnModel_AR, 
                                niter = 6500,
                                thin=2,
                                nburnin = 2500, 
                                nchains = 2) 

SummaryOutput(SamplesAR_War_43,
              params=c("betaJump","sigma_time",
                       "beta","drift","p","sdY","muY",
                       "a","N_t")) %>% print(n=400) #work

##### 2.2 Liu-Li ###############################################################
LiuLiModel_War <- nimbleModel(code=LiuLi_Model, 
                            constants = NimbleConst_UK, 
                            data=NimbleData_UK, 
                            buildDerivs = TRUE)

cLiuLiModel <- configureMCMC(LiuLiModel_War,monitors = c("sigma_eps",
                                                      "beta","betaJump",
                                                      "N_t","p",
                                                      "sigma_time","sdY","muY",
                                                      "drift",
                                                      "k",
                                                      "Y_t","J",
                                                      "mu"), 
                            enableWAIC = TRUE, useConjugacy = TRUE)


params_to_remove <-  c(
  paste0("b1[",1:10,"]"),
  paste0("b2[",1:10,"]"),
  "sdY","sigma_time","muY","p"
)

cLiuLiModel$removeSamplers(params_to_remove)

nimbleHMC::addHMC(cLiuLiModel, 
                  target = c("b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.82, maxTreeDepth = 11))

cLiuLiModel$addSampler(target = "muY",type="slice")
cLiuLiModel$addSampler(target = "sdY",type="slice")
cLiuLiModel$addSampler(target = "sigma_time",type="slice")
cLiuLiModel$addDefaultSampler(nodes = "p",print = TRUE)


bLiuLiModel <- buildMCMC(cLiuLiModel)
comLiuLiModel <- compileNimble(LiuLiModel_War,bLiuLiModel)


SamplesLiuLi_War_43 <- runMCMC(comLiuLiModel$bLiuLiModel, 
                                   niter = 8000,
                                   thin=4,
                                   nburnin = 4000, 
                                   nchains = 2) #set seed for reproducibility

SummaryOutput(SamplesLiuLi_War_43,
              params=c("betaJump","sigma_time",
                       "beta","drift","p","sdY","muY",
                       "a","N_t")) %>% print(n=400) #work

##### 2.3 Own Model MA #########################################################
OwnModel_MA_War <- nimbleModel(code=OwnModel_MA, 
                            constants = NimbleConst_UK, 
                            data=NimbleData_UK, 
                            buildDerivs = TRUE)

cOwnModel_MA <- configureMCMC(OwnModel_MA_War,monitors = c("sigma_eps",
                                                          "beta","betaJump",
                                                          "N_t","p",
                                                          "sigma_time",
                                                          "sdY","muY",
                                                          "drift",
                                                          "b",
                                                          "Y_t","J",
                                                          "b1","b2",
                                                          "mu","k"),
                              enableWAIC = TRUE, useConjugacy = FALSE)

#For Change of Samplers
params_to_remove <-  c(
  "muY","sdY","sigma_time",
  "b","p",
  paste0("b1[",1:xlength,"]"),
  paste0("b2[",1:xlength,"]"),
  paste0("Y_t[",1:(klength),"]")
)

cOwnModel_MA$removeSamplers(params_to_remove)
cOwnModel_MA$addSampler(target = "muY",type="slice")
cOwnModel_MA$addSampler(target = "b",type="slice")
cOwnModel_MA$addSampler(target = "sdY",type="slice")
cOwnModel_MA$addDefaultSampler(nodes="p", useConjugacy = TRUE)
cDiffModel_Own$addSampler(target = "sigma_time",type="slice")
# cDiffModel_Own$addDefaultSampler(nodes=paste0("k[",2:(klength -1),"]"),
#                                  useConjugacy = TRUE, print = FALSE)

for(j in 3:(klength)){
  cOwnModel_MA$addSampler(target=paste0("Y_t[",j,"]"),
                            type="slice")
}

nimbleHMC::addHMC(cOwnModel_MA, 
                  target = c("k[2:43]",
                             "b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.82, maxTreeDepth = 11))



bOwnModel_MA <- buildMCMC(cOwnModel_MA)
comOwnModel_MA <- compileNimble(OwnModel_MA_War,bOwnModel_MA)


SamplesMA_War_43 <- runMCMC(comOwnModel_MA$bOwnModel_MA, 
                                 niter = 6500,
                                 thin=2,
                                 nburnin = 2500, 
                                 nchains = 2) #set seed for reproducibility


SummaryOutput(Samples_OwnMod_War_MA,
              params=c("sigma_eps","betaJump","sigma_time",
                       "beta","drift","p","sdY","muY",
                       "b","N_t")) %>% print(n=400) #work


######### 2.4 Lee Carter Model (only for scenario II)) #########################
LC_War <-  nimbleModel(code = LC_Diff, 
                       constants = NimbleConst_UK, 
                       data=NimbleData_UK)

cLC_War <- configureMCMC(LC_War, 
                         print = TRUE, useConjugacy = TRUE,
                         monitors = c("k","sigma_eps",
                                      "beta",
                                      "sigma_time",
                                      "drift",#"b1",
                                      "mu"
                         ),
                         enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(
  paste0("b1[",1:10,"]")
)
cLC_War$removeSamplers(params_to_remove)
cLC_War$addSampler(target= c("b1[1:10]"),
                   type="AF_slice")

bLC_War <- buildMCMC(cLC_War)
comLC_War <- compileNimble(LC_War,
                           bLC_War)


SamplesLC_War_80<- runMCMC(comLC_War$bLC_War, 
                           niter = 4000,
                           thin=1,
                           nburnin = 2000, 
                           nchains = 2)


############# 3. Save Results ##################################################
save(SamplesAR_War_43,
     SamplesLiuLi_War_43,
     SamplesMA_War_43,
     file = file.path(getwd(),"Results/Samples_UKWar_43.RData"))

