pacman::p_load("nimble","nimbleHMC","tidyverse","abind","loo")

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 

load(file = file.path(getwd(),"Data/CovidData.RData")) # Load Data


######## 1. Data Preparation ###################################################
NAge <- nrow(ZMatSp)
NYear <- ncol(ZMatPl)
NCounty <- 3

NimbleConst <- list("N_AgeGroups"=NAge,
                    "N_Year"=NYear,
                    "N_Country"=3)

ZMatArray <- abind::abind(ZMatUS, ZMatSp, ZMatPl, along=3)

NimbleData<- list("ZMat"=ZMatArray) 

##### 2. Estimation of AR1 Model ###############################################
nimbleOptions(doADerrorTraps=FALSE)
OwnMod_Joint <-  nimbleModel(code=MultiPop_AR_GlobalN, 
                             constants = NimbleConst, 
                             data=NimbleData, buildDerivs = TRUE)

cOwnMod <- configureMCMC(OwnMod_Joint, 
                         print = TRUE, useConjugacy = FALSE,
                         monitors = c("k_Mat","sigma_eps",
                                      "beta","betaJump", 
                                      "N_t",
                                      "p",
                                      "sigma_time",
                                      "drift",
                                      "J_t_Mat","Y_t_Mat","a",
                                      "muY","sdY",
                                      "mu")) 
params_to_remove <-  c("sdY",
                       "muY", 
                       "p", 
                       "k_Mat", 
                       "a","b1","b2",
                       "rho12","rho13","rho23")

cOwnMod$removeSamplers(params_to_remove)


cOwnMod$addSampler(target = c("a[1:3]"), type="AF_slice")
cOwnMod$addSampler(target = c("rho12","rho13","rho23"), type="AF_slice")
cOwnMod$addSampler(target = c("sdY[1:3]"), type="AF_slice")
cOwnMod$addSampler(target = c("muY[1:3]"), type="AF_slice")

nimbleHMC::addHMC(cOwnMod, target = c("k_Mat[2:32, 1:3]",
                                      "b1[1:10, 1:3]",
                                      "b2[1:10, 1:3]","p"), replace = TRUE)


bOwn <- buildMCMC(cOwnMod)
comOwn <- compileNimble(OwnMod_Joint,
                        bOwn)

SamplesOwnJoint_AR <- runMCMC(comOwn$bOwn, 
                           niter  = 17500,
                           thin=10,
                           nburnin = 7500, 
                           nchains = 2)

SummaryOutput(SamplesOwnJoint_AR, 
              params=c("muY_Glob","sdY_Glob","mu_b_Glob",
                       "sd_b_Glob","drift","p","sigma_eps",
                       "muY","sdY","sigma_time","b","beta","betaJump",
                       "N_t")) %>% 
  print(n=400) 


###### Estimation of Joint MA model ###########################################
######## 3. MA1 Model ##########################################################
OwnMod_Joint <-  nimbleModel(code=MultiPop_MA_GlobalN, 
                             constants = NimbleConst, 
                             data=NimbleData, buildDerivs = TRUE)

cOwnMod <- configureMCMC(OwnMod_Joint, 
                         print = TRUE, useConjugacy = FALSE,
                         monitors = c("k_Mat","sigma_eps",
                                      "beta","betaJump", 
                                      "N_t",
                                      "p",
                                      "sigma_time",
                                      "drift",
                                      "J_t_Mat","Y_t_Mat","b",
                                      "muY","sdY",
                                      "mu",
                                      "rho12","rho13","rho23"
                         ),
                         enableWAIC = TRUE) #to compare Models


params_to_remove <-  c("sdY",
                       "muY", "p", 
                       "k_Mat", 
                       "b","b1","b2",
                       "rho12","rho13","rho23")

cOwnMod$removeSampler(params_to_remove)

cOwnMod$addSampler(target = c("b[1:3]"), type="AF_slice")
cOwnMod$addSampler(target = c("rho12","rho13","rho23"), type="AF_slice")
cOwnMod$addSampler(target = c("sdY[1:3]"), type="AF_slice")
cOwnMod$addSampler(target = c("muY[1:3]"), type="AF_slice")


nimbleHMC::addHMC(cOwnMod, target = c("k_Mat[2:32,1:3]",
                                      "b1[1:10,1:3]","b2[1:10, 1:3]","P"), replace = TRUE)


bOwn <- buildMCMC(cOwnMod)
comOwn <- compileNimble(OwnMod_Joint,
                        bOwn)

SamplesOwnJoint_MA <- runMCMC(comOwn$bOwn, 
                           niter  = 17500,
                           thin=10,
                           nburnin = 7500, 
                           nchains = 2)

SummaryOutput(SamplesOwnJoint_MA, 
              params=c("muY_Glob","sdY_Glob","mu_b_Glob",
                       "sd_b_Glob","drift","p","sigma_eps",
                       "muY","sdY","sigma_time","b","beta","betaJump",
                       "N_t")) %>% 
  print(n=400) 


save(SamplesOwnJoint_AR,
     SamplesOwnJoint_MA, 
     file = file.path(getwd(),"Results/SamplesMultiPop.RData"))



############ 4. Calculate In-Sample Metrics ####################################
load(file = file.path(getwd(),"Results/SamplesMultiPop.RData"))

LikeMatOwn_Joint_AR <- LikelihoodMatrixFunJoint(Samples = SamplesOwnJoint_AR,
                                                n = length(ZMatArray), 
                                                ZArray = ZMatArray,
                                                NCountry = 3)

WAIC_Joint <- loo::waic(log(LikeMatOwn_Joint_AR))
loo_Joint_AR <- loo::loo(log(LikeMatOwn_Joint_AR))


LikeMatOwn_Joint_MA <- LikelihoodMatrixFunJoint(Samples = SamplesOwnJoint_MA,
                                                n = length(ZMatArray), 
                                                ZArray = ZMatArray,
                                                NCountry = 3)

WAIC_Joint_MA <- loo::waic(log(LikeMatOwn_Joint_MA))
loo_Joint_MA <- loo::loo(log(LikeMatOwn_Joint_MA))
