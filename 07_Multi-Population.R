pacman::p_load("nimble","nimbleHMC","tidyverse")

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 

load(file = file.path(getwd(),"Data/CovidData_NewData_Starting1990.RData")) # Load Data


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
                                      "drift",#"b1","b2",
                                      "J_t_Mat","Y_t_Mat","a",
                                      "muY","sdY",
                                      "mu"),enableWAIC = TRUE) #to compare Models


params_to_remove <-  c("sdY",
                       "muY", 
                       "p", 
                       "k_Mat", 
                       "Y_t_Mat",
                       "a","b1","b2",
                       "rho12","rho13","rho23",
                       "muY_Glob", "sdY_Glob",
                       "a","b1","b2")

cOwnMod$removeSamplers(params_to_remove)


cOwnMod$addSampler(target = c("a[1:3]"), type="AF_slice")
cOwnMod$addSampler(target = c("rho12","rho13","rho23"), type="AF_slice")
cOwnMod$addSampler(target = c("sdY[1:3]"), type="AF_slice")
cOwnMod$addSampler(target = c("muY[1:3]"), type="AF_slice")

cOwnMod$addDefaultSampler(node = c("p"), useConjugacy = TRUE)

nimbleHMC::addHMC(cOwnMod, target = c("k_Mat[2:32,1:3]",
                                      "b1[1:10,1:3]","b2[1:10, 1:3]"), replace = TRUE)

for(c in 3:NYear){
  cOwnMod$addSampler(target = paste0("Y_t_Mat[",c,",1]"), type = "slice")
  cOwnMod$addSampler(target = paste0("Y_t_Mat[",c,",2]"), type = "slice")
  cOwnMod$addSampler(target = paste0("Y_t_Mat[",c,",3]"), type = "slice")
}

bOwn <- buildMCMC(cOwnMod)
comOwn <- compileNimble(OwnMod_Joint,
                        bOwn)

SamplesOwnJoint_AR <- runMCMC(comOwn$bOwn, 
                           niter  = 17500,
                           thin=10,
                           nburnin = 7500, 
                           nchains = 2)


###### Estimation of Joint MA model ###########################################
######## 2. MA1 Model ##########################################################
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
                                      "drift",#"b1","b2",
                                      "J_t_Mat","Y_t_Mat","b",
                                      "muY","sdY",
                                      #"muY_Glob","sdY_Glob",
                                      #"mu_b_Glob","sd_b_Glob",
                                      "mu",
                                      "rho12","rho13","rho23"
                         ),
                         enableWAIC = TRUE) #to compare Models


params_to_remove <-  c("sdY",
                       "muY", "p", 
                       "k_Mat", 
                       "Y_t_Mat",
                       "b","b1","b2",
                       "b1","b2",
                       "rho12","rho13","rho23"
)

cOwnMod$removeSamplers(params_to_remove)

cOwnMod$addSampler(target = c("sdY[1:3]"), type="AF_slice")
cOwnMod$addSampler(target = c("muY[1:3]"), type="AF_slice")
cOwnMod$addSampler(target = c("b[1:3]"), type="AF_slice")
cOwnMod$addSampler(target = c("rho12","rho13","rho23"), type="AF_slice")

cOwnMod$addDefaultSampler(node = c("p"), useConjugacy = TRUE)

nimbleHMC::addHMC(cOwnMod, target = c("k_Mat[2:32,1:3]",
                                      "b1[1:10,1:3]","b2[1:10, 1:3]"), replace = TRUE)


for(c in 3:NYear){
  cOwnMod$addSampler(target = paste0("Y_t_Mat[",c,",1]"), type = "slice")
  cOwnMod$addSampler(target = paste0("Y_t_Mat[",c,",2]"), type = "slice")
  cOwnMod$addSampler(target = paste0("Y_t_Mat[",c,",3]"), type = "slice")
}

bOwn <- buildMCMC(cOwnMod)
comOwn <- compileNimble(OwnMod_Joint,
                        bOwn)

SamplesOwnJoint_MA <- runMCMC(comOwn$bOwn, 
                           niter  = 17500,
                           thin=10,
                           nburnin = 7500, 
                           nchains = 2)

save(SamplesOwnJoint_AR,
     SamplesOwnJoint_MA, 
     file = file.path(getwd(),"Results/SamplesMultiPop.RData"))



############ 4. Calculate In-Sample Metrics ####################################
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
