## Estimation of Own and Liu,Li Model Covid Jumps #############################
pacman::p_load("nimble","nimbleHMC","tidyverse")

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 

load(file = file.path(getwd(),"Data/CovidData.RData")) # Load Data

############# 1. UNITED STATES #################################################
###### 1.1 Put Data into correct form for NIMBLE ###############################
#Create Death Matrix
NAge <- nrow(ZMatUS)
NYear <- ncol(ZMatUS)

NimbleConstUS <- list("N_AgeGroups"=NAge,
                       "N_Year"=NYear)

NimbleDataUS <- list("ZMat"=ZMatUS) 

############# 1.2 Liu-Li Model #################################################
LiuLi_US <-  nimbleModel(code=LiuLi_Model, 
                           constants = NimbleConstUS, 
                           data=NimbleDataUS)

cLiuLi_US <- configureMCMC(LiuLi_US, 
                             print = TRUE, useConjugacy = TRUE,
                             monitors = c("k",
                                          "sigma_eps",
                                          "beta","betaJump", "N_t","p",
                                          "sigma_time",
                                          "drift",#"b1","b2",
                                          "muY","sdY","mu","Y_t"
                                          
                             ),
                             enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
                       paste0("b2[",1:NAge,"]"),
                       paste0("Y_t[",1:(NYear),"]"),
                       paste0("k[",1:NYear,"]"),
                       "muY","sdY")

cLiuLi_US$removeSamplers(params_to_remove)
cLiuLi_US$addSampler(target = 
                       c(paste0("Y_t[3:",NYear+1,"]"),"muY","sdY"), 
                      type ="AF_slice")

nimbleHMC::addHMC(cLiuLi_US,
                  target = c("k[2:32]",
                             "b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.80, maxTreeDepth = 10))


bLiuLi_US <- buildMCMC(cLiuLi_US)
comLiuLi_US <- compileNimble(LiuLi_US,
                               bLiuLi_US)


SamplesLiuLi_US <- runMCMC(comLiuLi_US$bLiuLi_US, 
                          niter = 17500, #Increase Sample size a bit
                          thin=10,
                          nburnin = 7500, 
                          nchains = 2)


SummaryOutput(SamplesLiuLi_US, 
              params=c("beta","betaJump",
                       "drift","sigma_eps","sigma_time",
                       "p","a","muY","sdY")) %>% 
  print(n=250) 



######## 1.3 Own Model AR    ####################################################

##### 1.3.1 Estimate Model ####################################################
#Using a sum to zero constraint
OwnMod_US <-  nimbleModel(code=OwnModel_AR, 
                         constants = NimbleConstUS, 
                         data=NimbleDataUS, buildDerivs = TRUE)

cOwnMod_US <- configureMCMC(OwnMod_US, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("k","sigma_eps",
                                         "beta","betaJump", "N_t","p","a",
                                         "sigma_time",
                                         "drift",#"b1","b2",
                                         "J","Y_t",
                                         "muY","sdY","a","mu"
                                         #,"meanJ","sdJ"
                            ),
                            enableWAIC = TRUE) #to compare Models

############# 1.3.2 Adjust Samplers ###########################################
params_to_remove <-  c(
  "p","drift",
  paste0("k[",1:NYear,"]"),
  paste0("Y_t[",1:(NYear),"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "a","muY","sdY"
)

cOwnMod_US$removeSamplers(params_to_remove)
cOwnMod_US$addDefaultSampler(nodes=c("drift","p"), 
                             useConjugacy = TRUE)


cOwnMod_US$addSampler(target = "a", type="slice")

### HMC Sampler for betas and kappa
nimbleHMC::addHMC(cOwnMod_US,
                  target = c("k[2:32]",
                             "b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.80, maxTreeDepth = 10))

#AF slice sampler for Y, muY and sdY
cOwnMod_US$addSampler(target = c(paste0("Y_t[3:",NYear,"]"),"muY","sdY"),
                      type ="AF_slice")


bOwnMod_US <- buildMCMC(cOwnMod_US)
comOwnMod_US <- compileNimble(OwnMod_US,
                              bOwnMod_US)


SamplesOwn_US_AR <- runMCMC(comOwnMod_US$bOwnMod_US, 
                            niter = 27500,
                            thin=20,
                            nburnin = 7500, 
                            nchains = 2)

SummaryOutput(SamplesOwn_US_AR, 
              params=c("beta","betaJump","drift",
                       "sigma_time","sigma_eps",
                       "p","muY","a","sdY","Y_t")) %>% 
  print(n=250)  

#### 1.4 Own Model MA ##########################################################
OwnMod_US_MA <-  nimbleModel(code=OwnModel_MA, 
                          constants = NimbleConstUS, 
                          data=NimbleDataUS, buildDerivs = TRUE)

cOwnMod_US_MA <- configureMCMC(OwnMod_US_MA, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("k","sigma_eps",
                                         "beta","betaJump", "N_t","p",
                                         "sigma_time",
                                         "drift",#"b1","b2",
                                         "J","Y_t",
                                         "muY","sdY",
                                         "mu","b"),
                            enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(
  "p","drift",
  paste0("k[",1:NYear,"]"),
  paste0("Y_t[",1:(NYear),"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "b",#"a2",
  "muY","sdY"
)

cOwnMod_US_MA$removeSamplers(params_to_remove)
cOwnMod_US_MA$addDefaultSampler(nodes=c("drift","p"), 
                             useConjugacy = TRUE)

cOwnMod_US_MA$addSampler(target = "b", type="slice")

cOwnMod_US_MA$addSampler(target = c(paste0("Y_t[3:",NYear,"]"),"muY","sdY"),
                      type ="AF_slice")


nimbleHMC::addHMC(cOwnMod_US_MA,
                  target = c("k[2:32]",
                             "b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.80, maxTreeDepth = 10))


bOwnMod_US_MA <- buildMCMC(cOwnMod_US_MA)
comOwnMod_US_MA <- compileNimble(OwnMod_US_MA,
                              bOwnMod_US_MA)


SamplesOwn_US_MA1 <- runMCMC(comOwnMod_US_MA$bOwnMod_US_MA, 
                             niter = 27500,
                             thin=20,
                             nburnin = 7500, 
                             nchains = 2)

SummaryOutput(SamplesOwn_US_MA1, 
              params=c("beta","betaJump","drift",
                       "sigma_time","sigma_eps",
                       "p","muY","b","sdY","Y_t")) %>% 
  print(n=250)  

######### 1.4. Save Results ####################################################
save(SamplesLiuLi_US,
     SamplesOwn_US_AR,
     SamplesOwn_US_MA1,
     file = file.path(getwd(),"Results/SamplesUS.RData"))


####### 1.5 Calcualte WAIC AND LOO #############################################

#load data in case no new samples have been drawn
load(file.path(getwd(),"Results/SamplesUS.RData")) 

LikeMatOwnUS_AR <- LikelihoodMatrixFun(Samples = SamplesOwn_US_AR,
                                       n = length(ZMatUS),
                                       ZMat = ZMatUS)

LikeMatLiuUS <- LikelihoodMatrixFun(Samples = SamplesLiuLi_US,
                                    n = length(ZMatUS),
                                    ZMat = ZMatUS)

LikeMatOwnUS_MA <- LikelihoodMatrixFun(Samples = SamplesOwn_US_MA1,
                                       n = length(ZMatUS),
                                       ZMat = ZMatUS)

WAICOwnUS_AR <- loo::waic(log(LikeMatOwnUS_AR))
WAICOwnUS_MA <- loo::waic(log(LikeMatOwnUS_MA))
WAICLiuUS <- loo::waic(log(LikeMatLiuUS))

LOOOwnUS_AR <- loo::loo(log(LikeMatOwnUS_AR))
LOOOwnUS_MA <- loo::loo(log(LikeMatOwnUS_MA))
LOOLiuUS <- loo::loo(log(LikeMatLiuUS)) # loo IC = -2 elpd_loo

(CompDataFrameUS <- data.frame("Model"=c("Own_AR","Own_MA","Liu-Li"),
                               "LOOCV"= c(LOOOwnUS_AR$estimates[3,1],
                                          LOOOwnUS_MA$estimates[3,1],
                                          LOOLiuUS$estimates[3,1]),
                               "WAIC"=c(WAICOwnUS_AR$estimates[3,1],
                                        WAICOwnUS_MA$estimates[3,1],
                                        WAICLiuUS$estimates[3,1])))



# ####################### 2. Spain #############################################

#### 2.1 Data for NIMBLE #######################################################
NAge <- nrow(ZMatSp)
NYear <- ncol(ZMatSp)

NimbleConstSp <- list("N_AgeGroups"=NAge,
                      "N_Year"=NYear)

NimbleDataSp <- list("ZMat"=ZMatSp) 

############ 2.2 Liu,Li Model###################################################
LiuLi_Sp <-  nimbleModel(code=LiuLi_Model, 
                            constants = NimbleConstSp, 
                            data=NimbleDataSp)

cLiuLi_Sp <- configureMCMC(LiuLi_Sp, 
                           print = TRUE, useConjugacy = TRUE,
                           monitors = c("k",
                                        "sigma_eps",
                                        "beta","betaJump", "N_t","p",
                                        "sigma_time",
                                        "drift",#"b1","b2",
                                        "muY","sdY","mu","Y_t"
                                        
                           ),
                           enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
                       paste0("b2[",1:NAge,"]"),
                       paste0("Y_t[",1:(NYear),"]"),
                       paste0("k[",1:NYear,"]"),
                       "muY","sdY")

cLiuLi_Sp$removeSamplers(params_to_remove)
cLiuLi_Sp$addSampler(target = 
                       c(paste0("Y_t[3:",NYear+1,"]"),"muY","sdY"), 
                     type ="AF_slice")

nimbleHMC::addHMC(cLiuLi_Sp,
                  target = c("k[2:32]",
                             "b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.80, maxTreeDepth = 10))


bLiuLi_Sp <- buildMCMC(cLiuLi_Sp)
comLiuLi_Sp <- compileNimble(LiuLi_Sp,
                             bLiuLi_Sp)


SamplesLiuLi_Sp <- runMCMC(comLiuLi_Sp$bLiuLi_Sp, 
                           niter = 17500, #Increase Sample size a bit
                           thin=10,
                           nburnin = 7500, 
                           nchains = 2)


SummaryOutput(SamplesLiuLi_Sp, 
              params=c("beta","betaJump",
                       "drift","sigma_eps","sigma_time",
                       "p","muY","sdY")) %>% 
  print(n=250) 



#### 2.3 Own Model  ##########################################################
OwnMod_Sp <-  nimbleModel(code=OwnModel_AR, 
                          constants = NimbleConstSp, 
                          data=NimbleDataSp, buildDerivs = TRUE)

cOwnMod_Sp <- configureMCMC(OwnMod_Sp, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("k","sigma_eps",
                                         "beta","betaJump", "N_t","p","a",
                                         "sigma_time",
                                         "drift",#"b1","b2",
                                         "J","Y_t",
                                         "muY","sdY","a","mu"
                                         #,"meanJ","sdJ"
                            ),
                            enableWAIC = TRUE) #to compare Models

############# 1.3.2 Adjust Samplers ###########################################
params_to_remove <-  c(
  "p","drift",
  paste0("k[",1:NYear,"]"),
  paste0("Y_t[",1:(NYear),"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "a","muY","sdY"
)

cOwnMod_Sp$removeSamplers(params_to_remove)
cOwnMod_Sp$addDefaultSampler(nodes=c("drift","p"), 
                             useConjugacy = TRUE)


cOwnMod_Sp$addSampler(target = "a", type="slice")

### HMC Sampler for betas and kappa
nimbleHMC::addHMC(cOwnMod_Sp,
                  target = c("k[2:32]",
                             "b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.80, maxTreeDepth = 10))

#AF slice sampler for Y, muY and sdY
cOwnMod_Sp$addSampler(target = c(paste0("Y_t[3:",NYear,"]"),"muY","sdY"),
                      type ="AF_slice")


bOwnMod_Sp <- buildMCMC(cOwnMod_Sp)
comOwnMod_Sp <- compileNimble(OwnMod_Sp,
                              bOwnMod_Sp)


SamplesOwn_Sp_AR <- runMCMC(comOwnMod_Sp$bOwnMod_Sp, 
                            niter = 27500,
                            thin=20,
                            nburnin = 7500, 
                            nchains = 2)

SummaryOutput(SamplesOwn_Sp_AR, 
              params=c("beta","betaJump","drift",
                       "sigma_time","sigma_eps",
                       "p","muY","a","sdY","Y_t")) %>% 
  print(n=250)  


#### 2.4 Own Model MA ##########################################################
OwnMod_Sp_MA <-  nimbleModel(code=OwnModel_MA, 
                             constants = NimbleConstSp, 
                             data=NimbleDataSp, buildDerivs = TRUE)

cOwnMod_Sp_MA <- configureMCMC(OwnMod_Sp_MA, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("k","sigma_eps",
                                         "beta","betaJump", "N_t","p",
                                         "sigma_time",
                                         "drift",#"b1","b2",
                                         "J","Y_t",
                                         "muY","sdY",
                                         "mu","b"),
                            enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(
  "p","drift",
  paste0("k[",1:NYear,"]"),
  paste0("Y_t[",1:(NYear),"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "b",#"a2",
  "muY","sdY"
)

cOwnMod_Sp_MA$removeSamplers(params_to_remove)
cOwnMod_Sp_MA$addDefaultSampler(nodes=c("drift","p"), 
                             useConjugacy = TRUE)

cOwnMod_Sp_MA$addSampler(target = "b", type="slice")

cOwnMod_Sp_MA$addSampler(target = c(paste0("Y_t[3:",NYear,"]"),"muY","sdY"),
                      type ="AF_slice")


nimbleHMC::addHMC(cOwnMod_Sp_MA,
                  target = c("k[2:32]",
                             "b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.80, maxTreeDepth = 10))


bOwnMod_Sp_MA <- buildMCMC(cOwnMod_Sp_MA)
comOwnMod_Sp_MA <- compileNimble(OwnMod_Sp_MA,
                              bOwnMod_Sp_MA)


SamplesOwn_Sp_MA1 <- runMCMC(comOwnMod_Sp_MA$bOwnMod_Sp_MA, 
                             niter = 27500,
                             thin=20,
                             nburnin = 7500, 
                             nchains = 2)





#### 2.5 Save Results #########################################################
save(SamplesLiuLi_Sp, 
     SamplesOwn_Sp, file= file.path(getwd(),"Results/SamplesSp.RData"))


######### 2.6. Calculate WAIC and LOO-CV ######################################

#load data in case no new samples have been drawn
load(file.path(getwd(),"Results/SamplesSp.RData"))

#### Caclulation of WAIC and loo 
LikeMatOwnSp_AR <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesOwn_Sp_AR),
                                       n = length(ZMatSp),
                                       ZMat = ZMatSp)

LikeMatLiuSp <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesLiuLi_Sp),
                                    n = length(ZMatSp),
                                    ZMat = ZMatSp)

LikeMatOwnSp_MA <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesOwn_Sp_MA1),
                                       n = length(ZMatSp),
                                       ZMat = ZMatSp)

WAICOwnSp_AR <- loo::waic(log(LikeMatOwnSp_AR))
WAICOwnSp_MA <- loo::waic(log(LikeMatOwnSp_MA))
WAICLiuSp <- loo::waic(log(LikeMatLiuSp))

LOOOwnSp_AR <- loo::loo(log(LikeMatOwnSp_AR))
LOOOwnSp_MA <- loo::loo(log(LikeMatOwnSp_MA))
LOOLiuSp <- loo::loo(log(LikeMatLiuSp)) # loo IC = -2 elpd_loo

(CompDataFrameSp <- data.frame("Model"=c("Own_AR","Own_MA","Liu-Li"),
                               "LOOCV"= c(LOOOwnSp_AR$estimates[3,1],
                                          LOOOwnSp_MA$estimates[3,1],
                                          LOOLiuSp$estimates[3,1]),
                               "WAIC"=c(WAICOwnSp_AR$estimates[3,1],
                                        WAICOwnSp_MA$estimates[3,1],
                                        WAICLiuSp$estimates[3,1])))

# ############## 3. Poland #####################################################

#### 3.1 Put data into right form for NIMBLE ##################################
NAge <- nrow(ZMatPl)
NYear <- ncol(ZMatPl)

NimbleConstPl <- list("N_AgeGroups"=NAge,
                      "N_Year"=NYear)

NimbleDataPl <- list("ZMat"=ZMatPl) 


############ 2.2 Liu,Li Model###################################################
LiuLi_Pl <-  nimbleModel(code=LiuLi_Model, 
                         constants = NimbleConstPl, 
                         data=NimbleDataPl)

cLiuLi_Pl <- configureMCMC(LiuLi_Pl, 
                           print = TRUE, useConjugacy = TRUE,
                           monitors = c("k",
                                        "sigma_eps",
                                        "beta","betaJump", "N_t","p",
                                        "sigma_time",
                                        "drift",#"b1","b2",
                                        "muY","sdY","mu","Y_t"
                                        
                           ),
                           enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
                       paste0("b2[",1:NAge,"]"),
                       paste0("Y_t[",1:(NYear),"]"),
                       paste0("k[",1:NYear,"]"),
                       "muY","sdY")

cLiuLi_Pl$removeSamplers(params_to_remove)
cLiuLi_Pl$addSampler(target = 
                       c(paste0("Y_t[3:",NYear+1,"]"),"muY","sdY"), 
                     type ="AF_slice")

nimbleHMC::addHMC(cLiuLi_Pl,
                  target = c("k[2:32]",
                             "b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.80, maxTreeDepth = 10))


bLiuLi_Pl <- buildMCMC(cLiuLi_Pl)
comLiuLi_Pl <- compileNimble(LiuLi_Pl,
                             bLiuLi_Pl)


SamplesLiuLi_Pl <- runMCMC(comLiuLi_Pl$bLiuLi_Pl, 
                           niter = 17500, #Increase Sample size a bit
                           thin=10,
                           nburnin = 7500, 
                           nchains = 2)


SummaryOutput(SamplesLiuLi_Pl, 
              params=c("beta","betaJump",
                       "drift","sigma_eps","sigma_time",
                       "p","muY","sdY")) %>% 
  print(n=250) 



#### 2.3 Own Model  ##########################################################
OwnMod_Pl <-  nimbleModel(code=OwnModel_AR, 
                          constants = NimbleConstPl, 
                          data=NimbleDataPl, buildDerivs = TRUE)

cOwnMod_Pl <- configureMCMC(OwnMod_Pl, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("k","sigma_eps",
                                         "beta","betaJump", "N_t","p","a",
                                         "sigma_time",
                                         "drift",#"b1","b2",
                                         "J","Y_t",
                                         "muY","sdY","a","mu"
                                         #,"meanJ","sdJ"
                            ),
                            enableWAIC = TRUE) #to compare Models

############# 1.3.2 Adjust Samplers ###########################################
params_to_remove <-  c(
  "p","drift",
  paste0("k[",1:NYear,"]"),
  paste0("Y_t[",1:(NYear),"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "a","muY","sdY"
)

cOwnMod_Pl$removeSamplers(params_to_remove)
cOwnMod_Pl$addDefaultSampler(nodes=c("drift","p"), 
                             useConjugacy = TRUE)


cOwnMod_Pl$addSampler(target = "a", type="slice")

### HMC Sampler for betas and kappa
nimbleHMC::addHMC(cOwnMod_Pl,
                  target = c("k[2:32]",
                             "b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.80, maxTreeDepth = 10))

#AF slice sampler for Y, muY and sdY
cOwnMod_Pl$addSampler(target = c(paste0("Y_t[3:",NYear,"]"),"muY","sdY"),
                      type ="AF_slice")


bOwnMod_Pl <- buildMCMC(cOwnMod_Pl)
comOwnMod_Pl <- compileNimble(OwnMod_Pl,
                              bOwnMod_Pl)


SamplesOwn_Pl_AR <- runMCMC(comOwnMod_Pl$bOwnMod_Pl, 
                            niter = 27500,
                            thin=20,
                            nburnin = 7500, 
                            nchains = 2)

SummaryOutput(SamplesOwn_Pl_AR, 
              params=c("beta","betaJump","drift",
                       "sigma_time","sigma_eps",
                       "p","muY","a","sdY","Y_t")) %>% 
  print(n=250)  


#### 2.4 Own Model MA ##########################################################
OwnMod_Pl_MA <-  nimbleModel(code=OwnModel_MA, 
                             constants = NimbleConstPl, 
                             data=NimbleDataPl, buildDerivs = TRUE)

cOwnMod_Pl_MA <- configureMCMC(OwnMod_Pl_MA, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("k","sigma_eps",
                                         "beta","betaJump", "N_t","p",
                                         "sigma_time",
                                         "drift",#"b1","b2",
                                         "J","Y_t",
                                         "muY","sdY",
                                         "mu","b"),
                            enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(
  "p","drift",
  paste0("k[",1:NYear,"]"),
  paste0("Y_t[",1:(NYear),"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "b",#"a2",
  "muY","sdY"
)

cOwnMod_Pl_MA$removeSamplers(params_to_remove)
cOwnMod_Pl_MA$addDefaultSampler(nodes=c("drift","p"), 
                             useConjugacy = TRUE)

cOwnMod_Pl_MA$addSampler(target = "b", type="slice")

cOwnMod_US$addSampler(target = c(paste0("Y_t[3:",NYear,"]"),"muY","sdY"),
                      type ="AF_slice")


nimbleHMC::addHMC(cOwnMod_Pl_MA,
                  target = c("k[2:32]",
                             "b1[1:10]","b2[1:10]"),
                  type = "NUTS", replace = TRUE,
                  control = list(delta = 0.80, maxTreeDepth = 10))


bOwnMod_Pl_MA <- buildMCMC(cOwnMod_Pl_MA)
comOwnMod_Pl_MA <- compileNimble(OwnMod_Pl_MA,
                              bOwnMod_Pl_MA)


SamplesOwn_Pl_MA1 <- runMCMC(comOwnMod_Pl_MA$bOwnMod_Pl_MA, 
                             niter = 27500,
                             thin=20,
                             nburnin = 7500, 
                             nchains = 2)

############ 3.4 Save Results ##################################################

save(SamplesLiuLi_Pl, 
     SamplesOwn_Pl_AR,
     SamplesOwn_Pl_MA1, file= file.path(getwd(),"Results/SamplesPl.RData"))


########## 3.5. Calculation of WAIC AND LOO ####################################

#load data in case no new samples have been drawn
load(file.path(getwd(),"Results/SamplesPl.RData"))

LikeMatOwnPl_AR <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesOwn_Pl_AR),
                                       n = length(ZMatPl),
                                       ZMat = ZMatPl)

LikeMatLiuPl <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesLiuLi_Pl),
                                    n = length(ZMatPl),
                                    ZMat = ZMatPl)

LikeMatOwnPl_MA <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesOwn_Pl_MA1),
                                       n = length(ZMatPl),
                                       ZMat = ZMatPl)

WAICOwnPl_AR <- loo::waic(log(LikeMatOwnPl_AR))
WAICOwnPl_MA <- loo::waic(log(LikeMatOwnPl_MA))
WAICLiuPl <- loo::waic(log(LikeMatLiuPl))

LOOOwnPl_AR <- loo::loo(log(LikeMatOwnPl_AR))
LOOOwnPl_MA <- loo::loo(log(LikeMatOwnPl_MA))
LOOLiuPl <- loo::loo(log(LikeMatLiuPl)) # loo IC = -2 elpd_loo

(CompDataFramePl <- data.frame("Model"=c("Own_AR","Own_MA","Liu-Li"),
                               "LOOCV"= c(LOOOwnPl_AR$estimates[3,1],
                                          LOOOwnPl_MA$estimates[3,1],
                                          LOOLiuPl$estimates[3,1]),
                               "WAIC"=c(WAICOwnPl_AR$estimates[3,1],
                                        WAICOwnPl_MA$estimates[3,1],
                                        WAICLiuPl$estimates[3,1])))



