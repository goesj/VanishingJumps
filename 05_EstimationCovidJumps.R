## Estimation of Own and Liu,Li Model Covid Jumps #############################
library(nimble);library(tidyverse); library(loo);library(rstan)

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
                       "muY","sdY"
)

cLiuLi_US$removeSamplers(params_to_remove)
cLiuLi_US$addSampler(nodes = "muY",type="slice")
cLiuLi_US$addSampler(nodes = "sdY",type="slice")
cLiuLi_US$addSampler(target= paste0("b1[1:",NAge,"]"),
                       type="AF_slice") 
cLiuLi_US$addSampler(target= paste0("b2[1:",NAge,"]"),
                       type="AF_slice")


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



######## 1.3 Own Model     ####################################################

##### 1.3.1 Estimate Model ####################################################
#Using a sum to zero constraint
OwnMod_US <-  nimbleModel(code=OwnModel, 
                         constants = NimbleConstUS, 
                         data=NimbleDataUS)

cOwnMod_US <- configureMCMC(OwnMod_US, 
                           print = TRUE, 
                           useConjugacy = FALSE, #to save time configuring
                           monitors = c("k","sigma_eps",
                                        "beta","betaJump", "N_t","p","a",
                                        "sigma_time",
                                        "drift",
                                         "J","Y_t",
                                        "muY","sdY","a","mu"
                           ),
                           enableWAIC = TRUE) #to compare Models

############# 1.3.2 Adjust Samplers ###########################################
params_to_remove <-  c(
  "p","drift",
  paste0("k[",1:NYear,"]"),
  paste0("Y_t[",3:(NYear),"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "a","muY","sdY"
)

cOwnMod_US$removeSamplers(params_to_remove)

#Gibbs Sampler
cOwnMod_US$addDefaultSampler(nodes=c("drift","p"), 
                             useConjugacy = TRUE)
cOwnMod_US$addDefaultSampler(nodes=paste0("k[",2:NYear,"]"),
                             useConjugacy = TRUE, print = FALSE)

#Slice Sampler
cOwnMod_US$addSampler(nodes = "muY",type="slice")
cOwnMod_US$addSampler(nodes = "a",type="slice")
cOwnMod_US$addSampler(nodes = "sdY",type="slice")

# slice sampler for Y_t
for(j in 3:(NYear)){ 
  cOwnMod_US$addSampler(target=paste0("Y_t[",j,"]"),
                          type="slice")
}

#AF_Slice Sampler
cOwnMod_US$addSampler(target= paste0("b1[1:",NAge,"]"),
                          type="AF_slice") 
cOwnMod_US$addSampler(target= paste0("b2[1:",NAge,"]"),
                          type="AF_slice")


bOwnMod_US <- buildMCMC(cOwnMod_US)
comOwnMod_US <- compileNimble(OwnMod_US,
                            bOwnMod_US)


SamplesOwn_US <- runMCMC(comOwnMod_US$bOwnMod_US, 
                          niter = 15000,
                          thin=10,
                          nburnin = 5000, 
                          nchains = 2)

SummaryOutput(SamplesOwn_US, 
              params=c("beta","betaJump",
                       "drift","sigma_eps","sigma_time",
                       "p","a","muY","sdY")) %>% 
  print(n=250) 


######### 1.4. Save Results ####################################################
save(SamplesLiuLi_US,
     SamplesOwn_US,
     file = file.path(getwd(),"Results/SamplesUS.RData"))


####### 1.5 Calcualte WAIC AND LOO #############################################

#load data in case no new samples have been drawn
load(file.path(getwd(),"Results/SamplesUS.RData")) 

LikeMatOwn_US <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesOwn_US),
                                    n = length(ZMatUS),
                                    ZMat = ZMatUS)

LikeMatLiuLi_US <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesLiuLi_US),
                                    n = length(ZMatUS),
                                    ZMat = ZMatUS)

WAICOwnUS <- loo::waic(log(LikeMatOwn_US))
WAICLiuUS <- loo::waic(log(LikeMatLiuLi_US))

LOOOwnUS <- loo::loo(log(LikeMatOwn_US))
LOOLiuUS <- loo::loo(log(LikeMatLiuLi_US)) # loo IC = -2 elpd_loo

CompDataFrameUS <- data.frame("Model"=c("Own","Liu-Li"),
                              "LOOCV"= c(LOOOwnUS$estimates[3,1],
                                         LOOLiuUS$estimates[3,1]),
                              "WAIC"=c(WAICOwnUS$estimates[3,1],
                                       WAICLiuUS$estimates[3,1]))

save(CompDataFrameUS,
     file = file.path(getwd(),"Results/WAIC_US.RData"))



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
                                "drift","muY","sdY","mu","Y_t"
                              ),
                              enableWAIC = TRUE) 

params_to_remove <-  c(
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "muY","sdY"
)

cLiuLi_Sp$removeSamplers(params_to_remove)
cLiuLi_Sp$addSampler(nodes = "muY",type="slice")
cLiuLi_Sp$addSampler(nodes = "sdY",type="slice")

cLiuLi_Sp$addSampler(target= paste0("b1[1:",NAge,"]"),
                        type="AF_slice") 
cLiuLi_Sp$addSampler(target= paste0("b2[1:",NAge,"]"),
                        type="AF_slice")

bLiuLi_Sp <- buildMCMC(cLiuLi_Sp)
comLiuLi_Sp <- compileNimble(LiuLi_Sp,
                             bLiuLi_Sp)

SamplesLiuLi_Sp <- runMCMC(comLiuLi_Sp$bLiuLi_Sp, 
                           niter = 15000,
                           thin=10,
                           nburnin = 5000, 
                           nchains = 2)

SummaryOutput(SamplesLiuLi_Sp, 
              params=c("beta","betaJump",
                       "drift","sigma_eps","sigma_time",
                       "p","a","muY","sdY")) %>% 
  print(n=250) 


#### 2.3 Own Model  ##########################################################
OwnMod_Sp <-  nimbleModel(code=OwnModel, 
                          constants = NimbleConstSp, 
                          data=NimbleDataSp)

cOwnMod_Sp <- configureMCMC(OwnMod_Sp, 
                                print = TRUE, useConjugacy = FALSE,
                                monitors = c("k",
                                             "sigma_eps",
                                             "b1","b2",
                                             "sdY","muY",
                                             "N_t","p",
                                             "sigma_time",
                                             "drift","a",
                                             "beta","betaJump",
                                             "mu",#"sigma_squared" #for PPC and WAIC
                                             "Y_t","J"),
                                enableWAIC = TRUE) 

############# 2.2.1 Adjust Samplers ############################################
params_to_remove <-  c(
  "p","drift",
  paste0("k[",1:NYear,"]"),
  paste0("Y_t[",3:(NYear),"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "a","muY","sdY"
)

cOwnMod_Sp$removeSamplers(params_to_remove)

#Gibbs
cOwnMod_Sp$addDefaultSampler(nodes=c("drift","p"), 
                             useConjugacy = TRUE)
cOwnMod_Sp$addDefaultSampler(nodes=paste0("k[",2:NYear,"]"),
                             useConjugacy = TRUE, print = FALSE)

#Slice
cOwnMod_Sp$addSampler(nodes = "muY",type="slice")
cOwnMod_Sp$addSampler(nodes = "a",type="slice")
cOwnMod_Sp$addSampler(nodes = "sdY",type="slice")
# slice sampler for Y_t
for(j in 3:(NYear)){ #one time period longer due to differencing
  cOwnMod_Sp$addSampler(target=paste0("Y_t[",j,"]"),
                        type="slice")
}

#AF_Slice
cOwnMod_Sp$addSampler(target= paste0("b1[1:",NAge,"]"),
                      type="AF_slice") 
cOwnMod_Sp$addSampler(target= paste0("b2[1:",NAge,"]"),
                      type="AF_slice")

bOwnMod_Sp <- buildMCMC(cOwnMod_Sp)
comOwnMod_Sp <- compileNimble(OwnMod_Sp,
                              bOwnMod_Sp)

SamplesOwn_Sp <- runMCMC(comOwnMod_Sp$bOwnMod_Sp, 
                         niter = 15000,
                         thin=10,
                         nburnin = 5000, 
                         nchains = 2)

SummaryOutput(SamplesOwn_Sp, 
              params=c("beta","betaJump",
                       "drift","sigma_eps","sigma_time",
                       "p","a","muY","sdY")) %>% 
  print(n=450)  



#### 2.4 Save Results #########################################################
save(SamplesLiuLi_Sp, 
     SamplesOwn_Sp, file= file.path(getwd(),"Results/SamplesSp.RData"))


######### 2.5. Calculate WAIC and LOO-CV ######################################

#load data in case no new samples have been drawn
load(file.path(getwd(),"Results/SamplesSp.RData"))

#### Caclulation of WAIC and loo 
LikeMatOwn <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesOwn_Sp),
                                  n = length(ZMatSp),
                                  ZMat = ZMatSp)

LikeMatLiu <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesLiuLi_Sp),
                                  n = length(ZMatSp),
                                  ZMat = ZMatSp)

WAICOwn <- loo::waic(log(LikeMatOwn))
WAICLiu <- loo::waic(log(LikeMatLiu))

LOOLiu <- loo::loo(log(LikeMatLiu)) #loo IC = -2 elpd_loo
LOOOwn <- loo::loo(log(LikeMatOwn))


CompDataFrameSp <- data.frame("Model"=c("Own","Liu-Li"),
                              "LOOCV"= c(LOOOwn$estimates[3,1],
                                         LOOLiu$estimates[3,1]),
                              "WAIC"=c(WAICOwn$estimates[3,1],
                                       WAICLiu$estimates[3,1]))

save(CompDataFrameSp,
     file = file.path(getwd(),"Results/WAIC_Sp.RData"))

# ############## 3. ITALY #####################################################

#### 3.1 Put data into right form for NIMBLE ##################################
NAge <- nrow(ZMatIt)
NYear <- ncol(ZMatIt)

NimbleConstIt <- list("N_AgeGroups"=NAge,
                      "N_Year"=NYear)

NimbleDataIt <- list("ZMat"=ZMatIt) 

######### 3.2 Liu,Li Model #####################################################
LiuLi_It <-  nimbleModel(code=LiuLi_Model, 
                         constants = NimbleConstIt, 
                         data=NimbleDataIt)

cLiuLi_It <- configureMCMC(LiuLi_It, 
                           print = TRUE, useConjugacy = TRUE,
                           monitors = c("k",
                                        "sigma_eps",
                                        "beta","betaJump", "N_t","p",
                                        "sigma_time",
                                        "drift",
                                        "muY","sdY",
                                        "mu","Y_t"
                           ),
                           enableWAIC = TRUE) #to compare Models

cLiuLi_It$removeSamplers(paste0("b1[",1:NAge,"]"),
                         paste0("b2[",1:NAge,"]"),
                         "muY","sdY")

cLiuLi_It$addSampler(target= paste0("b1[1:",NAge,"]"),
                     type="AF_slice") 

cLiuLi_It$addSampler(target= paste0("b2[1:",NAge,"]"),
                     type="AF_slice")
cLiuLi_It$addSampler(nodes = "muY",type="slice")
cLiuLi_It$addSampler(nodes = "sdY",type="slice")

bLiuLi_It <- buildMCMC(cLiuLi_It)
comLiuLi_It <- compileNimble(LiuLi_It,
                             bLiuLi_It)


SamplesLiuLi_It <- runMCMC(comLiuLi_It$bLiuLi_It, 
                           niter = 15000,
                           thin=10,
                           nburnin = 5000, 
                           nchains = 2)


SummaryOutput(SamplesLiuLi_It, 
              params=c("beta","betaJump",
                       "drift","sigma_eps","sigma_time",
                       "p","a","muY","sdY")) %>% 
  print(n=250) 



####### 3.3. OWN MODEL #########################################################
OwnMod_It <-  nimbleModel(code=OwnModel, 
                          constants = NimbleConstIt, 
                          data=NimbleDataIt)

cOwnMod_It <- configureMCMC(OwnMod_It, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("k",
                                         "sigma_eps",
                                         "beta","betaJump","sdY","muY",
                                         "N_t","p",
                                         "sigma_time",
                                         "drift","a",
                                         "mu",
                                         "Y_t","J"
                            ),
                            enableWAIC = TRUE) 

############ 3.3.2 Adjust Samplers ############################################
params_to_remove <-  c(
  "p","drift",
  paste0("k[",1:NYear,"]"),
  paste0("Y_t[",3:(NYear),"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "a","muY","sdY"
)

cOwnMod_It$removeSamplers(params_to_remove)

#Gibbs Sampler
cOwnMod_It$addDefaultSampler(nodes=c("drift","p"), 
                             useConjugacy = TRUE)
cOwnMod_It$addDefaultSampler(nodes=paste0("k[",2:NYear,"]"),
                             useConjugacy = TRUE, print = FALSE)
#Slice
cOwnMod_It$addSampler(nodes = "muY",type="slice")
cOwnMod_It$addSampler(nodes = "a",type="slice")
cOwnMod_It$addSampler(nodes = "sdY",type="slice")

# slice sampler for Y_t
for(j in 3:(NYear)){ #one time period longer due to differencing
  cOwnMod_It$addSampler(target=paste0("Y_t[",j,"]"),
                        type="slice")
}

#AF Slice
cOwnMod_It$addSampler(target= paste0("b1[1:",NAge,"]"),
                      type="AF_slice") 
cOwnMod_It$addSampler(target= paste0("b2[1:",NAge,"]"),
                      type="AF_slice")



bOwnMod_It <- buildMCMC(cOwnMod_It)
comOwnMod_It <- compileNimble(OwnMod_It,
                              bOwnMod_It)

SamplesOwn_It <- runMCMC(comOwnMod_It$bOwnMod_It,
                         niter = 15000,
                         thin=10,
                         nburnin = 5000,
                         nchains = 2)

SummaryOutput(SamplesOwn_It, params=c("beta","betaJump",
                                      "drift","sigma_eps","sigma_time",
                                      "p","a","muY","sdY")) %>% print(n=450)  

############ 3.4 Save Results ##################################################

save(SamplesLiuLi_It, 
     SamplesOwn_It, file= file.path(getwd(),"Results/SamplesIt.RData"))


########## 3.5. Calculation of WAIC AND LOO ####################################

#load data in case no new samples have been drawn
load(file.path(getwd(),"Results/SamplesIt.RData"))

LikeMatOwn_It <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesOwn_It),
                                     n = length(ZMatIt),
                                     ZMat = ZMatIt)

LikeMatLiu_It <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesLiuLi_It),
                                     n = length(ZMatIt),
                                     ZMat = ZMatIt)

WAICOwnIt <- loo::waic(log(LikeMatOwn_It))
WAICLiuIt <- loo::waic(log(LikeMatLiu_It))

LOOOwnIt <- loo::loo(log(LikeMatOwn_It))
LOOLiuIt <- loo::loo(log(LikeMatLiu_It)) # loo IC = -2 elpd_loo


CompDataFrameIt <- data.frame("Model"=c("Own","Liu-Li"),
                              "LOOCV"= c(LOOOwnIt$estimates[3,1],
                                         LOOLiuIt$estimates[3,1]),
                              "WAIC"=c(WAICOwnIt$estimates[3,1],
                                       WAICLiuIt$estimates[3,1]))

save(CompDataFrameIt,
     file = file.path(getwd(),"Results/WAIC_It.RData"))


