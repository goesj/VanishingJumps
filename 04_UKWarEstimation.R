## Load necessary files and libraries ###
library(nimble);library(tidyverse); library(rstan);library(loo)

source("01_NimbleModels.R")
source("02_Functions.R") #needs rstan

############## 1. LOADING DATA #################################################
load(file = file.path(getwd(),"Data/UKWARData.RData")) # Load Data


############## 2. Estimating Parameters ########################################

#Estimating Parameters of both the Liu-Li model, as well as our Model 
#Route II Estimation Method in use (Mortality Improvements)

## Create Data and Constants for Nimble Models 
NimbleConst <- list("N_AgeGroups"=nrow(ZMatWar),
                    "N_Year"=ncol(ZMatWar))

NimbleData <- list("ZMat"=ZMatWar) #Mortality Improvement Rates

xlength <- nrow(ZMatWar) #Amount of Age Groups 
tlength <- ncol(ZMatWar) #Amount of Periods 

############ 2.1. Estimation of the Liu-Li Model ############################### 

## Liu, Li Model

#### ATTENTION ! CHANGE PRIORS IN NIMBLE CODE TO THE ONES SHOWN IN THE PAPER ###

#2.1.1 Initialize nimble Model
LiuLiModel_UK <- nimbleModel(code=LiuLi_Model, 
                        constants = NimbleConst, 
                        data=NimbleData, 
                        buildDerivs = FALSE)

#2.1.2 configure nimble Model, to change samplers
cLiuLiModel_UK <- configureMCMC(LiuLiModel_UK,monitors = c("sigma_eps",
                                                  "beta","betaJump",
                                                  "N_t","p",
                                                  "sigma_time","sdY","muY",
                                                  "drift",
                                                  "k",
                                                  "Y_t","J",
                                                  "mu"
                                                  ),
                               enableWAIC = TRUE, useConjugacy = TRUE)

#Define Parameters of which to remove samplers
params_to_remove <-  c(
  paste0("b1[",1:10,"]"),
  paste0("b2[",1:10,"]"),
  "sdY","sigma_time","muY"
)

#Remove Samplers of these parameters
cLiuLiModel_UK$removeSamplers(params_to_remove)

#Create New samplers for the removed parameters
cLiuLiModel_UK$addSampler(target= c("b1[1:10]"),
                  type="AF_slice")
cLiuLiModel_UK$addSampler(target= c("b2[1:10]"),
                  type="AF_slice")

cLiuLiModel_UK$addSampler(nodes = "muY",type="slice")
cLiuLiModel_UK$addSampler(nodes = "sdY",type="slice")
cLiuLiModel_UK$addSampler(nodes = "sigma_time",type="slice")


#2.1.3 Build nimble model
bLiuLiModel_UK <- buildMCMC(cLiuLiModel_UK)
#2.1.4 Compile nimble model 
comLiuLiModel_UK <- compileNimble(LiuLiModel_UK,bLiuLiModel_UK)


#2.1.5 Sample Parameters
SamplesLiuLi_War <- runMCMC(comLiuLiModel_UK$bLiuLiModel_UK, 
                       niter = 15000, 
                       thin=10,
                       nburnin = 5000, 
                       nchains = 2) 

#Create Summary Output of Paramters
SummaryOutput(SamplesLiuLi_War,
              params=c("beta","betaJump","drift","sigma_time",
              "sigma_eps","p","muY","sdY")) %>% print(n=150) 


######## 2.2 Estimation of Own Model  ########################################## 

# Estimation of Own Model with Vanishing Jumps

#2.2.1 Initialize nimble Model
OwnModel_UK <- nimbleModel(code=OwnModel, 
                           constants = NimbleConst, 
                           data=NimbleData, 
                           buildDerivs = FALSE)

#2.2.2 Configuration of nimble model 
cOwnModel_UK <- configureMCMC(OwnModel_UK,
                                monitors = c("sigma_eps","beta","betaJump",
                                             "N_t","p","sigma_time",
                                             "sdY","muY", "drift","a",
                                             "Y_t","J", "b1","b2",
                                             "mu","k"),
                            enableWAIC = TRUE, 
                            #set to false to speed up configuration
                            useConjugacy = FALSE) 

#Define Parameters of which to remove samplers
params_to_remove <-  c(
  "muY","sdY","sigma_time",
  "a","p",
  paste0("b1[",1:xlength,"]"),
  paste0("b2[",1:xlength,"]"),
  paste0("k[",2:(tlength),"]"),
  paste0("Y_t[",3:(tlength+1),"]")
)

cOwnModel_UK$removeSamplers(params_to_remove) #remove parameters
cOwnModel_UK$addSampler(nodes = "muY",type="slice")
cOwnModel_UK$addSampler(nodes = "a",type="slice")
cOwnModel_UK$addSampler(nodes = "sdY",type="slice")
#gibbs sampler for p
cOwnModel_UK$addDefaultSampler(nodes="p", useConjugacy = TRUE)
cOwnModel_UK$addSampler(nodes = "sigma_time",type="slice")
#gibbs sampler of k (attention, k1 is a constant)
cOwnModel_UK$addDefaultSampler(nodes=paste0("k[",2:tlength,"]"),
                                 useConjugacy = TRUE, print = FALSE)

# slice sampler for Y_t
for(j in 3:(tlength+1)){ #one time period longer due to differencing
  cOwnModel_UK$addSampler(target=paste0("Y_t[",j,"]"),
                            type="slice")
}
cOwnModel_UK$addSampler(target= c("b1[1:10]"),
                  type="AF_slice")
cOwnModel_UK$addSampler(target= c("b2[1:10]"),
                  type="AF_slice")


#2.2.3 build Model
bOwnModel_UK <- buildMCMC(cOwnModel_UK)
#2.2.4 compile Model
comOwnModel_UK <- compileNimble(OwnModel_UK,bOwnModel_UK)

#2.2.5 create samples 
Samples_OwnMod_UK <- runMCMC(comOwnModel_UK$bOwnModel_UK, 
                       niter = 15000,
                       thin=10,
                       nburnin = 5000, 
                       nchains = 2) #set seed for reproducibility

#Plot Summary Of samples, same column structure as in paper
SummaryOutput(Samples_OwnMod_UK,
              params=c("beta","betaJump","drift","sigma_time",
                       "sigma_eps","p","a","muY","sdY")) %>% print(n=400) #work


############# 3. Save Results ##################################################
save(Samples_OwnMod_UK,
     SamplesLiuLi_War,
     file = file.path(getwd(),"Results/Samples_UKWar.RData"))

#alternatively, load existing Samples
load(file.path(getwd(),"Results/Samples_UKWar.RData"))

########## 3.1 Calculate WAIC Values ###########################################
LikeMatOwn_UK <- LikelihoodMatrixFun(Samples = do.call(rbind, Samples_OwnMod_UK),
                                    n = length(ZMatWar),
                                    ZMat = ZMatWar)

LikeMatLiuLi_UK <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesLiuLi_War),
                                    n = length(ZMatWar),
                                    ZMat = ZMatWar)

WAICOwnUK <- loo::waic(log(LikeMatOwn_UK))
WAICLiuUK <- loo::waic(log(LikeMatLiuLi_UK))

LOOOwnUK <- loo::loo(log(LikeMatOwn_UK))
LOOLiuUK <- loo::loo(log(LikeMatLiuLi_UK)) # loo IC = -2 elpd_loo


CompDataFrameWar <- data.frame("Model"=c("Own","Liu-Li"),
                              "LOOCV"= c(LOOOwnUK$estimates[3,1],
                                         LOOLiuUK$estimates[3,1]),
                              "WAIC"=c(WAICOwnUK$estimates[3,1],
                                       WAICLiuUK$estimates[3,1]))

########## 3.2 Save Results ####################################################
save(CompDataFrameWar,
     file = file.path(getwd(),"Results/WAIC_UKWar.RData"))


