## RUNNING OF MODEL ###
library(nimble);library(tidyverse); library(nimbleHMC)

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 

load(file = file.path(getwd(),"Data/UKWARData.RData")) # Load Data




##### Liu, Li Model (with Mort Improvements)####################################
NewNimbleConst <- list("N_AgeGroups"=nrow(ZMatWar),
                       "N_Year"=ncol(ZMatWar))

NewNimbleData <- list("ZMat"=ZMatWar) 

xlength <- nrow(ZMatWar)

#Currently Works with the following informative priors:
#muy ~ N(1,2) (or N(0,5), only works sometimes))
#BetaJ ~ Dirichlet(1,1,1,3,3,2,,1,1,1)
DiffLiuModel <- nimbleModel(code=LiuLi_Diff_OC_CornerK, 
                               constants = NewNimbleConst, 
                               data=NewNimbleData, 
                               buildDerivs = FALSE)

cDiffModel <- configureMCMC(DiffLiuModel,monitors = c("sigma_eps",
                                                  "beta","betaJump",
                                                  "N_t","p",
                                                  "sigma_time","sdY","muY",
                                                  "drift",
                                                  "k",
                                                  "Y_t","J",
                                                  "mu"
                                                  #"b1","b2"
                                                  ),
                               enableWAIC = TRUE, useConjugacy = TRUE)


params_to_remove <-  c(
  paste0("b1[",1:10,"]"),
  paste0("b2[",1:10,"]"),
  "sdY","sigma_time","muY"
)

cDiffModel$removeSamplers(params_to_remove)

cDiffModel$addSampler(target= c("b1[1:10]"),
                  type="AF_slice")

cDiffModel$addSampler(target= c("b2[1:10]"),
                  type="AF_slice")

cDiffModel$addSampler(nodes = "muY",type="slice")
cDiffModel$addSampler(nodes = "sdY",type="slice")
cDiffModel$addSampler(nodes = "sigma_time",type="slice")


for(j in 1:108){
  cDiffModel$addSampler(target = paste0("k[",j,"]"),
                        type = 'slice')
}

for(j in 3:(ncol(ZMatWar)+1)){
  cDiffModel$addSampler(target = paste0("Y_t[",j,"]"),
                       type = 'slice')
}

bDiffModel <- buildMCMC(cDiffModel)
comDiff <- compileNimble(DiffLiuModel,bDiffModel)


SamplesLiuLi_OC_War <- runMCMC(comDiff$bDiffModel, 
                       niter = 4000,
                       thin=1,
                       nburnin = 2000, 
                       nchains = 2) #set seed for reproducibility

SummaryOutput(SamplesLiuLi_OC_War,
              params=c("beta","betaJump","a",
                       "drift","p","sigma_eps",
                       "muY","sdY","sigma_time"
                       )) %>% print(n=150) #work


######## OWN MODEL ############################################################# 
DiffOwnModel <- nimbleModel(code=LCJumpOwn_OC_CornerK, 
                            constants = NewNimbleConst, 
                            data=NewNimbleData, 
                            buildDerivs = FALSE)

cDiffModel_Own <- configureMCMC(DiffOwnModel,monitors = c("sigma_eps",
                                                      "beta","betaJump",
                                                      "N_t","p",
                                                      "sigma_time",
                                                      "sdY","muY",
                                                      "drift",
                                                      "a",
                                                      "Y_t","J",
                                                      "b1","b2",
                                                      "mu","k"
                                                      ),
                            enableWAIC = TRUE, useConjugacy = FALSE)

#For Change of Samplers
klength <- max(TotalDataWar$TInd)
xlength <- max(TotalDataWar$NewAgeInd)

params_to_remove <-  c(
  "muY","sdY","sigma_time",
  "a","p",
  paste0("b1[",1:xlength,"]"),
  paste0("b2[",1:xlength,"]"),
  paste0("k[",1:(klength -1),"]"),
  paste0("Y_t[",1:(klength),"]")
)

cDiffModel_Own$removeSamplers(params_to_remove)
cDiffModel_Own$addSampler(nodes = "muY",type="slice")
cDiffModel_Own$addSampler(nodes = "a",type="slice")
cDiffModel_Own$addSampler(nodes = "sdY",type="slice")
cDiffModel_Own$addDefaultSampler(nodes="p", useConjugacy = TRUE)
cDiffModel_Own$addSampler(nodes = "sigma_time",type="slice")
cDiffModel_Own$addDefaultSampler(nodes=paste0("k[",2:(klength -1),"]"),
                                 useConjugacy = TRUE, print = FALSE)

for(j in 3:(klength)){
  cDiffModel_Own$addSampler(target=paste0("Y_t[",j,"]"),
                            type="slice")
}
cDiffModel_Own$addSampler(target= c("b1[1:10]"),
                  type="AF_slice")

cDiffModel_Own$addSampler(target= c("b2[1:10]"),
                  type="AF_slice")

bDiffModel_Own <- buildMCMC(cDiffModel_Own)
comDiff_Own <- compileNimble(DiffOwnModel,bDiffModel_Own)


Samples_OwnMod_OC_War <- runMCMC(comDiff_Own$bDiffModel_Own, 
                       niter = 6500,
                       thin=2,
                       nburnin = 2500, 
                       nchains = 2) #set seed for reproducibility

SummaryOutput(Samples_OwnMod_OC_War,
              params=c("betaJump","sigma_time",
                       "beta","drift","p","sdY","muY","k",
                       "a","N_t","k")) %>% print(n=400) #work

comDiff_Own$bDiffModel_Own$getWAIC()

save(Samples_OwnMod_OC_War,
     SamplesLiuLi_OC_War,
     file = file.path(getwd(),"Results/Samples_UKWar.RData"))

SummaryOutput(Samples_OwnMod_OC,
              params=c("betaJump","sigma_time",
                       "beta","drift","p","sdY","muY",
                       "a")) %>% print(n=400)


LikeMatOwn_War <- LikelihoodMatrixFun(Samples = do.call(rbind, Samples_OwnMod_OC_War),
                                    n = length(ZMatWar),SingleVar = TRUE,
                                    ZMat = ZMatWar)

LikeMatLiu_War <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesLiuLi_OC_War),
                                    n = length(ZMatWar),SingleVar = TRUE,
                                    ZMat = ZMatWar)

WAICOwnWar <- loo::waic(log(LikeMatOwn_War))
WAICLiuWar <- loo::waic(log(LikeMatLiu_War))

LOOLiuWar <- loo::loo(log(LikeMatLiu_War)) # loo IC = -2 elpd_loo
LOOOwnWar <- loo::loo(log(LikeMatOwn_War))


CompDataFrameWar <- data.frame("Model"=c("Own","Liu-Li"),
                              "LOOCV"= c(LOOOwnWar$estimates[3,1],
                                         LOOLiuWar$estimates[3,1]),
                              "WAIC"=c(WAICOwnWar$estimates[3,1],
                                       WAICLiuWar$estimates[3,1]))

save(CompDataFrameWar,
     file = file.path(getwd(),"Results/WAIC_War.RData"))


############ Liu,Li Repara #####################################################
DiffLiuModel_Repa <- nimbleModel(code=LiuLi_Repara_Diff_OC_V2, 
                                 constants = NewNimbleConst, 
                                 data=NewNimbleData, 
                                 buildDerivs = FALSE)

cDiffModel_Repa <- configureMCMC(DiffLiuModel_Repa,monitors = c("sigma_eps",
                                                                "beta","betaJump",
                                                                "N_t","p",
                                                                "sigma_time","sdY","muY",
                                                                "drift",
                                                                #"k",
                                                                "mu","sigma_squared",
                                                                "VarY"),
                                 #"b1","b2"),
                                 enableWAIC = TRUE, useConjugacy = TRUE)
params_to_remove <-  c(
  paste0("b1[",1:10,"]"),
  paste0("b2[",1:10,"]"),
  "sigma_time","muY","VarY"
)
cDiffModel_Repa$removeSamplers(params_to_remove)
cDiffModel_Repa$addSampler(target= c("b1[1:10]"),
                           type="AF_slice")

cDiffModel_Repa$addSampler(target= c("b2[1:10]"),
                           type="AF_slice")

cDiffModel_Repa$addSampler(nodes = "muY",type="slice")
cDiffModel_Repa$addSampler(nodes = "VarY",type="slice")
cDiffModel_Repa$addSampler(nodes = "sigma_time",type="slice")


bDiffModel_Repa <- buildMCMC(cDiffModel_Repa)
comDiff_Repa <- compileNimble(DiffLiuModel_Repa,bDiffModel_Repa)

SamplesLiuLi_Repa <- runMCMC(comDiff_Repa$bDiffModel_Repa, 
                             niter = 4000,
                             thin=1,
                             nburnin = 2000, 
                             nchains = 2) #set seed for reproducibility


SummaryOutput(SamplesLiuLi_Repa,
              params=c("beta","betaJump","a",
                       "drift","p","sigma_eps",
                       "muY","sdY","sigma_time","N_t"
              )) %>% print(n=150) #work

comDiff_Repa$bDiffModel_Repa$getWAIC()

Samples_OwnMod_OC_War <- SamplesLiuLi_OC
Samples_OwnMod_OC_War_Rep <- SamplesLiuLi_Repa

save(Samples_OwnMod_OC_War,
     Samples_OwnMod_OC_War_Rep, file = file.path(getwd(),"Results/UKWarDatLiuMod.RData"))


##### LC MODEL #################################################################
LC_UKWar <-  nimbleModel(code=LC_Diff, 
                      constants = NewNimbleConst, 
                      data=NewNimbleData)

cLC_UKWar <- configureMCMC(LC_UKWar, 
                        print = TRUE, useConjugacy = TRUE,
                        monitors = c("k","sigma_eps",
                                     "beta",
                                     "sigma_time",
                                     "drift","b1","mu"
                        ),
                        enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(
  paste0("b1[",1:xlength,"]")
)

cLC_UKWar$removeSamplers(params_to_remove)
cLC_UKWar$addSampler(target= c("b1[1:10]"),
                  type="AF_slice")
bLC_UKWar <- buildMCMC(cLC_UKWar)
comLC_UKWar <- compileNimble(LC_UKWar,
                          bLC_UKWar)


Samples_LC_UKWar<- runMCMC(comLC_UKWar$bLC_UKWar, 
                        niter = 4000,
                        thin=1,
                        nburnin = 2000, 
                        nchains = 2)

SummaryOutput(Samples_LC_UKWar) %>% print(n=100)

save(Samples_LC_UKWar,
     file = file.path(getwd(),"Results/SamplesUKWar_LC.RData"))

LikeMatLC_War <- LikelihoodMatrixFun(Samples = do.call(rbind, Samples_LC_UKWar),
                                      n = length(ZMatWar),SingleVar = TRUE,
                                      ZMat = ZMatWar)

WAICLCWar <- loo::waic(log(LikeMatLC_War))

LOOLCWar <- loo::loo(log(LikeMatLC_War))


CompDataFrameWar_LC <- data.frame("Model"=c("Own"),
                               "LOOCV"= c(LOOLCWar$estimates[3,1]),
                               "WAIC"=c(WAICLCWar$estimates[3,1]))


save(CompDataFrameWar_LC,
     file = file.path(getwd(),"Results/WAIC_WAR_LC.RData"))

############# OWN MODEL REPARAMETERIZATION #####################################
DiffOwnModel <- nimbleModel(code=OwnMod_Repara_Diff_OC, 
                            constants = NewNimbleConst, 
                            data=NewNimbleData, 
                            buildDerivs = FALSE)

cDiffModel_Own <- configureMCMC(DiffOwnModel,monitors = c("sigma_eps",
                                                          "beta","betaJump",
                                                          "N_t","p",
                                                          "sigma_time",
                                                          "sdY","muY",
                                                          "drift",
                                                          "a","mu","sigma"),
                                enableWAIC = TRUE, useConjugacy = FALSE)


params_to_remove <-  c(
  "muY","sdY","a","p","sigma_eps"
)

cDiffModel_Own$removeSamplers(params_to_remove)
cDiffModel_Own$addSampler(nodes = "muY",type="slice")
cDiffModel_Own$addSampler(nodes = "a",type="slice")
cDiffModel_Own$addSampler(nodes = "sdY",type="slice")
cDiffModel_Own$addDefaultSampler(nodes="p", useConjugacy = TRUE)
cDiffModel_Own$addSampler(nodes = "sigma_time",type="slice")
cDiffModel_Own$addSampler(nodes = "sigma_eps",type="slice")



for(j in 1:(klength-1)){
  cDiffModel_Own$addSampler(target = paste0("k[",j,"]"),
                            type = 'slice')
}

bDiffModel_Own <- buildMCMC(cDiffModel_Own)
comDiff_Own <- compileNimble(DiffOwnModel,bDiffModel_Own)


SamplesDiff <- runMCMC(comDiff_Own$bDiffModel_Own, 
                       niter = 6000,
                       thin=2,
                       nburnin = 3000, 
                       nchains = 1) #set seed for reproducibility

SummaryOutput(SamplesDiff,
              params=c("betaJump","sigma_time",
                       "beta","drift","p","sdY","muY",
                       "N_t","a")) %>% print(n=400) #work

comDiff_Own$bDiffModel_Own$getWAIC()
comDiff_Own_Rep$bDiffModel_Own_Rep$getWAIC()


save(Samples_OwnMod_OC_War,
     Samples_OwnMod_OC_War_Rep,
     file = file.path(getwd(),"Results/UKWarDatOwnMod.RData"))

######## USING QR DECOMPOSITION ################################################
QRMod <- nimbleModel(code=LiuLi_Diff_QR,
                     constants = NewNimbleConst, 
                     data=NewNimbleData)

## test uncompiled function in uncompiled MCMC
QRModc <- configureMCMC(QRMod,monitors = c("k","sigma_eps",
                                           "QMat","BMat", "N_t","p",
                                           "muY","sdY","sigma_time",
                                           "drift","b1","b2",
                                           "Y_t","J","mu"
                                           ),
                        enableWAIC = TRUE)

params_to_remove <-  c(paste0("b1[",1:xlength,"]"),
                       paste0("b2[",1:xlength,"]"),
                       paste0("k[",1:(klength -1),"]"),
                       "sdY","muY")

QRModc$removeSamplers(params_to_remove)

QRModc$addSampler(target= c("b1[1:10]"),
                  type="AF_slice")

QRModc$addSampler(target= c("b2[1:10]"),
                  type="AF_slice")

for(j in 1:(klength -1)){
  QRModc$addSampler(target = paste0("k[",j,"]"),
                    type = 'slice'
  )
}

QRModc$addSampler(target= "sdY",
                  type="RW",log=TRUE)

QRModc$addSampler(target= "muY",
                  type="slice")

## test uncompiled function in uncompiled MCMC
buildQR <- buildMCMC(QRModc)
compQR <- compileNimble(QRMod,buildQR)


SamplesQR <- runMCMC(compQR$buildQR, 
                       niter = 5000,
                       nburnin = 2000,
                       thin=2,
                       nchains = 1) #set seed for reproducibility

SummaryOutput(SamplesQR,params=c("QMat","drift","p","sigma_eps",
                                   "muY","sdY","sigma_time",
                                 "N_t","k")) %>% print(n=150) 


compQR$buildQR$getWAIC()


save(SamplesQR,
     file = file.path(getwd(),"Results/LiuLiQR_WarUK.RData"))

########## OWN MODEL WITH QR ON MORT INPROVEMEN ################################
TestRunOwn_Diff <-  nimbleModel(code=OwnMod_Repara_Diff_QR, 
                                constants = NewNimbleConst, 
                                data=NewNimbleData)


cDiffModel_Own <- configureMCMC(TestRunOwn_Diff, 
                                print = TRUE, useConjugacy = FALSE,
                                monitors = c("sigma_eps",
                                             "QMat", "N_t","p",
                                             "k",
                                             "a",
                                             "muY","sdY","sigma_time",
                                             "drift","b1","b2","BMat"
                                ),
                                enableWAIC = TRUE)

params_to_remove <-  c(
  paste0("k[",1:(klength-1),"]"),
  paste0("b1[",1:xlength,"]"),
  paste0("b2[",1:xlength,"]"),
  "drift","sigma_time","muY","sdY")

cDiffModel_Own$removeSamplers(params_to_remove)

cDiffModel_Own$addSampler(target= c("b1[1:10]"),
                          type="AF_slice") 
cDiffModel_Own$addSampler(target= c("b2[1:10]"),
                          type="AF_slice")

cDiffModel_Own$addSampler(target="sdY", type="RW", log=TRUE)
cDiffModel_Own$addSampler(target="sigma_time", type="RW", log=TRUE)
cDiffModel_Own$addSampler(target="drift", type="slice")
cDiffModel_Own$addSampler(target="muY", type="slice")



for(j in 1:(klength-1)){
  cDiffModel_Own$addSampler(target = paste0("k[",j,"]"),
                            type = 'slice')
}


bDiffModel_Own <- buildMCMC(cDiffModel_Own)
comDiff_Own <- compileNimble(TestRunOwn_Diff,bDiffModel_Own)

SamplesDiff_Own <- runMCMC(comDiff_Own$bDiffModel_Own, 
                           niter = 5000,
                           thin=1,
                           nburnin = 2000, 
                           nchains = 1)


SummaryOutput(SamplesDiff_Own, params=c("QMat","drift","p","sigma_eps",
                                        "muY","sdY","sigma_time","a",
                                        "N_t","a","k")) %>% print(n=170) 

MCMCvis::MCMCsummary(SamplesDiff_Own,params=c("drift","sigma_time","p","a"))

##### AL OLD UNDERNEATH ########################################################
# ################## OLD PARAMETERIZATION ########################################
# TestRunOwn_Diff <-  nimbleModel(code=LCJumpOwn_QR, 
#                            constants = NewNimbleConst, 
#                            data=NewNimbleData)
# 
# 
# cDiffModel_Own <- configureMCMC(TestRunOwn_Diff, 
#                             print = TRUE, useConjugacy = FALSE,
#                             monitors = c("sigma_eps",
#                                          "QMat", "N_t","p",
#                                          "k",
#                                          "J","Y_t",
#                                          "a",
#                                          "muY","sdY","sigma_time",
#                                          "drift","b1","b2","BMat"
#                                          ),
#                             enableWAIC = TRUE)
# 
# params_to_remove <-  c(paste0("b1[",1:xlength,"]"),
#                        paste0("b2[",1:xlength,"]"),
#                        paste0("k[",1:(klength-1),"]"),
#                        paste0("Y_t[",1:(klength),"]"),
#                        "drift","sigma_time","muY","sdY")
# 
# cDiffModel_Own$removeSamplers(params_to_remove)
# 
# cDiffModel_Own$addSampler(target= c("b1[1:10]"),
#                        type="AF_slice") 
# cDiffModel_Own$addSampler(target= c("b2[1:10]"),
#                        type="AF_slice")
# 
# cDiffModel_Own$addSampler(target="sdY", type="RW", log=TRUE)
# cDiffModel_Own$addSampler(target="sigma_time", type="RW", log=TRUE)
# cDiffModel_Own$addSampler(target="drift", type="slice")
# cDiffModel_Own$addSampler(target="muY", type="slice")
# 
# 
# # cDiffSim_QR$addDefaultSampler(nodes = c("p","muY","drift"),
# #                               useConjugacy = TRUE)
# for(j in 1:(klength-1)){
#   cDiffModel_Own$addSampler(target = paste0("k[",j,"]"),
#                          type = 'slice')
# }
# for(j in 2:klength){
#   cDiffModel_Own$addSampler(target=paste0("Y_t[",j,"]"),
#                          type="slice")
# }
# 
# 
# bDiffModel_Own <- buildMCMC(cDiffModel_Own)
# comDiff_Own <- compileNimble(TestRunOwn_Diff,bDiffModel_Own)
# 
# SamplesDiff_Own <- runMCMC(comDiff_Own$bDiffModel_Own, 
#                        niter = 5000,
#                        thin=1,
#                        nburnin = 2000, 
#                        nchains = 1)
# 
# 
# SummaryOutput(SamplesDiff_Own, params=c("QMat","drift","p","sigma_eps",
#                                     "muY","sdY","sigma_time","a")) %>% print(n=170) 
# 
# MCMCvis::MCMCsummary(SamplesDiff_Own,params=c("drift","sigma_time","p","a"))
# 
# comDiff_Own$bDiffModel_Own$getWAIC()
# 
# 
# ts.plot(SamplesDiff$chain1[,"sdY"])
# lines(SamplesDiff$chain2[,"sdY"],col="red")
# 
# save(SamplesDiff_Own, 
#      file=file.path(getwd(),"SampsOwnModUK_War.Rdata"))
# ## Posterior predictive Checks
# Samples <- SamplesDiff_Own
# muPos <- grep("mu",colnames(Samples), 
#               fixed = TRUE)
# sdPos <- grep("sigma", colnames(Samples), 
#               fixed = TRUE)
# 
# muPos[1:(length(muPos)-1)] #without muY
# sdPos[1:(length(sdPos)-2)] #without simgaeps, sigmat, 
# 
# 
# Z_Rep <- matrix(data = 0,nrow = nrow(Samples),
#                 ncol = length(ZMat))
# 
# for(s in 1:nrow(Samples)){
#   Z_Rep[s,] <- rnorm(n = length(ZMat), 
#                      mean = Samples[s,muPos[1:(length(muPos)-1)]],
#                      sd = sqrt(Samples[s,sdPos[1:(length(sdPos)-2)]]))
# }
# 
# 
# LambdaArray_Rep <- array(data = 0, dim = c(nrow(LambdaMat),
#                                            ncol(LambdaMat),
#                                            nrow(Samples)))
# 
# LambdaArray_Rep[,1,] <- log(LambdaMat[,1]) #First year is the same
# for(s in 1:nrow(Samples)){
#   Z_Mat_Hat <- matrix(Z_Rep[s,],nrow = xlength, byrow = FALSE) #transformation of Z_Rep vector into matrix
#   LambdaArray_Rep[,2,s] <- log(LambdaMat[,1]) + Z_Mat_Hat[,1] #Second Year 
#   for(t in 3:ncol(LambdaMat)){ #starting for year 3, in recursive
#     #LambdaArray_US_Rep[,t,s] <- LambdaArray_US_Rep[,(t-1),s] + Z_Mat_Hat[,(t-1)] OOS FC
#     LambdaArray_Rep[,t,s] <- log(LambdaMat[,(t-1)]) + Z_Mat_Hat[,(t-1)] #In-Sample FC
#   }
# }
# 
# LowerPI <- LambdaArray_Rep %>% apply(., c(1,2),quantile, 0.1) %>% 
#   data.frame() %>% 
#   mutate("AgeInd"=1:xlength, .before = 1) %>% 
#   pivot_longer(., cols = 2:ncol(.),names_to = "Year", values_to = "lPI")
# 
# UpperPI <- LambdaArray_Rep %>% apply(., c(1,2),quantile, 0.9) %>% 
#   data.frame() %>% 
#   mutate("AgeInd"=1:xlength, .before = 1) %>% 
#   pivot_longer(., cols = 2:ncol(.),names_to = "Year", values_to = "uPI")
# 
# LowerPI %>% 
#   mutate("Year"=rep(1901:(1901+klength-1),
#                     xlength)) %>%
#   mutate("uPI"=UpperPI$uPI) %>% 
#   filter(AgeInd %in% c(3:5)) %>% 
#   ggplot(data=., aes(x=Year, group=AgeInd))+
#   geom_ribbon(aes(ymin=lPI, ymax=uPI, fill=AgeInd))+
#   ylab("Log Death Rate")+
#   geom_line(data=filter(LambdaVec,NewAgeInd %in% c(3:5)), 
#             aes(x=Year, y=log(Rate), group=NewAgeInd),
#             col="red")
# 
# ## Plotting ZMatrix
# LowerPI_Z <- 
#   Z_Rep %>% apply(., 2 , quantile, 0.1) %>% 
#   matrix(nrow = 7) %>% 
#   data.frame() %>% 
#   mutate("AgeInd"=1:7, .before = 1) %>% 
#   pivot_longer(., cols = 2:ncol(.),names_to = "Year", values_to = "lPI")
# 
# UpperPI_Z <- 
#   Z_Rep %>% apply(., 2 , quantile, 0.9) %>% 
#   matrix(nrow = 7) %>% 
#   data.frame() %>% 
#   mutate("AgeInd"=1:7, .before = 1) %>% 
#   pivot_longer(., cols = 2:ncol(.),names_to = "Year", values_to = "uPI")
# 
# LowerPI_Z %>% 
#   mutate("Year"=rep(2:69,7)) %>%
#   mutate("uPI"=UpperPI_Z$uPI) %>% 
#   filter(AgeInd %in% c(2:4)) %>% 
#   ggplot(data=., aes(x=Year, group=AgeInd))+
#   geom_ribbon(aes(ymin=lPI, ymax=uPI, fill=AgeInd))+
#   ylab("Log Death Rate")+
#   geom_line(data=filter(LambdaVec,NewAgeInd %in% c(2:4)), 
#             aes(x=TInd, y=ZVal, group=NewAgeInd),
#             col="red")
# 
# 
#   
#   
#   
# ### MA5 with Differenced log Rates #######################################
# Q <- 5
# NewNimbleConstMAQ <-  list("N_AgeGroups"=nrow(ZMat),
#                            "N_Year"=ncol(ZMat),
#                            "Q"=Q)
# 
# TestRunOwnMA_Diff <-  nimbleModel(code=LCJumpOwn_MAQ_QR, 
#                                 constants = NewNimbleConstMAQ, 
#                                 data= NewNimbleData)
# 
# cDiffModel_MA <- configureMCMC(TestRunOwnMA_Diff, 
#                             print = TRUE, useConjugacy = FALSE,
#                             monitors = c("k","sigma_eps",
#                                          "QMat", "N_t","p",
#                                          "J","muY","sdY","sigma_time",
#                                          "drift","b1","b2","BMat","Y_t","a"),
#                             enableWAIC = TRUE)
# 
# params_to_remove <-  c(paste0("b1[",1:10,"]"),
#                        paste0("b2[",1:10,"]"),
#                        paste0("k[",1:108,"]"),
#                        paste0("Y_t[",1:108,"]"),
#                        "muY","drift","sdY")
# 
# cDiffModel_MA$removeSamplers(params_to_remove)
# 
# cDiffModel_MA$addSampler(target= c("b1[1:10]"),
#                       type="AF_slice") 
# 
# cDiffModel_MA$addSampler(target= c("b2[1:10]"),
#                       type="AF_slice")
# 
# cDiffModel_MA$addSampler(target = "muY", type="slice")
# cDiffModel_MA$addSampler(target = "sdY", type="slice")
# cDiffModel_MA$addSampler(target = "drift", type="RW")
# 
# # cDiffSim_QR$addDefaultSampler(nodes = c("p","muY","drift"),
# #                               useConjugacy = TRUE)
# 
# for(j in 1:108){
#   cDiffModel_MA$addSampler(target = paste0("k[",j,"]"),
#                         type = 'slice')
#   
#   cDiffModel_MA$addSampler(target=paste0("Y_t[",j,"]"),
#                         type="slice")
# }
# 
# 
# bDiffModel_MA <- buildMCMC(cDiffModel_MA)
# comDiff_MA <- compileNimble(TestRunOwnMA_Diff,bDiffModel_MA)
# 
# 
# SamplesDiff_MA <- runMCMC(comDiff_MA$bDiffModel_MA, 
#                        niter = 5000,
#                        thin=2,
#                        nburnin = 2000, 
#                        nchains = 2)
# 
# SummaryOutput(SamplesDiff_MA, params=c("QMat","drift","p","sigma_eps",
#                                     "muY","sdY","sigma_time","a")) %>% print(n=150) 
# 
# comDiff_MA$bDiffModel_MA$getWAIC()
# 
# 
# ### NORMAL LC MODEL ON DIFFERENCED DATA #######################################
# TestRunLC_Diff <-  nimbleModel(code=LC_AltDirichBeta_Diff, 
#                                   constants = NewNimbleConst, 
#                                   data= NewNimbleData)
# 
# cDiffModel_LC <- configureMCMC(TestRunLC_Diff, 
#                                print = TRUE,
#                                monitors = c("k","sigma_eps","sigma_time",
#                                             "beta","b","drift"),
#                                enableWAIC = TRUE)
# cDiffModel_LC$removeSamplers(paste0("b[",1:10,"]"))
# 
# cDiffModel_LC$addSampler(target= c("b[1:10]"),
#                          type="AF_slice") 
# 
# bDiffModel_LC <- buildMCMC(cDiffModel_LC)
# comDiff_LC <- compileNimble(TestRunLC_Diff,
#                             bDiffModel_LC)
# 
# 
# SamplesDiff_LC <- runMCMC(comDiff_LC$bDiffModel_LC, 
#                           niter = 5000,
#                           thin=2,
#                           nburnin = 2000, 
#                           nchains = 2)
# 
# SummaryOutput(SamplesDiff_LC, 
#               params=c("beta","drift","sigma_eps","sigma_time")) %>% 
#   print(n=150) 
# 
# comOwnQR$bOwnQR$getWAIC()
# 
# 
# TestRunLC <- nimbleMCMC(code=LC_SumToZero,
#                          constants = JumpConst, data=JumpData, 
#                          monitors = c("kappa","alpha","beta","sigma_eps","drift"), 
#                          nburnin = 3000,
#                          niter = 8000, nchains = 2, thin = 5, WAIC = TRUE)
# 
# comOwnQR$bOwnQR$getWAIC()
# TestRunLC$WAIC
# 
# 
# ##### ESTIMATE MODEL ON DEATH RATES ###########################################
# TotalData <- TotalData %>% mutate("MortRate"=Y/Offset)
# 
# 
# RateConst <- list("N_Year"=max(TotalData$TInd),
#                  "N_AgeGroups"=max(TotalData$NewAgeInd),
#                  "N"=length(TotalData$Year),
#                  "age"=TotalData$NewAgeInd,
#                  "year"=TotalData$TInd)
# 
# #Data
# RateData <- list("y"=log(TotalData$MortRate))
# 
# RateInits <- list("drift"=0.1, #negative value
#                  "muY" = 2)  #positive value
# 
# 
# TestRunOwn_NP <- nimbleModel(code=LCJump_QR_NonPoisson,
#                      constants = RateConst, 
#                      data=RateData,
#                      inits = RateInits)
# 
# ## test uncompiled function in uncompiled MCMC
# NPc <- configureMCMC(TestRunOwn_NP,
#                      monitors = c("k","sigma_eps",
#                                   "QMat", "N_t","p",
#                                   "J","muY","sdY","sigma_time",
#                                   "drift","b1","b2","BMat","Y_t","kappa",
#                                   "alpha"),
#                         enableWAIC = TRUE, print=TRUE, useConjugacy = FALSE)
# 
# NYear <- max(TotalData$TInd)
# 
# params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
#                        paste0("b2[",1:NAge,"]"),
#                        paste0("alpha[",1:NAge,"]"),
#                        "sdY","muY")
# 
# NPc$removeSamplers(params_to_remove)
# 
# NPc$addSampler(target= c("b1[1:10]"),
#                   type="AF_slice")
# 
# NPc$addSampler(target= c("b2[1:10]"),
#                   type="AF_slice")
# 
# NPc$addSampler(target= "sdY",
#                   type="slice")
# 
# NPc$addSampler(target= "muY",
#                   type="slice")
# 
# NPc$addDefaultSampler(nodes= c("alpha[1:10]"),
#                             useConjugacy = TRUE)
# 
# # for(j in 1:NYear){
# #   NPc$addSampler(target = paste0("Y_t[",j,"]"),
# #                        type = 'slice'
# #   )
# # }
# 
# 
# 
# ## test uncompiled function in uncompiled MCMC
# buildNP <- buildMCMC(NPc)
# compNP <- compileNimble(TestRunOwn_NP,buildNP)
# 
# SamplesDiff <- runMCMC(compNP$buildNP, 
#                        niter = 5000,
#                        nburnin = 2000,
#                        thin=2,
#                        nchains = 1) #set seed for reproducibility
# 
# SummaryOutput(SamplesDiff,
#               params=c("QMat","drift","p","sigma_eps",
#                        "muY","sdY","sigma_time","alpha")) %>% print(n=150) 
# 
# 
# x11()
# ts.plot(SamplesDiff[,"drift"])
