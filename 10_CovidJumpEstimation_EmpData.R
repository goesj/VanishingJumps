## Estimation of Own and Liu,Li Model with Denmark Data ###
library(nimble);library(tidyverse); library(nimbleHMC);library(pacman)

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 

load(file = file.path(getwd(),"Data/CovidData.RData")) # Load Data
############# 1. UNITED STATES### #################################################

#Some plots...
LambdaVecUs %>% 
    ungroup() %>% 
    ggplot(data=., aes(x=Year, y=log(Rate), group=NewAgeInd))+
    geom_line(aes(col=NewAgeInd))+
    ylab("Log Death Rate")+
    facet_wrap(~NewAgeInd, scales = "free")
  
  
LambdaVecUs %>% 
  ggplot(data=., aes(x=Year, y=ZVal, group=factor(NewAgeInd)))+
  geom_line(aes(col=factor(NewAgeInd)))

LambdaVecUs %>% 
  ggplot(data=., aes(x=Year, y=log(Rate), group=NewAgeInd))+
  geom_line(aes(col=NewAgeInd))+
  ylab("Log Death Rate")+
  facet_wrap(~NewAgeInd,scales = "free")

###### 1.1 Put Data into correct form for NIMBLE ###############################
#Create Death Matrix
ZMatUS_LessAges <- ZMatUS[-c(1:8,19),]
NAge <- nrow(ZMatUS)
NYear <- ncol(ZMatUS)

NimbleConstUS <- list("N_AgeGroups"=NAge,
                       "N_Year"=NYear)

NimbleDataUS <- list("ZMat"=ZMatUS) 

###### 1.2 QR APPROACH #########################################################
##### 1.2.1 Estimate Model #####################################################
OwnQR_US <-  nimbleModel(code=OwnMod_Repara_Diff_QR, 
                          constants = NimbleConstUS, 
                          data=NimbleDataUS, buildDerivs = FALSE)

cOwnQR_US <- configureMCMC(OwnQR_US, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("k","sigma_eps",
                                         "QMat", "N_t","p",
                                         "sigma_time",
                                         "drift","b1","b2","BMat","a",
                                         #"Y_t"
                                         "muY","sdY"
                                         ),
                            enableWAIC = TRUE) #to compare Models

############# 1.2.2 Adjust Samplers ############################################
params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
                       paste0("b2[",1:NAge,"]"),
                       paste0("k[",1:NYear,"]"),
                       #paste0("Y_t[",1:(NYear+1),"]")
                       "sdY","muY","p","a"
                       )

cOwnQR_US$removeSamplers(params_to_remove)

cOwnQR_US$addSampler(target= paste0("b1[1:",NAge,"]"),
                      type="AF_slice") 

cOwnQR_US$addSampler(target= paste0("b2[1:",NAge,"]"),
                      type="AF_slice")


 cOwnQR_US$addSampler(target = "sdY", type="RW", log=TRUE)
 cOwnQR_US$addSampler(target = "muY", type="slice")
 cOwnQR_US$addSampler(target = "a", type="slice")
 cOwnQR_US$addSampler(target = "p", type="slice")
 

for(j in 1:NYear){
  cOwnQR_US$addSampler(target = paste0("k[",j,"]"),
                        type = 'slice')
}

for(j in 3:(NYear+1)){
  cOwnQR_US$addSampler(target = paste0("Y_t[",j,"]"),
                        type = 'slice')
}

bOwnQR_US <- buildMCMC(cOwnQR_US)
comOwnQR_US <- compileNimble(OwnQR_US,
                              bOwnQR_US)


SamplesOwnQR_US <- runMCMC(comOwnQR_US$bOwnQR_US, 
                            niter = 4000,
                            thin=1,
                            nburnin = 2000, 
                            nchains = 1)

SummaryOutput(SamplesOwnQR_US, 
              params=c("QMat","drift","p","sigma_eps",
                       "muY","sdY","sigma_time","a","N_t","Y_t"
                       )) %>% 
  print(n=250) 


Y_t_SampUS <- SamplesOwnQR_US[,grepl("Y_t",colnames(SamplesOwnQR_US))]
N_t_SampUS <- SamplesOwnQR_US[,grepl("N_t",colnames(SamplesOwnQR_US))]



#Only Y_t's with subsequent N_t's of non zero
Y_tN_t <- Y_t_SampUS*N_t_SampUS
Y_t_noZero <- as.vector(Y_tN_t)[which(as.vector(Y_tN_t)!=0)] 

## Estimates of Jump Compontent ?? 
mean(Y_t_noZero)
sd(Y_t_noZero)


### Stan Estimation of normal mixture (takes very! long)
Y_t_LessUS <- apply(Y_t_SampUS[,-c(1,2)],  #remove zeros
                    2, sample, size=500, 
                    replace=FALSE)
plot(density(Y_t_LessUS))

init_fun <- function(...) list(theta=c(0.05,0.95))

library(rstan)
set.seed(42)
dataStan <- list("K"=2,
                 "N"=length(as.vector(Y_t_LessUS)),
                 "y"=as.vector(Y_t_LessUS))

#Attention takes long !!!!
options(mc.cores = parallel::detectCores())
MixMod <- stan(file = file.path(getwd(),"Stan Code/NormalMixture.stan"),
               data=dataStan, chains = 2, init=init_fun, 
               iter = 4000, warmup = 2000,
               thin=2, control = list(adapt_delta = 0.81,
                                      max_treedepth=11))

summary(MixMod)$summary %>% round(.,4)

ParaY_t <- mixtools::normalmixEM(as.vector(Y_t_LessUS),
                                 lambda = c(0.05,0.95),
                                 k=2,
                                 arbmean = TRUE, 
                                 fast = TRUE,
                                 maxit = 2000)

summary(ParaY_t)

######## 1.3 SUM-TO-ONE (OC)##########################################
##### 1.3.1 Estimate Model ####################################################
#Using a sum to zero constraint
OwnMod_US <-  nimbleModel(code=LCJumpOwn_OC_CornerK, 
                         constants = NimbleConstUS, 
                         data=NimbleDataUS)

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
  paste0("Y_t[",1:(NYear+1),"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "a","muY","sdY"
)

cOwnMod_US$removeSamplers(params_to_remove)
cOwnMod_US$addDefaultSampler(nodes=c("drift","p"), 
                             useConjugacy = TRUE)

cOwnMod_US$addSampler(nodes = "muY",type="slice")
cOwnMod_US$addSampler(nodes = "a",type="slice")
cOwnMod_US$addSampler(nodes = "sdY",type="slice")


cOwnMod_US$addDefaultSampler(nodes=c("Y_t[3:41]"), 
                                 useConjugacy = TRUE)

cOwnMod_US$addDefaultSampler(nodes=c("k[2:40]"), 
                             useConjugacy = TRUE)

# for(j in 2:NYear){
#   cOwnMod_US$addSampler(target = paste0("k[",j,"]"),
#                         type = 'slice')
# }
# 
# for(j in 3:(NYear+1)){
#   cOwnMod_US$addSampler(target = paste0("Y_t[",j,"]"),
#                         type = 'slice')
# }

cOwnMod_US$addSampler(target= paste0("b1[1:",NAge,"]"),
                          type="AF_slice") 
cOwnMod_US$addSampler(target= paste0("b2[1:",NAge,"]"),
                          type="AF_slice")

bOwnMod_US <- buildMCMC(cOwnMod_US)
comOwnMod_US <- compileNimble(OwnMod_US,
                            bOwnMod_US)


SamplesOwn_US_NoRepara <- runMCMC(comOwnMod_US$bOwnMod_US, 
                          niter = 17500,
                          thin=10,
                          nburnin = 7500, 
                          nchains = 2)

SummaryOutput(SamplesOwn_US_NoRepara, 
              params=c("drift","p","sigma_eps",
                       "muY","sigma_time","beta","betaJump",
                       "a","sdY","N_t")) %>% 
  print(n=250) 

ts.plot(SamplesOwn_US_NoRepara$chain1[,"sdY"])
lines(SamplesOwn_US_NoRepara$chain2[,"sdY"],col="red")




####### 1.4. REPARA SUM-TO-ONE #######################################
#Compare that with The LIU,Li model using the same reparameterization
OwnMod_US_Rep <-  nimbleModel(code=OwnMod_Repara_Diff_OC_V2, 
                          constants = NimbleConstUS, 
                          data=NimbleDataUS)

cOwnMod_US_Rep <- configureMCMC(OwnMod_US_Rep, 
                            print = TRUE, useConjugacy = TRUE,
                            monitors = c(#"k",
                                         "sigma_eps",
                                         "beta","betaJump", "N_t","p",
                                         "sigma_time",
                                         "drift",#"b1","b2",
                                         #"J","Y_t",
                                         "muY","sdY","a",
                                         "mu","sigma_squared"
                            ),
                            enableWAIC = TRUE) #to compare Models

############# 1.4.1 Adjust Samplers ###########################################
params_to_remove <-  c(
  #paste0("k[",1:NYear,"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "a",
  "muY","sdY","sigma_time"
)

cOwnMod_US_Rep$removeSamplers(params_to_remove)
cOwnMod_US_Rep$addSampler(nodes = "muY",type="slice")
cOwnMod_US_Rep$addSampler(nodes = "a",type="slice")
cOwnMod_US_Rep$addSampler(nodes = "sdY",type="slice")
cOwnMod_US_Rep$addSampler(nodes = "sigma_time",type="slice")

for(j in 1:NYear){
  cOwnMod_US_Rep$addSampler(target = paste0("k[",j,"]"),
                        type = 'slice')
}

cOwnMod_US_Rep$addSampler(target= paste0("b1[1:",NAge,"]"),
                          type="AF_slice") 
cOwnMod_US_Rep$addSampler(target= paste0("b2[1:",NAge,"]"),
                          type="AF_slice")

bOwnMod_US_Rep <- buildMCMC(cOwnMod_US_Rep)
comOwnMod_US_Rep <- compileNimble(OwnMod_US_Rep,
                              bOwnMod_US_Rep)

SamplesOwn_US_Rep2 <- runMCMC(comOwnMod_US_Rep$bOwnMod_US_Rep, 
                            niter = 4000,
                            thin=1,
                            nburnin = 2000, 
                            nchains = 2)

SummaryOutput(SamplesOwn_US_Rep2, 
              params=c("drift","p","sigma_eps",
                       "muY","sigma_time","beta","betaJump",
                       "a","sdY")) %>% 
  print(n=300) 

comOwnMod_US_Rep$bOwnMod_US_Rep$getWAIC()

############# 1.5 Liu-Li Model ##############################################
LiLi_US_OC <-  nimbleModel(code=LiuLi_Diff_OC_CornerK, 
                           constants = NimbleConstUS, 
                           data=NimbleDataUS)

cLiLi_US_OC <- configureMCMC(LiLi_US_OC, 
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

cLiLi_US_OC$removeSamplers(params_to_remove)
cLiLi_US_OC$addSampler(nodes = "muY",type="slice")
cLiLi_US_OC$addSampler(nodes = "sdY",type="slice")
cLiLi_US_OC$addSampler(target= paste0("b1[1:",NAge,"]"),
                          type="AF_slice") 
cLiLi_US_OC$addSampler(target= paste0("b2[1:",NAge,"]"),
                          type="AF_slice")


for(j in 1:NYear){
  cLiLi_US_OC$addSampler(target = paste0("k[",j,"]"),
                      type = 'slice')
}



bLiLi_US_OC <- buildMCMC(cLiLi_US_OC)
comLiLi_US_OC <- compileNimble(LiLi_US_OC,
                               bLiLi_US_OC)


SamplesLiLi_US <- runMCMC(comLiLi_US_OC$bLiLi_US_OC, 
                          niter = 17500,
                          thin=10,
                          nburnin = 7500, 
                          nchains = 2)


SummaryOutput(SamplesLiLi_US, 
              params=c("QMat","drift","p","sigma_eps",
                       "muY","sigma_time","beta","betaJump","sdY")) %>% 
  print(n=250) 


save(SamplesLiLi_US,
     SamplesOwn_US_NoRepara,
     file = file.path(getwd(),"Results/SamplesUS_NoRep.RData"))



## Compare WAIC's
comOwnMod_US$bOwnMod_US$getWAIC()
comLiLi_US_OC$bLiLi_US_OC$getWAIC()


LikeMatOwnUS <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesOwn_US_NoRepara),
                                    n = length(ZMatUS),SingleVar = TRUE,
                                    ZMat = ZMatUS)

LikeMatLiuUS <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesLiLi_US),
                                    n = length(ZMatUS),SingleVar = TRUE,
                                    ZMat = ZMatUS)

WAICOwnUS <- loo::waic(log(LikeMatOwnUS))
WAICLiuUS <- loo::waic(log(LikeMatLiuUS))

LOOOwnUS <- loo::loo(log(LikeMatOwnUS))
LOOLiuUS <- loo::loo(log(LikeMatLiuUS)) # loo IC = -2 elpd_loo

CompDataFrameUS <- data.frame("Model"=c("Own","Liu-Li"),
                              "LOOCV"= c(LOOOwnUS$estimates[3,1],
                                         LOOLiuUS$estimates[3,1]),
                              "WAIC"=c(WAICOwnUS$estimates[3,1],
                                       WAICLiuUS$estimates[3,1]))

save(CompDataFrameUS,
     file = file.path(getwd(),"Results/WAIC_US_NoRep.RData"))



##### 1.6. Liu-Li  Repara ######################################################
LiLi_US_OC <-  nimbleModel(code=LiuLi_Repara_Diff_OC_V2, 
                           constants = NimbleConstUS, 
                           data=NimbleDataUS)

cLiLi_US_OC <- configureMCMC(LiLi_US_OC, 
                             print = TRUE, useConjugacy = TRUE,
                             monitors = c(#"k",
                                          "sigma_eps",
                                          "beta","betaJump", "N_t","p",
                                          "sigma_time",
                                          "drift",#"b1","b2",
                                          "muY","sdY","mu",#"Y_t"
                                          "sigma_squared"
                             ),
                             enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
                       paste0("b2[",1:NAge,"]"),
                       "muY","sdY"
)

cLiLi_US_OC$removeSamplers(params_to_remove)
cLiLi_US_OC$addSampler(nodes = "muY",type="slice")
cLiLi_US_OC$addSampler(nodes = "sdY",type="slice")
cLiLi_US_OC$addSampler(target= paste0("b1[1:",NAge,"]"),
                       type="AF_slice") 
cLiLi_US_OC$addSampler(target= paste0("b2[1:",NAge,"]"),
                       type="AF_slice")


for(j in 1:NYear){
  cLiLi_US_OC$addSampler(target = paste0("k[",j,"]"),
                         type = 'slice')
}


bLiLi_US_OC <- buildMCMC(cLiLi_US_OC)
comLiLi_US_OC <- compileNimble(LiLi_US_OC,
                               bLiLi_US_OC)

SamplesLiLi_US <- runMCMC(comLiLi_US_OC$bLiLi_US_OC, 
                          niter = 4000,
                          thin=1,
                          nburnin = 2000, 
                          nchains = 2)

SummaryOutput(SamplesLiLi_US, 
              params=c("drift","p","sigma_eps",
                       "muY","sigma_time","beta","betaJump","sdY")) %>% 
  print(n=250) 

save(SamplesLiLi_US,
     SamplesOwn_US_Rep,
     file = file.path(getwd(),"Results/SamplesUS_Rep2.RData"))



## Compare WAIC's
comOwnMod_US_Rep$bOwnMod_US_Rep$getWAIC()
comLiLi_US_OC$bLiLi_US_OC$getWAIC()


LikeMatOwnUS <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesOwn_US_Rep),
                                    n = length(ZMatUS),SingleVar = FALSE,
                                    ZMat = ZMatUS)

LikeMatLiuUS <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesLiLi_US),
                                    n = length(ZMatUS),SingleVar = FALSE,
                                    ZMat = ZMatUS)
WAICOwnUS <- loo::waic(log(LikeMatOwnUS))
WAICLiuUS <- loo::waic(log(LikeMatLiuUS))

LOOOwnUS <- loo::loo(log(LikeMatOwnUS))
LOOLiuUS <- loo::loo(log(LikeMatLiuUS)) # loo IC = -2 elpd_loo

CompDataFrameUS <- data.frame("Model"=c("Own","Liu-Li"),
                              "LOOCV"= c(LOOOwnUS$estimates[3,1],
                                         LOOLiuUS$estimates[3,1]),
                              "WAIC"=c(WAICOwnUS$estimates[3,1],
                                       WAICLiuUS$estimates[3,1]))

save(CompDataFrameUS,
     file = file.path(getwd(),"Results/WAIC_US_Rep2.RData"))


#### 1.7. LC Model Diff ########################################################
LC_US <-  nimbleModel(code=LC_Diff, 
                           constants = NimbleConstUS, 
                           data=NimbleDataUS)

cLC_US <- configureMCMC(LC_US, 
                             print = TRUE, useConjugacy = TRUE,
                             monitors = c("k","sigma_eps",
                                          "beta",
                                          "sigma_time",
                                          "drift","b1","mu"
                             ),
                             enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(
  paste0("b1[",1:NAge,"]")
)
cLC_US$removeSamplers(params_to_remove)
cLC_US$addSampler(target= c("b1[1:10]"),
                           type="AF_slice")
bLC_US <- buildMCMC(cLC_US)
comLC_US <- compileNimble(LC_US,
                          bLC_US)


Samples_LC_US<- runMCMC(comLC_US$bLC_US, 
                              niter = 4000,
                              thin=1,
                              nburnin = 2000, 
                              nchains = 2)

SummaryOutput(Samples_LC_US) %>% print(n=100)

save(Samples_LC_US,
     file = file.path(getwd(),"Results/SamplesUS_LC.RData"))



################ 2. ITALY #####################################################
### Some Plots
LambdaVecIt  %>% 
  mutate("DeathRate"=log(Rate)) %>% 
  filter(Year > 1980) %>%
  ggplot(aes(x=Year,y=DeathRate, group=factor(NewAgeInd)))+
  geom_line(aes(color=factor(NewAgeInd)))+
  facet_wrap(~NewAgeInd, scales = "free")

x11()
LambdaVecIt  %>% 
  mutate("YInd"=match(Year,unique(Year))) %>% 
  filter(NewAgeInd %in% c(6,7)) %>% 
  ggplot(aes(x=YInd,y=ZVal, group=factor(NewAgeInd)))+
  geom_line(aes(color=factor(NewAgeInd)))+
  scale_x_continuous(breaks=seq(1, 51, 2))+
  scale_color_brewer(palette="Dark2")

LambdaVecIt %>% 
  ggplot(data=., aes(x=Year, y=log(Rate), group=NewAgeInd))+
  geom_line(aes(col=NewAgeInd))+
  ylab("Log Death Rate")+
  facet_wrap(~NewAgeInd,scales = "free")


#### 2.1 Put data into right form for NIMBLE ##################################
NAge <- nrow(ZMatIt)
NYear <- ncol(ZMatIt)

NimbleConstIt <- list("N_AgeGroups"=NAge,
                       "N_Year"=NYear)

NimbleDataIt <- list("ZMat"=ZMatIt) 
##### 2.2 QR APPROACH #########################################################
##### 2.2.1 Estimate Model Parameters ##########################################
OwnQR_It <-  nimbleModel(code=LCJumpOwn_QR, 
                          constants = NimbleConstIt, 
                          data=NimbleDataIt)

cOwnQR_It <- configureMCMC(OwnQR_It, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("k","sigma_eps",
                                         "QMat","J","Y_t",
                                         "BMat","b1","b2",
                                        #"beta","betaJump","sdY","muY",
                                        "N_t","p",
                                         "sigma_time",
                                         "drift","a",
                                         "muY","sdY"
                                         ),
                            enableWAIC = TRUE) #to compare Models

###################### 2.2.2 Adjust Samplers ###################################
params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
                       paste0("b2[",1:NAge,"]"),
                       paste0("k[",1:NYear,"]"),
                       paste0("Y_t[",1:(NYear+1),"]"),
                       "p","muY","sdY")
                       

cOwnQR_It$removeSamplers(params_to_remove)

cOwnQR_It$addSampler(target= paste0("b1[1:",NAge,"]"),
                      type="AF_slice") 

cOwnQR_It$addSampler(target= paste0("b2[1:",NAge,"]"),
                      type="AF_slice")

for(j in 1:NYear){
  cOwnQR_It$addSampler(target = paste0("k[",j,"]"),
                        type = 'slice')
}

cOwnQR_It$addDefaultSampler(nodes = "p",useConjugacy = TRUE)
cOwnQR_It$addSampler(target= "muY",
                     type="slice")
cOwnQR_It$addSampler(target= "sdY",
                     type="slice")

for(j in 3:(NYear+1)){
  cOwnQR_It$addSampler(target = paste0("Y_t[",j,"]"),
                        type = 'slice')
}

bOwnQR_It <- buildMCMC(cOwnQR_It)
comOwnQR_It <- compileNimble(OwnQR_It,
                              bOwnQR_It)

SamplesOwnQR_It <- runMCMC(comOwnQR_It$bOwnQR_It, 
                            niter = 5000,
                            thin=2,
                            nburnin = 2000, 
                            nchains = 2)

SummaryOutput(SamplesOwnQR_It, 
              params=c("QMat","drift","p","sigma_eps",
                       "muY","sdY","sigma_time","a","N_t")) %>% 
  print(n=250) 


####### 2.3. SUM TO ONE APPROACH ###############################################
##### 2.3.1 Use Own Constraint Parametersization ################################
#Using a sum to zero constraint
OwnMod_It_Rep <-  nimbleModel(code=LCJumpOwn_OC_CornerK, 
                          constants = NimbleConstIt, 
                          data=NimbleDataIt)

cOwnMod_It_Rep <- configureMCMC(OwnMod_It_Rep, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("k",
                                         "sigma_eps",
                                        "beta","betaJump","sdY","muY",
                                        "N_t","p",
                                         "sigma_time",
                                         "drift","a",
                                        "mu",#"sigma_squared" #for PPC and WAIC
                                         "Y_t","J"
                                         ),
                            enableWAIC = TRUE) 

############# 2.3.2 Adjust Samplers ############################################
params_to_remove <-  c(paste0("k[",1:NYear,"]"),
                       paste0("b1[",1:NAge,"]"),
                       paste0("b2[",1:NAge,"]"),
                       paste0("Y_t[",1:(NYear+1),"]"),
                       "a","muY","sdY")
                                        
cOwnMod_It_Rep$removeSamplers(params_to_remove)
cOwnMod_It_Rep$addSampler(nodes = "muY",type="slice")
cOwnMod_It_Rep$addSampler(nodes = "a",type="slice")
cOwnMod_It_Rep$addSampler(nodes = "sdY",type="slice")
                                        
cOwnMod_It_Rep$addSampler(target= paste0("b1[1:",NAge,"]"),
                          type="AF_slice") 
                                        
cOwnMod_It_Rep$addSampler(target= paste0("b2[1:",NAge,"]"),
                          type="AF_slice")
                                        
for(j in 2:NYear){
  cOwnMod_It_Rep$addSampler(target = paste0("k[",j,"]"),
                            type = 'slice')
}

cOwnMod_It_Rep$addDefaultSampler(nodes=c("Y_t[3:42]"), 
                                 useConjugacy = TRUE)
                                        
bOwnMod_It_Rep <- buildMCMC(cOwnMod_It_Rep)
comOwnMod_It_Rep <- compileNimble(OwnMod_It_Rep,
                                  bOwnMod_It_Rep)
                                      
SamplesOwn_It <- runMCMC(comOwnMod_It_Rep$bOwnMod_It_Rep,
                         niter = 15000,
                         thin=10,
                         nburnin = 5000,
                         nchains = 2)
                                        
SummaryOutput(SamplesOwn_It, params=c("drift","p",
                                      "sigma_eps","beta","muY","sigma_time",
                                      "sdY","a")) %>% print(n=450)  

save(SamplesOwn_It,
     file = file.path(getwd(),"Results/SamplesIt_Own.RData"))


######### 2.4 Liu,Li Model #####################################################
LiLi_It_OC <-  nimbleModel(code=LiuLi_Repara_Diff_OC_V2, 
                           constants = NimbleConstIt, 
                           data=NimbleDataIt)

cLiLi_It_OC <- configureMCMC(LiLi_It_OC, 
                             print = TRUE, useConjugacy = TRUE,
                             monitors = c(#"k",
                                          "sigma_eps",
                                          "beta","betaJump", "N_t","p",
                                          "sigma_time",
                                          "drift",#"b1","b2",
                                          "muY","sdY","mu","sigma_squared"
                                          #"Y_t"
                             ),
                             enableWAIC = TRUE) #to compare Models

cLiLi_It_OC$removeSamplers(paste0("b1[",1:NAge,"]"),
                           paste0("b2[",1:NAge,"]"),
                           "muY","sdY")

cLiLi_It_OC$addSampler(target= paste0("b1[1:",NAge,"]"),
                          type="AF_slice") 

cLiLi_It_OC$addSampler(target= paste0("b2[1:",NAge,"]"),
                          type="AF_slice")
cLiLi_It_OC$addSampler(nodes = "muY",type="slice")
cLiLi_It_OC$addSampler(nodes = "sdY",type="slice")

bLiLi_It_OC <- buildMCMC(cLiLi_It_OC)
comLiLi_It_OC <- compileNimble(LiLi_It_OC,
                               bLiLi_It_OC)


SamplesLiLi_It <- runMCMC(comLiLi_It_OC$bLiLi_It_OC, 
                          niter = 4000,
                          thin=1,
                          nburnin = 2000, 
                          nchains = 2)


SummaryOutput(SamplesLiLi_It, 
              params=c("QMat","drift","p","sigma_eps",
                       "muY","sigma_time","beta","betaJump","sdY")) %>% 
  print(n=250) 


save(SamplesLiLi_It, 
     SamplesOwn_It, file= file.path(getwd(),"Results/SamplesIt_NoRep_2.RData"))


## Compare WAIC's
comOwnMod_It_Rep$bOwnMod_It_Rep$getWAIC()
comLiLi_It_OC$bLiLi_It_OC$getWAIC()


LikeMatOwnIt <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesOwn_It),
                                  n = length(ZMatIt),SingleVar = FALSE,
                                  ZMat = ZMatIt)

LikeMatLiuIt <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesLiLi_It),
                                  n = length(ZMatIt),SingleVar = FALSE,
                                  ZMat = ZMatIt)

WAICOwnIt <- loo::waic(log(LikeMatOwnIt))
WAICLiuIt <- loo::waic(log(LikeMatLiuIt))

LOOLiuIt <- loo::loo(log(LikeMatLiuIt)) # loo IC = -2 elpd_loo
LOOOwnIt <- loo::loo(log(LikeMatOwnIt))

CompDataFrameIt <- data.frame("Model"=c("Own","Liu-Li"),
                              "LOOCV"= c(LOOOwnIt$estimates[3,1],
                                         LOOLiuIt$estimates[3,1]),
                              "WAIC"=c(WAICOwnIt$estimates[3,1],
                                       WAICLiuIt$estimates[3,1]))

save(CompDataFrameIt,
     file = file.path(getwd(),"Results/WAIC_It_Rep2.RData"))


##### 2.5 LC Model Diff ########################################################
LC_It <-  nimbleModel(code=LC_Diff, 
                      constants = NimbleConstIt, 
                      data=NimbleDataIt)

cLC_It <- configureMCMC(LC_It, 
                        print = TRUE, useConjugacy = TRUE,
                        monitors = c("k","sigma_eps",
                                     "beta",
                                     "sigma_time",
                                     "drift","b1","mu"
                        ),
                        enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(
  paste0("b1[",1:NAge,"]")
)
cLC_It$removeSamplers(params_to_remove)
cLC_It$addSampler(target= c("b1[1:10]"),
                  type="AF_slice")
bLC_It <- buildMCMC(cLC_It)
comLC_It <- compileNimble(LC_It,
                          bLC_It)


Samples_LC_It<- runMCMC(comLC_It$bLC_It, 
                        niter = 4000,
                        thin=1,
                        nburnin = 2000, 
                        nchains = 2)

SummaryOutput(Samples_LC_It) %>% print(n=100)

save(Samples_LC_It,
     file = file.path(getwd(),"Results/SamplesIt_LC.RData"))




# ####################### 3. Spain #############################################
## Some Plots
LambdaVecSp  %>% 
  mutate("DeathRate"=log(Rate)) %>% 
  filter(Year > 1980) %>% 
  filter(NewAgeInd %in% c(9:18)) %>% 
  ggplot(aes(x=Year,y=DeathRate, group=factor(NewAgeInd)))+
  geom_line(aes(color=factor(NewAgeInd)))+
  facet_wrap(~NewAgeInd, scales="free")


LambdaVecSp  %>% 
  mutate("YInd"=match(Year,unique(Year))) %>% 
  ggplot(aes(x=YInd,y=ZVal, group=factor(NewAgeInd)))+
  geom_line(aes(color=factor(NewAgeInd)))+
  scale_x_continuous(breaks=seq(1, 51, 2))+
  #scale_color_brewer(palette="Dark2")
theme_set(theme_minimal(base_size = 10))


#### 3.1 Data for NIMBLE #######################################################
ZMatSp_LessAges <- ZMatSp[-c(1:8,19),]
ZMatSp
NAge <- nrow(ZMatSp)
NYear <- ncol(ZMatSp)

NimbleConstSp <- list("N_AgeGroups"=NAge,
                      "N_Year"=NYear)

NimbleDataSp <- list("ZMat"=ZMatSp) 
###### 3.2 SUM TO ZERO CONSTRAINT ##############################################
OwnMod_Sp_Rep <-  nimbleModel(code=LCJumpOwn_OC_CornerK, 
                              constants = NimbleConstSp, 
                              data=NimbleDataSp)

cOwnMod_Sp_Rep <- configureMCMC(OwnMod_Sp_Rep, 
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
                                             "Y_t","J"
                                ),
                                enableWAIC = TRUE) 

############# 3.2.1 Adjust Samplers ############################################
params_to_remove <-  c(
  paste0("k[",2:NYear,"]"),
  paste0("Y_t[",1:(NYear-1),"]"),
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  "a",
  "muY","sdY","sigma_time","sigma_eps"
)

cOwnMod_Sp_Rep$removeSamplers(params_to_remove)
cOwnMod_Sp_Rep$addSampler(nodes = "muY",type="slice")
cOwnMod_Sp_Rep$addSampler(nodes = "a",type="slice")
cOwnMod_Sp_Rep$addSampler(nodes = "sdY",type="slice")

cOwnMod_Sp_Rep$addSampler(target = c("sigma_time","sigma_eps"),
                          type="AF_slice")

cOwnMod_Sp_Rep$addSampler(target= paste0("b1[1:",NAge,"]"),
                          type="AF_slice") 

cOwnMod_Sp_Rep$addSampler(target= paste0("b2[1:",NAge,"]"),
                          type="AF_slice")

for(j in 2:NYear){
  cOwnMod_Sp_Rep$addSampler(target = paste0("k[",j,"]"),
                            type = 'slice')
}
for(j in 3:(NYear+1)){
  cOwnMod_Sp_Rep$addSampler(target = paste0("Y_t[",j,"]"),
                            type = 'slice')
}
# cOwnMod_Sp_Rep$addDefaultSampler(nodes=c("Y_t[3:42]"), 
#                                  useConjugacy = TRUE)

bOwnMod_Sp_Rep <- buildMCMC(cOwnMod_Sp_Rep)
comOwnMod_Sp_Rep <- compileNimble(OwnMod_Sp_Rep,
                                  bOwnMod_Sp_Rep)

SamplesOwn_Sp <- runMCMC(comOwnMod_Sp_Rep$bOwnMod_Sp_Rep, 
                         niter = 25000,
                         thin=20,
                         nburnin = 5000, 
                         nchains = 2)

SummaryOutput(SamplesOwn_Sp, 
              params=c("drift","p","sigma_eps","beta",
                       "muY","sigma_time","sdY","a","betaJump","N_t")) %>% 
  print(n=450)  

ts.plot(SamplesOwn_Sp$chain1[,"muY"])
lines(SamplesOwn_Sp$chain2[,"muY"],col="red")

SummaryOutput(SamplesOwn_Sp, 
              params = "k") %>% 
  select(mean) %>% summarise(sd(mean))



############ 3.3 Liu,Li Model###################################################
LiuLi_Sp_OC <-  nimbleModel(code=LiuLi_Repara_Diff_OC_V2, 
                            constants = NimbleConstSp, 
                            data=NimbleDataSp)

cLiuLi_Sp_OC <- configureMCMC(LiuLi_Sp_OC, 
                              print = TRUE, useConjugacy = TRUE,
                              monitors = c(#"k",
                                           "sigma_eps",
                                           "beta","betaJump", "N_t","p",
                                           "sigma_time",
                                           "drift",#"b1","b2",
                                           "muY","sdY",
                                           "mu","sigma_squared"
                                           #"Y_t"
                              ),
                              enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(
  paste0("b1[",1:NAge,"]"),
  paste0("b2[",1:NAge,"]"),
  #paste0("k[",1:NYear,"]"),
  "muY","sdY"
)

cLiuLi_Sp_OC$removeSamplers(params_to_remove)
cLiuLi_Sp_OC$addSampler(nodes = "muY",type="slice")
cLiuLi_Sp_OC$addSampler(nodes = "sdY",type="slice")

cLiuLi_Sp_OC$addSampler(target= paste0("b1[1:",NAge,"]"),
                        type="AF_slice") 

cLiuLi_Sp_OC$addSampler(target= paste0("b2[1:",NAge,"]"),
                        type="AF_slice")


for(j in 1:NYear){
  cLiuLi_Sp_OC$addSampler(target = paste0("k[",j,"]"),
                          type = 'slice')
}

bLiuLi_Sp_OC <- buildMCMC(cLiuLi_Sp_OC)
comLiuLi_Sp_OC <- compileNimble(LiuLi_Sp_OC,
                                bLiuLi_Sp_OC)

SamplesLiuLi_Sp <- runMCMC(comLiuLi_Sp_OC$bLiuLi_Sp_OC, 
                           niter = 4000,
                           thin=1,
                           nburnin = 2000, 
                           nchains = 2)

SummaryOutput(SamplesLiuLi_Sp, 
              params=c("QMat","drift","p","sigma_eps",
                       "muY","sigma_time","beta","betaJump","sdY","a")) %>% 
  print(n=250) 


save(SamplesLiuLi_Sp_Rep, 
     SamplesOwn_Sp_Rep, file= file.path(getwd(),"Results/SamplesSp_Rep2.RData"))

save(SamplesLiuLi_Sp, 
     SamplesOwn_Sp, file= file.path(getwd(),"Results/SamplesSp_NoRep.RData"))


#### Caclulation of WAIC and loo 
LikeMatOwn <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesOwn_Sp),
                                  n = length(ZMatSp),SingleVar = TRUE,
                                  ZMat = ZMatSp)
LikeMatLiu <- LikelihoodMatrixFun(Samples = do.call(rbind, SamplesLiuLi_Sp),
                                  n = length(ZMatSp),SingleVar = TRUE,
                                  ZMat = ZMatSp)

WAICOwn <- loo::waic(log(LikeMatOwn))
WAICLiu <- loo::waic(log(LikeMatLiu))

LOOLiu <- loo::loo(log(LikeMatLiu)) #loo IC = -2 elpd_loo
LOOOwn <- loo::loo(log(LikeMatOwn))


#Compare WAIC's
Details_LiuLi_Sp <- comLiuLi_Sp_OC$bLiuLi_Sp_OC$getWAIC()
Details_Own_Sp <- comOwnMod_Sp_Rep$bOwnMod_Sp_Rep$getWAIC()

CompDataFrameSp <- data.frame("Model"=c("Own","Liu-Li"),
                              "LOOCV"= c(LOOOwn$estimates[3,1],
                                         LOOLiu$estimates[3,1]),
                              "WAIC"=c(WAICOwn$estimates[3,1],
                                       WAICLiu$estimates[3,1]))

save(CompDataFrameSp,
     file = file.path(getwd(),"Results/WAIC_Sp_NoRep.RData"))

##### 3.4 LC Model Diff ########################################################
LC_Sp <-  nimbleModel(code=LC_Diff, 
                      constants = NimbleConstSp, 
                      data=NimbleDataSp)

cLC_Sp <- configureMCMC(LC_Sp, 
                        print = TRUE, useConjugacy = TRUE,
                        monitors = c("k","sigma_eps",
                                     "beta",
                                     "sigma_time",
                                     "drift",#"b1",
                                     "mu"
                        ),
                        enableWAIC = TRUE) #to compare Models

params_to_remove <-  c(
  paste0("b1[",1:NAge,"]")
)
cLC_Sp$removeSamplers(params_to_remove)
cLC_Sp$addSampler(target= c("b1[1:10]"),
                  type="AF_slice")

bLC_Sp <- buildMCMC(cLC_Sp)
comLC_Sp <- compileNimble(LC_Sp,
                          bLC_Sp)


Samples_LC_Sp<- runMCMC(comLC_Sp$bLC_Sp, 
                        niter = 4000,
                        thin=1,
                        nburnin = 2000, 
                        nchains = 2)

SummaryOutput(Samples_LC_Sp, 
              params = c("drift","sigma_time","beta")) %>% 
  print(n=100)

save(Samples_LC_Sp,
     file = file.path(getwd(),"Results/SamplesSp_LC.RData"))

LikeMatLC <- LikelihoodMatrixFun(Samples = do.call(rbind, Samples_LC_Sp),
                                  n = length(ZMatSp),SingleVar = TRUE,
                                  ZMat = ZMatSp)

##### 3.5 QR APPROACH #########################################################
##### 3.5.1 Estimate Model Parameters ##########################################
OwnQR_Sp <-  nimbleModel(code=LCJumpOwn_QR, 
                         constants = NimbleConstSp, 
                         data=NimbleDataSp)

cOwnQR_Sp <- configureMCMC(OwnQR_Sp, 
                           print = TRUE, useConjugacy = FALSE,
                           monitors = c("k","sigma_eps",
                                        "QMat","J","Y_t",
                                        "BMat","b1","b2",
                                        #"beta","betaJump","sdY","muY",
                                        "N_t","p",
                                        "sigma_time",
                                        "drift","a",
                                        "muY","sdY"
                           ),
                           enableWAIC = TRUE) #to compare Models

###################### 3.5.2 Adjust Samplers ###################################
params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
                       paste0("b2[",1:NAge,"]"),
                       paste0("k[",1:NYear,"]"),
                       paste0("Y_t[",1:(NYear+1),"]"),
                       "p","muY","sdY")


cOwnQR_Sp$removeSamplers(params_to_remove)

cOwnQR_Sp$addSampler(target= paste0("b1[1:",NAge,"]"),
                     type="AF_slice") 

cOwnQR_Sp$addSampler(target= paste0("b2[1:",NAge,"]"),
                     type="AF_slice")

for(j in 1:NYear){
  cOwnQR_Sp$addSampler(target = paste0("k[",j,"]"),
                       type = 'slice')
}

cOwnQR_Sp$addDefaultSampler(nodes = "p",useConjugacy = TRUE)
cOwnQR_Sp$addSampler(target= "muY",
                     type="slice")
cOwnQR_Sp$addSampler(target= "sdY",
                     type="slice")

for(j in 3:(NYear+1)){
  cOwnQR_Sp$addSampler(target = paste0("Y_t[",j,"]"),
                       type = 'slice')
}

bOwnQR_Sp <- buildMCMC(cOwnQR_Sp)
comOwnQR_Sp <- compileNimble(OwnQR_Sp,
                             bOwnQR_Sp)

SamplesOwnQR_Sp <- runMCMC(comOwnQR_Sp$bOwnQR_Sp, 
                           niter = 5000,
                           thin=2,
                           nburnin = 2000, 
                           nchains = 2)

SummaryOutput(SamplesOwnQR_Sp, 
              params=c("QMat","drift","p","sigma_eps",
                       "muY","sdY","sigma_time","a","N_t")) %>% 
  print(n=250) 



###### BELOW EVERYTHING OLD !!!! ###############################################
# # ########## 5. ENGLAND AND WALES ###############################################
# 
# #Get Death Data from 2021 and Transform into Long Format
# pacman::p_load(openxlsx)
# 
# #1. )Load Data
# ## Load Death Data from ONS 
# DeathEW_21 <- read.xlsx(xlsxFile = file.path(getwd(),
#                                              "Data/Deaths_EnglandWales_2021.xlsx"),
#                         colNames = TRUE, sheet="TotalDeathsAge")
# 
# 
# Death_EW_HMD <- read.table(file = file.path(getwd(),"Data/Deaths_5x1_EnglandWales.txt"),
#                            header = TRUE)
# 
# 
# #Load Death Data from HMD 
# Pop_EW_HMD <- read.table(file = file.path(getwd(),"Data/Exposures_5x1_EnglandWales.txt"),
#                          header = TRUE)
# 
# #Load Population from ONS
# PopEW_21_ONS <- read.xlsx(xlsxFile = file.path(getwd(),
#                                                "Data/UK_MidYearPop2021.xlsx"),
#                           colNames = TRUE, sheet="MYE1", 
#                           startRow = 35)
# 
# AgeW <- 10
# 
# #2.) combine Death Data
# AgeLabFun_EW <- function(AgeWidth=10,
#                          AgeLabels){
#   # Age Groups of 10. Starting from 0-9, 10-19,..., 90+
#   if(AgeWidth == 10){
#     HelperAgeLab <- data.frame("AgeOld" = unique(AgeLabels)[-1],
#                                "AgeNew" = c(rep(1:10,each=10),rep(10,6)))
#   }else {
#     #Age Groups of 5. Starting from 0-4, 5-9,..., 90+
#     HelperAgeLab <- data.frame("AgeOld" = unique(AgeLabels)[-1],
#                                "AgeNew" = c(rep(1:19,each=5),rep(19,11)))
#   }
# }
# 
# AgeHelperEW <- AgeLabFun_EW(AgeWidth = AgeW,
#                             AgeLabels = colnames(DeathEW_21)[-1])
# 
# DeathEW_21_Grouped <- 
#   DeathEW_21 %>% pivot_longer(., cols = 2:ncol(.),
#                               names_to = "Ages",
#                               values_to = "Deaths") %>% 
#   filter(Year.of.registration == 2021 & 
#            Ages != "All.ages") %>% 
#   mutate("AInd"= 1:length(unique(Ages))) %>% 
#   mutate("NewAgeInd"=AgeHelperEW$AgeNew[AInd])  %>% #add new Age ID
#   rename(Year = Year.of.registration) %>% 
#   group_by(NewAgeInd,Year) %>% #group by new Age ID
#   summarise(across(Deaths,sum))
# 
# 
# YearMin <- 1980
# 
# HelperHMD_EW <- AgeLabFun_HMD(AgeWidth = AgeW)
# 
# DeathsGrouped_HMD_EW <- Death_EW_HMD %>% filter(Year > YearMin) %>% 
#   select(Year,Age,Total) %>% 
#   mutate("AInd"=match(Age, unique(Age))) %>% #Get Age ID
#   mutate("NewAgeInd"=HelperHMD_EW$AgeNew[AInd]) %>% #add new Age ID
#   group_by(NewAgeInd,Year) %>% #group by new Age ID
#   summarise(across(Total,sum)) %>% 
#   rename(Deaths=Total) %>% 
#   select(NewAgeInd, Year, Deaths)
# 
# DeathsEW_Tot <- full_join(x=DeathsGrouped_HMD_EW,
#                           y = DeathEW_21_Grouped, 
#                           by=c("NewAgeInd","Year","Deaths"),
#                           na_matches="never")
# 
# #3.) Combine Population Data
# PopEW_21 <- PopEW_21_ONS %>% 
#   select(Age.Groups, England.and.Wales) %>% 
#   mutate("NewAgeInd"=if(AgeW==10){c(rep(1:9, each=2),10)}else{1:nrow(.)},
#          "Year"=2021) %>% 
#   rename(Pop = England.and.Wales) %>%  #rename column
#   group_by(NewAgeInd,Year) %>% #group by new Age ID
#   summarise(across(Pop,sum))
# 
# PopGrouped_HMD_EW <- 
#   Pop_EW_HMD %>% filter(Year > YearMin) %>% 
#   select(Year,Age,Total) %>% 
#   mutate("AInd"=match(Age, unique(Age))) %>% #Get Age ID
#   mutate("NewAgeInd"=HelperHMD_EW$AgeNew[AInd]) %>% #add new Age ID
#   group_by(NewAgeInd,Year) %>% #group by new Age ID
#   summarise(across(Total,sum)) %>% 
#   rename(Pop=Total)
# 
# 
# PopEW_Tot <- full_join(x=PopGrouped_HMD_EW,
#                        y = PopEW_21, 
#                        by=c("NewAgeInd","Year","Pop"),
#                        na_matches="never")
# 
# 
# TotDataEW <- dplyr::inner_join(x = DeathsEW_Tot,
#                                y = PopEW_Tot, 
#                                by =c("NewAgeInd","Year"),
#                                na_matches="never")
# 
# ################# 5.1. Estimate Model ##########################################
# ## Only take the oldest Age Groups 
# LambdaVecEW <- TotDataEW %>% 
#   mutate("Rate"=Deaths/Pop) %>% 
#   mutate("ZVal"=c(NA,diff(log(Rate)))) %>%
#   arrange(Year) %>% 
#   filter(NewAgeInd %in% c(4:9)) #only middle Age groups /4:9
# 
# 
# LambdaVecEW  %>%
#   mutate("DeathRate"=log(Rate)) %>%
#   mutate("AgeGroups"=factor(NewAgeInd)) %>%
#   ggplot(aes(x=Year,y=DeathRate, group=AgeGroups))+
#   geom_line(aes(color=AgeGroups))
# 
# 
# 
# LambdaVecEW  %>% 
#   mutate("YInd"=match(Year,unique(Year))) %>% 
#   ggplot(aes(x=YInd,y=ZVal, group=factor(NewAgeInd)))+
#   geom_line(aes(color=factor(NewAgeInd)))+
#   scale_x_continuous(breaks=seq(1, 60, 2))+
#   scale_color_brewer(palette="Paired")
# 
# LambdaMatEW <- matrix(LambdaVecEW$Rate,
#                       nrow = length(unique(LambdaVecEW$NewAgeInd)),
#                       byrow = FALSE)
# 
# #Differenced Log Death Rates
# ZMatEW <- apply(log(LambdaMatEW), 1, diff) %>% t()
# 
# 
# NAge <- nrow(ZMatEW)
# NYear <- ncol(ZMatEW)
# 
# NimbleConstEW <- list("N_AgeGroups"=NAge,
#                       "N_Year"=NYear)
# 
# NimbleDataEW <- list("ZMat"=ZMatEW) 
# 
# 
# ### Attention! only Age groups 4 to 10 will be in the Analysis
# ####### 5.1 QR APPROACH ########################################################
# #### Estimate Model (compare the Own Reparameterization with the Simple Model)
# 
# OwnQR_EW <-  nimbleModel(code=OwnMod_Repara_Diff_QR, 
#                          constants = NimbleConstEW, 
#                          data=NimbleDataEW, 
#                          buildDerivs = FALSE)
# 
# cOwnQR_EW <- configureMCMC(OwnQR_EW, 
#                            print = TRUE, useConjugacy = FALSE,
#                            monitors = c("k","sigma_eps",
#                                          "QMat","BMat","b1","b2",
#                                         #"beta","betaJump","sdY","muY",
#                                         "muY","sdY",
#                                         "N_t","p",
#                                         "sigma_time",
#                                         "drift","a"
#                            ),
#                            enableWAIC = TRUE) #to compare Models
# 
# ######### 5.1.2 Adjust Samplers###################################################
# params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
#                        paste0("b2[",1:NAge,"]"),
#                        paste0("k[",1:NYear,"]"),
#                        "muY","sdY","a","p"
#                        # "muY","sdY"
# )
# 
# cOwnQR_EW$removeSamplers(params_to_remove)
# cOwnQR_EW$addDefaultSampler(nodes=c("p"), useConjugacy = TRUE)
# 
# cOwnQR_EW$addSampler(target= paste0("b1[1:",NAge,"]"),
#                      type="AF_slice") 
# 
# cOwnQR_EW$addSampler(target= paste0("b2[1:",NAge,"]"),
#                      type="AF_slice")
# 
# 
# cOwnQR_EW$addSampler(target = "sdY", type="RW", log=TRUE)
# cOwnQR_EW$addSampler(target = "a", type="slice")
# cOwnQR_EW$addSampler(target = "muY", type="slice")
# 
# 
# for(j in 1:NYear){
#   cOwnQR_EW$addSampler(target = paste0("k[",j,"]"),
#                        type = 'slice')
# }
# 
# 
# 
# bOwnQR_EW <- buildMCMC(cOwnQR_EW)
# comOwnQR_EW <- compileNimble(OwnQR_EW,
#                              bOwnQR_EW)
# 
# 
# SamplesOwnQR_EW <- runMCMC(comOwnQR_EW$bOwnQR_EW, 
#                            niter = 4000,
#                            thin=1,
#                            nburnin = 2000, 
#                            nchains = 1)
# 
# SummaryOutput(SamplesOwnQR_EW, 
#               params=c("drift","p","sigma_eps",
#                        "muY","sdY","sigma_time","a",
#                        "QMat","k")) %>% 
#   print(n=200) 
# 
# ts.plot(SamplesOwnQR_EW[,"muY"])
# 
# comOwnQR_EW$bOwnQR_EW$getWAIC()
# 
# ##### SIMPLE MODEL #####
# OwnQR_EW_S <- nimbleModel(code=LCJumpOwn_QR_Simple, 
#                           constants = NimbleConstEW, 
#                           data=NimbleDataEW, 
#                           buildDerivs = TRUE)
# 
# 
# 
# cOwnQR_EW_S <- configureMCMC(OwnQR_EW_S, 
#                              print = TRUE, useConjugacy = FALSE,
#                              monitors = c("k","sigma_eps",
#                                           "QMat","BMat","b1","b2",
#                                           "N_t","p",
#                                           "sigma_time","a",
#                                           "drift"
#                              ),
#                              enableWAIC = TRUE) #to compare Models
# 
# 
# 
# 
# cOwnQR_EW_S$removeSamplers(c(paste0("b1[",1:NAge,"]"),
#                              paste0("b2[",1:NAge,"]"),
#                              paste0("k[",1:NYear,"]"),"a","p"
#                              # "muY","sdY"
# ))
# 
# cOwnQR_EW_S$addSampler(target= paste0("b1[1:",NAge,"]"),
#                        type="AF_slice") 
# 
# cOwnQR_EW_S$addSampler(target= paste0("b2[1:",NAge,"]"),
#                        type="AF_slice")
# 
# for(j in 1:NYear){
#   cOwnQR_EW_S$addSampler(target = paste0("k[",j,"]"),
#                          type = 'slice')
# }
# 
# cOwnQR_EW_S$addSampler(target = "a", type="slice")
# cOwnQR_EW_S$addDefaultSampler(nodes = "p", useConjugacy = TRUE)
# 
# 
# bOwnQR_EW_S <- buildMCMC(cOwnQR_EW_S)
# comOwnQR_EW_S <- compileNimble(OwnQR_EW_S,
#                                bOwnQR_EW_S)
# 
# SamplesOwnQR_EW_S <- runMCMC(comOwnQR_EW_S$bOwnQR_EW_S, 
#                              niter = 4000,
#                              thin=1,
#                              nburnin = 2000, 
#                              nchains = 1)
# 
# SummaryOutput(SamplesOwnQR_EW_S, 
#               params=c("QMat","drift","p","sigma_eps",
#                        "muY","sdY","sigma_time","a","N_t","Y_t",
#                        "J")) %>% 
#   print(n=300) 
# 
# 
# 
# comOwnQR_EW$bOwnQR_EW$getWAIC() 
# comOwnQR_EW_S$bOwnQR_EW_S$getWAIC()
# 
# 
# 
# ######### 5.2 SUM-TO-ONE-APPROACH  #############################################
# #### 5.2.1 Reparameterization ###################
# #Note: change of prior for muY for UK Model, tighter priors needed
# # muY ~ dchisq(df = 2)
# OwnMod_EW_Rep <-  nimbleModel(code=OwnMod_Repara_Diff_OC, 
#                          constants = NimbleConstEW, 
#                          data=NimbleDataEW, 
#                          buildDerivs = FALSE)
# 
# 
# cOwnMod_EW_Rep <- configureMCMC(OwnMod_EW_Rep, 
#                            print = TRUE, useConjugacy = FALSE,
#                            monitors = c("k","sigma_eps",
#                                         "beta","betaJump",
#                                         "muY","sdY",
#                                         "N_t","p",
#                                         "sigma_time",
#                                         "drift","a"
#                            ),
#                            enableWAIC = TRUE) #to compare Models
# 
# params_to_remove <-  c(
#   #"p","drift",
#   paste0("b1[",1:6,"]"),
#   paste0("b2[",1:6,"]"),
#   paste0("k[",1:40,"]"),
#   "a","muY","sdY","sigma_time"
# )
# 
# cOwnMod_EW_Rep$removeSamplers(params_to_remove)
# 
# cOwnMod_EW_Rep$addSampler(nodes = "muY",type="slice")
# cOwnMod_EW_Rep$addSampler(nodes = "a",type="slice")
# cOwnMod_EW_Rep$addSampler(nodes = "sdY",type="slice")
# cOwnMod_EW_Rep$addSampler(nodes = "sigma_time",type="slice")
# cOwnMod_EW_Rep$addSampler(target= c("b1[1:6]"),
#                           type="AF_slice")
# 
# cOwnMod_EW_Rep$addSampler(target= c("b2[1:6]"),
#                           type="AF_slice")
# for(j in 1:NYear){
#   cOwnMod_EW_Rep$addSampler(target = paste0("k[",j,"]"),
#                         type = 'slice')
# }
# 
# bOwnMod_EW_Rep <- buildMCMC(cOwnMod_EW_Rep)
# comOwnMod_EW_Rep <- compileNimble(OwnMod_EW_Rep,
#                               bOwnMod_EW_Rep)
# 
# 
# SamplesOwn_EW_Rep <- runMCMC(comOwnMod_EW_Rep$bOwnMod_EW_Rep, 
#                             niter = 4000,
#                             thin=1,
#                             nburnin = 2000, 
#                             nchains = 1)
# 
# SummaryOutput(SamplesOwn_EW_Rep, 
#               params=c("drift","p","sigma_eps",
#                        "muY","sigma_time","beta","betaJump","sdY","a","k")) %>% 
#   print(n=250)  
# 
# save(SamplesOwn_EW_Rep,
#      file = file.path(getwd(),"Results/SamplesEW_Own.RData"))
# 
# comOwnMod_EW_Rep$bOwnMod_EW_Rep$getWAIC()
# 
# ##### 5.2.2 Sum to One standard Repara #######################################
# OwnMod_EW_Rep <-  nimbleModel(code=LCJumpOwn_OC, 
#                               constants = NimbleConstEW, 
#                               data=NimbleDataEW, 
#                               buildDerivs = FALSE)
# 
# 
# cOwnMod_EW_Rep <- configureMCMC(OwnMod_EW_Rep, 
#                                 print = TRUE, useConjugacy = TRUE,
#                                 monitors = c("k","sigma_eps",
#                                              "beta","betaJump",
#                                              "muY","sdY",
#                                              "N_t","p",
#                                              "sigma_time",
#                                              "drift","a",
#                                              "Y_t"
#                                 ),
#                                 enableWAIC = TRUE) #to compare Models
# 
# params_to_remove <-  c(
#   paste0("b1[",1:6,"]"),
#   paste0("b2[",1:6,"]"),
#   "a","muY","sdY","sigma_time"
# )
# 
# 
# cOwnMod_EW_Rep$removeSamplers(params_to_remove)
# 
# cOwnMod_EW_Rep$addSampler(nodes = "muY",type="slice")
# cOwnMod_EW_Rep$addSampler(nodes = "a",type="slice")
# cOwnMod_EW_Rep$addSampler(nodes = "sdY",type="slice")
# cOwnMod_EW_Rep$addSampler(nodes = "sigma_time",type="slice")
# cOwnMod_EW_Rep$addSampler(target= c("b1[1:6]"),
#                       type="AF_slice")
# 
# cOwnMod_EW_Rep$addSampler(target= c("b2[1:6]"),
#                       type="AF_slice")
# 
# bOwnMod_EW_Rep <- buildMCMC(cOwnMod_EW_Rep)
# comOwnMod_EW_Rep <- compileNimble(OwnMod_EW_Rep,
#                                   bOwnMod_EW_Rep)
# 
# 
# SamplesOwn_EW_Rep <- runMCMC(comOwnMod_EW_Rep$bOwnMod_EW_Rep, 
#                              niter = 4000,
#                              thin=1,
#                              nburnin = 2000, 
#                              nchains = 2)
# SummaryOutput(SamplesOwn_EW_Rep, 
#               params=c("drift","p","sigma_eps",
#                        "muY","sigma_time","beta","betaJump","sdY","a")) %>% 
#   print(n=250)  
# 
# #### 5.3 Liu,Li Model #########################################################
# # Note, change prior to p ~ dbern(1,25), a ~ dbern(1,10)
# LiLi_EW_OC <-  nimbleModel(code=LiuLi_Repara_Diff_OC, 
#                            constants = NimbleConstEW, 
#                            data=NimbleDataEW)
# 
# cLiLi_EW_OC <- configureMCMC(LiLi_EW_OC, 
#                              print = TRUE, useConjugacy = TRUE,
#                              monitors = c("k","sigma_eps",
#                                           "beta","betaJump", "N_t","p",
#                                           "sigma_time",
#                                           "drift",#"b1","b2",
#                                           "muY","sdY"#,"meanJ","sdJ"
#                              ),
#                              enableWAIC = TRUE) #to compare Models
# 
# bLiLi_EW_OC <- buildMCMC(cLiLi_EW_OC)
# comLiLi_EW_OC <- compileNimble(LiLi_EW_OC,
#                                bLiLi_EW_OC)
# 
# 
# SamplesLiLi_EW <- runMCMC(comLiLi_EW_OC$bLiLi_EW_OC, 
#                           niter = 4000,
#                           thin=2,
#                           nburnin = 2000, 
#                           nchains = 2)
# 
# 
# SummaryOutput(SamplesLiLi_EW, 
#               params=c("QMat","drift","p","sigma_eps",
#                        "muY","sigma_time","beta","betaJump","sdY")) %>% 
#   print(n=250) 
# 
# ## Compare WAIC's
# comOwnMod_EW_Rep$bOwnMod_EW_Rep$getWAIC()
# comLiLi_EW_OC$bLiLi_EW_OC$getWAIC()
# 
# 
# # #### NOT IN USE



# ####################### 3. FRANCE ##############################################
# library(openxlsx)
# #Eurostat Death Data starts from 1998
# #Eurostat Population Data starts from 1991
# 
# #HMD Data from 1960 - 1991 respectively 1960 - 1998 for Pop/Deaths
# 
# DeathsFr_Eu <- read.xlsx(xlsxFile = 
#                           file.path(getwd(),"Data/DeathsFrance_Eurostat.xlsx"),
#                         startRow = 8, sheet = 3)
# 
# DeathsFr_Eu22 <- read.xlsx(xlsxFile = 
#                              file.path(getwd(),"Data/WeeklyDeaths2022_It_Fr.xlsx"),
#                            startRow = 25, sheet = "FranceTotal")
# 
# AgeW <- 10
# HelperAgeLab <- AgeLabFun_EU(AgeWidth = AgeW, 
#                              AgeLabels = DeathsFr_Eu$`AGE.(Labels)`)
# 
# DeathsFr_EU21 <- DeathsFr_Eu %>% 
#   filter(`AGE.(Labels)`!= "Total") %>% #remove total age group 
#   pivot_longer(.,                      #transform into long format 
#                cols = 3:ncol(DeathsFr_Eu), 
#                names_to = "Year", 
#                values_to = "Deaths") %>% 
#   mutate("AInd"=match(`AGE.(Labels)`, unique(`AGE.(Labels)`))) %>% #Get Age ID
#   mutate("NewAgeInd"=HelperAgeLab$AgeNew[AInd]) %>% #add new Age ID
#   group_by(NewAgeInd,Year) %>% #group by new Age ID
#   summarise(across(Deaths,sum)) %>% 
#   filter(Year >1997) %>% # data is available from 1998 onward
#   mutate(Year = as.numeric(Year)) # transform year into numeric
# 
# DeathsFr_Eu22 <- DeathsFr_Eu22 %>%
#   slice(-20) %>% #remove Unknown Deaths (Unknown age group)
#   mutate("NewAgeInd"= if(AgeW==5){c(1:19)}else{c(rep(1:9, each=2),10)}) %>% 
#   group_by(NewAgeInd) %>% 
#   summarise("Deaths"=sum(Total)) %>% 
#   mutate("Year" = 2022)
# 
# 
# DeathsGroupedFr_EU <- full_join(x=DeathsFr_EU21,
#                                 y = DeathsFr_Eu22, 
#                                 by=c("NewAgeInd","Year","Deaths"),
#                                 na_matches="never")
# 
# 
# 
# ## Get HMD Data of deaths and combine 
# # HMD DATA starts from 1960
# DeathsFr_HMD <- read.table(file = file.path(getwd(),"Data/Deaths_5x1_France.txt"),
#                            header = TRUE)
# YearMin <- 1980
# 
# HelperHMD <- AgeLabFun_HMD(AgeWidth = AgeW)
# 
# 
# DeathsGroupedFr_HMD <- DeathsFr_HMD %>% filter(Year > YearMin) %>% 
#   select(Year,Age,Total) %>% 
#   mutate("AInd"=match(Age, unique(Age))) %>% #Get Age ID
#   mutate("NewAgeInd"=HelperHMD$AgeNew[AInd]) %>% #add new Age ID
#   group_by(NewAgeInd,Year) %>% #group by new Age ID
#   summarise(across(Total,sum)) %>% 
#   filter(Year < 1998) %>% 
#   mutate(Deaths = Total) %>% 
#   select(NewAgeInd, Year, Deaths)
# 
# 
# DeathsFr_Tot <- full_join(x= DeathsGroupedFr_HMD, 
#                           y= DeathsGroupedFr_EU, 
#                           by=c("NewAgeInd","Year","Deaths"),
#                           na_matches="never")
# 
# 
# ## Get Exposure. Value is based on annual (January 1st) population estimates
# ## Population, that is Exposure estimates starts from 1991
# PopFr_EU <- read.xlsx(xlsxFile = 
#                         file.path(getwd(),"Data/PopulationFrance_Eurostat.xlsx"),
#                       startRow = 8, sheet = 3)
# PopFr_HMD <- read.table(file = file.path(getwd(),"Data/Exposures_5x1_France.txt"),
#                         header = TRUE)
# 
# PopGrouped_EU_Fr <- PopFr_EU %>% 
#   filter(`AGE.(Labels)`!= "Total") %>% #remove total age group 
#   pivot_longer(.,                      #transform into long format 
#                cols = 3:ncol(PopFr_EU), 
#                names_to = "Year", 
#                values_to = "Pop") %>% 
#   filter(Year > 1991) %>% 
#   mutate("AInd"=match(`AGE.(Labels)`, unique(`AGE.(Labels)`))) %>% #Get Age ID
#   mutate("NewAgeInd"=HelperAgeLab$AgeNew[AInd]) %>% #add new Age ID
#   group_by(NewAgeInd,Year) %>% #group by new Age ID
#   summarise(across(Pop,sum)) %>% 
#   mutate("Year"=as.numeric(Year)) # transform year into integer
# 
# PopGrouped_EU_HMD <- 
#   PopFr_HMD %>% filter(Year > YearMin & Year < 1992) %>% 
#   select(Year,Age,Total) %>% 
#   mutate("AInd"=match(Age, unique(Age))) %>% #Get Age ID
#   mutate("NewAgeInd"=HelperHMD$AgeNew[AInd]) %>% #add new Age ID
#   group_by(NewAgeInd,Year) %>% #group by new Age ID
#   summarise(across(Total,sum)) %>% 
#   mutate(Pop = Total) %>% 
#   select(-3)
# 
# PopFr_Tot <- rbind(PopGrouped_EU_HMD, 
#                    PopGrouped_EU_Fr)
# 
# #### Total Data
# TotDataFr <- dplyr::inner_join(x = DeathsFr_Tot,
#                               y = PopFr_Tot, 
#                               by =c("NewAgeInd","Year"),
#                               na_matches="never")
# 
# 
# ####### 3.1 Estimate Model #####################################################
# #Only Use Ages of 50 + 
# LambdaVec <- TotDataFr %>% 
#   mutate("Rate"=Deaths/Pop) %>% 
#   arrange(Year)
# 
# LambdaMatFr <- matrix(LambdaVec$Rate,
#                       nrow = length(unique(LambdaVec$NewAgeInd)), 
#                       ncol = length(unique(LambdaVec$Year)),
#                       byrow = FALSE)
# ## plot death rates
# LambdaVec %>% 
#   ggplot(data=., aes(x=Year, y=log(Rate), group=NewAgeInd))+
#   geom_line(aes(col=NewAgeInd))
# 
# 
# 
# #Differenced Log Death Rates
# ZMatFr <- apply(log(LambdaMatFr), 1, diff) %>% t()
# 
# NAge <- nrow(ZMatFr)
# NYear <- ncol(ZMatFr)
# 
# NimbleConstFr <- list("N_AgeGroups"=NAge,
#                       "N_Year"=NYear)
# 
# NimbleDataFr <- list("ZMat"=ZMatFr) 
# 
# #### Estimate Model 
# OwnQR_Fr <-  nimbleModel(code=OwnMod_Repara_Diff_OC, 
#                          constants = NimbleConstFr, 
#                          data=NimbleDataFr, buildDerivs = FALSE)
# 
# cOwnQR_Fr <- configureMCMC(OwnQR_Fr, 
#                            print = TRUE, useConjugacy = FALSE,
#                            monitors = c("k","sigma_eps",
#                                        "N_t","p",
#                                        "beta","betaJump",
#                                         "muY","sdY","sigma_time",
#                                         "drift",
#                                        #"b1","b2",
#                                         #"BMat","QMat","Y_t","J",
#                                         "a"),
#                            enableWAIC = TRUE) #to compare Models
# 
# ####################### 3.2 Adjust Samplers ####################################
# params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
#                        paste0("b2[",1:NAge,"]"),
#                        paste0("k[",1:NYear,"]"),
#                        paste0("Y_t[",1:(NYear+1),"]"),
#                        "sdY","muY")
# 
# cOwnQR_Fr$removeSamplers(params_to_remove)
# 
# cOwnQR_Fr$addSampler(target= paste0("b1[1:",NAge,"]"),
#                      type="AF_slice") 
# 
# cOwnQR_Fr$addSampler(target= paste0("b2[1:",NAge,"]"),
#                      type="AF_slice")
# 
# 
# cOwnQR_Fr$addSampler(target= c("sdY","muY"),
#                      type="HMC", 
#                      control = list(nwarmup=2000,  #increase delta
#                                     delta=0.82,
#                                     maxTreeDepth = 11)) #increase tree depth
# 
# 
# for(j in 1:NYear){
#   cOwnQR_Fr$addSampler(target = paste0("k[",j,"]"),
#                        type = 'slice')
# }
# 
# for(j in 3:(NYear+1)){
#   cOwnQR_Fr$addSampler(target = paste0("Y_t[",j,"]"),
#                        type = 'slice')
# }
# 
# bOwnQR_Fr <- buildMCMC(cOwnQR_Fr)
# comOwnQR_Fr <- compileNimble(OwnQR_Fr,
#                              bOwnQR_Fr)
# 
# 
# 
# params_to_remove <-  c(
#   "p","drift",
#   paste0("k[",1:NYear,"]"),
#   "a","muY","sdY"
# )
# 
# cOwnQR_Fr$removeSamplers(params_to_remove)
# cOwnQR_Fr$addDefaultSampler(nodes=c("drift","p"), 
#                                  useConjugacy = TRUE)
# 
# cOwnQR_Fr$addSampler(nodes = "muY",type="slice")
# cOwnQR_Fr$addSampler(nodes = "a",type="slice")
# cOwnQR_Fr$addSampler(nodes = "sdY",type="slice")
# 
# for(j in 1:NYear){
#   cOwnQR_Fr$addSampler(target = paste0("k[",j,"]"),
#                             type = 'slice')
# }
# 
# bOwnQR_Fr <- buildMCMC(cOwnQR_Fr)
# comOwnQR_Fr <- compileNimble(OwnQR_Fr,
#                              bOwnQR_Fr)
# 
# 
# SamplesOwnQR_Fr <- runMCMC(comOwnQR_Fr$bOwnQR_Fr, 
#                            niter  = 5000,
#                            thin=1,
#                            nburnin = 2000, 
#                            nchains = 1)
# 
# SummaryOutput(SamplesOwnQR_Fr, 
#               params=c("QMat","drift","p","sigma_eps",
#                        "muY","sdY","sigma_time","a","N_t")) %>% 
#   print(n=150) 

################### 4. GERMANY #################################################
pacman::p_load(openxlsx)

## Load Population Data from Eurostat 
PopGer_Eu <- read.xlsx(xlsxFile = file.path(getwd(),
                                            "Data/PopulationGermany_Eurostat.xlsx"),
                       colNames = TRUE, sheet="Total", 
                       startRow = 8, na.strings = ":")


## Load Death Data from Eurostat as well as provisional Deaths
DeathGer_EU <- read.xlsx(xlsxFile = file.path(getwd(),
                                              "Data/DeathsGermany_Eurostat.xlsx"),
                         colNames = TRUE, sheet="Total", 
                         startRow = 8, na.strings = ":")

DeathGer_22_weekly <- read.xlsx(xlsxFile = file.path(getwd(),
                                              "Data/ProvisionalCoronaDeaths_Germany.xlsx"),
                         colNames = TRUE, sheet="DeathsSummarised")

# Helper Function for Deaths. 
#Provisional Death data is available for age groups 0-15, 15-30, then in 5 years till 95+
AgeLabFunGer <- function(AgeWidth=10, AgeLabels){
  if(AgeWidth==10){
    HelperAgeLab <- data.frame("AgeOld"=AgeLabels[-1],
                               "AgeNew"= c(rep(1,15),rep(2,15),
                                           rep(3:9,each=10),9,9))
  } else {
    HelperAgeLab <- data.frame("AgeOld"=AgeLabels[-1],
                               "AgeNew"= c(rep(1,15),rep(2,15),
                                           rep(3:14,each=5), 
                                           rep(15,10),15,15))
  }
}
AgeW <- 10
HelperAgeLabGer <- AgeLabFunGer(AgeWidth = AgeW,
                                AgeLabels = unique(DeathGer_EU$`AGE.(Labels)`))


#1. Combine Death Data
DeathsGer_21_Grouped <- 
  DeathGer_EU %>% pivot_longer(cols = 3:ncol(.),
                            values_to = "Deaths", 
                            names_to="Year") %>% 
  filter(`AGE.(Labels)`!="Total") %>%  # remove Total
  mutate("AInd"=match(`AGE.(Labels)`, unique(`AGE.(Labels)`))) %>% 
  mutate("NewAgeInd"=HelperAgeLabGer$AgeNew[AInd]) %>% 
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise(across(Deaths,\(x) sum(x, na.rm = TRUE))) %>% #remove NA's
  mutate("Year"=as.numeric(Year)) #Death as numeric

DeathGer_22 <- DeathGer_22_weekly %>% 
  slice(-nrow(.)) %>% #remove last row
  select(Age.Group,Total) %>% 
  mutate("NewAgeInd"=if(AgeW==10){c(1,2,rep(3:9,each=2))}else{c(1:15,15)}) %>% 
  group_by(NewAgeInd) %>% 
  summarise("Deaths"=sum(Total)) %>% 
  mutate("Year"=2022)


#Add both in one
DeathsGer_Tot <- full_join(x= DeathsGer_21_Grouped, 
                          y= DeathGer_22, 
                          by=c("NewAgeInd","Year","Deaths"),
                          na_matches="never")

#2. Population Data
PopGer <- PopGer_Eu %>% pivot_longer(cols = 3:ncol(.),
                             values_to = "Pop", 
                             names_to="Year") %>% 
  filter(`AGE.(Labels)`!="Total") %>%  # remove Total
  mutate("AInd"=match(`AGE.(Labels)`, unique(`AGE.(Labels)`))) %>% 
  mutate("NewAgeInd"=HelperAgeLabGer$AgeNew[AInd]) %>% 
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise(across(Pop,\(x) sum(x, na.rm = TRUE))) %>% #remove NA's
  mutate("Year"=as.numeric(Year)) #Death as numeric

#3. Combine Both
TotDataGer <- dplyr::inner_join(x = DeathsGer_Tot,
                               y = PopGer, 
                               by =c("NewAgeInd","Year"),
                               na_matches="never")



####### 4.1 Estimate Model #####################################################
#Only Use Ages of 50 + 
LambdaVec <- TotDataGer %>% 
  mutate("Rate"=Deaths/Pop) %>% 
  arrange(Year)

LambdaVec %>% 
  ggplot(data=., aes(x=Year, y=log(Rate), group=NewAgeInd))+
  geom_line(aes(col=NewAgeInd))+
  ylab("Log Death Rate")+
  facet_wrap(~NewAgeInd,scales = "free")

LambdaMatGer <- matrix(LambdaVec$Rate,
                      nrow = length(unique(LambdaVec$NewAgeInd)), 
                      ncol = length(unique(LambdaVec$Year)),
                      byrow = FALSE)
## plot death rates
LambdaVec %>% 
  ggplot(data=., aes(x=Year, y=log(Rate), group=NewAgeInd))+
  geom_line(aes(col=NewAgeInd))


#Differenced Log Death Rates
ZMatGer <- apply(log(LambdaMatGer), 1, diff) %>% t()

NAge <- nrow(ZMatGer)
NYear <- ncol(ZMatGer)

NimbleConstGer <- list("N_AgeGroups"=NAge,
                      "N_Year"=NYear)

NimbleDataGer <- list("ZMat"=ZMatGer) 

#### Estimate Model 
OwnQR_Ger <-  nimbleModel(code=LCJumpOwn_QR, 
                         constants = NimbleConstGer, 
                         data=NimbleDataGer, buildDerivs = TRUE)

cOwnQR_Ger <- configureMCMC(OwnQR_Ger, 
                           print = TRUE, useConjugacy = FALSE,
                           monitors = c("k","sigma_eps",
                                        "QMat", "N_t","p",
                                        "J","muY","sdY","sigma_time",
                                        "drift","b1","b2","BMat","Y_t","a"),
                           enableWAIC = TRUE) #to compare Models

############### 4.2 Adjust Samplers ############################################
params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
                       paste0("b2[",1:NAge,"]"),
                       paste0("k[",1:NYear,"]"),
                       paste0("Y_t[",1:(NYear+1),"]"),
                       "sdY","muY")

cOwnQR_Ger$removeSamplers(params_to_remove)

cOwnQR_Ger$addSampler(target= paste0("b1[1:",NAge,"]"),
                     type="AF_slice") 

cOwnQR_Ger$addSampler(target= paste0("b2[1:",NAge,"]"),
                     type="AF_slice")


cOwnQR_Ger$addSampler(target= c("sdY","muY"),
                     type="HMC", 
                     control = list(nwarmup=2000,  #increase delta
                                    delta=0.82,
                                    maxTreeDepth = 11)) #increase tree depth


for(j in 1:NYear){
  cOwnQR_Ger$addSampler(target = paste0("k[",j,"]"),
                       type = 'slice')
}

for(j in 3:(NYear+1)){
  cOwnQR_Ger$addSampler(target = paste0("Y_t[",j,"]"),
                       type = 'slice')
}

bOwnQR_Ger <- buildMCMC(cOwnQR_Ger)
comOwnQR_Ger <- compileNimble(OwnQR_Ger,
                             bOwnQR_Ger)


SamplesOwnQR_Ger <- runMCMC(comOwnQR_Ger$bOwnQR_Ger, 
                           niter  = 5000,
                           thin=1,
                           nburnin = 2000, 
                           nchains = 1)

SummaryOutput(SamplesOwnQR_Ger, 
              params=c("QMat","drift","p","sigma_eps",
                       "muY","sdY","sigma_time","a","N_t")) %>% 
  print(n=150) 



# ################################################################################
# ############### DIFFERENT COUNTRIES ############################################
# ###### SWEDEN ##################################################################
# SwedenDeaths <- read.table(file = file.path(getwd(),"Data/Deaths_5x1_Sweden.txt"),
#                            header = TRUE)
# SwedenExp <- read.table(file = file.path(getwd(),"Data/Exposures_5x1_Sweden.txt"),
#                         header = TRUE)
# YearMin <-  1960
# 
# ###### 2. Put Data in correct Format ###########################################
# #Create Data Frame
# ## Problem: Too few people die in the old ages, differencing does not work
# #Cant be zero in sucessive years. Change highest age group to 95+
# # and do 0-5
# Helper <- data.frame("AgeOld"=1:24,
#                      "AgeNew"=c(1,1,2:19,rep(20,4)))
# 
# 
# # Problem 2: Population of Denmark is too small. Random Noise in Data is very big
# # - try and group ages together even further
# Helper <- data.frame("AgeOld"=1:24,
#                      "AgeNew"=c(1,2,rep(seq(3,10,by=1),each=2),rep(11,6)))
# 
# TotalData <- data.frame("Y"=SwedenDeaths$Total,
#                         "Offset"=SwedenExp$Total,
#                         "Year"=SwedenDeaths$Year,
#                         "AgeGroup"=SwedenDeaths$Age) %>% 
#   filter(Year > YearMin) %>% #only 1901 onward
#   mutate("AInd"=match(AgeGroup, unique(AgeGroup))) %>% #Get Age ID
#   mutate("NewAgeInd"=Helper$AgeNew[AInd]) %>% #add new Age ID
#   group_by(NewAgeInd,Year) %>% #group by new Age ID
#   summarise(across(c(Y,Offset),sum)) %>% #sum over Age ID
#   mutate("TInd"=match(Year, unique(Year))) %>% 
#   mutate("Y"=round(Y,0), #why are there decimal points in data?
#          "Offset"=round(Offset,0)) %>% 
#   arrange(Year)
# 
# 
# 
# 
# #Create Matrix of Differenced log Death Rates
# LambdaMatSw <- matrix(TotalData$Y/TotalData$Offset,
#                       nrow = max(TotalData$NewAgeInd), 
#                       byrow = FALSE)
# 
# ts.plot(LambdaMatSw[4,])
# 
# #Differenced Log Death Rates
# ZMatSw <- apply(log(LambdaMatSw), 1, diff) %>% t()
# 
# NAge <- nrow(ZMatSw)
# NYear <- ncol(ZMatSw)
# 
# NimbleConstSw <- list("N_AgeGroups"=NAge,
#                       "N_Year"=NYear)
# 
# NimbleDataSw <- list("ZMat"=ZMatSw) 
# 
# 
# 
# ### 3. Run Liu, Li Model ######################################################
# QRMod_Sw <- nimbleModel(code=LCJumpOwn_QR,
#                         constants = NimbleConstSw, 
#                         data=NimbleDataSw)
# 
# QRModc <- configureMCMC(QRMod_Sw,monitors = c("k","sigma_eps",
#                                               "QMat", "N_t","p",
#                                               "J","muY","sdY","sigma_time",
#                                               "drift","b1","b2","BMat","Y_t","a"),
#                         enableWAIC = TRUE, useConjugacy = FALSE)
# 
# params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
#                        paste0("b2[",1:NAge,"]"),
#                        paste0("k[",1:NYear,"]"))
# 
# QRModc$removeSamplers(params_to_remove)
# 
# QRModc$addSampler(target= paste0("b1[1:",NAge,"]"),
#                   type="AF_slice")
# 
# QRModc$addSampler(target= paste0("b2[1:",NAge,"]"),
#                   type="AF_slice")
# 
# for(j in 1:NYear){
#   QRModc$addSampler(target = paste0("k[",j,"]"),
#                     type = 'slice'
#   )
# }
# 
# 
# ## test uncompiled function in uncompiled MCMC
# buildQR <- buildMCMC(QRModc)
# compQR <- compileNimble(QRMod_Sw,
#                         buildQR)
# 
# SamplesLiuLi_Sw <- runMCMC(compQR$buildQR, 
#                            niter = 5000,
#                            nburnin = 2000,
#                            thin=2,
#                            nchains = 1) #set seed for reproducibility
# 
# SummaryOutput(SamplesLiuLi_Sw, 
#               params=c("QMat","drift","p","sigma_eps",
#                        "muY","sdY","sigma_time","a")) %>% 
#   print(n=150) 
# 
# ############# JAPAN ##########################################################
# JapanDeaths <- read.table(file = file.path(getwd(),"Data/Deaths_5x1_Japan.txt"),
#                           header = TRUE)
# JapanExp <- read.table(file = file.path(getwd(),"Data/Exposures_5x1_Japan.txt"),
#                        header = TRUE)
# YearMin <- 1950
# ###### 2. Put Data in correct Format ###########################################
# #Create Data Frame
# ## Problem: Too few people die in the old ages, differencing does not work
# #Cant be zero in sucessive years. Change highest age group to 95+
# # and do 0-5
# Helper <- data.frame("AgeOld"=1:24,
#                      "AgeNew"=c(1,1,2:19,rep(20,4)))
# 
# 
# # Problem 2: Population of Denmark is too small. Random Noise in Data is very big
# # - try and group ages together even further
# Helper <- data.frame("AgeOld"=1:24,
#                      "AgeNew"=c(1,2,rep(seq(3,10,by=1),each=2),rep(11,6)))
# 
# TotalDataJap <- data.frame("Y"=JapanDeaths$Total,
#                            "Offset"=JapanExp$Total,
#                            "Year"=JapanDeaths$Year,
#                            "AgeGroup"=JapanDeaths$Age) %>% 
#   filter(Year > YearMin) %>% #only 1901 onward
#   mutate("AInd"=match(AgeGroup, unique(AgeGroup))) %>% #Get Age ID
#   mutate("NewAgeInd"=Helper$AgeNew[AInd]) %>% #add new Age ID
#   group_by(NewAgeInd,Year) %>% #group by new Age ID
#   summarise(across(c(Y,Offset),sum)) %>% #sum over Age ID
#   mutate("TInd"=match(Year, unique(Year))) %>% 
#   mutate("Y"=round(Y,0), #why are there decimal points in data?
#          "Offset"=round(Offset,0)) %>% 
#   arrange(Year)
# 
# LambdaMatJap <- matrix(TotalDataJap$Y/TotalDataJap$Offset,
#                        nrow = max(TotalDataJap$NewAgeInd), 
#                        byrow = FALSE)
# 
# 
# #Differenced Log Death Rates
# ZMatJap <- apply(log(LambdaMatJap), 1, diff) %>% t()
# 
# NAge <- nrow(ZMatJap)
# NYear <- ncol(ZMatJap)
# 
# NimbleConstJap <- list("N_AgeGroups"=NAge,
#                        "N_Year"=NYear)
# 
# NimbleDataJap <- list("ZMat"=ZMatJap) 
# 
# ## Estimate Model ###########
# OwnQR_Jap <-  nimbleModel(code=LCJumpOwn_QR, 
#                           constants = NimbleConstJap, 
#                           data=NimbleDataJap)
# 
# cOwnQR_Jap <- configureMCMC(OwnQR_Jap, 
#                             print = TRUE, useConjugacy = FALSE,
#                             monitors = c("k","sigma_eps",
#                                          "QMat", "N_t","p",
#                                          "J","muY","sdY","sigma_time",
#                                          "drift","b1","b2","BMat","Y_t","a"),
#                             enableWAIC = TRUE) #to compare Models
# 
# #### Adjust Samplers ########
# params_to_remove <-  c(paste0("b1[",1:NAge,"]"),
#                        paste0("b2[",1:NAge,"]"),
#                        paste0("k[",1:NYear,"]"),
#                        paste0("Y_t[",1:(NYear+1),"]"),
#                        "sdY","muY")
# 
# cOwnQR_Jap$removeSamplers(params_to_remove)
# 
# cOwnQR_Jap$addSampler(target= paste0("b1[1:",NAge,"]"),
#                       type="AF_slice") 
# 
# cOwnQR_Jap$addSampler(target= paste0("b2[1:",NAge,"]"),
#                       type="AF_slice")
# 
# 
# cOwnQR_Jap$addSampler(target = "sdY", type="slice")
# cOwnQR_Jap$addSampler(target = "muY", type="slice")
# 
# 
# for(j in 1:NYear){
#   cOwnQR_Jap$addSampler(target = paste0("k[",j,"]"),
#                         type = 'slice')
# }
# 
# for(j in 1:(NYear+1)){
#   cOwnQR_Jap$addSampler(target = paste0("Y_t[",j,"]"),
#                         type = 'slice')
# }
# 
# bOwnQR_Jap <- buildMCMC(cOwnQR_Jap)
# comOwnQR_Jap <- compileNimble(OwnQR_Jap,
#                               bOwnQR_Jap)
# 
# 
# SamplesOwnQR_Jap <- runMCMC(comOwnQR_Jap$bOwnQR_Jap, 
#                             niter = 5000,
#                             thin=1,
#                             nburnin = 4000, 
#                             nchains = 2)
# 
# SummaryOutput(SamplesOwnQR_Jap, 
#               params=c("QMat","drift","p","sigma_eps",
#                        "muY","sdY","sigma_time","a","N_t")) %>% 
#   print(n=150) 
# 
# ts.plot(SamplesOwnQR_Jap$chain2[,"p"])
# 
# 
# 
# 
