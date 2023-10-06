### Simulation Study###
library(nimble);library(tidyverse); library(nimbleHMC)

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 

############ SIMULATION COMPARISON ###########################################
## Compare Own Approach with that or reparameterization
library(nimbleHMC)
## Compare Liu,Li Model using QR and sum-to-one constraint for
# both the simple and the complex model (estimation of jump parameter)


#1. First Simulate some Data,
set.seed(42)
NAge <- 10
# beta <- c(0.181,0.172,0.149,0.078,0.06,0.07,0.07,0.08,0.08,0.06)
# betaJ <- c(0.02,0.07,0.07,0.0337,0.008,0.17,0.1083,0.28,0.21,0.03)

#US estimates
beta <- c(0.117,0.172,0.145,0.108,0.1,0.0996,0.0693,0.069,0.0651,0.055)
betaJ <- c(0.02,0.03,0.04,0.127,0.158,0.155,0.132,0.117,0.111,0.11)


QMat <- qr.Q(qr(cbind(beta,betaJ))) #QR decomposition of beta
QMat <- QMat*(-1)

#From England Model
QMat2 <- SummaryOutput(SamplesOwnQR_EW_S, 
                       params=c("QMat")) %>% 
  select(2) %>% pull() %>% 
  matrix(data=., ncol=2, byrow = FALSE)
#Make sure that QMat has norm of 1, mean estimates have a norm of slightly less,
#also add 1 more Age group with small effect
QMat3 <- rbind(QMat2,
               c(0.1,0.1)) %>% 
  apply(., 2, function(x){x/norm(x, type="2")})

# 1. Type of constraints 
#(drift = -0.4, sigma_time = 0.05)
NYear <-300
# vt <- rnorm(n = NYear,mean = 0,sd = 0.05)
# drift <- -0.4
drift <- -0.2
vt <- c(rnorm(n = 1, mean = 0,sd = 1),0,rnorm(n = NYear-2,mean = 0,sd = 0.1))

RWDrift <- cumsum(drift+vt)
#RWDrift <- RWDrift-mean(RWDrift) #sum to zero constraint

p <- 0.15
N_t <- rbinom(n= NYear,size = 1,prob = p)
Y_t <- rnorm(n = NYear,mean = 1.5,sd = 0.4)
# N_t <- c(rep(0,10),1,1,rep(0,20),1,1,1,rep(0,40),1,1,rep(0,23))
# Y_t <- rnorm(n = NYear,mean = 2,sd = 0.5)

# p <- 0.2
# N_t <- rbinom(n= NYear,size = 1,prob = p)
# Y_t <- rnorm(n = NYear,mean = 1,sd = 0.3)
#p = 0.2, muY = 1, sdY = 0.3


#Step 2, let J be a RV
J_t <- c(0,0,N_t[-c(1:2)]*Y_t[-c(1:2)])


ax <- seq(-5,0.5,0.5)

## Simulation of log death rates 
LambdaMatLiu <- matrix(0, nrow = NAge, ncol = NYear)
for(x in 1:NAge){
  for (t in 1:NYear) {
    LambdaMatLiu[x,t] <- ax[x] + 
      beta[x]*RWDrift[t] + betaJ[x]*J_t[t]+
      rnorm(n=1, sd=0.03)
  }
}

LambdaMatLiuQR <- matrix(0, nrow = NAge, ncol = NYear)
for(x in 1:NAge){
  for (t in 1:NYear) {
    LambdaMatLiuQR[x,t] <- #ax[x] + 
      QMat[x,1]*RWDrift[t] + QMat[x,2]*J_t[t]+
      rnorm(n=1, sd=0.3)
  }
}



LambdaMatLiu %>% as.data.frame() %>% 
  mutate("Age"=as.factor(1:10), .before=1) %>% 
  pivot_longer(., cols=2:101, names_to = "Year", 
               values_to = "LogRate") %>% 
  mutate("Time"=rep(1:100,10)) %>% 
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()

LambdaMatLiuQR %>% as.data.frame() %>% 
  mutate("Age"=as.factor(1:10), .before=1) %>% 
  pivot_longer(., cols=2:101, names_to = "Year", 
               values_to = "LogRate") %>% 
  mutate("Time"=rep(1:100,10)) %>% 
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()


ZMatSimLiu <- apply(LambdaMatLiu, 1, diff) %>% t() 
ZMatSimQR <- apply(LambdaMatLiuQR, 1, diff) %>% t() 

ZMatSimLiu %>% as.data.frame() %>% 
  mutate("Age"=as.factor(1:10), .before=1) %>% 
  pivot_longer(., cols=2:100, names_to = "Year", 
               values_to = "LogRate") %>% 
  mutate("Time"=rep(1:99,10)) %>% 
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()


ZMatSimQR %>% as.data.frame() %>% 
  mutate("Age"=as.factor(1:10), .before=1) %>% 
  pivot_longer(., cols=2:100, names_to = "Year", 
               values_to = "LogRate") %>% 
  mutate("Time"=rep(1:99,10)) %>% 
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()


#sigma_epsZ 
sqrt(2*0.03^2) #sigma_epsZ (Standard deviation of 2*sigma_eps^2)
mean(N_t) # p
mean(J_t[which(J_t!=0)]) #muY
sd(J_t[which(J_t!=0)]) #sdY
mean(Y_t[which(N_t!=0)]) #muY
sd(Y_t[which(N_t!=0)]) #muY

sd(diff(RWDrift))
mean(diff(RWDrift))



#################### LIU, LI ModelCompare estimates of own model with that of simpler model 
SimLiuConst <- list("N_Year"=ncol(ZMatSimLiu),
                    "N_AgeGroups"=nrow(ZMatSimLiu))

SimLiuDat <- list("ZMat"=ZMatSimLiu) 

LiuLiSimMod <- nimbleModel(code=LiuLi_Repara_Diff_OC_Test, 
                           constants = SimLiuConst, 
                           data=SimLiuDat, buildDerivs = FALSE)

cLiuLiSimMod <- configureMCMC(LiuLiSimMod, 
                              print = TRUE, useConjugacy = TRUE,
                              monitors = c("N_t","muY",
                                           "sdY","p","sigma_eps",
                                           "beta","betaJump",
                                           "sigma_time",
                                           "drift",
                                           "k","VarY"
                              ),
                              enableWAIC = TRUE)

cLiuLiSimMod$removeSamplers(c("VarY", 
                              "muY",
                               #paste0("b1[",1:10,"]"),
                              # paste0("b2[",1:10,"]"),
                              #paste0("Y_std[",3:100,"]"),
                              "sigma_time",
                              "sigma_eps"
                            ))


cLiuLiSimMod$addSampler(target = "muY",type="slice")
cLiuLiSimMod$addSampler(target = "VarY",type="slice")

cLiuLiSimMod$addSampler(target = c("sigma_time",
                                   "sigma_eps"),type="AF_slice")


cLiuLiSimMod$addSampler(target = "sigma_time", type="RW", log=TRUE)
cLiuLiSimMod$addSampler(target = "sigma_eps", type="slice")

cLiuLiSimMod$addSampler(target= c("b1[1:10]"),
                        type="AF_slice")

cLiuLiSimMod$addSampler(target= paste0("b2[1:",NAge,"]"),
                        type="AF_slice")



bLiuLiSimMod <- buildMCMC(cLiuLiSimMod)
comLiuLiSimMod <- compileNimble(LiuLiSimMod,
                                bLiuLiSimMod)

SamplesLiuLiSimMod <- runMCMC(comLiuLiSimMod$bLiuLiSimMod, 
                              niter = 4000,
                              thin=1,
                              nburnin = 2000, 
                              nchains = 1)

SummaryOutput(SamplesLiuLiSimMod, 
              params=c("sigma_eps",
                       "p","muY","sdY","sigma_time","beta","betaJump",
                       "drift","VarY")) %>% print(n=250) 

par(mfrow=c(2,1))
plot(SamplesLiuLiSimMod[,"muY"],type="l")
plot(SamplesLiuLiSimMod[,"sigma_time"],col="red",type="l")

cor(SamplesLiuLiSimMod[,"muY"],
    SamplesLiuLiSimMod[,"sdY"])

comLiuLiSimMod$bLiuLiSimMod$getWAIC()
SamplesLiuLiSimMod %>%
  bayesplot::mcmc_pairs(pars = c("sdY", "sigma_time","sigma_eps"))


save(SamplesLiuLiSimMod,
     SamplesLiuLiSimMod_Rep,
     file = file.path(getwd(),"Results/SimulationStandardVsRepara.RData"))

########### Test of Parameterization ######################################
SimLiuConst <- list("N_Year"=ncol(ZMatSimLiu),
                    "N_AgeGroups"=nrow(ZMatSimLiu))

SimLiuDat <- list("ZMat"=ZMatSimLiu) 

LiuLiSimMod <- nimbleModel(code=LiuLi_Repara_Diff_OC_V2, 
                           constants = SimLiuConst, 
                           data=SimLiuDat, buildDerivs = FALSE)

cLiuLiSimMod <- configureMCMC(LiuLiSimMod, 
                              print = TRUE, useConjugacy = TRUE,
                              monitors = c("N_t","muY",
                                           "sdY","p","sigma_eps",
                                           "beta","betaJump",
                                           "sigma_time",
                                           "drift",
                                           #"k",
                                           "VarY"
                              ),
                              enableWAIC = TRUE)

cLiuLiSimMod$removeSamplers(c("VarY", 
                              "muY",
                              #paste0("b1[",1:10,"]"),
                              # paste0("b2[",1:10,"]"),
                              #paste0("Y_std[",3:100,"]"),
                              "sigma_time", "sigma_eps"
))


cLiuLiSimMod$addSampler(target = "muY",type="slice")
cLiuLiSimMod$addSampler(target = "VarY",type="slice")

cLiuLiSimMod$addSampler(target = c("sigma_time",
                                   "sigma_eps"),type="AF_slice")


LiuLiSimMod2 <- nimbleModel(code=LiuLi_Repara_Diff_OC_V2_Test, 
                           constants = SimLiuConst, 
                           data=SimLiuDat, buildDerivs = FALSE)

cLiuLiSimMod2 <- configureMCMC(LiuLiSimMod2, 
                              print = TRUE, useConjugacy = TRUE,
                              monitors = c("N_t","muY",
                                           "sdY","p","sigma_eps",
                                           "beta","betaJump",
                                           "sigma_time",
                                           "drift",
                                           #"k",
                                           "sdY"
                              ),
                              enableWAIC = TRUE)

cLiuLiSimMod2$removeSamplers(c("sdY", 
                              "muY",
                              #paste0("b1[",1:10,"]"),
                              # paste0("b2[",1:10,"]"),
                              #paste0("Y_std[",3:100,"]"),
                              "sigma_time","sigma_eps"
))


cLiuLiSimMod2$addSampler(target = "muY",type="slice")
cLiuLiSimMod2$addSampler(target = "sdY",type="slice")

cLiuLiSimMod2$addSampler(target = c("sigma_time",
                                   "sigma_eps"),type="AF_slice")

bLiuLiSimMod2 <- buildMCMC(cLiuLiSimMod2)
comLiuLiSimMod2 <- compileNimble(LiuLiSimMod2,
                                bLiuLiSimMod2)

SamplesLiuLiSimMod2 <- runMCMC(comLiuLiSimMod2$bLiuLiSimMod2, 
                              niter = 4000,
                              thin=1,
                              nburnin = 2000, 
                              nchains = 1)

bLiuLiSimMod <- buildMCMC(cLiuLiSimMod)
comLiuLiSimMod <- compileNimble(LiuLiSimMod,
                                bLiuLiSimMod)

SamplesLiuLiSimMod <- runMCMC(comLiuLiSimMod$bLiuLiSimMod, 
                              niter = 4000,
                              thin=1,
                              nburnin = 2000, 
                              nchains = 1)

SummaryOutput(SamplesLiuLiSimMod, 
              params=c("sigma_eps",
                       "p","muY","sdY","sigma_time","beta","betaJump",
                       "drift","VarY")) %>% print(n=250) 

SummaryOutput(SamplesLiuLiSimMod2, 
              params=c("sigma_eps",
                       "p","muY","sdY","sigma_time","beta","betaJump",
                       "drift","VarY")) %>% print(n=250) 

comLiuLiSimMod$bLiuLiSimMod$getWAIC()
comLiuLiSimMod2$bLiuLiSimMod2$getWAIC()

SummaryOutput(SamplesLiuLiSimMod, params = c("k")) %>% print(n=150)
SummaryOutput(SamplesLiuLiSimMod2, params = c("k")) %>% print(n=150)
######## LIU, LI Repara ######################################################
SimLiuConst <- list("N_Year"=ncol(ZMatSimLiu),
                    "N_AgeGroups"=nrow(ZMatSimLiu))

SimLiuDat <- list("ZMat"=ZMatSimLiu) 

LiuLiSimMod_Rep <- nimbleModel(code=LiuLi_Repara_Diff_OC_V2, 
                           constants = SimLiuConst, 
                           data=SimLiuDat, buildDerivs = FALSE)

cLiuLiSimMod_Rep <- configureMCMC(LiuLiSimMod_Rep, 
                              print = TRUE, useConjugacy = TRUE,
                              monitors = c("N_t","muY",
                                           "sdY","p","sigma_eps",
                                           "beta","betaJump",
                                           "sigma_time","drift",
                                           "VarY"
                                           #"k",
                                           #"mu","sigma_squared"
                                           #"Y_t","J"
                              ),
                              enableWAIC = TRUE)

cLiuLiSimMod_Rep$removeSamplers(c(#"sdY", 
                                  "muY",
                                  "VarY",
                              paste0("b1[",1:10,"]"),
                              paste0("b2[",1:10,"]"),
                              "sigma_time"))

#cLiuLiSimMod_Rep$addSampler(target = "sdY",type="slice")
cLiuLiSimMod_Rep$addSampler(target = "VarY",type="slice")
cLiuLiSimMod_Rep$addSampler(target = "muY",type="slice")
cLiuLiSimMod_Rep$addSampler(target = "sigma_time",type="slice")

cLiuLiSimMod_Rep$addSampler(target= c("b1[1:10]"),
                        type="AF_slice")

cLiuLiSimMod_Rep$addSampler(target= c("b2[1:10]"),
                        type="AF_slice")

#moderate Correlation between muY and sdY
# cLiuLiSimMod$addSampler(target = c("muY","sdY"),
#                         type="AF_slice")


bLiuLiSimMod_Rep <- buildMCMC(cLiuLiSimMod_Rep)
comLiuLiSimMod_Rep <- compileNimble(LiuLiSimMod_Rep,
                                bLiuLiSimMod_Rep)

SamplesLiuLiSimMod_Rep <- runMCMC(comLiuLiSimMod_Rep$bLiuLiSimMod_Rep, 
                              niter = 4000,
                              thin=1,
                              nburnin = 2000, 
                              nchains = 1)

SummaryOutput(SamplesLiuLiSimMod_Rep,
              params=c("sigma_eps","p","muY","sdY",
                       "sigma_time","beta","betaJump","drift")) %>% 
  print(n=250) 




## Posterior Predictive Checks#
Z_RepSim <- ReplicatedDataFun(Samples = SamplesLiuLiSimMod_Rep,
                             n = length(ZMatSimLiu), SingleVar = FALSE)

Z_RepSimOC <- ReplicatedDataFun(Samples = SamplesLiuLiSimMod,
                                n = length(ZMatSimLiu), SingleVar = TRUE)


LambdaSimLiu <- ReplicatedLambdaFun(Samples = SamplesLiuLiSimMod_Rep,
                                    LambdaMat = exp(LambdaMatLiu), 
                                    Z_Rep = Z_RepSim)

LambdaSimLiuOC <- ReplicatedLambdaFun(Samples = SamplesLiuLiSimMod,
                                    LambdaMat = exp(LambdaMatLiu), 
                                    Z_Rep = Z_RepSimOC)


AgeTimeGrid <- expand.grid("NewAgeInd"=1:nrow(LambdaMatLiu), "Year"=1:ncol(LambdaMatLiu))

LambdaVecLiu <- data.frame("Rate"=exp(as.numeric(LambdaMatLiu))) %>% 
  mutate(AgeTimeGrid)

library(patchwork)

p1 <- LambdaPlot(LambdaArray = LambdaSimLiu, 
           LambdaMat = exp(LambdaMatLiu), 
           LambdaVec = LambdaVecLiu, 
           startYear = 1, AgesUsed = c(4))

p2 <- LambdaPlot(LambdaArray = LambdaSimLiuOC, 
           LambdaMat = exp(LambdaMatLiu), 
           LambdaVec = LambdaVecLiu, 
           startYear = 1, AgesUsed = c(4))

#Not seen, but width of jump is slightly wider in Reparameterized model, than in one with parameter

LambdaSimLiu %>% apply(., c(1,2),mean )
LambdaSimLiuOC %>% apply(., c(1,2),mean )

comLiuLiSimMod$bLiuLiSimMod$getWAIC() #higher WAIC, due to higher likelihood?
comLiuLiSimMod_Rep$bLiuLiSimMod_Rep$getWAIC()

SummaryOutput(SamplesLiuLiSimMod, params = "J") %>% 
  filter(mean >0.1) %>% 
  summarise(ME = mean(mean),
            SE = sd(mean))

###################### QR Model ################################################
SimLiuConstQR <- list("N_Year"=ncol(ZMatSimQR),
                      "N_AgeGroups"=nrow(ZMatSimQR))

SimLiuDatQR <- list("ZMat"=ZMatSimQR) 

LiuLiSimModQR <- nimbleModel(code=LiuLi_Repara_QR, 
                             constants = SimLiuConstQR, 
                             data=SimLiuDatQR, buildDerivs = FALSE)

cLiuLiSimModQR <- configureMCMC(LiuLiSimModQR, 
                                print = TRUE, useConjugacy = TRUE,
                                monitors = c("N_t",
                                             #"Y_t",
                                             "muY","sdY",
                                             "p","sigma_eps",
                                             "QMat","k",
                                             "sigma_time","drift", "b1","b2"
                                ),
                                enableWAIC = TRUE)

cLiuLiSimModQR$removeSamplers(c(paste0("b1[",1:10,"]"),
                                paste0("b2[",1:10,"]")),
                              "sigma_time","sdY","muY")

cLiuLiSimModQR$addSampler(target = "sdY",type="RW", log=TRUE)
cLiuLiSimModQR$addSampler(target = "sigma_time",type="RW", log=TRUE)
cLiuLiSimModQR$addSampler(target = "muY",type="slice")


cLiuLiSimModQR$addSampler(target= c("b1[1:10]"),
                          type="AF_slice")

cLiuLiSimModQR$addSampler(target= c("b2[1:10]"),
                          type="AF_slice")

bLiuLiSimModQR <- buildMCMC(cLiuLiSimModQR)
comLiuLiSimModQR <- compileNimble(LiuLiSimModQR,
                                  bLiuLiSimModQR)


SamplesLiuLiSimModQR <- runMCMC(comLiuLiSimModQR$bLiuLiSimModQR, 
                                niter = 4000,
                                thin=1,
                                nburnin = 2000, 
                                nchains = 1)#, 

SummaryOutput(SamplesLiuLiSimModQR, 
              params=c("sigma_eps","muY","sdY",
                       "p","sigma_time","QMat","muY","sdY",
                       "drift")) %>% print(n=250) 

comLiuLiSimModQR$bLiuLiSimModQR$getWAIC()
ts.plot(SamplesLiuLiSimModQR[,"QMat[2, 1]"])


########### OWN MODEL SIMULATION ###############################################
J_t <- numeric(NYear)
a <- 0.2
J_t[1:2] <- 0
for(j in 3:(NYear)){
  J_t[j] <- a*J_t[j-1]+N_t[j]*Y_t[j]
}

EpsMat <- sapply(rep(10,NYear), function(x) rnorm(n=x, mean=0, sd=0.03))

## Simulation of log death rates 
LambdaMatOwn <- matrix(0, nrow = NAge, ncol = NYear)
for(x in 1:NAge){
  for (t in 1:NYear) {
    LambdaMatOwn[x,t] <- #ax[x] + 
      beta[x]*RWDrift[t] + betaJ[x]*J_t[t]+
      EpsMat[x,t]
  }
}

LambdaMatOwnQR <- matrix(0, nrow = NAge, ncol = NYear)
for(x in 1:NAge){
  for (t in 1:NYear) {
    LambdaMatOwnQR[x,t] <- #ax[x] + 
      QMat[x,1]*RWDrift[t] + QMat[x,2]*J_t[t]+
      #QMat3[x,1]*RWDrift[t] + QMat3[x,2]*J_t[t]+
      EpsMat[x,t]
  }
}

ZMatOwnTest <- matrix(0, nrow = NAge, ncol = NYear-1)
for (t in 1:(NYear-1)) {
  for(x in 1:NAge){
    ZMatOwnTest[x,t] <- beta[x]*(RWDrift[t+1]-RWDrift[t])+
      betaJ[x]*(J_t[t+1]-J_t[t])+
      EpsMat[x,t+1]-EpsMat[x,t]
  }
}

#is the same as lambda mat and then differentiating
ZMatOwnQRTest <- matrix(0, nrow = NAge, ncol = NYear)
for (t in 1:(NYear-1)) {
  for(x in 1:NAge){
    ZMatOwnQRTest[x,t] <- QMat3[x,1]*(RWDrift[t+1]-RWDrift[t])+
      QMat3[x,2]*(J_t[t+1]-J_t[t])+
      EpsMat[x,t+1]-EpsMat[x,t]
  }
}



LambdaMatOwn %>% as.data.frame() %>% 
  mutate("Age"=as.factor(1:10), .before=1) %>% 
  pivot_longer(., cols=2:101, names_to = "Year", 
               values_to = "LogRate") %>% 
  mutate("Time"=rep(1:100,10)) %>% 
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()



LambdaMatOwnQR %>% as.data.frame() %>% 
  mutate("Age"=as.factor(1:10), .before=1) %>% 
  pivot_longer(., cols=2:101, names_to = "Year", 
               values_to = "LogRate") %>% 
  mutate("Time"=rep(1:100,10)) %>% 
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()

## Own Model

ZMatSimOwn <- apply(LambdaMatOwn, 1, diff) %>% t() 
ZMatSimOwnQR <- apply(LambdaMatOwnQR, 1, diff) %>% t() 

all(round(ZMatSimOwn,4) == round(ZMatOwnTest,4))

ZMatSimOwn %>% as.data.frame() %>% 
  mutate("Age"=as.factor(1:10), .before=1) %>% 
  pivot_longer(., cols=2:100, names_to = "Year", 
               values_to = "LogRate") %>% 
  mutate("Time"=rep(1:99,10)) %>% 
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()


ZMatSimOwnQR %>% as.data.frame() %>%
  mutate("Age"=as.factor(1:10), .before=1) %>%
  pivot_longer(., cols=2:100, names_to = "Year",
               values_to = "LogRate") %>%
  mutate("Time"=rep(1:99,10)) %>%
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()

sqrt(2*0.03^2) #sigma_epsZ (Standard deviation of 2*sigma_eps^2)
mean(N_t) # p
mean(J_t[which(J_t!=0)]) #muY
sd(J_t[which(J_t!=0)]) #sdY
mean(Y_t[which(N_t!=0)]) #muY
sd(Y_t[which(N_t!=0)]) #muY

## Compare estimates of own model with that of simpler model 
SimOwnConst <- list("N_Year"=ncol(ZMatSimOwn),
                    "N_AgeGroups"=nrow(ZMatSimOwn))

SimOwnDat <- list("ZMat"=ZMatSimOwn) 

OwnSimMod_Rep <- nimbleModel(code=OwnMod_Repara_Diff_OC_V2, 
                         constants = SimOwnConst, 
                         data=SimOwnDat, buildDerivs = FALSE)

cOwnSimMod_Rep <- configureMCMC(OwnSimMod_Rep, 
                            print = TRUE, useConjugacy = TRUE,
                            monitors = c("N_t","p","sigma_eps",
                                         "beta","betaJump",
                                         #"k",
                                         "sigma_time",
                                         "drift", "beta","betaJump",
                                         "a","CovVec","muY","sdY"
                            ),
                            enableWAIC = TRUE)

cOwnSimMod_Rep$removeSamplers(c( 
  paste0("b1[",1:10,"]"),
  paste0("b2[",1:10,"]"),
  #paste0("k[",1:(NYear-1),"]"),
  "sigma_time","sigma_eps",
  "sdY",
  "a","muY"))

# cOwnSimMod$removeSamplers("sdY",
#                           "sigma_time")


cOwnSimMod_Rep$addSampler(target = "sdY",type="slice")
cOwnSimMod_Rep$addSampler(target = c("sigma_time","sigma_eps"),
                          type="AF_slice")
cOwnSimMod_Rep$addSampler(target = "a",type="slice")
cOwnSimMod_Rep$addSampler(target = "muY",type="slice")


cOwnSimMod_Rep$addSampler(target= c("b1[1:10]"),
                            type="AF_slice")

cOwnSimMod_Rep$addSampler(target= c("b2[1:10]"),
                            type="AF_slice")

bOwnSimMod_Rep <- buildMCMC(cOwnSimMod_Rep)
comOwnSimMod_Rep <- compileNimble(OwnSimMod_Rep,
                              bOwnSimMod_Rep)

SamplesOwnMod_Rep <- runMCMC(comOwnSimMod_Rep$bOwnSimMod_Rep, 
                         niter = 3000,
                         thin=1,
                         nburnin = 1500, 
                         nchains = 1)

SummaryOutput(SamplesOwnMod_Rep, params = c("beta","betaJump",
                                            "muY","sdY","drift","sigma_time",
                                            "a","sigma_eps")) %>% 
  print(n=500)

comOwnSimMod$bOwnSimMod$getWAIC()
SamplesOwnMod_Rep %>%
  bayesplot::mcmc_pairs(pars = c("sdY", "sigma_time","sigma_eps"))

###### Standard PARAMETERIZATION ###############################################
SimOwnConst <- list("N_Year"=ncol(ZMatSimOwn),
                    "N_AgeGroups"=nrow(ZMatSimOwn))

SimOwnDat <- list("ZMat"=ZMatSimOwn) 

OwnSimMod <- nimbleModel(code=LCJumpOwn_OC_CornerK, 
                         constants = SimOwnConst, 
                         data=SimOwnDat, buildDerivs = FALSE)

cOwnSimMod <- configureMCMC(OwnSimMod, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("N_t","p","sigma_eps",
                                         "beta","betaJump",
                                         "k",
                                         "sigma_time","drift", "beta","betaJump",
                                         "a","muY","sdY",
                                         "J","Y_t"
                            ),
                            enableWAIC = TRUE)

cOwnSimMod$removeSamplers(c( 
  paste0("b1[",1:10,"]"),
  paste0("b2[",1:10,"]"),
  paste0("k[",1:(NYear-1),"]"),
  paste0("Y_t[",1:(NYear),"]"),
  "sigma_time","sdY",
  "a","muY"))



cOwnSimMod$addSampler(target = "sdY",type="slice")
cOwnSimMod$addSampler(target = "sigma_time",type="slice")
cOwnSimMod$addSampler(target = "a",type="slice")
cOwnSimMod$addSampler(target = "muY",type="slice")
cOwnSimMod$addDefaultSampler(nodes = "k[2:299]",useConjugacy = TRUE)

# for(j in 2:(NYear-1)){
#   cOwnSimMod$addSampler(target = paste0("k[",j,"]"),
#                         type = 'slice')
# }

for(j in 3:(NYear)){
  cOwnSimMod$addSampler(target = paste0("Y_t[",j,"]"),
                        type = 'slice')
}


cOwnSimMod$addSampler(target= c("b1[1:10]"),
                      type="AF_slice")

cOwnSimMod$addSampler(target= c("b2[1:10]"),
                      type="AF_slice")


bOwnSimMod <- buildMCMC(cOwnSimMod)
comOwnSimMod <- compileNimble(OwnSimMod,
                              bOwnSimMod)

SamplesOwnMod <- runMCMC(comOwnSimMod$bOwnSimMod, 
                         niter = 4000,
                         thin=1,
                         nburnin = 2000, 
                         nchains = 1)
SummaryOutput(SamplesOwnMod, 
              params=c("sigma_eps",
                       "p","muY","sdY","sigma_time","beta","betaJump",
                       "drift","a")) %>% print(n=250) 
 
ts.plot(SamplesOwnMod[,"sdY"])
ts.plot(SamplesOwnMod[,"sigma_time"])

cor(SamplesOwnMod[,"sdY"],
    SamplesOwnMod[,"sigma_time"])

#Number of Parameters
OwnSimMod$getNodeNames(stochOnly = TRUE) %>% grep("ZMat",x = ., invert = TRUE) %>% 
  length()

#estimated number of effective parameters (way to low ..)
comOwnSimMod$bOwnSimMod$getWAIC() %>% .$pWAIC

ZMatSimOwn %>% apply(.,1, acf,lag.max = 5, plot=FALSE)

## using QR ##
SimOwnConstQR <- list("N_Year"=ncol(ZMatSimOwnQR),
                      "N_AgeGroups"=nrow(ZMatSimOwnQR))

SimOwnDatQR <- list("ZMat"=ZMatSimOwnQR) 

OwnSimModQR <- nimbleModel(code=OwnMod_Repara_Diff_QR, 
                           constants = SimOwnConstQR, 
                           data=SimOwnDatQR, 
                           buildDerivs = FALSE)

cOwnSimModQR <- configureMCMC(OwnSimModQR, 
                              print = TRUE, useConjugacy = FALSE,
                              monitors = c("N_t","p","sigma_eps","muY","sdY",
                                           "QMat","k",
                                           "sigma_time","drift", 
                                           "b1","b2",
                                           "a","CovVec"
                              ),
                              enableWAIC = TRUE)
cOwnSimModQR$removeSamplers(c( 
  paste0("b1[",1:10,"]"),
  paste0("b2[",1:10,"]")),
  "sigma_time")

cOwnSimModQR$addSampler(target = "sdY",type="RW", log=TRUE)
cOwnSimModQR$addSampler(target = "sigma_time",type="RW", log=TRUE)

cOwnSimModQR$addSampler(target= c("b1[1:10]"),
                        type="AF_slice")

cOwnSimModQR$addSampler(target= c("b2[1:10]"),
                        type="AF_slice")


cOwnSimModQR$addSampler(target= c("b1[1:10]",
                                  "b2[1:10]"),
                        type="HMC",
                        control = list(nwarmup=2000,  #increase delta
                                       delta=0.81,
                                       maxTreeDepth = 10)) #increase tree depth

cOwnSimModQR$addSampler(target= c("b2[1:10]"),
                        type="HMC",
                        control = list(nwarmup=2000,  #increase delta
                                       delta=0.81,
                                       maxTreeDepth = 10)) #increase tree depth
# 
bOwnSimModQR <- buildMCMC(cOwnSimModQR)
comOwnSimModQR <- compileNimble(OwnSimModQR,
                                bOwnSimModQR)

InitsQRRest <- list("a"=0.3,
                    "sigma_time"=0.025,
                    "sigma_eps"=0.025,
                    "muY"=1,
                    "sdY"=0.3)
SamplesOwnModQR <- runMCMC(comOwnSimModQR$bOwnSimModQR, 
                           niter = 3000,
                           thin=1,
                           nburnin = 1500, 
                           nchains = 1)

SummaryOutput(SamplesOwnModQR, 
              params=c("sigma_eps",
                       "p","muY","sdY","sigma_time","beta","betaJump",
                       "drift","N_t","a")) %>% print(n=250) 

comLiuLiSimMod$bLiuLiSimMod$getWAIC()

ts.plot(SamplesOwnModQR[,"N_t[10]"])



##### Successive Jumps Comparison ##############################################
set.seed(42)
NAge <- 10

#US estimates
beta <- c(0.117,0.172,0.145,0.108,0.1,0.0996,0.0693,0.069,0.0651,0.055)
betaJ <- c(0.02,0.03,0.04,0.127,0.158,0.155,0.132,0.117,0.111,0.11)

QMat <- qr.Q(qr(cbind(beta,betaJ))) #QR decomposition of beta
QMat <- QMat*(-1)

NYear <-100
vt <- rnorm(n = NYear,mean = 0,sd = 0.1)
drift <- -0.15

RWDrift <- cumsum(drift+vt)
RWDrift <- RWDrift-mean(RWDrift) #sum to zero constraint


N_t <- c(rep(0,10),1,1,rep(0,20),1,1,1,rep(0,40),1,1,rep(0,23))
Y_t <- rnorm(n = NYear,mean = 1.5,sd = 0.4)


#Step 2, let J be a RV
J_t <- c(0,N_t[-1]*Y_t[-1])

LambdaMatLiu <- matrix(0, nrow = NAge, ncol = NYear)
for(x in 1:NAge){
  for (t in 1:NYear) {
    LambdaMatLiu[x,t] <- #ax[x] + 
      beta[x]*RWDrift[t] + betaJ[x]*J_t[t]+
      rnorm(n=1, sd=0.025)
  }
}

LambdaMatLiu %>% as.data.frame() %>% 
  mutate("Age"=as.factor(1:10), .before=1) %>% 
  pivot_longer(., cols=2:101, names_to = "Year", 
               values_to = "LogRate") %>% 
  mutate("Time"=rep(1:100,10)) %>% 
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()+
  geom_point()

ZMatSimLiu <- apply(LambdaMatLiu, 1, diff) %>% t() 
x11()
ZMatSimLiu %>% as.data.frame() %>% 
  mutate("Age"=as.factor(1:10), .before=1) %>% 
  pivot_longer(., cols=2:100, names_to = "Year", 
               values_to = "LogRate") %>% 
  mutate("Time"=rep(1:99,10)) %>% 
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()+
  geom_point()

J_t_AR <- numeric(NYear)
a <- 0.3
J_t_AR[1:2] <- 0
for(j in 3:(NYear)){
  J_t_AR[j] <- a*J_t_AR[j-1]+N_t[j]*Y_t[j]
}

EpsMat <- sapply(rep(10,100), function(x) rnorm(n=x, mean=0, sd=0.025))

## Simulation of log death rates 
LambdaMatOwn <- matrix(0, nrow = NAge, ncol = NYear)
for(x in 1:NAge){
  for (t in 1:NYear) {
    LambdaMatOwn[x,t] <- #ax[x] + 
      beta[x]*RWDrift[t] + betaJ[x]*J_t_AR[t]+
      EpsMat[x,t]
  }
}

LambdaMatOwn %>% as.data.frame() %>% 
  mutate("Age"=as.factor(1:10), .before=1) %>% 
  pivot_longer(., cols=2:101, names_to = "Year", 
               values_to = "LogRate") %>% 
  mutate("Time"=rep(1:100,10)) %>% 
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()

ZMatSimOwn <- apply(LambdaMatOwn, 1, diff) %>% t() 
x11()
ZMatSimOwn %>% as.data.frame() %>% 
  mutate("Age"=as.factor(1:10), .before=1) %>% 
  pivot_longer(., cols=2:100, names_to = "Year", 
               values_to = "LogRate") %>% 
  mutate("Time"=rep(1:99,10)) %>% 
  ggplot(data=., aes(x=Time, y=LogRate, group=Age, col=Age))+
  geom_line()



SimLiuConst <- list("N_Year"=ncol(ZMatSimLiu),
                    "N_AgeGroups"=nrow(ZMatSimLiu))

SimLiuDat <- list("ZMat"=ZMatSimLiu) 

LiuLiSimMod <- nimbleModel(code=LiuLi_Repara_Diff_OC, 
                           constants = SimLiuConst, 
                           data=SimLiuDat, buildDerivs = FALSE)

cLiuLiSimMod <- configureMCMC(LiuLiSimMod, 
                              print = TRUE, useConjugacy = TRUE,
                              monitors = c("N_t","muY",
                                           "sdY","p","sigma_eps",
                                           "beta","betaJump","k",
                                           "sigma_time","drift"
                              ),
                              enableWAIC = TRUE)
bLiuLiSimMod <- buildMCMC(cLiuLiSimMod)
comLiuLiSimMod <- compileNimble(LiuLiSimMod,
                                bLiuLiSimMod)

SamplesLiuLiSimMod <- runMCMC(comLiuLiSimMod$bLiuLiSimMod, 
                              niter = 4000,
                              thin=1,
                              nburnin = 2000, 
                              nchains = 1)

SummaryOutput(SamplesLiuLiSimMod, 
              params=c("sigma_eps",
                       "p","muY","sdY","sigma_time","beta","betaJump",
                       "drift","N_t")) %>% print(n=250) 

### OWN MODEL
SimOwnConst <- list("N_Year"=ncol(ZMatSimOwn),
                    "N_AgeGroups"=nrow(ZMatSimOwn))

SimOwnDat <- list("ZMat"=ZMatSimOwn) 

OwnSimMod <- nimbleModel(code=OwnMod_Repara_Diff_OC, 
                         constants = SimOwnConst, 
                         data=SimOwnDat, buildDerivs = FALSE)

cOwnSimMod <- configureMCMC(OwnSimMod, 
                            print = TRUE, useConjugacy = TRUE,
                            monitors = c("N_t","p","sigma_eps",
                                         "beta","betaJump","k",
                                         "sigma_time","drift", "beta","betaJump",
                                         "a","CovVec","muY","sdY"
                            ),
                            enableWAIC = TRUE)


cOwnSimMod$removeSamplers(c(
  paste0("b1[",1:10,"]"),
  paste0("b2[",1:10,"]"),
  "sigma_time",
  "a","muY"))


cOwnSimMod$addSampler(target = "sigma_time",type="RW", log=TRUE)
cOwnSimMod$addSampler(target = "a",type="slice")
cOwnSimMod$addSampler(target = "muY",type="slice")
cOwnSimMod$addSampler(target= c("b1[1:10]"),
                        type="AF_slice")

cOwnSimMod$addSampler(target= c("b2[1:10]"),
                        type="AF_slice")
for(j in 1:(NYear-1)){
  cOwnSimMod$addSampler(target = paste0("k[",j,"]"),
                        type = 'slice')
}

bOwnSimMod <- buildMCMC(cOwnSimMod)
comOwnSimMod <- compileNimble(OwnSimMod,
                              bOwnSimMod)

SamplesOwnMod <- runMCMC(comOwnSimMod$bOwnSimMod, 
                         niter = 4000,
                         thin=1,
                         nburnin = 200, 
                         nchains = 1)

SummaryOutput(SamplesOwnMod, 
              params=c("sigma_eps",
                       "p","muY","sdY","sigma_time","beta","betaJump",
                       "drift","N_t","a","muY","sdY")) %>% print(n=250) 


#### EVERYTHING BELOW IS NOT IN USE !!!!########################################

### Trying to estimate mode on death rates instead of diff death rates 
set.seed(42)
NAge <- 10

#US estimates
beta <- c(0.117,0.172,0.145,0.108,0.1,0.0996,0.0693,0.069,0.0651,0.055)
betaJ <- c(0.02,0.03,0.04,0.127,0.158,0.155,0.132,0.117,0.111,0.11)


QMat <- qr.Q(qr(cbind(beta,betaJ))) #QR decomposition of beta
QMat <- QMat*(-1)


# 1. Type of constraints 
#(drift = -0.4, sigma_time = 0.05)
NYear <-100
# vt <- rnorm(n = NYear,mean = 0,sd = 0.05)
# drift <- -0.4
vt <- rnorm(n = NYear,mean = 0,sd = 0.1)
drift <- -0.15


RWDrift <- cumsum(drift+vt)
RWDrift <- RWDrift-mean(RWDrift) #sum to zero constraint

p <- 0.1
N_t <- rbinom(n= NYear,size = 1,prob = p)
Y_t <- rnorm(n = NYear,mean = 1.5,sd = 0.4)


J_t <- numeric(NYear)
a <- 0.2
J_t[1:2] <- 0
for(j in 3:(NYear)){
  J_t[j] <- a*J_t[j-1]+N_t[j]*Y_t[j]
}

EpsMat <- sapply(rep(10,100), function(x) rnorm(n=x, mean=0, sd=0.025))

ax <- seq(-5,-0.5,0.5)

LambdaMatOwnQR <- matrix(0, nrow = NAge, ncol = NYear)
for(x in 1:NAge){
  for (t in 1:NYear) {
    LambdaMatOwnQR[x,t] <- ax[x] + 
      QMat[x,1]*RWDrift[t] + QMat[x,2]*J_t[t]+EpsMat[x,t]
  }
}

SimOwnConst <- list("N_Year"=ncol(LambdaMatOwnQR),
                    "N_AgeGroups"=nrow(LambdaMatOwnQR))

SimOwnDat <- list("ZMat"=LambdaMatOwnQR) 

OwnSimMod <- nimbleModel(code=OwnMod_Repara_QR, 
                         constants = SimOwnConst, 
                         data=SimOwnDat)

cOwnSimMod <- configureMCMC(OwnSimMod, 
                            print = TRUE, useConjugacy = FALSE,
                            monitors = c("N_t","p","sigma_eps",
                                         "k","kappa","alpha",
                                         "sigma_time","drift", "QMat","BMat",
                                         "a","muY","sdY"
                            ),
                            enableWAIC = TRUE)

cOwnSimMod$removeSamplers(c( 
  paste0("b1[",1:10,"]"),
  paste0("b2[",1:10,"]"),
  "sigma_time",
  "a","muY","alpha","p")
  )

# cOwnSimMod$removeSamplers("sdY",
#                           "sigma_time")



cOwnSimMod$addDefaultSampler(nodes=c("alpha[1:10]"), 
                         useConjugacy = TRUE)

cOwnSimMod$addSampler(target = "sigma_time",type="RW", log=TRUE)
cOwnSimMod$addSampler(target = "a",type="slice")
cOwnSimMod$addSampler(target = "p",type="slice")
cOwnSimMod$addSampler(target = "muY",type="slice")

cOwnSimMod$addSampler(target= c("b1[1:10]"),
                      type="AF_slice")

cOwnSimMod$addSampler(target= c("b2[1:10]"),
                      type="AF_slice")


bOwnSimMod <- buildMCMC(cOwnSimMod)

comOwnSimMod <- compileNimble(OwnSimMod,
                              bOwnSimMod)

SamplesOwnMod <- runMCMC(comOwnSimMod$bOwnSimMod, 
                         niter = 4000,
                         thin=1,
                         nburnin = 2000, 
                         nchains = 1)

SummaryOutput(SamplesOwnMod, 
              params=c("sigma_eps",
                       "p","muY","sdY","sigma_time","beta","betaJump",
                       "drift","N_t","a","muY","sdY")) %>% print(n=250) 


################ Simulated Liu, Li Model #######################################
set.seed(42)
NAge <- 10
betaJ <- c(0.01,0.08,0.08,0.38,0.31,0.07,0.0183,0.0137,0.008,0.03)
beta <- c(0.12,0.262,0.158,0.06,0.07,0.08,0.07,0.06,0.06,0.06)

# betaJ <- extraDistr::rdirichlet(n=1, alpha = rep(1/15,1))
# beta <- extraDistr::rdirichlet(n=1, alpha = rep(1/15,1))

NYear <-400
vt <- rnorm(n = NYear,mean = 0,sd = 0.1)
drift <- -0.2

RWDrift <- cumsum(drift+vt)
RWDrift <- RWDrift-mean(RWDrift)

p <- 0.1
N_t <- rbinom(n= NYear,size = 1,prob = p)
Y_t <- rnorm(n = NYear,mean = 2,sd = 1)

J_t <- c(0,N_t*Y_t)

ax <- seq(-5,0.5,0.5)

Indicies <- expand.grid("Age"=1:NAge,"Time"=1:NYear)

## Simlaution of Poisson Data
LambdaMat <- matrix(0, nrow = NAge, ncol=NYear)
for (t in 1:NYear) {
  for(x in 1:NAge){
    LambdaMat[x,t] <- ax[x] + beta[x]*RWDrift[t] + betaJ[x]*J_t[t+1]+rnorm(n=1, sd=0.1)
  }
}
#Poisson Distributed Deaths
phi <- 1000000
YMat <- apply(exp(LambdaMat)*100, 1:2, 
              function(gxt) stats::rnbinom(1,size = phi, prob = phi / (phi + gxt)))

ZMatSim <- matrix(0, nrow = NAge, ncol = NYear)
for (t in 1:NYear) {
  for(x in 1:NAge){
    ZMatSim[x,t] <- beta[x]*(drift+vt[t])+
      betaJ[x]*(J_t[t+1]-J_t[t])+rnorm(n=1, sd=0.1)
  }
}

SimConst <- list("N_AgeGroups"=nrow(ZMatSim),
                       "N_Year"=ncol(ZMatSim))

SimData <- list("ZMat"=ZMatSim) 

DiffLiuModel_Sim <- nimbleModel(code=LiuLi_Repara_Diff_OC, 
                            constants = SimConst, 
                            data=SimData)

cDiffModel <- configureMCMC(DiffLiuModel_Sim, 
                            print = TRUE)
#try HMC sampling
cDiffModel$addMonitors("k","beta","betaJump","N_t","p",
                       "sigma_time","sigma_eps")

params_to_remove <-  c("sdY","sigma_time",
                       paste0("b2[",1:NAge,"]"),
                       paste0("b1[",1:NAge,"]"))

cDiffModel$removeSamplers(params_to_remove)

cDiffModel$addSampler(target= c("b1[1:10]"),
                      type="AF_slice" ) #increase tree depth

cDiffModel$addSampler(target= c("b2[1:10]"),
                      type="AF_slice")  #increase tree depth

cDiffModel$addSampler(target= "sdY",
                      type="slice")  
cDiffModel$addSampler(target= "sigma_time",
                      type="slice")  


bDiffModel <- buildMCMC(cDiffModel)
comDiff <- compileNimble(DiffLiuModel_Sim,bDiffModel)


SamplesDiff <- runMCMC(comDiff$bDiffModel, 
                       niter = 4000,
                       thin=1,
                       nburnin = 2000, 
                       nchains = 2)

SummaryOutput(SamplesDiff,
              params=c("betaJump",
                       "beta","drift","p","sigma_eps",
                       "muY","sdY","sigma_time")) %>% print(n=150) #works

x11()
ts.plot(SamplesDiff$chain1[,"muY"])


########## QR SIMULATION #######################################################
#Use Real Beta Parameters as Starting point
betaJ <- c(0.01,0.08,0.08,0.38,0.31,0.07,0.0183,0.0137,0.008,0.03)
beta <- c(0.12,0.262,0.158,0.06,0.07,0.08,0.07,0.06,0.06,0.06)

QMat <- qr.Q(qr(cbind(beta,betaJ))) #QR decomposition of beta

NAge <- length(beta)
set.seed(420)
NYear <-200
vt <- rnorm(n = NYear,mean = 0,sd = 0.2)
drift <- -0.2

p <- 0.3
N_t <- rbinom(n= NYear+1,size = 1,prob = p)
Y_t <- rnorm(n = NYear+1,mean = 4,sd = 1)

J_t <- c(0,N_t*Y_t) #alt (but why?)
J_t <- c(0,0,N_t[-c(1:2)]*Y_t[-c(1:2)])

ZMatSim <- matrix(0, nrow = NAge, ncol = NYear)
for (t in 1:NYear) {
  for(x in 1:NAge){
    ZMatSim[x,t] <- QMat[x,1]*(drift+vt[t])+
                    QMat[x,2]*(J_t[t+1]-J_t[t])+rnorm(n=1, sd=0.1)
  }
}

SimConst <- list("N_AgeGroups"=nrow(ZMatSim),
                 "N_Year"=ncol(ZMatSim))

SimData <- list("ZMat"=ZMatSim) 

QRMod_Sim <- nimbleModel(code=LiuLi_Diff_QR,
                     constants = SimConst, 
                     data=SimData)

## test uncompiled function in uncompiled MCMC
QRModc <- configureMCMC(QRMod_Sim,monitors = c("k","sigma_eps",
                                           "QMat", "N_t","p",
                                           "J","muY","sdY","sigma_time",
                                           "drift","b1","b2","BMat","Y_t"),
                        enableWAIC = TRUE, print=TRUE)

params_to_remove <-  c(paste0("b1[",1:10,"]"),
                       paste0("b2[",1:10,"]"),
                       paste0("k[",1:NYear,"]"))

QRModc$removeSamplers(params_to_remove)

QRModc$addSampler(target= c("b1[1:10]"),
                  type="AF_slice")

QRModc$addSampler(target= c("b2[1:10]"),
                  type="AF_slice")

for(j in 1:NYear){
  QRModc$addSampler(target = paste0("k[",j,"]"),
                    type = 'slice'
  )
}

bDiffModelQR <- buildMCMC(QRModc)

comDiffQR <- compileNimble(QRMod_Sim,
                            bDiffModelQR)

SamplesQR <- runMCMC(comDiffQR$bDiffModelQR, 
                       niter = 5000,
                       nburnin = 2000,
                       thin=2,
                       nchains = 1) #set seed for reproducibility

SummaryOutput(SamplesQR,params=c("QMat","drift","p","sigma_eps",
                                   "muY","sdY","sigma_time")) %>% print(n=150) 


######## OWN MODEL SIMULATION #############################################
set.seed(42)
NAge <- 10
# beta <- c(0.181,0.172,0.149,0.078,0.06,0.07,0.07,0.08,0.08,0.06)
# betaJ <- c(0.02,0.07,0.07,0.0337,0.008,0.17,0.1083,0.28,0.21,0.03)

#US estimates
beta <- c(0.117,0.172,0.145,0.108,0.1,0.0996,0.0693,0.069,0.0651,0.055)
betaJ <- c(0.02,0.03,0.04,0.127,0.158,0.155,0.132,0.117,0.111,0.11)

# 1. Type of constraints 
NYear <-100
# vt <- rnorm(n = NYear,mean = 0,sd = 0.05)
# drift <- -0.4
vt <- rnorm(n = NYear,mean = 0,sd = 0.1)
drift <- -0.15


RWDrift <- cumsum(drift+vt)
RWDrift <- RWDrift-mean(RWDrift) #sum to zero constraint

p <- 0.1
N_t <- rbinom(n= NYear+1,size = 1,prob = p)
Y_t <- rnorm(n = NYear+1,mean = 2,sd = 0.4)

#Step 2, let J be a RV
J_t <- c(0,N_t[-1]*Y_t[-1])


#Step 3: AR MODEL
J_t <- numeric(NYear+1)
a <- 0.3
J_t[1:2] <- 0
for(j in 3:length(J_t)){
  J_t[j] <- a*J_t[j-1]+N_t[j]*Y_t[j]
}

ZMatSim <- matrix(0, nrow = NAge, ncol = NYear)
for (t in 1:NYear) {
  for(x in 1:NAge){
    ZMatSim[x,t] <- beta[x]*(drift+vt[t])+
                    betaJ[x]*(J_t[t+1]-J_t[t])+
                     rnorm(n=1, sd=0.05)
  }
}

#ages first, then years
Indicies <- expand.grid("Age"=1:NAge,
                        "Year"=1:NYear)

ZMatVec <- data.frame("ZVal"=as.vector(ZMatSim)) %>% 
  mutate("Age"=Indicies$Age,
         "Year"=Indicies$Year)

ZMatVec %>% 
  ggplot(data=., aes(x=Year, y=ZVal, group=Age))+
  geom_line(aes(col=factor(Age)))+
  scale_x_continuous(breaks=seq(1, 200, 5))


SimConst <- list("N_AgeGroups"=nrow(ZMatSim),
                 "N_Year"=ncol(ZMatSim))

SimData <- list("ZMat"=ZMatSim) 

DiffModel_Sim_QR <- nimbleModel(code=OwnMod_Repara_Diff_OC, 
                                constants = SimConst, 
                                data=SimData)

cDiffSim_QR <- configureMCMC(DiffModel_Sim_QR,
                             print = TRUE, useConjugacy = TRUE,
                             monitors = c("k","sigma_eps",
                                          "beta", "betaJump","p",
                                          "sigma_time",
                                          "drift","b1","b2","N_t","a",
                                          "muY","sdY"
                             ),
                             enableWAIC = TRUE)

params_to_remove <-  c(paste0("b1[",1:10,"]"),
                       paste0("b2[",1:10,"]"),
                       #paste0("k[",1:NYear,"]"),
                       #"drift","sigma_time","p",
                       #paste0("Y_t[",1:NYear,"]"),
                       "sdY", "muY"
                       )

cDiffSim_QR$removeSamplers(params_to_remove)

cDiffSim_QR$addSampler(target= c("b1[1:10]"),
                       type="AF_slice") 

cDiffSim_QR$addSampler(target= c("b2[1:10]"),
                       type="AF_slice")


cDiffSim_QR$addSampler(target = "sdY", type="slice")
cDiffSim_QR$addSampler(target = "muY", type="slice")

cDiffSim_QR$addSampler(target = "sigma_time", type="slice")
cDiffSim_QR$addDefaultSampler(nodes= c("drift","p"),
                            useConjugacy = TRUE)


for(j in 1:NYear){
  cDiffSim_QR$addSampler(target = paste0("k[",j,"]"),
                         type = 'slice')
  
  # cDiffSim_QR$addSampler(target=paste0("Y_t[",j,"]"),
  #                        type="slice")
}



bDiffSim_QR <- buildMCMC(cDiffSim_QR)
comDiffSim_QR <- compileNimble(DiffModel_Sim_QR,
                               bDiffSim_QR)

SamplesSim_QR <- runMCMC(comDiffSim_QR$bDiffSim_QR, 
                       niter = 4000,
                       thin=1,
                       nburnin = 2000, 
                       nchains = 1)

SummaryOutput(SamplesSim_QR,
              params=c("drift","p","sigma_eps",
                       "muY","sdY","sigma_time","a","beta","betaJump")) %>% 
  print(n=250) #works

x11()
ts.plot(SamplesSim_QR[,"sigma_time"])



Yts <- which(grepl("Y_t",colnames(SamplesSim_QR)))*N_t
YtNZ <- Yts[Yts!=0]
SamplesSim_QR[,YtNZ] %>% apply(., 2, mean) %>% mean()
SamplesSim_QR[,YtNZ] %>% apply(., 2, sd) %>% mean()


Y_t_Samp <- SamplesSim_QR[,grepl("Y_t",colnames(SamplesSim_QR))]

Y_t_Less <- apply(Y_t_Samp, 2, sample, size=150, replace=FALSE)
as.vector(Y_t_Less)

init_fun <- function(...) list(theta=c(0.2,0.8))

library(rstan)
set.seed(42)
dataStan <- list("K"=2,
                 "N"=length(as.vector(Y_t_Less)),
                 "y"=as.vector(Y_t_Less))

#attention takes long
options(mc.cores = parallel::detectCores())
MixMod <- stan(file = file.path(getwd(),"Stan Code/NormalMixture.stan"),
               data=dataStan, chains = 2, init=init_fun, 
               iter = 4000, warmup = 2000,
               thin=2, control = list(adapt_delta = 0.81,
                                      max_treedepth=11))

summary(MixMod)$summary %>% round(.,4)
##### MAQ Model ###############################################################
set.seed(42)
NAge <- 10
betaJ <- c(0.01,0.08,0.08,0.38,0.31,0.07,0.0183,0.0137,0.008,0.03)
beta <- c(0.12,0.262,0.158,0.06,0.07,0.08,0.07,0.06,0.06,0.06)

QMat <- qr.Q(qr(cbind(beta,betaJ))) #QR decomposition of beta

# betaJ <- extraDistr::rdirichlet(n=1, alpha = rep(1/15,1))
# beta <- extraDistr::rdirichlet(n=1, alpha = rep(1/15,1))

NYear <-100
vt <- rnorm(n = NYear,mean = 0,sd = 0.2)
drift <- -0.2

RWDrift <- cumsum(drift+vt)
RWDrift <- RWDrift-mean(RWDrift)
Q <- 5
p <- 0.3
N_t <- rbinom(n= NYear+Q,size = 1,prob = p)
Y_t <- rnorm(n = NYear+Q,mean = 2,sd = 1)

J_t <- numeric(NYear+Q)
a <- 0.4
J_t[1:Q] <- 0

for(j in (Q+1):length(J_t)){
  J_t[j] <- N_t[j]*Y_t[j]
  for(q in 1:min(j-1,Q)){
    J_t[j] <- J_t[j]+(a^q)*(N_t[j-q]*Y_t[j-q])
  }
}

ZMatSim <- matrix(0, nrow = NAge, ncol = NYear)
for (t in 1:NYear) {
  for(x in 1:NAge){
    ZMatSim[x,t] <- QMat[x,1]*(drift+vt[t])+QMat[x,2]*(J_t[t+1]-J_t[t])+rnorm(n=1, sd=0.1)
  }
}

SimConst <- list("N_AgeGroups"=nrow(ZMatSim),
                 "N_Year"=ncol(ZMatSim),#starting at t=2,
                 "Q"=Q)

SimData <- list("ZMat"=ZMatSim) 

Inits <- list("N_t"=rep(0, NYear+Q),
              "Y_t"=rnorm(n=NYear+Q))

DiffLiuModel_Sim <- nimbleModel(code=LCJumpOwn_MAQ_QR, 
                                constants = SimConst, 
                                data=SimData,
                                inits = Inits)


cDiffModel <- configureMCMC(DiffLiuModel_Sim, 
                            print = TRUE, 
                            useConjugacy = FALSE,
                            monitors = c("k","sigma_eps",
                              "QMat", "N_t","p",
                              "J","muY","sdY","sigma_time",
                              "drift","b1","b2","BMat","Y_t","a"))

RemoveParam <- c("a","sigma_time","sigma_eps","drift","muY","sdY",
                 paste0("b1[",1:10,"]"),
                 paste0("b2[",1:10,"]"))

cDiffModel$removeSamplers(RemoveParam)

cDiffModel$addSampler(target= c("b1[1:10]"),
                       type="AF_slice") 

cDiffModel$addSampler(target= c("b2[1:10]"),
                       type="AF_slice")

cDiffModel$addSampler(target="sdY",
                      type="slice")

cDiffModel$addSampler(target="drift",
                      type="slice")

cDiffModel$addSampler(target="muY",
                      type="slice")

cDiffModel$addSampler(target="a", #geht nicht
                      type="RW")

cDiffModel$addSampler(target="sigma_time",
                      type="slice")
cDiffModel$addSampler(target="sigma_eps",
                      type="slice")


bDiffModel <- buildMCMC(cDiffModel)
comDiff <- compileNimble(DiffLiuModel_Sim,bDiffModel)


SamplesDiff <- runMCMC(comDiff$bDiffModel, 
                       niter = 4000,
                       thin=1,
                       nburnin = 1500, 
                       nchains = 2)

SummaryOutput(SamplesDiff,
              params=c("QMat","drift","p","sigma_eps",
                       "muY","sdY","sigma_time","a")) %>% print(n=150) #works

x11()
ts.plot(SamplesDiff$chain1[,"sdY"])




#### Change of Y_t; now RV instead of Paramter      ############################
set.seed(42)
NAge <- 10
betaJ <- c(0.01,0.08,0.08,0.38,0.31,0.07,0.0183,0.0137,0.008,0.03)
beta <- c(0.12,0.262,0.158,0.06,0.07,0.08,0.07,0.06,0.06,0.06)
QMat <- qr.Q(qr(cbind(beta,betaJ))) #QR decomposition of beta

NYear <-100
vt <- rnorm(n = NYear,mean = 0,sd = 0.2)
drift <- -0.1

RWDrift <- cumsum(drift+vt)
RWDrift <- RWDrift-mean(RWDrift) #sum to zero constraint

p <- 0.2
N_t <- rbinom(n= NYear,size = 1,prob = p)
Y_t <- rnorm(n = NYear,mean = 0.5,sd = 0.3)

#J <- N_t

#Step 2, let J be a RV
J_t <- c(0,N_t[-1]*Y_t[-1])

#Step 3: Let J be an AR Model 
J_t <- numeric(NYear)
a <- 0.4
J_t[1] <- 0
for(j in 2:(NYear)){
  J_t[j] <- a*J_t[j-1]+N_t[j]*Y_t[j]
}

ax <- seq(-5,0.5,0.5)+3


## Simlaution of Poisson Data
LambdaMat <- matrix(0, nrow = NAge, ncol=NYear)
for(x in 1:NAge){
  for (t in 1:NYear) {
    LambdaMat[x,t] <- #ax[x] + 
      QMat[x,1]*RWDrift[t] + QMat[x,2]*J_t[t]+
      rnorm(n=1, sd=0.1)
  }
}

IndiciesSim <- expand.grid("Age"=1:NAge,
                        "Time"=1:NYear)

SimConst <- list("N_Year"=max(IndiciesSim$Time),
                  "N_AgeGroups"=max(IndiciesSim$Age))

SimData <- list("ZMat"=LambdaMat)



LiuLi_II_Sim <- nimbleModel(code=LiuLi_Repara_QR,
                          constants = SimConst, 
                          data=SimData)

## test uncompiled function in uncompiled MCMC
LiuLiModc <- configureMCMC(LiuLi_II_Sim,
                           monitors = c("k","sigma_eps",
                                        "QMat", "N_t","p",
                                        "muY","sdY","sigma_time",
                                        "drift","b1","b2"
                                                   ),
                           enableWAIC = TRUE, print=TRUE, 
                           useConjugacy = TRUE)

params_to_remove <- c(paste0("b1[",1:10,"]"),
                      paste0("b2[",1:10,"]"))

LiuLiModc$removeSamplers(params_to_remove)

LiuLiModc$addSampler(target= c("b1[1:10]"),
                     type="AF_slice")

LiuLiModc$addSampler(target= c("b2[1:10]"),
                     type="AF_slice")

LiuLiModc$addDefaultSampler(nodes= c("alpha[1:10]"),
                            useConjugacy = TRUE)


for(j in 1:NYear){
  LiuLiModc$addSampler(target = paste0("k[",j,"]"),
                       type = 'slice'
  )
}

bDiffLiuLiMod <- buildMCMC(LiuLiModc)

comLiuLi <- compileNimble(LiuLi_II_Sim,
                          bDiffLiuLiMod)

SamplesLiuLi <- runMCMC(comLiuLi$bDiffLiuLiMod, 
                        niter = 4000,
                        nburnin = 2000,
                        thin=1,
                        nchains = 1) #set seed for reproducibility

SummaryOutput(SamplesLiuLi,params=c("sigma_eps",
                                    "QMat","p",
                                    "muY","sdY","sigma_time",
                                    "drift","alpha","N_t")) %>% print(n=250) 

Nts <- SummaryOutput(SamplesLiuLi,
                     params=c("N_t")) %>% select(mean)

pull(Nts)-N_t

Y_t[70]

x11()
ts.plot(SamplesLiuLi[,"sdY"])
### Liu,Li Route II Estimation via Differencing ###
set.seed(42)
NAge <- 10
betaJ <- c(0.01,0.08,0.08,0.38,0.31,0.07,0.0183,0.0137,0.008,0.03)
beta <- c(0.12,0.262,0.158,0.06,0.07,0.08,0.07,0.06,0.06,0.06)
QMat <- qr.Q(qr(cbind(beta,betaJ))) #QR decomposition of beta

NYear <-400
vt <- rnorm(n = NYear,mean = 0,sd = 0.2)
drift <- -0.1

RWDrift <- cumsum(drift+vt)
RWDrift <- RWDrift-mean(RWDrift) #sum to zero constraint

p <- 0.2
N_t <- rbinom(n= NYear,size = 1,prob = p)
Y_t <- rnorm(n = NYear,mean = 4,sd = 1)

#J <- N_t

#Step 2, let J be a RV
J <- N_t*Y_t

LambdaMat <- matrix(0, nrow = NAge, ncol=NYear)
for(x in 1:NAge){
  for (t in 1:NYear) {
    LambdaMat[x,t] <- ax[x] + QMat[x,1]*RWDrift[t] + QMat[x,2]*J[t]+
                      rnorm(n=1, sd=0.1)
  }
}

ZMatSim <- apply(LambdaMat, 1, diff) %>% t()


SimConst <- list("N_AgeGroups"=nrow(ZMatSim),
                 "N_Year"=ncol(ZMatSim))

SimData <- list("ZMat"=ZMatSim) 

QRMod_Sim <- nimbleModel(code=LiuLi_Diff_QR,
                         constants = SimConst, 
                         data=SimData)

## test uncompiled function in uncompiled MCMC
LiuLiModc <- configureMCMC(QRMod_Sim,monitors = c("k","sigma_eps","sigma_time",
                                                   "QMat","p","b1","b2","Y_t","N_t",
                                                   "drift","muY","sdY","J"),
                           enableWAIC = TRUE, print=TRUE)

params_to_remove <- c(paste0("b1[",1:10,"]"),
                      paste0("b2[",1:10,"]"))

LiuLiModc$removeSamplers(params_to_remove)

LiuLiModc$addSampler(target= c("b1[1:10]"),
                     type="AF_slice")

LiuLiModc$addSampler(target= c("b2[1:10]"),
                     type="AF_slice")

bDiffLiuLiMod <- buildMCMC(LiuLiModc)

comLiuLi <- compileNimble(QRMod_Sim,
                          bDiffLiuLiMod)

SamplesLiuLi <- runMCMC(comLiuLi$bDiffLiuLiMod, 
                        niter = 5000,
                        nburnin = 2000,
                        thin=1,
                        nchains = 1) #set seed for reproducibility

SummaryOutput(SamplesLiuLi,params=c("QMat","drift","muY","sigma_eps",
                                    "sigma_time","p","sdY")) %>% print(n=190) 

x11()
ts.plot(SamplesLiuLi[,"QMat[4, 2]"])




######## SIMULATION OF LIU,LI Route I Approach #################################
set.seed(42)
NAge <- 10
betaJ <- c(0.01,0.08,0.08,0.38,0.31,0.07,0.0183,0.0137,0.008,0.03)
beta <- c(0.12,0.262,0.158,0.06,0.07,0.08,0.07,0.06,0.06,0.06)

#use QR decomposition for betas
QMat <- qr.Q(qr(cbind(beta,betaJ))) #QR decomposition of beta

QMat <- QMat*(-1) #if all are negative


NYear <-150
vt <- rnorm(n = NYear,mean = 0,sd = 0.2)
drift <- -0.1

RWDrift <- cumsum(drift+vt)
RWDrift <- RWDrift-mean(RWDrift) #sum to zero constraint

RWDrift <- RWDrift - tail(RWDrift,1) #corner Constraint


p <- 0.2
N_t <- rbinom(n= NYear,size = 1,prob = p)
Y_t <- rnorm(n = NYear,mean = 4,sd = 1)

Y_t <- rep(4, NYear)

J_t <- c(N_t*Y_t,1*4) #Let there be a Jump in the end of size mean

ax <- seq(-5,0.5,0.5)

Indicies <- expand.grid("Age"=1:NAge,"Time"=1:NYear)

## Simlaution of Poisson Data
LambdaMat <- matrix(0, nrow = NAge, ncol=NYear)
for (t in 1:NYear) {
  for(x in 1:NAge){
    LambdaMat[x,t] <- ax[x] + QMat[x,1]*RWDrift[t] + QMat[x,2]*J_t[t+1]+rnorm(n=1, sd=0.1)
  }
}
#Poisson Distributed Deaths
phi <- 1000000
YMat <- apply(exp(LambdaMat)*50, 1:2, 
              function(gxt) stats::rnbinom(1,size = phi, 
                                           prob = phi / (phi + gxt)))


IndiciesSim <- expand.grid("Age"=1:NAge,
                           "Time"=1:NYear)

SimConst <- list("N_Year"=max(IndiciesSim$Time),
                 "N_AgeGroups"=max(IndiciesSim$Age),
                 "N"=length(YMat),
                 "age"=Indicies$Age,
                 "year"=Indicies$Time,
                 "Offset"=rep(30,length(YMat)))
#Data
SimData <- list("y"=as.vector(YMat))




LiuLiI_Sim <- nimbleModel(code=LiuLiModel_Corner_QR,
                          constants = SimConst, 
                          data=SimData)

## test uncompiled function in uncompiled MCMC
LiuLiModc <- configureMCMC(LiuLiI_Sim,monitors = c("k","sigma_eps",
                                                   "QMat", "N_t","p",
                                                   "J","muY","sdY","sigma_time",
                                                   "drift","b1","b2","Y_t"),
                           enableWAIC = TRUE, print=TRUE, useConjugacy = FALSE)

params_to_remove <-  c(paste0("b1[",1:10,"]"),
                       paste0("b2[",1:10,"]"))

LiuLiModc$removeSamplers(params_to_remove)

LiuLiModc$addSampler(target= c("b1[1:10]"),
                     type="AF_slice")

LiuLiModc$addSampler(target= c("b2[1:10]"),
                     type="AF_slice")

for(j in 1:NYear){
  LiuLiModc$addSampler(target = paste0("k[",j,"]"),
                       type = 'slice'
  )
}

bDiffLiuLiMod <- buildMCMC(LiuLiModc)

comLiuLi <- compileNimble(LiuLiI_Sim,
                          bDiffLiuLiMod)

SamplesLiuLi <- runMCMC(comLiuLi$bDiffLiuLiMod, 
                        niter = 12000,
                        nburnin = 10000,
                        thin=1,
                        nchains = 1) #set seed for reproducibility

SummaryOutput(SamplesLiuLi,params=c("QMat","drift","p","sigma_eps",
                                    "muY","sdY","sigma_time","k")) %>% print(n=190) 

x11()
ts.plot(SamplesLiuLi[,"k[8]"])

########### New Para with Poisson Dist###########################
set.seed(42)
NAge <- 10
betaJ <- c(0.01,0.08,0.08,0.38,0.31,0.07,0.0183,0.0137,0.008,0.03)
beta <- c(0.12,0.262,0.158,0.06,0.07,0.08,0.07,0.06,0.06,0.06)
QMat <- qr.Q(qr(cbind(beta,betaJ))) #QR decomposition of beta

NYear <-100
vt <- rnorm(n = NYear,mean = 0,sd = 0.2)
drift <- -0.1

RWDrift <- cumsum(drift+vt)
RWDrift <- RWDrift-mean(RWDrift) #sum to zero constraint

p <- 0.2
N_t <- rbinom(n= NYear,size = 1,prob = p)
Y_t <- rnorm(n = NYear,mean = 0.5,sd = 0.2)

#J <- N_t

#Step 2, let J be a RV
J <- c(0,N_t[-1]*Y_t[-1])


Indicies <- expand.grid("Age"=1:NAge,
                        "Time"=1:NYear)

## Simlaution of Poisson Data
LambdaMat <- matrix(0, nrow = NAge, ncol=NYear)
for(x in 1:NAge){
  for (t in 1:NYear) {
    LambdaMat[x,t] <- #ax[x] + 
      QMat[x,1]*RWDrift[t] + QMat[x,2]*J[t]+
      rnorm(n=1, sd=0.1)
  }
}

#Sample from Dist
YMat <- apply(exp(LambdaMat)*10, 1:2, 
              function(gxt) stats::rpois(1,lambda = gxt))

SimConst2 <- list("N_Year"=max(IndiciesSim$Time),
                  "N_AgeGroups"=max(IndiciesSim$Age),
                  "Offset"=matrix(data = 10,nrow = NAge, 
                                  ncol = NYear))

SimData2 <- list("YMat"=YMat)

LiuLi_Sim_Pois <- nimbleModel(code=LiuLi_Repara_2_Poi,
                              constants = SimConst2, 
                              data=SimData2)

## test uncompiled function in uncompiled MCMC
LiuLiModc_Pois <- configureMCMC(LiuLi_Sim_Pois,
                                monitors = c("k","sigma_eps",
                                             "QMat", "N_t","p",
                                             "muY","sdY","sigma_time",
                                             "drift","b1","b2",
                                             "ZMat"
                                ),
                                enableWAIC = TRUE, print=TRUE, useConjugacy = TRUE)

params_to_remove <- c(paste0("b1[",1:10,"]"),
                      paste0("b2[",1:10,"]"))

LiuLiModc_Pois$removeSamplers(params_to_remove)

LiuLiModc_Pois$addSampler(target= c("b1[1:10]"),
                          type="AF_slice")

LiuLiModc_Pois$addSampler(target= c("b2[1:10]"),
                          type="AF_slice")

bDiffLiuLiMod_Pois <- buildMCMC(LiuLiModc_Pois)

comLiuLi_Pois <- compileNimble(LiuLi_Sim_Pois,
                               bDiffLiuLiMod_Pois)

SamplesLiuLi_Pois <- runMCMC(comLiuLi_Pois$bDiffLiuLiMod_Pois, 
                             niter = 4000,
                             nburnin = 2000,
                             thin=1,
                             nchains = 1) #set seed for reproducibility

SummaryOutput(SamplesLiuLi,params=c("sigma_eps",
                                    "QMat","p",
                                    "muY","sdY","sigma_time",
                                    "drift","alpha","N_t")) %>% print(n=250) 

#### Nimble Jump Repara Test ############
NTime <- 200
a <- 0.3
p <- 0.1
N_t <- rbinom(n= NTime,size = 1,prob = p)
Y_t <- rnorm(n = NTime,mean = 4,sd = 1.5)

J_t <- N_t*Y_t
J_t <- numeric(NTime)
a <- 0.3
J_t[1] <- 0
for(j in 3:length(J_t)){
  J_t[j] <- a*J_t[j-1]+N_t[j]*Y_t[j]
}

Y <- rnorm(n = NTime,mean = J_t,sd=0.2)
ts.plot(Y[1:10])

SimConstJTest <- list("NTime"=NTime)

SimDataJTes <- list("y"=Y) 

JumpTestMod <- nimbleModel(code=NimbleJumpTest, 
                           constants = SimConstJTest, 
                           data=SimDataJTes)

cJumpTestMod <- configureMCMC(JumpTestMod, 
                              print = TRUE, useConjugacy = FALSE,
                              monitors = c("a","y","N_t","muY",
                                           "sdY","p","sigma_eps"
                              ),
                              enableWAIC = TRUE)
cJumpTestMod$removeSamplers(c("p","muY","sdY"))
cJumpTestMod$addDefaultSampler(nodes=c("muY","p"), 
                               useConjugacy = TRUE)

cJumpTestMod$addSampler(target = "sdY",type="RW", log=TRUE)

bJumpTestMod <- buildMCMC(cJumpTestMod)
comJumpTestMod <- compileNimble(JumpTestMod,
                                bJumpTestMod)

SamplesJumpTestMod <- runMCMC(comJumpTestMod$bJumpTestMod, 
                              niter = 4000,
                              thin=4,
                              nburnin = 2000, 
                              nchains = 1)
SamplesJumpTestMod

SummaryOutput(SamplesJumpTestMod, 
              params=c("a","p","muY","sdY","sigma_eps")) %>% 
  print(n=150) 

ts.plot(SamplesJumpTestMod[,"p"])

## Difference model 
NTime <- 200
a <- 0.3
p <- 0.15
N_t <- rbinom(n= NTime,size = 1,prob = p)
Y_t <- rnorm(n = NTime,mean = 2,sd = 0.5)

#äJ_t <- N_t*Y_t
J_t <- numeric(NTime)
a <- 0.3
J_t[1] <- 0
for(j in 2:length(J_t)){
  J_t[j] <- a*J_t[j-1]+N_t[j]*Y_t[j]
}

eps <- rnorm(n = NTime, mean = 0, sd=0.1)
Y <- J_t+eps
Ytilde <- c(eps[1],diff(Y))
Ytilde2 <- diff(Y)

SimConstJTest <- list("NTime"=length(Ytilde))
SimDataJTes <- list("y"=Ytilde) 

SimConstJTest2 <- list("NTime"=length(Ytilde2))
SimDataTes2 <- list("y"=Ytilde2)


JumpTestMod <- nimbleModel(code=NimbleJumpTestDiff, 
                           constants = SimConstJTest, 
                           data=SimDataJTes)


cJumpTestMod <- configureMCMC(JumpTestMod, 
                              print = TRUE, useConjugacy = FALSE,
                              monitors = c("a","y","N_t","muY",
                                           "sdY","p","sigma_eps"
                              ),
                              enableWAIC = TRUE)

cJumpTestMod$removeSamplers(c("p","muY","sdY"))

cJumpTestMod$addDefaultSampler(nodes=c("muY","p"), 
                               useConjugacy = TRUE)

cJumpTestMod$addSampler(target = "sdY",type="RW", log=TRUE)

bJumpTestMod <- buildMCMC(cJumpTestMod)
comJumpTestMod <- compileNimble(JumpTestMod,
                                bJumpTestMod)

SamplesJumpTestMod <- runMCMC(comJumpTestMod$bJumpTestMod, 
                              niter = 4000,
                              thin=1,
                              nburnin = 2000, 
                              nchains = 1)


SummaryOutput(SamplesJumpTestMod, 
              params=c("a","p","muY","sdY","sigma_eps")) %>% 
  print(n=150) 


JumpTestMod2 <- nimbleModel(code=NimbleJumpTestDiff_2, 
                            constants = SimConstJTest2, 
                            data=SimDataTes2)


cJumpTestMod2 <- configureMCMC(JumpTestMod2, 
                               print = TRUE, useConjugacy = FALSE,
                               monitors = c("a","y","N_t","muY",
                                            "sdY","p","sigma_eps"
                               ),
                               enableWAIC = TRUE)

cJumpTestMod2$removeSamplers(c("p","muY","sdY"))

cJumpTestMod2$addDefaultSampler(nodes=c("muY","p"), 
                                useConjugacy = TRUE)

cJumpTestMod2$addSampler(target = "sdY",type="RW", log=TRUE)
bJumpTestMod2 <-  buildMCMC(cJumpTestMod2)
comJumpTestMod2 <- compileNimble(JumpTestMod2,
                                 bJumpTestMod2)

SamplesJumpTestMod2 <- runMCMC(comJumpTestMod2$bJumpTestMod2, 
                               niter = 4000,
                               thin=1,
                               nburnin = 2000, 
                               nchains = 1)

SummaryOutput(SamplesJumpTestMod2, 
              params=c("a","p","muY","sdY","sigma_eps")) %>% 
  print(n=150) 



#Find the Error of Function##


#Test Function in R (works)
debugonce(CovZMat)
CovZMat(N_t = N_t,a = 0.2,t = 3,sigma2 = 2)

# Test uncompiled function in uncompiled model
InitsTest <- list("muY"=3,
                  "sdY"=1, "a"=0.3,
                  "N_t"=rbinom(n= NTime,size = 1,prob = p),
                  "p"=0.3,"sigma_eps"=0.2)

JumpTestMod2 <- nimbleModel(code=NimbleJumpTestDiff_2, 
                            constants = SimConstJTest2, 
                            data=SimDataTes2)


JumpTestMod2$x
JumpTestMod2$simulate("x")
?buildMCMC
# test uncompiled function in uncompiled MCMC
BTest <- buildMCMC(JumpTestMod2, monitors=c("x","a","sdY","muY"))
BTest$run(niter=3) #works
BTest$mvSamples$x



#Test compiled Function on its own
CCovMat <- compileNimble(CovZMat)
CCovMat(N_t = N_t, a = 0.9,t = 18,sigma2 = 2)

#Test compiled function in compiled model and in MCMC
CTest <- compileNimble(JumpTestMod2)

CTest2 <- compileNimble(BTest, project = JumpTestMod2)

