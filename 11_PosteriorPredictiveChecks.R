#### POSTERIOR PREDICTIVE CHEKCS ##############
library(nimble);library(tidyverse); library(nimbleHMC)

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 


load(file = file.path(getwd(),"Results/SamplesUS_Own.RData"))
load(file = file.path(getwd(),"Results/SamplesIt_Own.RData"))
load(file = file.path(getwd(),"Results/SamplesEW_Own.RData"))
load(file = file.path(getwd(),"Data/CovidData.RData")) #load Data

SamplesItSingle <- do.call(rbind, 
                           SamplesOwnNew_It) %>% 
  data.frame() %>% 
  mutate(Country="It")

SamplesUSSingle <- do.call(rbind, 
                           SamplesOwn_US_Rep) %>% 
  data.frame() %>% 
  mutate(Country="US")

SamplesEWSingle <- do.call(rbind, 
                           SamplesOwn_EW_Rep) %>% 
  data.frame() %>% 
  mutate(Country = "UK")


#Function for Creation of Replicated Data of Z Matrix as well as Lambda (log Death rates)
ReplicatedDataFun <- function(Samples, n, SingleVar = FALSE){
## Find positions of mu and sigma in Samples
  muPos <- grep("mu",colnames(Samples), 
                fixed = TRUE)
  
  if(SingleVar==TRUE){ #if sigma_eps is only variance
    sdPos <-grep("sigma_eps", colnames(Samples), 
                 fixed = TRUE) 
  } else {
    sdPos <- grep("sigma_squared", colnames(Samples), 
                  fixed = TRUE)
  }
  
  muPosReal <- muPos[1:(length(muPos)-1)] #without muY
  
  
  Z_Rep <- matrix(data = 0,
                  nrow = nrow(Samples),
                  ncol = n)
  if(SingleVar == TRUE){
    for(s in 1:nrow(Samples)){
      Z_Rep[s,] <- rnorm(n = n, 
                         mean = as.numeric(Samples[s,muPosReal]),
                         sd = Samples[s, sdPos])
    } 
  } else {
    for(s in 1:nrow(Samples)){
      Z_Rep[s,] <- rnorm(n = n, 
                         mean = as.numeric(Samples[s,muPosReal]),
                         sd = as.numeric(sqrt(Samples[s,sdPos])))
    }
  }
  
  return(Z_Rep)
}
ReplicatedLambdaFun <- function(Samples,LambdaMat, Z_Rep){
  #Creation of New log Death Rates
  LambdaArray_Rep <- array(data = 0, dim = c(nrow(LambdaMat), #age 
                                             ncol(LambdaMat), #time
                                             nrow(Samples)))  #Iteration
  
  LambdaArray_Rep[,1,] <- log(LambdaMat[,1]) #First year is the same
  for(s in 1:nrow(Samples)){
    Z_Mat_Hat <- matrix(Z_Rep[s,],nrow = nrow(LambdaMat), byrow = FALSE) #transformation of Z_Rep vector into matrix
    LambdaArray_Rep[,2,s] <- log(LambdaMat[,1]) + Z_Mat_Hat[,1] #Second Year 
    for(t in 3:ncol(LambdaMat)){ #starting for year 3, in recursive
      #LambdaArray_US_Rep[,t,s] <- LambdaArray_US_Rep[,(t-1),s] + Z_Mat_Hat[,(t-1)] OOS FC
      LambdaArray_Rep[,t,s] <- log(LambdaMat[,(t-1)]) + Z_Mat_Hat[,(t-1)] #In-Sample FC
    }
  }
  return(LambdaArray_Rep)
}

Z_RepUS <- ReplicatedDataFun(Samples = SamplesOwn_US_Rep,
                             n = length(ZMatUS))

Z_RepWar <- ReplicatedDataFun(Samples = SamplesLiuLi_OC,
                             n = length(ZMatWar), SingleVar = TRUE)

Z_RepWar_Repa <- ReplicatedDataFun(Samples = SamplesLiuLi_Repa,
                              n = length(ZMatWar), SingleVar = FALSE)

Z_RepWar_Own <- ReplicatedDataFun(Samples = Samples_OwnMod_OC,
                                  n = length(ZMatWar), SingleVar = TRUE)

LambdaArray_US <- ReplicatedLambdaFun(Samples = SamplesOwn_US_Rep,
                                      LambdaMat = LambdaMatUS, Z_Rep = Z_RepUS)

LambdaArray_War <- ReplicatedLambdaFun(Samples = SamplesLiuLi_OC,
                                       LambdaMat = LambdaMatWar, 
                                       Z_Rep = Z_RepWar)

LambdaArray_War_Repa <- ReplicatedLambdaFun(Samples = SamplesLiuLi_Repa,
                                       LambdaMat = LambdaMatWar, 
                                       Z_Rep = Z_RepWar_Repa)

LambdaArray_War_Own <- ReplicatedLambdaFun(Samples = Samples_OwnMod_OC,
                                           LambdaMat = LambdaMatWar,
                                           Z_Rep = Z_RepWar_Own)



LambdaPlot <- function(LambdaArray, Quants=c(0.1,0.9),
                       LambdaMat, LambdaVec,
                       AgesUsed = 1:10,
                       startYear){
  LowerPI <- LambdaArray %>% apply(., c(1,2),quantile, Quants[1]) %>% 
    data.frame() %>% 
    mutate("AgeInd"=1:nrow(LambdaMat), .before = 1) %>% 
    pivot_longer(., cols = 2:ncol(.),names_to = "Year", values_to = "lPI")
  
  UpperPI <- LambdaArray %>% apply(., c(1,2),quantile, Quants[2]) %>% 
    data.frame() %>% 
    mutate("AgeInd"=1:nrow(LambdaMat), .before = 1) %>% 
    pivot_longer(., cols = 2:ncol(.),names_to = "Year", values_to = "uPI")
  Plot <- LowerPI %>% 
    mutate("Year"=rep(startYear:(startYear+ncol(LambdaMat)-1),
                      nrow(LambdaMat))) %>%
    mutate("uPI"=UpperPI$uPI) %>% 
    filter(AgeInd %in% AgesUsed) %>% 
    ggplot(data=., aes(x=Year, group=AgeInd))+
    geom_ribbon(aes(ymin=lPI, ymax=uPI, fill=AgeInd))+
    ylab("Log Death Rate")+
    geom_line(data=filter(LambdaVec,NewAgeInd %in% AgesUsed), 
              aes(x=Year, y=log(Rate), group=NewAgeInd),
              col="red")
  return(Plot)
}

LambdaPlot(LambdaArray = LambdaArray_War_Repa, 
           LambdaMat = LambdaMatWar, 
           LambdaVec = LambdaVec, startYear = 1901, AgesUsed = c(4))


ZMatPlot <- function(ZMat, Quants = c(0.1,0.9), Z_Rep, AgesUsed = 1:10){
  AgeTimeGrid <- expand.grid("Age"=1:nrow(ZMat), "Time"=1:ncol(ZMat))
  ZMaLong <- ZMat %>% 
    as.vector() %>% data.frame("Val"=.) %>% 
    mutate(AgeTimeGrid)
  
  ### Plot Of ZMat ###
  LowerPIZ <- Z_Rep %>% apply(., 2,quantile, min(Quants)) %>% 
    data.frame("LPI"=.) %>% 
    mutate(AgeTimeGrid) 
  
  HigherPIZ <- Z_Rep %>% apply(., 2,quantile, max(Quants)) %>% 
    data.frame("HPI"=.) %>% 
    mutate(AgeTimeGrid) %>% 
    mutate("Mean"=apply(Z_Rep,2,mean))
  
  Plot <- LowerPIZ %>% 
    mutate("HPI"=HigherPIZ$HPI,
           "Mean"=HigherPIZ$Mean) %>% 
    filter(Age %in% AgesUsed) %>% 
    ggplot(aes(x=Time, group=Age))+
    geom_ribbon(aes(ymin=LPI, ymax=HPI, fill=Age))+
    geom_line(aes(y = Mean),col="red")+
    geom_line(data=filter(ZMaLong, Age %in% AgesUsed), 
              aes(x=Time, y=Val, group=Age),
              col="black")
  return(Plot)
}

ZMatPlot(ZMat = ZMatWar, 
         Quants = c(0.05,0.95),Z_Rep = Z_RepWar_Repa,AgesUsed = 6)

ZMatPlot(ZMat = ZMatWar, 
         Quants = c(0.05,0.95),Z_Rep = Z_RepWar,AgesUsed = 6)

#Creation of Likelihood Matrix
LikeMat <- matrix(0,nrow = nrow(Samples),
                  ncol = length(ZMat))
for(s in 1:nrow(Samples)){
  LikeMat[s,] <- dnorm(x=as.vector(ZMat),
                       mean = as.numeric(Samples[s,muPosReal]),
                       sd = as.numeric(sqrt(Samples[s,sdPosReal])))
}



#WAIC calculation by HAND (different to Nimble but same to loo package)
lppd <- LikeMat %>% apply(., 2, mean) %>% log() %>% sum()
pWaic <- LikeMat %>% log() %>% apply(., 2, var) %>% sum

lppd-pWaic
loo::waic(log(LikeMat))


################## loo and WAIC of actual data##################################
LikeMat_Spain <- LikelihoodMatrixFun(Samples = SamplesOwn_Sp$chain2,
                                   n = length(ZMatSp), ZMat = ZMatSp,
                                   SingleVar = FALSE)


loo::waic(log(LikeMat_Spain))
loo::loo(log(LikeMat_Spain))


library(bayesplot)
color_scheme_set("purple")
ppc_dens_overlay(y = as.vector(ZMatWar), 
                 yrep = Z_RepWar[1:300,])

x11()
ppc_dens_overlay_grouped(y = as.vector(ZMatWar), 
                 yrep = Z_RepWar[1:300,],
                 group = rep(1:10, ncol(ZMatWar)))

ppc_stat_grouped(y = as.vector(ZMatWar),
                 yrep = Z_RepWar_Repa[1:300,],
                 stat="max", binwidth = 0.005,
                 group = rep(1:10, ncol(ZMatWar)))+
  ggtitle("Maximum Value by Age Group")

ppc_dens_overlay(y = as.vector(log(LambdaMatWar)), 
                 yrep = t(apply(LambdaArray_War, 3, c))[1:300,])

ppc_stat_grouped(y = as.vector(log(LambdaMatWar)),
                 yrep =  t(apply(LambdaArray_War_Own, 3, c))[1:300,],
                 stat="max", binwidth = 0.005,
                 group = rep(1:10, ncol(LambdaMatWar)))+
  ggtitle("Maximum Value by Age Group")



#Summed Pearson Residuals
muPos <- grep("mu",colnames(SamplesLiuLi_Repa), 
              fixed = TRUE)
sdPos <- grep("sigma", colnames(SamplesLiuLi_Repa), 
              fixed = TRUE)

#ri <-  (yi - E(yi)) / sd(yi)
MeanVals <-SamplesLiuLi_Repa[,muPos[1:(length(muPos)-1)]] 
sdVals <- SamplesLiuLi_Repa[,sdPos[1:(length(sdPos)-2)]]

#Chi-Squared Discrepency 
ChiSq <- sum((as.numeric(ZMatWar) - apply(MeanVals,2,mean))^2/
  apply(sqrt(sdVals),2,mean))

## ChiSqTestValue
ChiSq_Rep <- numeric(2000)
for(s in 1:nrow(SamplesLiuLi_OC)){
  ChiSq_Rep[s] <- sum((Z_RepWar[s,] - MeanVals[s,])^2/
                        sqrt(sdVals[s,]))
}
hist(ChiSq_Rep,breaks = 100,probability = TRUE)
abline(v = ChiSq)



#Pearson Residuals
xSeq <- seq(-4,4,0.05)
png(filename = file.path(getwd(),"Bilder/PearsonResiduals.png"))
hist(RiSum,breaks = 100,probability = TRUE)
lines(xSeq,dnorm(xSeq,0,1),col="red")
dev.off()

#Sum of Residuals
RiSumVal <- sum(RiSum)

Ri_Sample <- numeric(nrow(SamplesOwnNew_US))
for(s in 1:nrow(SamplesOwnNew_US)){
  Ri_Sample[s] <- sum((ppRep[s,]-MeanVals[s,])/sdVals[s,])
}
hist(Ri_Sample,breaks = 100)
abline(v=RiSumVal,col="red")

#Log Prob
comOwnMod_US$OwnMod_US$getLogProb()

#Calculate log probability by hand
ppLog <- numeric(length(ZMatUSVec))
for(i in 1:length(ZMatUSVec)){
  LogS <- numeric(nrow(SamplesOwnNew_US))
  for(s in 1:nrow(SamplesOwnNew_US)){
    LogS[s] <- 
      dnorm(ZMatUSVec[i], 
          mean = SamplesOwnNew_US[s,
                                  muPos[1:(length(muPos)-1)][i]
          ],
          sd = sqrt(SamplesOwnNew_US[s,
                                     sdPos[1:(length(sdPos)-2)][i]
          ]
          ))  
  }
  ppLog[i] <- log(mean(LogS))
}
#lppd = sum(ppLog)

#do it for every draw from posterior predictive
ppLogMat <- matrix(data = 0,nrow = nrow(SamplesOwnNew_US),
                   ncol = length(ZMatUSVec))

for(i in 1:length(ZMatUSVec)){
  for(s in 1:nrow(SamplesOwnNew_US)){
    ppLogMat[s,i] <- 
      dnorm(Y_Rep[s,i], 
            mean = SamplesOwnNew_US[s,
                                    muPos[1:(length(muPos)-1)][i]
            ],
            sd = sqrt(SamplesOwnNew_US[s,
                                       sdPos[1:(length(sdPos)-2)][i]
            ]
            ),log = TRUE)  
  }
}

ppc_stat(y = ppLog, 
         yrep=ppLogMat, 
         stat="sum", binwidth = 5)+ title("lppd")




dim(ZMatUS)


Z_RepUS <- ReplicatedDataFun(Samples = SamplesItSingle,
                             n = length(ZMatIt))

### Test of PPC's ##########
set.seed(42)
YModel <- rnorm(n = 100,mean = 5,sd = 2)

NimMod <- nimbleCode({
  mu ~ dnorm(mean = 5,sd = 10)
  sd ~ dgamma(1,1)
  for(i in 1:N){
    y[i] ~ dnorm(mean=mu, sd = sd)
  }
})

Samples <- nimbleMCMC(code=NimMod, constants = list("N"=100),
                      data=list("y"=YModel), niter = 5000,nburnin = 2000)

Y_Model_Rep <- matrix(data = 0,
                      nrow = 3000,
                      ncol = 100)

for(s in 1:3000){
  Y_Model_Rep[s,] <- rnorm(n = 100, 
                     mean = Samples[s,1],
                     sd = Samples[s,2])
}
dim(Y_Model_Rep)

ChiSq <- sum((YModel - mean(Samples[,1]))^2/
               mean(Samples[,2]))

## ChiSqTestValue
ChiSq_Rep <- numeric(3000)
for(s in 1:3000){
  ChiSq_Rep[s] <- sum((Y_Model_Rep[s,] - Samples[s,1])^2/
                        Samples[s,2])
}
hist(ChiSq_Rep,breaks = 100,probability = TRUE)
abline(v = ChiSq)




## Alternative, write Own Function for posterior predictive checks (without Nimble)##
ppParameters_OwnMod_OC <- function(Samples, NAges, NYear){
  NMat <- Samples[,grep("N_t", colnames(Samples))]
  KMat <- Samples[,grep("k", colnames(Samples))]
  betaMat <- Samples[,grep("beta", colnames(Samples))[1:11]]
  betaJumpMat <- Samples[,grep("betaJump", colnames(Samples))]
  
  aSamps <- Samples[,match("a",colnames(Samples))]
  pSamps <- Samples[,match("p",colnames(Samples))]
  
  muYSamps <- Samples[,match("muY",colnames(Samples))]
  sdYSamps <- Samples[,match("sdY",colnames(Samples))]
  
  sigma_epsSamps <- Samples[,match("sigma_eps",colnames(Samples))]
  
  muArray <- array(data = 0, 
                   dim = c(NAges,NYear,nrow(Samples)))
  
  sdArray <- array(data = 0, 
                   dim = c(NAges,NYear,nrow(Samples))) 
  
  for(s in 1:nrow(Samples)){ #for each draw
    #1. Define R Matrix
    RMat_s <- RTildeMat(Ntime = ncol(KMat),
                        a = aSamps[s])
    
    
    #2. Loop over time and age
    for(x in 1:NAges){
      for(t in 1:NYear){
        muArray[x,t,s] <- #alpha[x]+
          betaMat[s,x]*KMat[s,t]+
          betaJumpMat[s,x]*muYSamps[s]*(
            inprod(RMat_s[t+1,],NMat[s,])- #J_t
              inprod(RMat_s[t,],  NMat[s,])    #J_t-1
          )
        
        
        CovVec <- CovZMat(N_t = NMat[s,],
                          a =  aSamps[s],
                          t =  t,
                          sigma2 = pow(sdYSamps[s],2))
        
        sdArray[x,t,s] <- pow(sigma_epsSamps[s],2)  +
          pow(betaJumpMat[s,x]*sdYSamps[s]*
                inprod(RMat_s[t+1,],
                       NMat[s,]),2)+
          pow(betaJumpMat[s,x]*sdYSamps[s]*
                inprod(RMat_s[t,],
                       NMat[s,]),2) -  #Var second part
          2*pow(betaJumpMat[s,x],2)*CovVec #-Covariance
      }
    }
  }
  return(list(muArr = muArray,
              sdArr = sdArray))
}



# #Employ posterior predictive checking 
# dataNodes <- OwnMod_US$getNodeNames(dataOnly = TRUE)
# parentNodes <- OwnMod_US$getParents(dataNodes, stochOnly = TRUE)  # `getParents` is new in nimble 0.11.0
# 
# simNodes <- OwnMod_US$getDependencies(parentNodes, self = FALSE)  
# 
# 
# ppSamples <- matrix(0, nrow = 2000, 
#                     ncol = length(OwnMod_US$expandNodeNames(dataNodes, returnScalarComponents = TRUE)))
# 
# ## Determine ordering of variables in `mvSamples` modelValues and therefore in `samples`
# vars <- comOwnMod_US$bOwnMod_US$mvSamples$getVarNames()
# 
# ## Quick check of variable ordering
# colnames(SamplesOwnNew_US$chain1)
# 
# 
# set.seed(1)
# system.time({
#   for(i in seq_len(2000)) {
#     values(comOwnMod_US$OwnMod_US,vars) <- SamplesOwnNew_US$chain1[i, ] # assign 'flattened' values
#     comOwnMod_US$OwnMod_US$simulate(simNodes, includeData = TRUE)
#     ppSamples[i, ] <- values(comOwnMod_US$OwnMod_US, dataNodes)
#   }
# })
# 
# 
# ### Posterior Predictive Sampling in 1 Function
# ppSampler <- nimbleFunction(
#   setup = function(model, mcmc) {
#     dataNodes <- model$getNodeNames(dataOnly = TRUE)
#     parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
#     cat("Stochastic parents of data are: ", paste(parentNodes, sep = ','), ".\n")
#     simNodes <- model$getDependencies(parentNodes, self = FALSE)
#     vars <- mcmc$mvSamples$getVarNames()  # need ordering of variables in mvSamples / samples matrix
#     cat("Using posterior samples of: ", paste(vars, sep = ','), ".\n")
#     nData <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
#   },
#   run = function(samples = double(2)) {
#     niter <- dim(samples)[1]
#     ppSamples <- matrix(nrow = niter, ncol = nData)   
#     for(i in 1:niter) {
#       values(model, vars) <<- samples[i, ]
#       model$simulate(simNodes, includeData = TRUE)
#       ppSamples[i, ] <- values(model, dataNodes)
#     }
#     return(ppSamples)
#     returnType(double(2))
# })
# 
# ppSampler_basic <- ppSampler(model = OwnMod_US, 
#                              mcmc = bOwnMod_US)
# 
# cppSampler_basic <- compileNimble(ppSampler_basic, 
#                                   project = OwnMod_US)
# set.seed(1)
# ppSamples_via_nf <- cppSampler_basic$run(SamplesOwnNew_US$chain1)
# 
