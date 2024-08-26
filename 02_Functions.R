########## FUNCTIONS ###########################################################
pacman::p_load("rstan","tidyverse","bayestestR","MCMCvis","truncnorm")
library(rstan);library(tidyverse);library(bayestestR)

#### FUNCTION FOR CREATION OF SUMMARY OUTPUT FROM NIMBLE SAMPLES #############
SummaryOutput <- function(MCMCSampler, params="all", quantiles=c(0.1,0.9)){
  
  if(is.matrix(MCMCSampler)){ #n_chains is equal to one
    n_chains <- 1
    n_samples <- nrow(MCMCSampler)
    n_params <- ncol(MCMCSampler)
    names_params <- colnames(MCMCSampler) #get Parameter names
    chain_output <- list(MCMCSampler) #transform chain output into list
    
  } else { #multiple chains
    n_chains   <- length(MCMCSampler) #get Number of Chains
    n_samples  <- nrow(MCMCSampler[[1]]) #get Number of samples
    n_params   <- ncol(MCMCSampler[[1]]) #get Number of Parameters
    names_params <- colnames(MCMCSampler[[1]]) #get Parameternames
  }
  
  
  #Transform into Format for Stan functions
  MR <- reorganise_mortality(chain_output = MCMCSampler, 
                             params_vec =  1:n_params, 
                             n_samples = n_samples,
                             n_chains = n_chains,
                             names_params = names_params)
  
  #Get names of all parameters without brakets
  names <- vapply(strsplit(names(MR), 
                           split = "[", 
                           fixed = TRUE), 
                  `[`, 1, FUN.VALUE = character(1))
  
  if(params[1]=="all"){ #check if 
    Summary <- sapply(MR, function(x){c(mean(x),
                                        bayestestR::map_estimate(x),
                                        sd(x),
                                        quantile(x, probs=quantiles, na.rm = TRUE),
                                        rstan::Rhat(x),
                                        rstan::ess_bulk(x),
                                        rstan::ess_tail(x))}) %>% t() %>%  
      as_tibble(., .name_repair = "unique") %>% 
      rename_with(.fn = function(x) {c("mean","MAP", "sd",quantiles, "Rhat", "bulk_ess", "tail_ess")}, 
                  .cols = everything()) %>% 
      mutate("Param"=names_params, .before=1)
  } else {
    if(length(params)==1){
      Ind <- grep(paste0("\\<",params), names)
    } else { #Parameters have more than one length
      Ind <- sapply(params, function(y) which(names %in% y)) %>% unlist()
    }
    Summary <- sapply(MR[Ind], function(x){c(mean(x),
                                             #bayestestR::map_estimate(as.vector(x))$MAP_Estimate,
                                             sd(x),
                                             quantile(x, probs=quantiles, na.rm = TRUE),
                                             rstan::Rhat(x),
                                             rstan::ess_bulk(x),
                                             rstan::ess_tail(x))}) %>% t() %>%  
      as_tibble(., .name_repair = "unique") %>% 
      rename_with(.fn = function(x) {c("mean",
                                       #"MAP", 
                                       "sd",quantiles, "Rhat", "bulk_ess", "tail_ess")}, 
                  .cols = everything()) %>% 
      mutate("Param"=names(MR)[Ind], .before=1)
  }
  
  return(Summary)
  
}

#Helper Function for Summary Like Output (taken from Theo Rashid)
reorganise_mortality <- function(
    chain_output, params_vec, n_samples, n_chains, names_params
) {
  
  mr <- list() #empty list
  for (i in params_vec) {
    par_mat <- matrix(data = 0, nrow = n_samples, ncol = n_chains)
    for (c in 1:n_chains) {
      par_mat[, c] <- chain_output[[c]][, i]
    }
    name <- names_params[i]
    mr[[name]] <- par_mat
  }
  return(mr)
}

####### HELPER FUNCTIONS FOR DATA MANIPULATION #################################
AgeLabFun_EU <- function(AgeLabels){
    #Age Groups of 10. Starting from <5,5-14,15-25,...,85+ 
    HelperAgeLab <- data.frame("AgeOld" = unique(AgeLabels)[-1],
                               "AgeNew" = c(rep(1,5),rep(2:10,each=10),rep(10,6)))
  return(HelperAgeLab)
}

#Function for EU22 Age Groups
AgeLabFun_EU23 <- function(){
  #Age Groups of 10. Starting from <5,5-14,15-25,...,85+ 
  val <- val <-  c(1,rep(2:9, each = 2), rep(10,2))
  return(val)
}

AgeLabFun_HMD <- function(){
    #Age Groups of 10. Starting from <5,5-14,15-25,...,85+ 
    HelperHMD <- data.frame("AgeOld"=1:24,
                            "AgeNew"=c(1,1,rep(2:9,each=2),rep(10,6)))
  return(HelperHMD)
}

##### FUNCTION FOR CALCULATION OF Likelihood MATRIX (FOR WAIC/LOOCV)############ 
LikelihoodMatrixFun <- function(Samples, n, ZMat){
  
  muSamp <- MCMCvis::MCMCchains(Samples, params = "mu")
  sdSamp <- MCMCvis::MCMCchains(Samples, params = "sigma_eps")
  
  LikeMat <- matrix(0,
                    nrow = nrow(muSamp),
                    ncol = n)
  
  ZMatVec <- as.vector(ZMat) #Transfrom Matrix in Vector (row wise )
  
  #if(SingleVar == TRUE){
  for(s in 1:nrow(muSamp)){
    LikeMat[s,] <- dnorm(ZMatVec, 
                         mean = muSamp[s,],
                         sd = sdSamp[s])
  } 
  return(LikeMat)
}

LikelihoodMatrixFunJoint <- function(Samples, n, ZArray, NCountry){
  
  muSamp <- MCMCvis::MCMCchains(Samples, params = "mu")
  sdSamp <- MCMCvis::MCMCchains(Samples, params = "sigma_eps")
  
  LikeMat <- matrix(0,
                    nrow = nrow(muSamp),
                    ncol = n)
  
  #Helper Values For Loop
  NObsCount <- ncol(muSamp)/NCountry #number of observations by country
  DevPointsLower <- seq(1, ncol(muSamp), NObsCount) #Find n's beloning to each country
  DevPointsUpper <- seq(NObsCount, ncol(muSamp), NObsCount)
  
  for(c in 1:NCountry){
    ZMatVec <- as.vector(ZArray[,,c]) #Transfrom Matrix in Vector (row wise )
    
    for(s in 1:nrow(muSamp)){
      LikeMat[s,DevPointsLower[c]:
                DevPointsUpper[c]] <- dnorm(ZMatVec, 
                                            mean = muSamp[s,DevPointsLower[c]:
                                                            DevPointsUpper[c]],
                                            sd = sdSamp[s,c])
    } 
    
  }
  return(LikeMat)
}


######### FUNCTION TO CREATE FORECASTS FROM MORTALITY IMPROVEMENTRATES #########
FutureZ <- function(Samples, S, H, Mod = "AR", LC = FALSE){
  
  # @H: Forecast of H time periods ahead
  # @S: Amount of draws from the posterior
  # @NAge: Amount of age groups in the data
  # @Samples: Posterior Samples from nimble
  # @Mod: Choose Either AR, MA or Liu
  # @LC: Binary if LC model should be forecasted
  
  if(!Mod %in% c("AR","MA","Liu")){
    print("Mod must either be AR, MA or Liu")
    stop()
  }
  
  #Get Random indicies from posterior. All should have the same draw 
  PostDraw <- sample(1:nrow(MCMCchains(Samples)), 
                     size = S, 
                     replace = TRUE)
  
  #1) Generation of future kt's
  driftVec <- MCMCchains(Samples, params = "drift")[PostDraw,]
  
  sigmaTime <- MCMCchains(Samples, params = "sigma_time")[PostDraw,]
  
  # generation of future kt's, for each posterior draw, generate H future kt's
  #Row's are posterior draws, columns are time periods
  FutureKtMat <- sapply(1:H, function(x){
    rnorm(n = S, mean = driftVec, sd = sigmaTime)})
  
  if(LC == FALSE){ ## Lee Carter Model Forecasts
    #2) Generation of future Jt's (be very careful about indices..)
    
    #2.1) First Generation of future N'ts
    # Note that Z1 = J2 - J1, thus ZT = JT+1 - JT, hence we have T+1 J's and N's
    pVec <- MCMCchains(Samples, params = "p")[PostDraw,]
    
    #Row's are posterior draws, columns are time periods,
    #for each posterior draw, generate new N'ts
    FutureNtMat <- sapply(1:H, function(x)rbinom(n = S,size = 1, prob = pVec))
    
    #Generation of N_t and Y_t Matrix 
    NtMatrix <- MCMCchains(Samples, params = "N_t")[PostDraw,]
    YtMatrix <- MCMCchains(Samples, params = "Y_t")[PostDraw,]
    
    #Total N't (observed and future values of N_t)
    NtTotMat <- cbind(NtMatrix,FutureNtMat)
    
    
    #2.2) Generate new Values of J't
    muYVec <- MCMCchains(Samples, params = "muY")[PostDraw,]
    sdYVec <- MCMCchains(Samples, params = "sdY")[PostDraw,]
    
    JtMat <- MCMCchains(Samples, params = "J")[PostDraw,]
    
    #Creation of mean and sd of future J's
    #one more column, since Zx,T = JT+1 - JT. Thus first column is JT, second JT+1 and so on
    #First column is last "observed" JT
    
    FutureYt <- sapply(1:H, function(x){ 
      truncnorm::rtruncnorm(n = S, mean = muYVec, sd = sdYVec, a = 0)})
    
    #Total N't (observed and future values of N_t)
    YtTotMat <- cbind(YtMatrix, FutureYt)
    
    FutureJt <- matrix(data = 0,nrow = S,ncol = H+1) #one column more
    FutureJt[,1] <- JtMat[,ncol(JtMat)] 
    
    if(Mod =="AR"){
      aVec <- MCMCchains(Samples, params = "a")[PostDraw,]
      for(j in 2:(H+1)){ #Var(Jt)=Rt*N*sigma*N*Rt
        h <- j-1  
        FutureJt[,j] <- aVec*FutureJt[,j-1]+FutureNtMat[,h]*FutureYt[,h]
      }
    } else if(Mod =="MA"){
      bVec <- MCMCchains(Samples, params = "b")[PostDraw,]
      TMax <- ncol(NtMatrix)
      for(j in 2:(H+1)){ 
        h <- j-1  
        FutureJt[,j] <- bVec*NtTotMat[,(TMax+h-1)]*YtTotMat[,(TMax+h-1)]+FutureNtMat[,h]*FutureYt[,h]
      }
    } else { #Liu Li Model
      aVec <- 0 
      for(j in 2:(H+1)){ #Var(Jt)=Rt*N*sigma*N*Rt
        h <- j-1  
        FutureJt[,j] <- aVec*FutureJt[,j-1]+FutureNtMat[,h]*FutureYt[,h]
      }
    }
    
    
    #3.) Plug all together to generate new Zx,t+h
    betaMat <- MCMCchains(Samples, params = "beta")[PostDraw,]
    betaJumpMat <- MCMCchains(Samples, params = "betaJump")[PostDraw,]
    
    sigma_epsVec <- MCMCchains(Samples, params = "sigma_eps")[PostDraw,]
    
    NAge <- ncol(betaMat)
    FutureZArray <- array(data = 0,dim = c(NAge,H,S),
                          dimnames = list("Age"=1:NAge, "Time"=1:H, "It"=1:S))
    
    for(s in 1:S){
      #via Matrix multiplication (outer product)
      FutureZArray[,,s] <- 
        betaMat[s,]%*%t(FutureKtMat[s,])+ #Normal effect
        betaJumpMat[s,]%*%t(FutureJt[s,2:(H+1)]) - #JT+1
        betaJumpMat[s,]%*%t(FutureJt[s,1:H]) +#JT
        sapply(1:H, function(x){rnorm(n = NAge, 
                                      mean = 0, sd = sigma_epsVec[s])}) #error Term
    }
    
    ReturnList <- list("Rates"=FutureZArray,
                       "Jt"=FutureJt,
                       "betaJump"=betaJumpMat)
    
  } else {
    
    betaMat <- MCMCchains(Samples, params = "beta")[PostDraw,]
    
    sigma_epsVec <- MCMCchains(Samples, params = "sigma_eps")[PostDraw,]
    
    NAge <- ncol(betaMat)
    
    FutureZArray <- array(data = 0,dim = c(NAge,H,S),
                          dimnames = list("Age"=1:NAge, "Time"=1:H, "It"=1:S))
    for(s in 1:S){
      #via Matrix multiplication (outer product)
      FutureZArray[,,s] <- 
        betaMat[s,]%*%t(FutureKtMat[s,])+ #Normal effect
        sapply(1:H, function(x){rnorm(n = NAge, 
                                      mean = 0, sd = sigma_epsVec[s])}) #error Term
    }
    
    ReturnList <- list("Rates"=FutureZArray)
  }
  
  return(ReturnList)
}
