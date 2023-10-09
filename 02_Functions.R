########## FUNCTIONS ###########################################################
library(rstan);library(tidyverse)

#### FUNCTION FOR CREATION OF SUMMARY OUTPUT FROM NIMBLE SAMPLES #############
SummaryOutput <- function(MCMCSampler, params="all", quantiles=c(0.1,0.9)){
  
  chain_output <- MCMCSampler
  if(is.matrix(chain_output)){ #n_chains is equal to one
    n_chains <- 1
    n_samples <- nrow(chain_output)
    n_params <- ncol(chain_output)
    names_params <- colnames(chain_output) #get Parameter names
    chain_output <- list(chain_output) #transform chain output into list
    
  } else { #multiple chains
    n_chains   <- length(chain_output) #get Number of Chains
    n_samples  <- nrow(chain_output[[1]]) #get Number of samples
    n_params   <- ncol(chain_output[[1]]) #get Number of Parameters
    names_params <- colnames(chain_output[[1]]) #get Parameternames
  }
  
  
  #Transform into Format for Stan functions
  MR <- reorganise_mortality(chain_output = chain_output, 
                             params_vec =  1:n_params, 
                             n_samples = n_samples,
                             n_chains = n_chains,
                             names_params = names_params)
  
  if(params[1]=="all"){ #check if 
    Summary <- sapply(MR, function(x){c(mean(x),
                                        sd(x),
                                        quantile(x, probs=quantiles, na.rm = TRUE),
                                        rstan::Rhat(x),
                                        rstan::ess_bulk(x),
                                        rstan::ess_tail(x))}) %>% t() %>%  
      as_tibble(., .name_repair = "unique") %>% 
      rename_with(.fn = function(x) {c("mean", "sd",quantiles, "Rhat", "bulk_ess", "tail_ess")}, 
                  .cols = everything()) %>% 
      mutate("Param"=names_params, .before=1)
  } else {
    if(length(params)==1){
      Ind <- grep(paste0("\\<",params), names(MR))
    } else { #Parameters have more than one length
      #Match each parameter on its own, to avoid situations where a single letter i.e. a is matched to multiple words
      Ind <- sapply(params, 
                    function(y){grep(paste0("\\<",y), #match beginning of word
                            x = names(MR))}) %>% unlist() %>% unique() #only unique names
    }
    Summary <- sapply(MR[Ind], function(x){c(mean(x),
                                             sd(x),
                                             quantile(x, probs=quantiles, na.rm = TRUE),
                                             rstan::Rhat(x),
                                             rstan::ess_bulk(x),
                                             rstan::ess_tail(x))}) %>% t() %>%  
      as_tibble(., .name_repair = "unique") %>% 
      rename_with(.fn = function(x) {c("mean", "sd",quantiles, "Rhat", "bulk_ess", "tail_ess")}, 
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
  # Age Groups of 10. Starting from 0-9, 10-19,..., 90+
    HelperAgeLab <- data.frame("AgeOld" = unique(AgeLabels)[-1],
                               "AgeNew" = c(rep(1:10,each=10),10))
  return(HelperAgeLab)
}

#Function for EU22 Age Groups
AgeLabFun_EU22 <- function(){
    # Age Groups of 10. Starting from 0-9, 10-19,..., 90+
    val <- c(rep(1:9, each=2),10,10)
  return(val)
}

AgeLabFun_HMD <- function(Type = 1){
  if(Type == 1){
    # Age Groups of 10. Starting from 0-9, 10-19,..., 90+
    HelperHMD <- data.frame("AgeOld"=1:24,
                            "AgeNew"=c(1,1,1,rep(2:10,each=2),rep(10,3)))
  }else { #Only for UK Data 
    #Age Groups of <1, 1-4, 5-15, 15-25,...
    HelperHMD <- data.frame("AgeOld"=1:24,
                            "AgeNew"=c(1,2,rep(3:12,each=2),13,13))
                            "AgeNew"=c(1,1,rep(2:18),rep(19,5))
  }
  return(HelperHMD)
}

##### FUNCTION FOR CALCULATION OF Likelihood MATRIX (FOR WAIC/LOOCV)############ 
LikelihoodMatrixFun <- function(Samples, n,ZMat){
  
  muPos <- grep("mu[",colnames(Samples), 
                fixed = TRUE)
  
  sdPos <-grep("sigma_eps", colnames(Samples), 
                 fixed = TRUE) 
  
  LikeMat <- matrix(0,
                    nrow = nrow(Samples),
                    ncol = n)
  
  ZMatVec <- as.vector(ZMat)
  for(s in 1:nrow(Samples)){
      LikeMat[s,] <- dnorm(ZMatVec, 
                           mean = as.numeric(Samples[s,muPos]),
                           sd = Samples[s, sdPos])
    } 
  
  return(LikeMat)
}

######### FUNCTION TO CREATE FORECASTS FROM MORTALITY IMPROVEMENTRATES #########
FutureZ <- function(Samples,S ,H ,OwnMod = TRUE,NAge){
  
  # @H: Forecast of H time periods ahead
  # @S: Amount of draws from the posterior
  # @NAge: Amount of age groups in the data
  # @Samples: Posterior Samples from nimble
  # @OwnMod: Binary Variable to choose either Liu,Li or Own model
  
  #Get Random indicies from posterior. All should have the same draw 
  PostDraw <- sample(1:nrow(Samples), 
                     size = S, 
                     replace = FALSE)
  
  #1) Generation of future kt's
  driftVec <- Samples[PostDraw, #Select S posterior draws
                      grep("drift",colnames(Samples))] #extract drift
  
  sigmaTime <- Samples[PostDraw,
                       grep("sigma_time",colnames(Samples))]
  
  # generation of future kt's, for each posterior draw, generate H future kt's
  #Row's are posterior draws, columns are time periods
  FutureKtMat <- sapply(1:H, function(x){
    rnorm(n = S, mean = driftVec, sd = sigmaTime)})
  
    #2) Generation of future Jt's 
    
    #2.1) First Generation of future N'ts
    # Note that Z1 = J2 - J1, thus ZT = JT+1 - JT, hence we have T+1 J's and N's
    pVec <- Samples[PostDraw, #Select S posterior draws
                    grep(paste0("\\<","p"), #starts with p to select "p" only 
                         colnames(Samples))] 
    
    #Row's are posterior draws, columns are time periods,
    #for each posterior draw, generate new N'ts
    FutureNtMat <- sapply(1:H, function(x)rbinom(n = S,size = 1, prob = pVec))
    
    #Generation of R Matrix 
    NtMatrix <- Samples[PostDraw, #Select S posterior draws
                        grep(paste0("\\<","N"), #starts with N to select "N_t" only 
                             colnames(Samples))] 
    
    
    #Total N't Matrix (observed and future values of N_t)
    NtTotMat <- cbind(NtMatrix,FutureNtMat)
    
    #2.2) Generate new Values of J't
    muYVec <- Samples[PostDraw, #Select S posterior draws
                      grep(paste0("\\<","muY"), #select muY 
                           colnames(Samples))]
    sdYVec <- Samples[PostDraw, #Select S posterior draws
                      grep(paste0("\\<","sdY"), #select sdY 
                           colnames(Samples))]
    if(OwnMod != TRUE){
      aVec <- numeric(S) #aVec filled with zeros
    } else {
      aVec <- Samples[PostDraw, #Select S posterior draws
                      grep(paste0("\\<","a"), #starts with a to select "a" only 
                           colnames(Samples))]
    }
    
    #Creation of mean and sd of future J's
    #one more column, since Zx,T = JT+1 - JT. Thus first column is JT, second JT+1 and so on
    #First column is last "observed" JT
    JtMat <- Samples[PostDraw, #Select S posterior draws
                     grep(paste0("\\<","J"), #starts with N to select "N_t" only 
                          colnames(Samples))]
    
    
    FutureYt <- sapply(1:H, function(x){ 
      rnorm(n = S, mean = muYVec, sd = sdYVec)})
    
    FutureJt <- matrix(data = 0,nrow = S,ncol = H+1) #one column more
    FutureJt[,1] <- JtMat[,ncol(JtMat)]
    
    for(j in 2:(H+1)){ #Var(Jt)=Rt*N*sigma*N*Rt
      h <- j-1  
      FutureJt[,j] <- aVec*FutureJt[,j-1]+FutureNtMat[,h]*FutureYt[,h]
    }
    
    #3.) Plug all together to generate new Zx,t+h
    betaMat <- Samples[PostDraw, #Select S posterior draws
                       grep("beta", #starts with p to select "p" only 
                            colnames(Samples))][,1:NAge]
    
    betaJumpMat <- Samples[PostDraw, #Select S posterior draws
                           grep(paste0("\\<","betaJump"), #starts with p to select "p" only 
                                colnames(Samples))]
    
    sigma_epsVec <- Samples[PostDraw, #Select S posterior draws
                            grep(paste0("\\<","sigma_eps"), #starts with p to select "p" only 
                                 colnames(Samples))]
    
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
    
  return(FutureZArray)
}
