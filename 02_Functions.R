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
  
  muPos <- grep("mu",colnames(Samples), 
                fixed = TRUE)
  
  sdPos <-grep("sigma_eps", colnames(Samples), 
                 fixed = TRUE) 

  muPosReal <- muPos[1:(length(muPos)-1)] #without muY
  
  LikeMat <- matrix(0,
                    nrow = nrow(Samples),
                    ncol = n)
  
  ZMatVec <- as.vector(ZMat)
  for(s in 1:nrow(Samples)){
      LikeMat[s,] <- dnorm(ZMatVec, 
                           mean = as.numeric(Samples[s,muPosReal]),
                           sd = Samples[s, sdPos])
    } 
  
  return(LikeMat)
}