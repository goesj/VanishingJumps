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
LikelihoodMatrixFun <- function(Samples, n,ZMat, SingleVar = FALSE){
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
  
  LikeMat <- matrix(0,
                    nrow = nrow(Samples),
                    ncol = n)
  
  ZMatVec <- as.vector(ZMat)
  if(SingleVar == TRUE){
    for(s in 1:nrow(Samples)){
      LikeMat[s,] <- dnorm(ZMatVec, 
                           mean = as.numeric(Samples[s,muPosReal]),
                           sd = Samples[s, sdPos])
    } 
  } else {
    for(s in 1:nrow(Samples)){
      LikeMat[s,] <- dnorm(ZMatVec, 
                           mean = as.numeric(Samples[s,muPosReal]),
                           sd = as.numeric(sqrt(Samples[s,sdPos])))
    }
  }
  
  return(LikeMat)
}

########### NIMBLE FUNCTIONS ###################################################
CovZMat <- nimbleFunction(
  run=function(N_t=double(1), #vector 
               a = double(0), #scalar
               t = double(0), #scalar
               sigma2=double(0)){ #scalar
    returnType(double(0)) #scalar
    #Must be written using a loop
    Cov <- nimNumeric(length=t)
    for(j in 1:t){ #calculate according to formula
      Cov[j] <- N_t[t+1-j]*pow_int(a,(2*j)-1)*sigma2 #scalar*scalar*scalar
    }
    CovSum <- sum(Cov) #sum over all indicies
    return(CovSum)
  }, buildDerivs = list(run = list(ignore = 'j'))
)

##########EVERYTHING BELOW NOT IN USE ! ########################################

### Age Lab Functions Including other types of Age grouping
AgeLabFun_EU_MoreAgeGroups <- function(Type = 1,
                         AgeLabels){
  # Age Groups of 10. Starting from 0-9, 10-19,..., 90+
  if(Type == 1){
    HelperAgeLab <- data.frame("AgeOld" = unique(AgeLabels)[-1],
                               "AgeNew" = c(rep(1:10,each=10),10))
  }else if(Type == 2){ 
    #Age Groups of 10. Starting from <5,5-14,15-25,...,85+ 
    HelperAgeLab <- data.frame("AgeOld" = unique(AgeLabels)[-1],
                               "AgeNew" = c(rep(1,5),rep(2:10,each=10),rep(10,6)))
    
  }else{ #Age Groups of 5. Starting from 0-4, 5-9,..., 90
    HelperAgeLab <- data.frame("AgeOld" = unique(AgeLabels)[-1],
                               "AgeNew" = c(rep(1:19,each=5),rep(19,6)))
  }
  return(HelperAgeLab)
}

#Function for EU22 Age Groups
AgeLabFun_EU22_MoreAgeGroups <- function(Type=1){
  if(Type==1){
    # Age Groups of 10. Starting from 0-9, 10-19,..., 90+
    val <- c(rep(1:9, each=2),10,10)
  } else if(Type == 2){
    #Age Groups of 10. Starting from <5,5-14,15-25,...,85+ 
    val <-  c(1,rep(2:9, each = 2), rep(10,3))
  } else {
    #Age Groups of 5. Starting from 0-4, 5-9,..., 90+ 
    val <-c(1:19,19) 
  }
  return(val)
}

AgeLabFun_HMD_MoreAgeGroups <- function(Type = 1){
  if(Type == 1){
    # Age Groups of 10. Starting from 0-9, 10-19,..., 90+
    HelperHMD <- data.frame("AgeOld"=1:24,
                            "AgeNew"=c(1,1,1,rep(2:10,each=2),rep(10,3)))
  }else if(Type == 2){ 
    #Age Groups of 10. Starting from <5,5-14,15-25,...,85+ 
    HelperHMD <- data.frame("AgeOld"=1:24,
                            "AgeNew"=c(1,1,rep(2:9,each=2),rep(10,6)))
    
  }else {
    # #Age Groups of 5. Starting from 0-4, 5-9,..., 90+ 
    HelperHMD <- data.frame("AgeOld"=1:24,
                            "AgeNew"=c(1,1,rep(2:18),rep(19,5)))
  }
  return(HelperHMD)
}


########### Posterior Predictive Checks Functions###############################
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
ZMatPlot <- function(ZMat, Quants = c(0.1,0.9), Z_Rep, AgesUsed = 1:10, facetWrap=FALSE){
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
  
  
  if(facetWrap==TRUE){
    Plot <- LowerPIZ %>% 
      mutate("HPI"=HigherPIZ$HPI,
             "Mean"=HigherPIZ$Mean) %>% 
      filter(Age %in% AgesUsed) %>% 
      ggplot(aes(x=Time, group=Age))+
      geom_ribbon(aes(ymin=LPI, ymax=HPI, fill=Age))+
      geom_line(aes(y = Mean),col="red")+
      geom_line(data=filter(ZMaLong, Age %in% AgesUsed), 
                aes(x=Time, y=Val, group=Age),
                col="black")+
      facet_wrap(~Age)
  } else {
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
  }
  return(Plot)
}


########## NIMBLE FUNCTIONS ####################################################
SumToZero <- nimbleFunction( #Function to Sum Vector to Zero
  run = function(x = double(1)){ #input Vector
    returnType(double(1))
    xtilde <- x - mean(x) #vectorized version
    return(xtilde)
  },
  buildDerivs = "run"
)

SumToOneNorm <- nimbleFunction( #Function to Sum Vector to Zero
  run = function(x = double(1)){ #input Vector
    returnType(double(1)) #output Vector
    sum_x <- sum(x) #scalar
    xtilde <- x/sum_x #vectorized version
    return(xtilde)
  },
  buildDerivs = "run"
)

gramschmidt_Nimble <- nimble::nimbleFunction(
  run = function(x = double(1), #vector
                 y = double(1)) #vector
  {
    returnType(double(1)) #vector as return  
    
    #Check if only zeros in Model
    # allzero <- all(y==0)
    # if(allzero==TRUE)
    # ynew <- numeric(length(y)) #zero vector with length of y
    # else
      #Step 1: Find non-zeros
      NonZeroInd <- which(y!=0) #vector
      NZero <- length(y)-length(NonZeroInd) #amout of zeros
      
      #Step 2: Do Gram Schmidt w/o the Zero Inds
      v2 <- numeric(length(NonZeroInd),init=FALSE)
      
      v1 <- x[NonZeroInd]
      v2 <- y[NonZeroInd] - inprod(v1,y[NonZeroInd])/inprod(v1,v1)*v1 #Gram Schmidt Orthogonalization
      
      #Step 3: Place zeros in correct position
      ynew <- numeric(length(y))
      ynew[NonZeroInd] <- v2 #overwrite zeros with new vector in positions of non-zeros
      
    
    
    return(ynew) #only return ynew, as x stays unchanged
  },
  buildDerivs = "run"
)
#Function to save some nodes for faster compilation and Building
calculateLambda <- nimbleFunction(
  run = function(eta = double(0), #scalar
                 Offset = double(0)){ #scalar
    returnType(double(0)) #scalar
    lambda <-  exp(eta)*Offset
    return(lambda)
  },
  buildDerivs = "run"
)


gramschmidt_Nimble2 <- nimbleFunction(
  run = function(x = double(1), #vector
                 y = double(1)) #vector
  {
    returnType(double(1)) #vector as return  
      ynew <- y - (inprod(x,y)/inprod(x,x)*x) #Gram Schmidt Orthogonalization (vector - scalar * vector = vector)
      
      return(ynew) #only return ynew, as x stays unchanged
  },
  buildDerivs = "run"
)

### QR composition in NIMBLE ############
qr_Nimble <- nimbleFunction(
  run=function(A = double(2)) #matrix
  {
    returnType(double(2))
    n <- dim(A)[1] #Number of rows 
    m <- dim(A)[2] # number of Pazrameters
    Q <- matrix(0, nrow=n, ncol=m)
    R <- matrix(0, nrow=m, ncol=m)
    for (j in 1:m) { # seq along columns
      v <- A[,j] # value (value of current column)
      if(j == 1){
        R[j,j] <- sqrt(inprod(v,v)) #Norm
        Q[,j] <- v / R[j,j]
      } else{
        for (i in 1:(j-1)) {
          R[i,j] <- inprod(Q[,i],v)  #returns a scalar
          v <- v - (R[i,j] * Q[,i]) #scalar times vector
        }
        R[j,j] <- sqrt(inprod(v,v)) #Norm
        Q[,j] <- v / R[j,j]  
      }
    }
    return(Q)
  }, buildDerivs = list(run = list(ignore = c('i','j')))
)

CovMat <- nimbleFunction(
  run = function(A = double(0), 
                 sigma = double(0)) {
    returnType(double(2))
    IMat <- diag(A-1)
    JMat <- matrix(value=1, nrow = (A-1), ncol=(A-1))
    CovMat <- sigma*(IMat - (1/A)*JMat)
    return(CovMat)
  },
  buildDerivs = list(run = list(ignore = "A"))
)


#Adjusted QR with Resulting Matrix having length 1 instead of Norm 1
qr_Nimble_Adj <- nimbleFunction(
  run=function(A = double(2)) #matrix
  {
    returnType(double(2))
    n <- dim(A)[1] #Number of rows 
    m <- dim(A)[2] # number of Pazrameters
    Q <- matrix(0, nrow=n, ncol=m)
    R <- matrix(0, nrow=m, ncol=m)
    for (j in 1:m) { # seq along columns
      v <- A[,j] # value (value of current column)
      if(j == 1){
        R[j,j] <- sqrt(inprod(v,v)) #Norm
        Q[,j] <- v / R[j,j]
      } else{
        for (i in 1:(j-1)) {
          R[i,j] <- inprod(Q[,i],v)  #returns a scalar
          v <- v - (R[i,j] * Q[,i]) #scalar times vector
        }
        R[j,j] <- sqrt(inprod(v,v)) #Norm
        Q[,j] <- v / R[j,j]  
      }
    }
    CVec <- numeric(m)
    for(k in 1:m){
      CVec[k] <- 1/sum(Q[,k])
    }
    QTilde <- Q%*%diag(CVec) #Scale Vectors of Q that they have length 1 
    return(QTilde)
  }, buildDerivs = "run"
)

RTildeMat <- nimbleFunction(
  run = function(Ntime = double(0),
                 a = double(0)){
    returnType(double(2)) # matrix
    #Create Matrix, so that J = RN with N = (N1,...,NT+1)
    # Note that Z1 = J2 - J1, thus ZT = JT+1 - JT
    RMat <- matrix(0, nrow=(Ntime+1), ncol=(Ntime+1))
    for(t in 1:(Ntime+1)){ #over all rows
      for(j in 1:t){ #over columns till t
        RMat[t,j] <- pow_int(a,t-j)
      } 
    }
    return(RMat)
  }, buildDerivs = list(run = list(ignore = c('t','j')))
)

### Creation of R Matrix ###############################
RMatrix <- nimbleFunction(
  run = function(Ntime = double(0),
                 a = double(0)){
    returnType(double(2)) # matrix
    
    RMat <- matrix(0, nrow=Ntime, ncol=Ntime)
    for(t in 1:Ntime){ #over all rows
      for(j in 1:t){ #over columns till t
        RMat[t,j] <- pow_int(a,t-j)
      } 
    }
    return(RMat)
  }, buildDerivs = list(run = list(ignore = c('t','j')))
)






# ## Does not work!!!
# MAQFun <- nimbleFunction(
#   run = function(n = double(1), #Vector Of Jump Occurence
#                  y = double(1), #Vector of Jump Intensity
#                  Q = double(0), #Upper Limit of MAQ
#                  t = double(0), #current time index
#                  a = double(0)) #MAQ Parameter
#   {
#     returnType(double(0)) #returns a scalar
#     Error <- 0 #Start with zero
#     for(q in 1:min(t,Q)){ #add to current Error 
#       Error <- Error + pow(a,q)*(n[q]*y[q])
#     }
#     return(Error)
#   }
# )


CovMatMult <- nimbleFunction(
  run= function(A = double(0), 
                sigma=double(0)){
    returnType(double(2))
    CovMatTot <- matrix(value=0, nrow = (A-1)*2, ncol=(A-1)*2)
    CovMatTot[1:(A-1), 1:(A-1)] <-  CovMat(A,sigma) #Result of Nimble Function
    CovMatTot[A:(2*A-2), A:(2*A-2)] <- CovMat(A,sigma) #Result of Nimble Function
  return(CovMatTot)
  },
  buildDerivs = "run"
)




