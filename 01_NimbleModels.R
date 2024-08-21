##### NIMBLE FUNCTIONS ######
pacman::p_load("nimble")


## SINGLE POPULATION MODELS ####################################################

#Liu, Li Model with Sum to One on Beta and Corner Constraint on kappa 1
LiuLi_Model <- nimbleCode({
  
  # Liu,Li Model with Jumps
  #
  # Time Effect (iid Normal Prior)+
  # Age Effect (Dirichtlet Prior) +
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  
  #Priors on Time Parameters
  for(t in 2:(N_Year)){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #differenced Random Walk
  }
  
  k[1] <- drift
  sigma_time ~ T(dnorm(mean=0, sd=1), min = 0, ) #truncated normal
  drift ~ dnorm(mean=0, sd=5) #normal
  
  #other time effects
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0 #Corner Constraint on Delta J1 = 0
  for (t in 3:(N_Year+1)){ #one year longer for differencing
    N_t[t] ~ dbern(p)
    Y_t[t] ~ T(dnorm(mean=muY, sd=sdY), min = 0, max = ) 
    J[t] <- N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape=1, shape2 = 20)
  muY ~ T(dnorm(mean=0, sd=4), min = 0, max = ) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, max = ) #truncated normal
  
  #Prior on Age Effects
  #Dirichlet as Normalized Gamma
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
    b2[x]~ dgamma(shape = 1, rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  sigma_eps ~ T(dnorm(mean=0, sd=2), min = 0, )
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mean = mu[x,t],
                        sd=sigma_eps)
    }
  }
})

# AR Model with Sum to Zero constraint on beta and corner constraint on kappa1
OwnModel_AR <- nimbleCode({
  
  # Own Model with autoregressive jump component 

  # Time Effect (iid Normal Prior)+
  # Age Effect (Dirichtlet Prior) +
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  
  #Priors on Time Parameters
  for(t in 2:(N_Year)){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #Differenced Random Walk
  }
  
  k[1] <- drift #Corner constraint
  sigma_time ~ T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  drift ~ dnorm(mean=0, sd=5) 
  
  #Other Time Effect
  N_t[1:2] <- 0 #Corner Constraint on Delta J_1
  Y_t[1:2] <- 0
  J[1:2] <- 0
  
  N_t[N_Year+1] <- 0 #Corner Constraint on N_T 
  Y_t[N_Year+1] <- 0
  
  for (t in 3:(N_Year)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ T(dnorm(mean = muY, sd = sdY), min = 0, max = )
    J[t] <- a*J[t-1]+N_t[t]*Y_t[t]
  }
  
  J[N_Year+1] <- a*J[N_Year] 
  
  #low values are preferred
  p ~ dbeta(shape=1, shape2 = 20)
  
  muY ~ T(dnorm(mean=0, sd=4), min = 0, max = ) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, max = ) #truncated normal
  
  #Vanishing Parameter
  a ~ T(dnorm(mean = 0, sd = 0.4), min = 0, max = 1) # low values are preferred 
  
  # Age Effects
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
    b2[x]~ dgamma(shape = 1, rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  
  sigma_eps ~ T(dnorm(mean=0, sd=2), min = 0, max = ) #truncated normal
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
                 betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
    
  }
})

# MA Model with Sum to Zero constraint on beta and corner constraint on kappa1
OwnModel_MA <- nimbleCode({
  
  # Own Model with moving average jump component 
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Sum to zero Decomposition) +
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  for(t in 2:(N_Year)){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formulation
  }
  sigma_time ~  T(dnorm(mean=0, sd=1), min = 0, max = ) #truncated normal
  k[1] <- drift
  
  #alternatively drift must be negative
  drift ~ dnorm(mean=0, sd=5) 
  
  #other time effects (with Corner Constraint)
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0
  
  N_t[N_Year+1] <- 0
  Y_t[N_Year+1] <- 0
  J[N_Year+1] <- b*N_t[N_Year]*Y_t[N_Year] #No Jump Last Year
  
  for (t in 3:(N_Year)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ T(dnorm(mean = muY, sd = sdY), min = 0, max = )
    error[t] <- b*N_t[(t-1)]*Y_t[(t-1)]
    # #Multiple assignments not allowed
    J[t] <-  N_t[t]*Y_t[t]+error[t]
  }
  # 
  
  
  p ~ dbeta(shape1 = 1,shape2 = 20) #low values are preferred
  
  muY ~ T(dnorm(mean=0, sd=4), min = 0, max = ) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, max = ) #truncated normal

  
  b ~  T(dnorm(mean=0, sd=0.4), min = 0, max = 1)
  
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
    b2[x]~ dgamma(shape = 1, rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  sigma_eps ~ T(dnorm(mean=0, sd=2), min = 0, max = ) #truncated normal
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
                 betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
    
  }
})

############## MULTI POPULATION MODELS #########################################


MultiPop_AR_GlobalN <- nimbleCode({
  
  # Multi-population Version of AR model with global N_t parameter 
  #Prior Parameterisation
  
  ## Trivariate Normal modelled in terms of conditionals and marginals due to better performance
  rho12 ~ T(dnorm(mean = 0, sd = 0.5), min = -1, max = 1)
  rho13 ~ T(dnorm(mean = 0, sd = 0.5), min = -1, max = 1)
  rho23 ~ T(dnorm(mean = 0, sd = 0.5), min = -1, max = 1)
  
  
  for(t in 2:N_Year){
    k_Mat[t, 1] ~ dnorm(mean = drift[1], sd = sigma_time[1])
    
    meanMult1[t] <- drift[2] + rho12 * sigma_time[2] / sigma_time[1] *
      (k_Mat[t, 1] - drift[1])
    
    SdMult1[t] <- sigma_time[2] * sqrt(1 - pow(rho12, 2))
    
    k_Mat[t, 2] ~ dnorm(mean = meanMult1[t], sd = SdMult1[t])
    
    meanMult2[t] <- drift[3] + 
      (rho13 * sigma_time[3] / sigma_time[1]) * (k_Mat[t, 1] - drift[1]) +
      (rho23 * sigma_time[3] / sigma_time[2]) * (k_Mat[t, 2] - drift[2])
    
    SdMult2[t] <- sigma_time[3] * sqrt(1 - pow(rho13, 2) - pow(rho23, 2))
    
    k_Mat[t, 3] ~ dnorm(mean = meanMult2[t], sd =  SdMult2[t])
  }
  
  k_Mat[1,1:N_Country] <- drift[1:N_Country]
  
  for(c in 1:N_Country){
    sigma_time[c] ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
    drift[c] ~ dnorm(mean = 0, sd = 5)
    a[c] ~ T(dnorm(mean=0, sd=0.4), min = 0, max = 1)
  }
 
  p ~ dbeta(shape = 1, shape2 = 20)

  #Jump Parameters
  N_t[1:2] <- 0
  Y_t_Mat[1:2,1:N_Country] <- 0
  J_t_Mat[1:2,1:N_Country] <- 0
  
  for (t in 3:(N_Year)){ #one year longer
    N_t[t] ~ dbern(p)
    for(c in 1:N_Country){
      Y_t_Mat[t,c] ~ T(dnorm(mean = muY[c], sd = sdY[c]), min = 0, max = )
      J_t_Mat[t,c] <- a[c]*J_t_Mat[t-1,c]+N_t[t]*Y_t_Mat[t,c]
    }
  }
  N_t[N_Year+1] <- 0 #Corner Constraint 
  Y_t_Mat[N_Year+1,1:N_Country] <- 0
  
  for(c in 1:N_Country){
    J_t_Mat[N_Year+1,c] <- a[c]*J_t_Mat[N_Year,c]
  }
  
  #Jump Intensity Parameters
  for(c in 1:N_Country){
    muY[c] ~ T(dnorm(mean=0, sd=4), min = 0,)
    sdY[c] ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
  }
  
  #other time effects
  # Dirichlet Distribution is standardized Gamma Dist
  for(c in 1:N_Country){
    for(x in 1:N_AgeGroups){
      b1[x,c]~ dgamma(shape = 1, rate = 1)
      b2[x,c]~ dgamma(shape = 1, rate = 1)
    }
    
    # Dirichlet Distribution is standardized Gamma Dist
    beta[1:(N_AgeGroups-1),c] <- b1[1:(N_AgeGroups-1),c]/sum(b1[1:N_AgeGroups,c])
    beta[N_AgeGroups,c] <- 1 - sum(beta[1:(N_AgeGroups-1),c])
    
    betaJump[1:(N_AgeGroups-1),c] <- b2[1:(N_AgeGroups-1),c]/sum(b2[1:N_AgeGroups,c])
    betaJump[N_AgeGroups,c] <- 1 - sum(betaJump[1:(N_AgeGroups-1),c])
  }
  
  for(c in 1:N_Country){
    sigma_eps[c] ~ T(dnorm(mean=0, sd=2), min = 0,)
  }
  
  
  #Putting all Parameters together
  for(c in 1:N_Country){
    for(x in 1:N_AgeGroups){
      for(t in 1:N_Year){
        mu[x,t,c] <- beta[x,c]*k_Mat[t,c]+betaJump[x,c]*(J_t_Mat[t+1,c]-J_t_Mat[t,c])
        ZMat[x,t,c] ~ dnorm(mean = mu[x,t,c],
                            sd=sigma_eps[c])
      }
    }
  }
})

MultiPop_MA_GlobalN <- nimbleCode({
  
  # Multi-population Version of AR model with global N_t parameter 
  
  #Prior Parameterisation
  rho12 ~ T(dnorm(mean = 0, sd = 0.5), min = -1, max = 1)
  rho13 ~ T(dnorm(mean = 0, sd = 0.5), min = -1, max = 1)
  rho23 ~ T(dnorm(mean = 0, sd = 0.5), min = -1, max = 1)
  
  ## Trivariate Normal modelled in terms of conditionals and marginals
  for(t in 2:N_Year){
    k_Mat[t,1] ~ dnorm(mean = drift[1], sd = sigma_time[1])
    
    meanMult1[t] <- drift[2] + rho12 * sigma_time[2] / sigma_time[1] * (k_Mat[t, 1] -
                                                                          drift[1])
    SdMult1[t] <- sigma_time[2] * sqrt(1 - rho12 ^ 2)
    
    k_Mat[t, 2] ~ dnorm(mean = meanMult1[t], sd = SdMult1[t])
    
    meanMult2[t] <- drift[3] + (rho13 * sigma_time[3] / sigma_time[1]) *
      (k_Mat[t, 1] - drift[1]) +
      (rho23 * sigma_time[3] / sigma_time[2]) * (k_Mat[t, 2] - drift[2])
    
    SdMult2[t] <- sigma_time[3] * sqrt(1 - rho13^2 - rho23^2)
    
    k_Mat[t, 3] ~ dnorm(mean = meanMult2[t], sd =  SdMult2[t])
  }
  
  k_Mat[1,1:N_Country] <- drift[1:N_Country]
  
  for(c in 1:N_Country){
    sigma_time[c] ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
    drift[c] ~ dnorm(mean = 0, sd = 5)
    b[c] ~ T(dnorm(mean=0, sd=0.4), min = 0, max = 1)
  }
  
  p ~ dbeta(shape = 1, shape2 = 20)

  N_t[1:2] <- 0
  Y_t_Mat[1:2,1:N_Country] <- 0
  J_t_Mat[1:2,1:N_Country] <- 0
  
  
  
  for (t in 3:(N_Year)){ #one year longer
    N_t[t] ~ dbern(p)
    for(c in 1:N_Country){
      Y_t_Mat[t,c] ~ T(dnorm(mean = muY[c], sd = sdY[c]),min = 0, max = )
      error[t,c] <- b[c]*N_t[(t-1)]*Y_t_Mat[(t-1),c]
      J_t_Mat[t,c] <- error[t,c]+N_t[t]*Y_t_Mat[t,c]
    }
  }
  
  N_t[N_Year+1] <- 0 #Corner Constraint 
  Y_t_Mat[N_Year+1,1:N_Country] <- 0
  
  for(c in 1:N_Country){
    J_t_Mat[N_Year+1,c] <-  b[c]*N_t[N_Year]*Y_t_Mat[N_Year,c] #No Jump Last Year 
  }
  
  
  for(c in 1:N_Country){
    muY[c] ~ T(dnorm(mean=0, sd=4), min = 0,)
    sdY[c] ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  for(c in 1:N_Country){
    for(x in 1:N_AgeGroups){
      b1[x,c]~ dgamma(shape = 1, rate = 1)
      b2[x,c]~ dgamma(shape = 1, rate = 1)
    }
    
    # Dirichlet Distribution is standardized Gamma Dist
    beta[1:(N_AgeGroups-1),c] <- b1[1:(N_AgeGroups-1),c]/sum(b1[1:N_AgeGroups,c])
    beta[N_AgeGroups,c] <- 1 - sum(beta[1:(N_AgeGroups-1),c])
    
    betaJump[1:(N_AgeGroups-1),c] <- b2[1:(N_AgeGroups-1),c]/sum(b2[1:N_AgeGroups,c])
    betaJump[N_AgeGroups,c] <- 1 - sum(betaJump[1:(N_AgeGroups-1),c])
  }
  
  
  for(c in 1:N_Country){
    sigma_eps[c] ~ T(dnorm(mean=0, sd=2), min = 0,)
  }
  
  
  #Putting all Parameters together
  for(c in 1:N_Country){
    for(x in 1:N_AgeGroups){
      for(t in 1:N_Year){
        mu[x,t,c] <- beta[x,c]*k_Mat[t,c]+betaJump[x,c]*(J_t_Mat[t+1,c]-J_t_Mat[t,c])
        ZMat[x,t,c] ~ dnorm(mean = mu[x,t,c],
                            sd = sigma_eps[c])
      }
    }
  }
})
