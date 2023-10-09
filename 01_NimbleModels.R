##### NIMBLE FUNCTIONS ######
library(nimble)


## ALL MODELS BELOW ARE ON DIFFERENCED LOG DEATH RATES

#Liu, Li Model with Sum to One on Beta and Corner Constraint on kappa 1
LiuLi_Model <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
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
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #Differenced Random Walk
  }
  
  k[1] <- drift #Corner constraint
  sigma_time ~ T(dnorm(mean=0, sd=1), min = 0, ) #truncated normal
  drift ~ dnorm(mean=0, sd=2) #normal
  
  #other time effects
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0 #Corner Constraint on Delta J1 = 0
  for (t in 3:(N_Year+1)){ #one year longer for differencing
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape=1, shape2 = 25)
  muY ~ T(dnorm(mean=0, sd=5), min = 0, max = 1000) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, max = 1000) #truncated normal
  
  #Prior on Age Effects
  #Dirichlet as Normalized Gamma
  #alphaJump[1:N_AgeGroups] <- c(rep(1,3),5,5,5,5,rep(1,3))
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
    #b2[x]~ dgamma(shape = alphaJump[x], rate = 1)
    b2[x]~ dgamma(shape = 1, rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  sigma_eps ~ T(dnorm(mean=0, sd=1), min = 0, )
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mean = mu[x,t],
                        sd=sigma_eps)
    }
  }
})

# Own Model with Sum to Zero constraint on beta and corner constraint on kappa1
OwnModel <- nimbleCode({
  
  # Own Model on Mortality Improvement Rates (Route II)

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
  sigma_time ~ T(dnorm(mean=0, sd=1), min = 0, ) #truncated normal
  drift ~ dnorm(mean=0, sd=2) #normal
  
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0 #Corner Constraint on Delta J1 = 0
  for (t in 3:(N_Year+1)){ #one year longer for differencing
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- a*J[t-1]+N_t[t]*Y_t[t]
  }

  #low values are preferred
  p ~ dbeta(shape=1, shape2 = 25)
  
  muY ~ T(dnorm(mean=1, sd=2), min = 0, max = ) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, max = ) #truncated normal
  
  #Vanishing Parameter
  a ~ dbeta(shape1 = 1, shape2 = 5) # low values are preferred 
  
  # Age Effects
  #alphaJump[1:N_AgeGroups] <- c(rep(1,3),5,5,5,5,rep(1,3))
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
    b2[x]~ dgamma(shape = 1, rate = 1)
    #b2[x]~ dgamma(shape = alphaJump[x], rate = 1)
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

## Lee Carter Model differenced Data
LC_Diff <- nimbleCode({
  
  # Lee Carter Model differenced
  #
  # Time Effect (iid Normal Prior)+
  # Age Effect (Dirichtlet Prior) +
  # Overdispersion Parameter (normal prior)+
  #
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:(N_Year)){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  
  sigma_time ~  T(dnorm(mean=0, sd=2), min = 0, ) 
  drift ~ dnorm(mean = 0,sd = 5)
  
  #Dirichlet Prior
  #alphaJump[1:N_AgeGroups] <- c(rep(1,3),3,3,3,rep(1,4))
  # a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  # beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
  }

  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
   
  
  sigma_eps ~ T(dnorm(mean=0, sd=2), min = 0, ) 
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})


