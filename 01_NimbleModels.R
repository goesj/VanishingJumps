##### NIMBLE FUNCTIONS ######
library(nimble)


#### MODELS IN USE FOR ANALYSIS#################################################

## ALL MODELS BELOW ARE ON DIFFERENCED LOG DEATH RATES

## Liu,Li Model with QR for Beta
#Liu, Li Model with Sum to One on Beta and Corner Constraint on kappa 1
LiuLi_Diff_OC_CornerK <- nimbleCode({
  
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
  #differenced random Walk 
  for(t in 2:(N_Year)){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  
  k[1] <- drift
  sigma_time ~ T(dnorm(mean=0, sd=0.5), min = 0, ) #truncated normal
  drift ~ dnorm(mean=0, sd=2) #truncated normal
  
  #other time effects
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape=1, shape2 = 20)
  muY ~ T(dnorm(mean=0, sd=5), min = 0, max = 1000) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, max = 1000) #truncated normal
  

  alphaJump[1:N_AgeGroups] <- c(rep(1,3),3,3,3,3,rep(1,3))
  #alphaN[1:N_AgeGroups] <- c(3,3,3,1,1,1,1,rep(1,3))
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
    #b1[x]~ dgamma(shape = alphaN[x], rate = 1)
    b2[x]~ dgamma(shape = alphaJump[x], rate = 1)
    #b2[x]~ dgamma(shape = 1, rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  sigma_eps ~ dunif(0,2)
  
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
LCJumpOwn_OC_CornerK <- nimbleCode({
  
  # Own Model on Differenced Rates
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
  #differenced random Walk 
  for(t in 2:(N_Year)){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formulation
  }
  sigma_time ~  T(dnorm(mean=0, sd=1), min = 0, max = ) #truncated normal
  k[1] <- drift
  
  #alternatively drift must be negative
  drift ~ T(dnorm(mean=0, sd=2), min = -1000, max = ) #truncated normal
  
  #other time effects (with Corner Constraint)
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0
  
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- a*J[t-1]+N_t[t]*Y_t[t]
  }

  #p ~ dbeta(shape1 = 1,shape2 = 5) #low values are preferred
  p ~ dbeta(shape=1, shape2 = 20)
  
  muY ~ T(dnorm(mean=1, sd=2), min = 0, max = ) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, max = ) #truncated normal
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1, shape2 = 5) # low values are preferred 
  
  # Age Effects
  # a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  # 
  # beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # 
  
  alphaJump[1:N_AgeGroups] <- c(rep(1,3),5,5,5,5,rep(1,3))
  #alphaJump[1:N_AgeGroups] <- c(rep(1,N_AgeGroups))
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
    #b2[x]~ dgamma(shape = 1, rate = 1)
    b2[x]~ dgamma(shape = alphaJump[x], rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  sigma_eps ~ T(dnorm(mean=0, sd=1), min = 0, max = ) #truncated normal
  #sigma_eps ~ dunif(min=0, max=5)
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
        betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
    
  }
})

## Lee Carter Model Differenced Data
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


#### NOT IN USE !!!!!###########################################################
LiuLi_Diff_QR <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Normal Prior Prior) +
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formulation
  }
  sigma_time ~  dunif(0,2) 
  #drift ~ T(dnorm(mean=0, sd=3), min = -1000, max = 0) #truncated normal
  drift ~ dnorm(mean=0, sd=3)
  #other time effects
  N_t[1] <- 0
  Y_t[1] <- 0
  J[1] <- 0
  for (t in 2:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape1 = 1, shape2 = 20)
  muY ~ dchisq(df = 2)
  #muY ~ T(dnorm(mean=1, sd=4), min = 0, max = 1000) #truncated normal
  #muY ~ dgamma(shape =1.4,scale = 4)
  sdY ~ dunif(min=0, max=5)
  #sdY ~ dgamma(shape = 1,rate = 1)
  
  #Normal Prior
  for(x in 1:N_AgeGroups){
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- QMat[x,1]*k[t]+
        QMat[x,2]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

#Liu Li Model, Simple Version (no hyperpriors on the Jump Parameter Y), QR Beta
LiuLi_Diff_QR_S <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Normal Prior Prior) +
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formulation
  }
  sigma_time ~  dunif(0,2) 
  #drift ~ T(dnorm(mean=0, sd=3), min = -1000, max = 0) #truncated normal
  drift ~ T(dnorm(mean=0, sd=3), min = -1000, max = 0) #truncated normal
  
  #other time effects
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0
  
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = 0, sd = 0.5)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape1 = 1,shape2 = 4) #low values are preferred
  #muY ~ T(dnorm(mean=1, sd=4), min = 0, max = 1000) #truncated normal
  # muY ~ dgamma(shape =1.4,scale = 4)
  # sdY ~ dunif(min=0, max=5)
  #sdY ~ dgamma(shape = 1,rate = 1)
  
  #Normal Prior
  for(x in 1:N_AgeGroups){
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- QMat[x,1]*k[t]+
        QMat[x,2]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})
#Liu, Li Model with Sum to One on Beta
LiuLi_Diff_OC <- nimbleCode({
  
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
  #differenced random Walk 
  for(t in 1:(N_Year)){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  
  sigma_time ~ T(dnorm(mean=0, sd=0.5), min = 0, ) #truncated normal
  drift ~ dnorm(mean=0, sd=2) #truncated normal
  
  #other time effects
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape1 = 1,shape2 = 25) #low values are preferred
  muY ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
  #muY ~ dnorm(mean=0, sd=5)
  #sdY ~ dunif(min=0, max=5)
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  #Dirichlet Prior
  
  # for(x in 1:N_AgeGroups){
  #   b1[x]~ dnorm(mean = 0, sd=5)
  #   b2[x]~ dnorm(mean = 0, sd=5)
  # }
  # 
  #Sum equal to one
  # beta[1:N_AgeGroups] <- b1[1:N_AgeGroups]/sum(b1[1:N_AgeGroups])
  # betaJump[1:N_AgeGroups] <- b2[1:N_AgeGroups]/sum(b2[1:N_AgeGroups])
  # 
  # a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  
  #
  # beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # betaJump[1:N_AgeGroups] ~ ddirch(alpha=alphaJump[1:N_AgeGroups])
  #alphaJump[1:N_AgeGroups] <- c(rep(1,3),3,3,1,1,rep(1,3))
  #alphaN[1:N_AgeGroups] <- c(3,3,3,1,1,1,1,rep(1,3))
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
    #b1[x]~ dgamma(shape = alphaN[x], rate = 1)
    #b2[x]~ dgamma(shape = alphaJump[x], rate = 1)
    b2[x]~ dgamma(shape = 1, rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  # # #
  #Mult Normal Prior for Betas
  # meanVec[1:(N_AgeGroups-1)] <- rep((1/N_AgeGroups),N_AgeGroups-1)
  # #meanVecJ[1:(N_AgeGroups-1)] <- c(rep(1/9,2),rep(3/9,4),rep(1,3))
  # 
  # NormMat[1:(N_AgeGroups-1),
  #         1:(N_AgeGroups-1)] <- CovMat(A=N_AgeGroups,
  #                                      sigma = sigma_age)
  # 
  # b1[1:(N_AgeGroups-1)] ~ dmnorm(mean= meanVec[1:(N_AgeGroups-1)],
  #                                cov = NormMat[1:(N_AgeGroups-1),
  #                                              1:(N_AgeGroups-1)])
  # 
  # b2[1:(N_AgeGroups-1)] ~ dmnorm(mean= meanVec[1:(N_AgeGroups-1)],
  #                                cov = NormMat[1:(N_AgeGroups-1),
  #                                              1:(N_AgeGroups-1)])
  # 
  # beta[1] <-  1 - sum(b1[1:(N_AgeGroups-1)])
  # beta[2:N_AgeGroups] <-  b1[1:(N_AgeGroups-1)]
  # 
  # betaJump[1] <- 1 - sum(b2[1:(N_AgeGroups-1)])
  # betaJump[2:N_AgeGroups] <- b2[1:(N_AgeGroups-1)]
  # # 
  # sigma_age ~ dunif(0,2)
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mean = mu[x,t],
                        sd=sigma_eps)
    }
  }
})

#LiuLi_Repara_Joint (Model for joint estimation of Countries (hierarchical))
LiuLi_Diff_OC_Joint <- nimbleCode({
  
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
  #differenced random Walk 
  for(t in 1:N_Year){
    for(c in 1:N_Country){ #number of countries
      k_Mat[t,c] ~ dnorm(mean = drift[c],sd=sigma_time[c]) #State space model formualtion
    }
  }
  for(c in 1:N_Country){
    sigma_time[c] ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
    drift[c] ~ dnorm(mean = 0, sd = 2)
    p[c] ~ dbeta(shape1 = 1,shape2 = 20) #joint prior
    #a[c] ~ dbeta(shape=a_Glob, shape2 =b_Glob)
  }
  # a_Glob ~ T(dnorm(mean=1, sd=5), min = 0,) #truncated normal
  # b_Glob ~ T(dnorm(mean=5, sd=5), min = 0,) #truncated normal
  
  N_t_Mat[1:2,1:N_Country] <- 0 #Corner Constraint N_t[1:2] <- 0
  Y_t_Mat[1:2,1:N_Country] <- 0
  J_t_Mat[1:2,1:N_Country] <- 0
  for (t in 3:(N_Year+1)){ #one year longer
    for(c in 1:N_Country){
      N_t_Mat[t,c] ~ dbern(p[c])
      Y_t_Mat[t,c] ~ dnorm(mean = muY[c], sd = sdY[c])
      J_t_Mat[t,c] <- N_t_Mat[t,c]*Y_t_Mat[t,c]
    }
  }
  
  #Global prior
  muY_Glob ~ T(dnorm(mean=0, sd=2), min = 0,)
  sdY_Glob ~ T(dnorm(mean=0, sd=1), min = 0,) #truncated normal
  # Tau_Glob ~ T(dnorm(mean=0, sd=0.5), min = 0,)
  # TauSD_Glob ~T(dnorm(mean=0, sd=0.5), min = 0,)
  
  for(c in 1:N_Country){
    muY[c] ~ T(dnorm(mean=muY_Glob, sd=0.5), min = 0,)
    sdY[c] ~ T(dnorm(mean=sdY_Glob, sd=0.5), min = 0,) #truncated normal
  }
  
  #other time effects
  alphaDirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  
  for(c in 1:N_Country){
    beta[1:N_AgeGroups,c] ~ ddirch(alpha=alphaDirich[1:N_AgeGroups])
    betaJump[1:N_AgeGroups,c] ~ ddirch(alpha=alphaDirich[1:N_AgeGroups])
  }
  
  
  for(c in 1:N_Country){
    sigma_eps[c] ~ dunif(0,2)
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

#Own Model QR Beta
LCJumpOwn_QR <- nimbleCode({
  
  # Own Model on Differenced Rates
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (QR Decomposition) +
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formulation
  }
  
  # sigma_time <- exp(sigma_time_log)
  #sigma_time ~  dunif(0,2)
  sigma_time ~ T(dnorm(mean=0, sd=3), min = 0,) 
  
  
  #alternatively drift must be negative
  #drift ~ dunif(min = -5, max = 0) #making the sign of QR identifiable
  drift ~ dnorm(0,sd=2) 
  #drift ~ T(dnorm(mean=0, sd=3), min = -1000, max = 0) #truncated normal
  
  #other time effects (with Corner Constraint)
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0
  
  # N_t[1] <- 0
  # Y_t[1] <- 0
  # J[1] <- 0
  
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- a*J[t-1]+N_t[t]*Y_t[t]
  }
  
  #p ~ dbeta(shape1 = 1,shape2 = 5) #low values are preferred
  p ~ dbeta(shape=1, shape2 = 20)
  
  muY ~ T(dnorm(mean=1, sd=1), min = 0, ) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1, shape2 = 5) # low values are preferred 
  
  
  for(x in 1:N_AgeGroups){
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(min= 0, max = 2) #weakly informative 
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- QMat[x,1]*k[t]+
        QMat[x,2]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
    
  }
})
#Own Model QR Beta, simple version (see Info above)
LCJumpOwn_QR_Simple <- nimbleCode({
  
  # Own Model on Differenced Rates
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (QR Decomposition) +
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  
  # sigma_time_log ~ dLogGamma(shape=1, rate=1)
  # sigma_time <- exp(sigma_time_log)
  sigma_time ~  dunif(0,2) 
  
  
  #alternatively drift must be negative
  drift ~ dunif(min = -5, max = 0) #making the sign of QR identifiable
  #drift ~ T(dnorm(mean=0, sd=3), min = -1000, max = 0) #truncated normal
  #drift ~ dnorm(0,sd=2) 
  #drift ~ T(dnorm(mean=0, sd=3), min = -1000, max = 0) #truncated normal
  
  #other time effects (with Corner Constraint)
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0
  
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = 0, sd = 0.5)
    J[t] <- a*J[t-1]+N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape1 = 1,shape2 = 20) #low values are preferred
  #p ~ dbeta(shape1 = 1,shape2 = 5) #low values are preferred
  
  #muY ~ T(dnorm(mean=1, sd=4), min = 0, max = 1000) #truncated normal
  #muY ~ dgamma(shape =1.4,scale = 4) #similar in shape to above truncated normal
  #muY ~ dnorm(mean = 0, sd = 5)
  #sdY ~ dgamma(shape = 2,rate = 1) #weakly informative
  
  #sdY ~ dunif(min=0, max=5)
  # muY ~ dnorm(mean = 0, sd = 5)
  # sdY ~ dunif(min = 0,  max = 2)
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1, shape2 = 5) # low values are preferred 
  
  
  for(x in 1:N_AgeGroups){
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(min= 0, max = 2) #weakly informative 
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- QMat[x,1]*k[t]+
        QMat[x,2]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

# Own Model with Sum to Zero constraint
LCJumpOwn_OC <- nimbleCode({
  
  # Own Model on Differenced Rates
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
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formulation
  }
  #sigma_time ~ dunif(min=0, max=5)
  sigma_time ~  T(dnorm(mean=0, sd=2), min = 0, max = 1000) #truncated normal
  
  #alternatively drift must be negative
  drift ~ T(dnorm(mean=0, sd=2), min = -1000, max = 0) #truncated normal
  
  #other time effects (with Corner Constraint)
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0
  
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- a*J[t-1]+N_t[t]*Y_t[t]
  }
  
  #p ~ dbeta(shape1 = 1,shape2 = 5) #low values are preferred
  p ~ dbeta(shape=1, shape2 = 20)
  
  muY ~ T(dnorm(mean=0, sd=2), min = 0, max = 1000) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=1), min = 0, max = 1000) #truncated normal
  
  
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1, shape2 = 5) # low values are preferred 
  
  # Age Effects
  # a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  # 
  # beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # 
  
  #alphaJump[1:N_AgeGroups] <- c(rep(1,3),3,3,1,1,rep(1,3))
  #alphaJump[1:N_AgeGroups] <- c(rep(1,N_AgeGroups))
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
  
  sigma_eps ~ T(dnorm(mean=0, sd=5), min = 0, max = 1000) #truncated normal
  #sigma_eps ~ dunif(min=0, max=5)
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
        betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
    
  }
})

#Own Model for all Countries jointly (Model for joint estimation of Countries (hierarchical))
LCJumpOwn_OC_Joint <- nimbleCode({
  
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
  #differenced random Walk 
  for(t in 1:N_Year){
    for(c in 1:N_Country){ #number of countries
      k_Mat[t,c] ~ dnorm(mean = drift[c],sd=sigma_time[c]) #State space model formualtion
    }
  }
  for(c in 1:N_Country){
    sigma_time[c] ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
    drift[c] ~ dnorm(mean = 0, sd = 2)
    p[c] ~ dbeta(shape1 = 1,shape2 = 20) #joint prior
    #a[c] ~ dbeta(shape=a_Glob, shape2 =b_Glob)
    a[c] ~ dbeta(shape1 = 1, shape2 = 5)
  }
  
  
  # a_Glob ~ T(dnorm(mean=1, sd=5), min = 0,) #truncated normal
  # b_Glob ~ T(dnorm(mean=5, sd=5), min = 0,) #truncated normal
  
  N_t_Mat[1:2,1:N_Country] <- 0 #Corner Constraint N_t[1:2] <- 0
  Y_t_Mat[1:2,1:N_Country] <- 0
  J_t_Mat[1:2,1:N_Country] <- 0
  for (t in 3:(N_Year+1)){ #one year longer
    for(c in 1:N_Country){
      N_t_Mat[t,c] ~ dbern(p[c])
      Y_t_Mat[t,c] ~ dnorm(mean = muY[c], sd = sdY[c])
      J_t_Mat[t,c] <- a[c]*J_t_Mat[t-1,c]+N_t_Mat[t,c]*Y_t_Mat[t,c]
    }
  }
  
  #Global prior
  muY_Glob ~ T(dnorm(mean=0, sd=2), min = 0,)
  sdY_Glob ~ T(dnorm(mean=0, sd=1), min = 0,) #truncated normal
  # Tau_Glob ~ T(dnorm(mean=0, sd=0.5), min = 0,)
  # TauSD_Glob ~T(dnorm(mean=0, sd=0.5), min = 0,)
  
  for(c in 1:N_Country){
    muY[c] ~ T(dnorm(mean=muY_Glob, sd=0.5), min = 0,)
    sdY[c] ~ T(dnorm(mean=sdY_Glob, sd=0.5), min = 0,) #truncated normal
  }
  # muY[1]~T(dnorm(mean=0, sd=3), min = 0,)
  # muY[2]~T(dnorm(mean=0, sd=3), min = 0,)
  
  
  #other time effects
  alphaDirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  
  for(c in 1:N_Country){
    beta[1:N_AgeGroups,c] ~ ddirch(alpha=alphaDirich[1:N_AgeGroups])
    betaJump[1:N_AgeGroups,c] ~ ddirch(alpha=alphaDirich[1:N_AgeGroups])
  }
  
  
  
  #alphaJump[1:N_AgeGroups] <- c(rep(1,3),3,3,1,1,rep(1,3))
  #alphaN[1:N_AgeGroups] <- c(3,3,3,1,1,1,1,rep(1,3))
  # for(x in 1:N_AgeGroups){
  #   b1[x]~ dgamma(shape = 1, rate = 1)
  #   #b1[x]~ dgamma(shape = alphaN[x], rate = 1)
  #   #b2[x]~ dgamma(shape = alphaJump[x], rate = 1)
  #   b2[x]~ dgamma(shape = 1, rate = 1)
  # }
  # 
  # # Dirichlet Distribution is standardized Gamma Dist
  # beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  # beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  # 
  # betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  # betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  # 
  # 
  for(c in 1:N_Country){
    sigma_eps[c] ~ dunif(0,2)
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

############ Reparameterization of ZMatrix (see Pdf in Theory folder)###########
#Liu, Li Model QR Beta
LiuLi_Repara_QR <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1)
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  
  N_t[1:2] <- 0 #Corner Constraint 
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
  }
  p ~ dbeta(shape1 = 1,shape2 = 20) 
  muY ~ dnorm(mean=0, sd=4) #weakly informative
  sdY ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
  
  #Works with both parameterization on real data
  #a ~ dbeta(shape1 = 1,shape2 = 4) # low values are preferred 
  
  ###### AGE PARAMETERS########
  for(x in 1:N_AgeGroups){
    #alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(0,2)
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      
      mu[x,t] <- #alpha[x]+
        QMat[x,1]*k[t]+
        QMat[x,2]*muY*(N_t[t+1]-N_t[t])
      
      sigma[x,t] <- pow(sigma_eps,2)  +
        pow(QMat[x,2]*N_t[t+1]*sdY,2) + #Var first part
        pow(QMat[x,2]*N_t[t]*sdY,2)     #Var second part
      
      
      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma[x,t])
    }
  }
  
})

#Liu Li Model, sum to one beta (here kappa is still a parameter)
LiuLi_Repara_Diff_OC <- nimbleCode({
  
  # Differenced Liu,Li Model with Own Constraints (no QR)
  #
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
  #drift ~ dunif(min = -2,max = 0) 
  drift ~ dnorm(mean = 0, sd = 2)
  
  N_t[1:2] <- 0 #Corner Constraint 
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
  }
  #p ~ dbeta(shape1 = 1,shape2 = 4) 
  p ~ dbeta(shape1 = 1,shape2 = 20) 
  #muY ~ dnorm(mean=0, sd=4) #weakly informative
  #muY ~ dchisq(df = 2)
  muY ~ T(dnorm(mean=0, sd=5), min = 0,)
  #muY ~ T(dnorm(mean=2, sd=1), min = 0,) #truncated normal
  #sdY ~ dgamma(shape = 1,rate = 1) #weakly informative
  sdY ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
  # muY~ dnorm(mean=2, sd=0.05)
  # sdY ~ T(dnorm(mean=1.3, sd=0.05), min = 0, ) #truncated normal
  # 
  #Works with both parameterization on real data
  #a ~ dbeta(shape1 = 1,shape2 = 4) # low values are preferred 
  
  ###### AGE PARAMETERS########
  #alphaJump[1:N_AgeGroups] <- c(rep(1,3),3,3,1,1,rep(1,3))
  #alphaN[1:N_AgeGroups] <- c(3,4,3,1,1,1,1,rep(1,3))
  
  for(x in 1:N_AgeGroups){
    #alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dgamma(shape = 1, rate = 1)
    #b2[x] ~dgamma(shape = alphaJump[x], rate = 1)
    b2[x]~ dgamma(shape = 1, rate = 1)
  }
  
  #Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  
  # a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  # 
  # #
  # beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  sigma_eps ~  T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      
      mu[x,t] <- #alpha[x]+
        beta[x]*k[t]+
        betaJump[x]*muY*(N_t[t+1]-N_t[t])
      
      sigma_squared[x,t] <- pow(sigma_eps,2)  +
        pow_int(betaJump[x]*N_t[t+1]*sdY,2) + #Var first part
        pow_int(betaJump[x]*N_t[t]*sdY,2)     #Var second part
      
      
      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma_squared[x,t])
    }
  }
  
})

#Liu, Li with all Parameters as Random Variables
LiuLi_Repara_Diff_OC_V2<- nimbleCode({
  
  # Differenced Liu,Li Model with Own Constraints (no QR)
  #
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  # for(t in 1:N_Year){
  #   k[t] ~ dnorm(mean = 0,sd=1) #State space model formualtion
  # }
  sigma_time ~  T(dnorm(mean=0, sd=2), min = 0,)
  drift ~ dnorm(0,sd=2) 
  
  #Corner Constraint
  N_t[1:2] <- 0
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
  }
  p ~ dbeta(shape1 = 1,shape2 = 20) 
  muY ~ T(dnorm(mean=0, sd=4), min = 0,) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
  #sdY <- abs(sqrt(VarY))
  # 
  #Works with both parameterization on real data
  #a ~ dbeta(shape1 = 1,shape2 = 4) # low values are preferred 
  
  ###### AGE PARAMETERS########
  for(x in 1:N_AgeGroups){
    #alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dgamma(shape = 1, rate = 1)
    b2[x]~ dgamma(shape = 1, rate = 1)
  }
  
  
  #Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  
  # a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  # 
  # beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  #QR Decomposition to get orthogonal betas
  # BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  # BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  # 
  # QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  # 
  sigma_eps ~ dunif(0,2)
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      
      mu[x,t] <- #alpha[x]+
        beta[x]*drift+ #non centered
        betaJump[x]*muY*(N_t[t+1]-N_t[t])
      
      sigma_squared[x,t] <- pow(sigma_eps,2)+
        pow(beta[x]*sigma_time,2)+
        pow(betaJump[x]*sdY,2)*(N_t[t+1]+N_t[t])
      #pow(betaJump[x]*N_t[t+1]*sdY,2) + #Var first part
      #pow(betaJump[x]*N_t[t]*sdY,2)     #Var second part
      
      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma_squared[x,t])
    }
  }
  
})

#Own Model, sum to one beta (here kappa is still a parameter)
OwnMod_Repara_Diff_OC <- nimbleCode({
  
  # Differenced Liu,Li Model with Own Constraints (no QR)
  #
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  for(t in 2:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  #drift ~ dunif(min = -2,max = 0) 
  k[1] <- drift
  drift ~ dnorm(mean = 0 ,sd = 2)
  
  N_t[1:2] <- 0 #Corner Constraint 
  R[1:(N_Year+1),
    1:(N_Year+1)] <- RTildeMat(Ntime = N_Year,a = a)
  
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
  }
  #p ~ dbeta(shape1 = 1,shape2 = 20) 
  p ~ dbeta(shape1 = 1,shape2 = 20) 
  #muY ~ dchisq(df = 1)
  #muY ~ T(dnorm(mean=0, sd=4), min = 0, ) #truncated normal
  muY ~ T(dnorm(mean=0, sd=4), min = 0, ) #truncated normal
  #sdY ~ dunif(min=0, max=5)
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 =1, shape2 = 5)
  #a ~ dbeta(shape1 = 1,shape2 = 1) # low values are preferred 
  
  ###### AGE PARAMETERS########
  # alphaN[1:N_AgeGroups] <- c(3,3,3,1,1,1,1,rep(1,3))
  # alphaJump[1:N_AgeGroups] <- c(rep(1,3),1,1,1,1,rep(1,3))
  # 
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
  #
  # a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  # 
  # beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  
  #QR Decomposition to get orthogonal betas
  # BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  # BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  # 
  # QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  # 
  sigma_eps ~  T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  
  for(t in 1:N_Year){
    CovVec[t] <- CovZMat(N_t = N_t[1:(N_Year+1)],
                         a =  a,
                         t =  t,
                         sigma2 = pow(sdY,2))
  }
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      
      mu[x,t] <- #alpha[x]+
        beta[x]*k[t]+
        betaJump[x]*muY*(
          inprod(R[t+1,1:(N_Year+1)],N_t[1:(N_Year+1)])- #J_t
            inprod(R[t,1:(N_Year+1)],N_t[1:(N_Year+1)])  #J_t-1
        )
      
      sigma_squared[x,t] <- pow(sigma_eps,2)  + #sigma_eps
        pow(betaJump[x]*sdY*            #betaJ*sdY*J_t
              inprod(R[t+1,1:(N_Year+1)],
                     N_t[1:(N_Year+1)]),2)+  
        pow(betaJump[x]*sdY*                  #betaJ*sdY*J_t-1
              inprod(R[t,1:(N_Year+1)],
                     N_t[1:(N_Year+1)]),2) -  
        2*pow(betaJump[x],2)*CovVec[t] #-Covariance
      
      
      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma_squared[x,t])
    }
  }
  
})

#Own Model, sum to one beta all Parameters as Random Variables
OwnMod_Repara_Diff_OC_V2 <- nimbleCode({
  
  # Differenced Liu,Li Model with Own Constraints (no QR)
  #
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  # for(t in 1:N_Year){
  #   k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  # }
  sigma_time ~  T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  drift ~ dunif(min = -2,max = 0) 
  
  N_t[1:2] <- 0 #Corner Constraint 
  R[1:(N_Year+1),
    1:(N_Year+1)] <- RTildeMat(Ntime = N_Year,a = a)
  
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
  }
  p ~ dbeta(shape1 = 1,shape2 = 25) 
  #muY ~ dchisq(df = 1)
  muY ~ T(dnorm(mean=1, sd=1), min = 0, ) #truncated normal
  sdY ~ T(dnorm(mean=0.5, sd=0.5), min = 0, ) #truncated normal
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1,shape2 = 5) # low values are preferred 
  
  
  ###### AGE PARAMETERS########
  for(x in 1:N_AgeGroups){
    #alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dgamma(shape = 1, rate = 1)
    b2[x]~ dgamma(shape = 1, rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  #
  # 
  sigma_eps ~  T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  
  for(t in 1:N_Year){
    CovVec[t] <- CovZMat(N_t = N_t[1:(N_Year+1)],
                         a =  a,
                         t =  t,
                         sigma2 = pow(sdY,2))
  }
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      
      mu[x,t] <- #alpha[x]+
        beta[x]*drift+
        betaJump[x]*muY*(
          inprod(R[t+1,1:(N_Year+1)],N_t[1:(N_Year+1)])- #J_t
            inprod(R[t,1:(N_Year+1)],N_t[1:(N_Year+1)])    #J_t-1
        )
      
      sigma_squared[x,t] <- pow(sigma_eps,2)  + #sigma_eps
        pow(sigma_time*beta[x],2)+ #sigma_time
        pow(betaJump[x]*sdY*            #betaJ*sdY*J_t
              inprod(R[t+1,1:(N_Year+1)],
                     N_t[1:(N_Year+1)]),2)+  
        pow(betaJump[x]*sdY*                  #betaJ*sdY*J_t-1
              inprod(R[t,1:(N_Year+1)],
                     N_t[1:(N_Year+1)]),2) -  
        2*pow(betaJump[x],2)*CovVec[t] #-Covariance
      
      
      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma_squared[x,t])
    }
  }
  
})

OwnMod_Repara_Diff_OC_V2_Joint <- nimbleCode({
  
  # Differenced Liu,Li Model with Own Constraints (no QR)
  #
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  for(c in 1:N_Country){
    sigma_time[c] ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
    drift[c] ~ dnorm(mean = 0, sd = 2)
    p[c] ~ dbeta(shape1 = 1,shape2 = 20) #joint prior
    a[c] ~ dbeta(shape1 = 1, shape2 = 5) 
  }
  
  
  ##### JUMP PARAMETERS #################
  for(c in 1:N_Country){
    R_Array[1:(N_Year+1),
            1:(N_Year+1),c] <- RTildeMat(Ntime = N_Year,a = a[c])
  }
  
  N_t_Mat[1:2,1:N_Country] <- 0 #Corner Constraint N_t[1:2] <- 0
  for (t in 3:(N_Year+1)){ #one year longer
    for(c in 1:N_Country){
      N_t_Mat[t,c] ~ dbern(p[c])
    }
  }
  
  #Global prior
  muY_Glob ~ T(dnorm(mean=0, sd=2), min = 0,)
  sdY_Glob ~ T(dnorm(mean=0, sd=1), min = 0,) #truncated normal
  # Tau_Glob ~ T(dnorm(mean=0, sd=0.5), min = 0,)
  # TauSD_Glob ~T(dnorm(mean=0, sd=0.5), min = 0,)
  
  for(c in 1:N_Country){
    muY[c] ~ T(dnorm(mean=muY_Glob, sd=0.5), min = 0,)
    sdY[c] ~ T(dnorm(mean=sdY_Glob, sd=0.5), min = 0,) #truncated normal
  }
  
  #other time effects
  alphaDirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  
  for(c in 1:N_Country){
    beta[1:N_AgeGroups,c] ~ ddirch(alpha=alphaDirich[1:N_AgeGroups])
    betaJump[1:N_AgeGroups,c] ~ ddirch(alpha=alphaDirich[1:N_AgeGroups])
  }
  
  
  # ###### AGE PARAMETERS########
  # for(x in 1:N_AgeGroups){
  #   #alpha[x]~ dnorm(0, sd=5)
  #   b1[x]~ dgamma(shape = 1, rate = 1)
  #   b2[x]~ dgamma(shape = 1, rate = 1)
  # }
  # 
  # # Dirichlet Distribution is standardized Gamma Dist
  # beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  # beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  # 
  # betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  # betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  #
  # 
  for(c in 1:N_Country){
    sigma_eps[c] ~ dunif(0,2)
  }
  for(c in 1:N_Country){
    for(t in 1:N_Year){
      CovMat[t,c] <- CovZMat(N_t = N_t_Mat[1:(N_Year+1),c],
                             a =  a[c],
                             t =  t,
                             sigma2 = pow(sdY[c],2))
    }
  }
  for(c in 1:N_Country){
    for(x in 1:N_AgeGroups){
      for(t in 1:N_Year){
        mu[x,t,c] <- #alpha[x]+
          beta[x,c]*drift[c]+
          betaJump[x,c]*muY[c]*(
            inprod(R_Array[t+1,1:(N_Year+1),c],N_t_Mat[1:(N_Year+1),c])- #J_t
              inprod(R_Array[t,1:(N_Year+1),c],N_t_Mat[1:(N_Year+1),c])    #J_t-1
          )
        
        sigma_squared[x,t,c] <- pow(sigma_eps[c],2)  + #sigma_eps
          pow(sigma_time[c]*beta[x,c],2)+ #sigma_time
          pow(betaJump[x,c]*sdY*            #betaJ*sdY*J_t
                inprod(R_Array[t+1,1:(N_Year+1),c],
                       N_t_Mat[1:(N_Year+1)],c),2)+  
          pow(betaJump[x,c]*sdY[c]*                  #betaJ*sdY*J_t-1
                inprod(R_Array[t,1:(N_Year+1),c],
                       N_t_Mat[1:(N_Year+1)],c),2) -  
          2*pow(betaJump[x,c],2)*CovMat[t,c] #-Covariance
        
        
        ZMat[x,t,c] ~ dnorm(mu[x,t,c],var=sigma_squared[x,t,c])
      }
    }
  }
})
#own Model QR beta
OwnMod_Repara_Diff_QR <- nimbleCode({
  
  # Differenced Liu,Li Model with Own Constraints (no QR)
  #
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  drift ~ dunif(min = -2,max = 0) 
  
  N_t[1] <- 0 #Corner Constraint 
  R[1:(N_Year+1),
    1:(N_Year+1)] <- RTildeMat(Ntime = N_Year,a = a)
  
  for (t in 2:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
  }
  p ~ dbeta(shape1 = 1,shape2 = 25) 
  #muY ~ dnorm(mean=1, sd=1) #weakly informative
  #muY ~ dnorm(mean = 0, sd = 4)
  muY ~ T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1,shape2 = 5) # low values are preferred 
  
  ###### AGE PARAMETERS########
  # for(x in 1:N_AgeGroups){
  #   #alpha[x]~ dnorm(0, sd=5)
  #   b1[x]~ dgamma(shape = 1, rate = 1)
  #   b2[x]~ dgamma(shape = 1, rate = 1)
  # }
  # 
  # # Dirichlet Distribution is standardized Gamma Dist
  # beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  # beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  # 
  # betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  # betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  for(x in 1:N_AgeGroups){
    #alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dnorm(0, sd=5)
    b2[x]~ dnorm(0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  # 
  sigma_eps ~ dunif(0,2)
  
  for(t in 1:N_Year){
    CovVec[t] <- CovZMat(N_t = N_t[1:(N_Year+1)],
                         a =  a,
                         t =  t,
                         sigma2 = pow(sdY,2))
  }
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      
      mu[x,t] <- #alpha[x]+
        QMat[x,1]*k[t]+
        QMat[x,2]*muY*(
          inprod(R[t+1,1:(N_Year+1)],N_t[1:(N_Year+1)])- #J_t+1
            inprod(R[t,1:(N_Year+1)],N_t[1:(N_Year+1)])    #J_t
        )
      
      sigma[x,t] <- pow(sigma_eps,2)  +
        pow(QMat[x,2]*sdY*
              inprod(R[t+1,1:(N_Year+1)],
                     N_t[1:(N_Year+1)]),2)+
        pow(QMat[x,2]*sdY*
              inprod(R[t,1:(N_Year+1)],
                     N_t[1:(N_Year+1)]),2)-  #Var second part 
        2*pow(QMat[x,2],2)*CovVec[t] #-2 Covariance
      
      
      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma[x,t])
    }
  }
  
})

OwnMod_Repara_Diff_QR_V2 <- nimbleCode({
  
  # Differenced Liu,Li Model with Own Constraints (no QR)
  #
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  # for(t in 1:N_Year){
  #   k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  # }
  sigma_time ~  T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  drift ~ dunif(min = -2,max = 0) 
  
  N_t[1:2] <- 0 #Corner Constraint 
  R[1:(N_Year+1),
    1:(N_Year+1)] <- RTildeMat(Ntime = N_Year,a = a)
  
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
  }
  p ~ dbeta(shape1 = 1,shape2 = 20) 
  
  muY ~ dnorm(mean=1, sd=1) #weakly informative
  #muY ~ dnorm(mean = 0, sd = 4)
  #muY ~ dchisq(df=2)
  sdY ~ T(dnorm(mean=0, sd=5), min = 0, ) #truncated normal
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1,shape2 = 10) # low values are preferred 
  
  ###### AGE PARAMETERS########
  # for(x in 1:N_AgeGroups){
  #   #alpha[x]~ dnorm(0, sd=5)
  #   b1[x]~ dgamma(shape = 1, rate = 1)
  #   b2[x]~ dgamma(shape = 1, rate = 1)
  # }
  # 
  # # Dirichlet Distribution is standardized Gamma Dist
  # beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  # beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  # 
  # betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  # betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  for(x in 1:N_AgeGroups){
    #alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dnorm(0, sd=5)
    b2[x]~ dnorm(0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  # 
  sigma_eps ~ dunif(0,2)
  
  for(t in 1:N_Year){
    CovVec[t] <- CovZMat(N_t = N_t[1:(N_Year+1)],
                         a =  a,
                         t =  t,
                         sigma2 = pow(sdY,2))
  }
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      
      mu[x,t] <- #alpha[x]+
        QMat[x,1]*drift+
        QMat[x,2]*muY*(
          inprod(R[t+1,1:(N_Year+1)],N_t[1:(N_Year+1)])- #J_t+1
            inprod(R[t,1:(N_Year+1)],N_t[1:(N_Year+1)])    #J_t
        )
      
      sigma[x,t] <- pow(sigma_eps,2)  + 
        pow(sigma_time*QMat[x,2],2)+
        pow(QMat[x,2]*sdY*
              inprod(R[t+1,1:(N_Year+1)],
                     N_t[1:(N_Year+1)]),2)+
        pow(QMat[x,2]*sdY*
              inprod(R[t,1:(N_Year+1)],
                     N_t[1:(N_Year+1)]),2)-  #Var second part 
        2*pow(QMat[x,2],2)*CovVec[t] #-2 Covariance
      
      
      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma[x,t])
    }
  }
  
})



###### ALL MODELS ##############################################################
LiuLi_Repara_Diff_OC_Test <- nimbleCode({
  
  # Differenced Liu,Li Model with Own Constraints (no QR)
  #
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  #sigma_time ~ T(dnorm(mean=0, sd=0.5), min = 0,) #truncated normal
  sigma_time <- 0.1
  #drift ~ dunif(min = -2,max = 0) 
  drift ~ dnorm(mean = 0, sd = 2)
  
  N_t[1:2] <- 0 #Corner Constraint 
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
  }
  
  p ~ dbeta(shape1 = 1,shape2 = 20) 
  muY ~ T(dnorm(mean=0, sd=5), min = 0,)
  VarY ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
  sdY <- abs(sqrt(VarY))
  #VarY <- sdY^2
  #sdY <-  0.5
  #Works with both parameterization on real data
  #a ~ dbeta(shape1 = 1,shape2 = 4) # low values are preferred 
  
  ###### AGE PARAMETERS########
  alphaN[1:N_AgeGroups] <- c(8,8,8,3,3,2,1,rep(1,3))
  alphaJump[1:N_AgeGroups] <- c(rep(1,3),3,7,7,3,rep(3,3))
  
  # for(x in 1:N_AgeGroups){
  #   #alpha[x]~ dnorm(0, sd=5)
  #   b1[x] ~ dgamma(shape = alphaN[x], rate = 1)
  #   b2[x] ~ dgamma(shape = alphaJump[x], rate = 1)
  #   #b2[x]~ dgamma(shape = 1, rate = 1)
  # }
  # 
  #Dirichlet Distribution is standardized Gamma Dist
  # beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  # beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  # 
  # betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  # betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  # 
  
  #a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  
  #
  beta[1:N_AgeGroups] ~ ddirch(alpha=alphaN[1:N_AgeGroups])
  betaJump[1:N_AgeGroups] ~ ddirch(alpha=alphaJump[1:N_AgeGroups])
  
  sigma_eps ~  T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      
      mu[x,t] <- #alpha[x]+
        beta[x]*k[t]+
        betaJump[x]*muY*(N_t[t+1]-N_t[t])
      
      sigma_squared[x,t] <- pow(sigma_eps,2)  +
        pow_int(betaJump[x],2)*VarY*(N_t[t+1]+N_t[t])
      
      
      
      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma_squared[x,t])
    }
  }
  
})
LiuLi_Repara_Diff_OC_V2_Test<- nimbleCode({
  
  # Differenced Liu,Li Model with Own Constraints (no QR)
  #
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  # for(t in 1:N_Year){
  #   k[t] ~ dnorm(mean = 0,sd=1) #State space model formualtion
  # }
  sigma_time ~  T(dnorm(mean=0, sd=2), min = 0,)
  drift ~ dnorm(0,sd=2) 
  
  #Corner Constraint
  N_t[1:2] <- 0
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
  }
  p ~ dbeta(shape1 = 1,shape2 = 20) 
  muY ~ T(dnorm(mean=0, sd=4), min = 0,) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0,) #truncated normal
  #sdY <- abs(sqrt(VarY))
  
  #Works with both parameterization on real data
  
  ###### AGE PARAMETERS########
  for(x in 1:N_AgeGroups){
    #alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dgamma(shape = 1, rate = 1)
    b2[x]~ dgamma(shape = 1, rate = 1)
  }
  
  
  #Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  
  # a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  # 
  # beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  #QR Decomposition to get orthogonal betas
  # BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  # BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  # 
  # QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  # 
  sigma_eps ~ dunif(0,2)
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      
      mu[x,t] <- #alpha[x]+
        beta[x]*drift+ #non centered
        betaJump[x]*muY*(N_t[t+1]-N_t[t])
      
      sigma_squared[x,t] <- pow(sigma_eps,2)+
        pow(beta[x]*sigma_time,2)+
        pow(betaJump[x]*sdY,2)*(N_t[t+1]+N_t[t])
      #pow(betaJump[x]*N_t[t+1]*sdY,2) + #Var first part
      #pow(betaJump[x]*N_t[t]*sdY,2)     #Var second part
      
      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma_squared[x,t])
    }
  }
  
})
OwnMod_Repara_Diff_OC_Test <- nimbleCode({
  
  # Differenced Liu,Li Model with Own Constraints (no QR)
  #
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  #drift ~ dunif(min = -2,max = 0) 
  drift ~ dnorm(mean = 0 ,sd = 2)
  
  N_t[1:2] <- 0 #Corner Constraint 
  R[1:(N_Year+1),
    1:(N_Year+1)] <- RTildeMat(Ntime = N_Year,a = a)
  
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
  }
  p ~ dbeta(shape1 = 1,shape2 = 20) 
  muY ~ T(dnorm(mean=0, sd=5), min = 0, ) #truncated normal
  sdY ~ T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 =1, shape2 = 5)
  #a ~ dbeta(shape1 = 1,shape2 = 1) # low values are preferred 
  
  ###### AGE PARAMETERS########
  alphaN[1:N_AgeGroups] <- c(3,3,3,1,1,1,1,rep(1,3))
  alphaJump[1:N_AgeGroups] <- c(rep(1,3),1,3,3,1,rep(1,3))
  # 
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = alphaN[x], rate = 1)
    b2[x]~ dgamma(shape = alphaJump[x], rate = 1)
    #b2[x]~ dgamma(shape = alphaJump[x], rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  #
  # a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  # 
  # beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  sigma_eps ~  T(dnorm(mean=0, sd=2), min = 0, ) #truncated normal
  
  for(t in 1:N_Year){
    CovVec[t] <- CovZMat(N_t = N_t[1:(N_Year+1)],
                         a =  a,
                         t =  t,
                         sigma2 = pow(sdY,2))
  }
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      
      mu[x,t] <- #alpha[x]+
        beta[x]*k[t]+
        betaJump[x]*muY*(
          inprod(R[t+1,1:(N_Year+1)],N_t[1:(N_Year+1)])- #J_t
            inprod(R[t,1:(N_Year+1)],N_t[1:(N_Year+1)])  #J_t-1
        )
      
      sigma_squared[x,t] <- pow(sigma_eps,2)  + #sigma_eps
        pow(betaJump[x]*sdY*            #betaJ*sdY*J_t
              inprod(R[t+1,1:(N_Year+1)],
                     N_t[1:(N_Year+1)]),2)+  
        pow(betaJump[x]*sdY*                  #betaJ*sdY*J_t-1
              inprod(R[t,1:(N_Year+1)],
                     N_t[1:(N_Year+1)]),2) -  
        2*pow(betaJump[x],2)*CovVec[t] #-Covariance
      
      
      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma_squared[x,t])
    }
  }
  
})
##Lee Carter Models
#LC Model with Sum to Zero constraint on Kappa
LC_SumToZero <- nimbleCode({
  
  # Lee Carter Model
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  
  k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  kappa[1:N_Year] <- SumToZero(k[1:N_Year]) #SumToZero Constraint 
  
  sigma_time ~  dunif(0,2) #prior on time sd
  drift ~ dnorm(0,sd=2) #prior on drift term
  
  #Paramterisation for Age Effects (age-specific Intercept)
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
  }
  
  #Dirichlet Prior
  a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  #Prior on Overdispersion Paramter
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+beta[age[i]]*kappa[year[i]]
    log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) # Effect
    lambda[i] <- eta[i]*Offset[i]
    y[i] ~ dpois(lambda[i])
  }
})

LC_SumToZero_AltDirichBeta <- nimbleCode({
  
  # Lee Carter Model
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  
  k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  kappa[1:N_Year] <- SumToZero(k[1:N_Year]) #SumToZero Constraint 
  
  sigma_time ~  dunif(0,2) #prior on time sd
  drift ~ dnorm(0,sd=2) #prior on drift term
  
  #Paramterisation for Age Effects (age-specific Intercept)
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
  }
  
  #Alternative Representation of Dirichlet prior
  for(x in 1:N_AgeGroups){
    b[x] ~ dgamma(shape = 1, rate = 1)
  }
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b[1:(N_AgeGroups-1)]/sum(b[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  #Prior on Overdispersion Paramter
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+beta[age[i]]*kappa[year[i]]
    log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) # Effect
    lambda[i] <- eta[i]*Offset[i]
    y[i] ~ dpois(lambda[i])
  }
})

LC_SumToZero_MultNormBeta <- nimbleCode({
  
  # Lee Carter Model
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  
  k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  kappa[1:N_Year] <- SumToZero(k[1:N_Year]) #SumToZero Constraint 
  
  sigma_time ~  dunif(0,2) #prior on time sd
  drift ~ dnorm(0,sd=2) #prior on drift term
  
  #Paramterisation for Age Effects (age-specific Intercept)
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
  }
  
  #Beta Distributed as Multivariate Normal (see Wong et. al )
  meanVec[1:(N_AgeGroups-1)] <- rep((1/N_AgeGroups),N_AgeGroups-1)
  NormMat[1:(N_AgeGroups-1),
          1:(N_AgeGroups-1)] <- CovMat(A=N_AgeGroups,
                                       sigma = sigma_age)
  
  b[1:(N_AgeGroups-1)] ~ dmnorm(mean= meanVec[1:(N_AgeGroups-1)],
                                cov = NormMat[1:(N_AgeGroups-1),
                                              1:(N_AgeGroups-1)])
  
  beta[1] <-  1 - sum(b[1:(N_AgeGroups-1)])
  beta[2:N_AgeGroups] <-  b[1:(N_AgeGroups-1)]

  
  #Prior on Overdispersion Paramter
  sigma_eps ~ dunif(0,2)
  sigma_age ~ dunif(0,4)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+beta[age[i]]*kappa[year[i]]
    log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) # Effect
    lambda[i] <- eta[i]*Offset[i]
    y[i] ~ dpois(lambda[i])
  }
})

LC_AltDirichBeta_Diff <- nimbleCode({
  
  # Lee Carter Model
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  for(t in 1:N_Year){
    k[t] ~ dnorm(drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) #prior on time sd
  drift ~ dnorm(0,sd=2) #prior on drift term
  
  #Alternative Representation of Dirichlet prior
  for(x in 1:N_AgeGroups){
    b[x] ~ dgamma(shape = 1, rate = 1)
  }
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b[1:(N_AgeGroups-1)]/sum(b[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  #Prior on Overdispersion Paramter
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

#LC Model with Corner Constraint on Kappa
LC_Corner <- nimbleCode({
  
  # Lee Carter Model
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  kappa[1] <- 0 #Corner Constraint
  for(t in 2:N_Year){
    kappa[t] ~ dnorm(kappa[t-1]+drift,sd=sigma_time) #Random Walk with Drift
  }
  
  sigma_time ~  dunif(0,2) #Uniform Prior on Variance Parameters
  drift ~ dnorm(0,sd=2) 
  
  #Paramterisation for Age Effects
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
  }
  #Dirichlet Prior
  a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+beta[age[i]]*kappa[year[i]]
    eta[i] ~ dnorm(mu[i],sd=sigma_eps) # Effect
    lambda[i] <- calculateLambda(eta=eta[i], #may decrease building Time
                                 Offset = Offset[i])
    y[i] ~ dpois(lambda[i])
  }
})

## Liu Li Model
LiuLiModel_SumZero <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1)
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  kappa[1:N_Year] <- SumToZero(k[1:N_Year]) #SumToZero Constraint 
  N_t[1] <- 0
  Y_t[1] <- 0
  #other time effects
  for (t in 2:N_Year){
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t] #Jump effect (deterministic)
  }
  
  p ~ dbeta(shape1 = 2, shape2 = 5)
  muY ~ dnorm(mean = 0, sd = 10)
  sdY ~ dunif(min = 0,  max = 2)
  
  #Paramterisation for Age Effects
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
  }
  #Dirichlet Prior
  a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+beta[age[i]]*kappa[year[i]]+ #LC Model
              betaJump[age[i]]*J[year[i]] #Jump Part
    log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) # Effect
    lambda[i] <- eta[i]*Offset[i]
    y[i] ~ dpois(lambda[i])
  }
})

LiuLiModel_SumZero_QR <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1)
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dunif(min=-5, max=0) 
  k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  kbar <-mean(k[1:N_Year])  #mean of k vector
  # Sum to zero constraint
  for(t in 1:N_Year){
    kappa[t] <- k[t]-kbar
  }
  
  #other time effects
  for (t in 1:N_Year){
    J[t] <- N_t[t]*Y_t[t] #Jump effect (deterministic)
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
  }
  
  p ~ dbeta(shape1 = 1,shape2 = 6) 
  #muY ~ dgamma(shape = 2,rate = 1) #weakly informative
  muY ~ dnorm(mean=0, sd=4) #weakly informative
  sdY ~ dgamma(shape = 1,rate = 1) #weakly informative
  
  #Paramterisation for Age Effects
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+QMat[age[i],1]*kappa[year[i]]+ #LC Model
                          QMat[age[i],2]*J[year[i]] #Jump Effect
    log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) 
    lambda[i] <- eta[i]*Offset[i]
    y[i] ~ dpois(lambda[i])
  }
})

LiuLiModel_Corner_QR <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1)
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Priors on Time Parameters
  k[1] ~ dnorm(k0,sd=sigma_time)
  k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  for(t in 2:(N_Year-1)){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  k[N_Year] <- 0 # last year effect is zero 
  
  sigma_time ~  dunif(0,2) 
  drift ~ dunif(min=-5, max=0) 
  
  # Jump Effect
  J[1] <- 0 #Corner Constraint for Jump Part
  Y_t[1] <- 0 # First Jump effect is zero 
  N_t[1] <- 0 # First Jump effect is zero 
  
  #other time effects
  for (t in 2:N_Year){
    J[t] <- N_t[t]*Y_t[t] #Jump effect (deterministic)
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
  }
  
  p ~ dbeta(shape1 = 1,shape2 = 6) 
  #muY ~ dgamma(shape = 2,rate = 1) #weakly informative
  
  muY ~ dnorm(mean=0,   sd=5) #weakly informative
  sdY ~ dunif(min=0, max=3) #weakly informative
  
  #Paramterisation for Age Effects
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+QMat[age[i],1]*k[year[i]]+ #LC Model
      QMat[age[i],2]*J[year[i]] #Jump Effect
    log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) 
    lambda[i] <- eta[i]*Offset[i]
    y[i] ~ dpois(lambda[i])
  }
})

#### DIfferenced Liu,Li Models##################################################

LiuLiModel_SumZero_Diff <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  
  #alternative for kappa
  # k[1]~dnorm(k0, sd=sigma_time)
  # for(t in 2:N_Year){
  #   k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  # }
  # sigma_time ~  dunif(0,2) 
  # drift ~ dnorm(0,sd=2) 
  # k0 ~ dnorm(0,5)
  # 
  # kappa[1:N_Year] <- SumToZero(k[1:N_Year])
  # kDiff[1:(N_Year-1)] <- kappa[2:N_Year]-kappa[1:(N_Year-1)]
  # 
  #other time effects
  N_t[1] <- 0
  Y_t[1] <- 0
  for (t in 2:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
  }

  p ~ dbeta(shape1 = 2, shape2 = 5)
  muY ~ dnorm(mean = 0, sd = 10)
  sdY ~ dunif(min = 0,  max = 2)

 
  #Dirichlet Prior
  a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  sigma_eps ~ dunif(0,2)
  

for(x in 1:N_AgeGroups){
    for(l in 1:(N_Year+1)){ #one year longer due to N0Y0 in the Model
      J[x,l] <- betaJump[x]*Y_t[l]
    }
  }
  #Putting all Parameters together
for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
        N_t[(t+1)]*J[x,(t+1)]-
        N_t[t]*J[x,t]
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

LiuLiModel_SumZero_Diff_Restated <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  
  #other time effects
  N_t[1] <- 0
  Y_t[1] <- 0
  J[1] <- 0
  for (t in 2:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape1 = 2, shape2 = 5)
  muY ~ dnorm(mean = 0, sd = 5)
  sdY ~ dunif(min = 0,  max = 2)
  
  #Dirichlet Prior
  a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  a_dirichJump[1:N_AgeGroups] <- c(rep(1,3),3,3,3,rep(1,4)) # more informative Hyperprior
  
  beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirichJump[1:N_AgeGroups])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
                 betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

LiuLiModel_SumZero_Diff_Restated_HMC <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  #other time effects
  N_t[1] <- 0
  Y_t[1] <- 0
  J[1] <- 0
  for (t in 2:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape1 = 2, shape2 = 5)
  muY ~ dnorm(mean = 0, sd = 5)
  sdY ~ dunif(min = 0,  max = 2)
  
  #Dirichlet Prior
  # a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  # beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # #betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  for(i in 1:N_AgeGroups){
    b1[i] ~ dnorm(mean=0,sd=5)
    b2[i] ~ dnorm(mean=0,sd=5)
  }
 
  beta[1:N_AgeGroups] <- b1[1:N_AgeGroups]/sum(b1[1:N_AgeGroups])
  betaJump[1:N_AgeGroups] <- b2[1:N_AgeGroups]/sum(b2[1:N_AgeGroups])
  
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
        betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

LiuLiModel_SumZero_Diff_Restated_NonDirichletBeta <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  
  #other time effects
  N_t[1] <- 0
  Y_t[1] <- 0
  J[1] <- 0
  for (t in 2:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape1 = 2, shape2 = 5)
  muY ~ dnorm(mean = 0, sd = 5)
  sdY ~ dunif(min = 0,  max = 2)
  
  #Dirichlet Prior
  alphaJump[1:N_AgeGroups] <- c(rep(1,3),3,3,3,rep(1,4))
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
    b2[x]~ dgamma(shape = alphaJump[x], rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
                betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

LiuLiModel_SumZero_Diff_Restated_MultNormBeta <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  
  #other time effects
  N_t[1] <- 0
  Y_t[1] <- 0
  J[1] <- 0
  for (t in 2:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  
  p ~ dbeta(shape1 = 2, shape2 = 5)
  muY ~ dnorm(mean = 0, sd = 5)
  sdY ~ dunif(min = 0,  max = 2)
  
  
  #Mult Normal Prior for Betas
  meanVec[1:((N_AgeGroups-1)*2)] <- rep((1/N_AgeGroups),(N_AgeGroups-1)*2)
  
  NormMat[1:((N_AgeGroups-1)*2), 
          1:((N_AgeGroups-1)*2)] <- CovMatMult(A=N_AgeGroups,
                                               sigma = sigma_age)
  
  bTot[1:((N_AgeGroups-1)*2)] ~ dmnorm(mean = meanVec[1:((N_AgeGroups-1)*2)],
                                       cov = NormMat[1:((N_AgeGroups-1)*2), 
                                                     1:((N_AgeGroups-1)*2)])
  
  beta[1] <-  1 - sum(bTot[1:(N_AgeGroups-1)])
  beta[2:N_AgeGroups] <-  bTot[1:(N_AgeGroups-1)]
  
  betaJump[1] <- 1 - sum(bTot[N_AgeGroups:(2*N_AgeGroups-2)])
  betaJump[2:N_AgeGroups] <- bTot[N_AgeGroups:(2*N_AgeGroups-2)]
  
  sigma_eps ~ dunif(0,2)
  sigma_age ~ dunif(0,2)
  
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
                 betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

LiuLiModel_SumZero_Diff_Restated_MultNormBetaS <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  
  #other time effects
  N_t[1] <- 0
  Y_t[1] <- 0
  J[1] <- 0
  for (t in 2:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = 1, sd = sdY)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  
  p ~ dbeta(shape1 = 2, shape2 = 5)
  #muY ~ dnorm(mean = 0, sd = 5)
  sdY ~ dunif(min = 0,  max = 2)
  
  
  #Mult Normal Prior for Betas
  meanVec[1:(N_AgeGroups-1)] <- rep((1/N_AgeGroups),N_AgeGroups-1)
  NormMat[1:(N_AgeGroups-1),
          1:(N_AgeGroups-1)] <- CovMat(A=N_AgeGroups,
                                       sigma = sigma_age)
  
  b1[1:(N_AgeGroups-1)] ~ dmnorm(mean= meanVec[1:(N_AgeGroups-1)],
                                cov = NormMat[1:(N_AgeGroups-1),
                                              1:(N_AgeGroups-1)])
  
  b2[1:(N_AgeGroups-1)] ~ dmnorm(mean= meanVec[1:(N_AgeGroups-1)],
                                 cov = NormMat[1:(N_AgeGroups-1),
                                               1:(N_AgeGroups-1)])
  
  beta[1] <-  1 - sum(b1[1:(N_AgeGroups-1)])
  beta[2:N_AgeGroups] <-  b1[1:(N_AgeGroups-1)]
  
  betaJump[1] <- 1 - sum(b2[1:(N_AgeGroups-1)])
  betaJump[2:N_AgeGroups] <- b2[1:(N_AgeGroups-1)]
  
  sigma_eps ~ dunif(0,2)
  sigma_age ~ dunif(0,2)
  
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
                 betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

LiuLiModel_SumZero_Diff_Restated_OrthBeta <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Normal Prior Prior) +
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formulation
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  
  #other time effects
  N_t[1] <- 0
  Y_t[1] <- 0
  J[1] <- 0
  for (t in 2:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = 1, sd = sdY)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  Jtilde[1:N_Year] <- J[2:(N_Year+1)]-J[1:N_Year]
  
  p ~ dbeta(shape1 = 2, shape2 = 5)
  #muY ~ dnorm(mean = 0, sd = 5)
  sdY ~ dunif(min = 0,  max = 2)
  
  #Normal Prior
  for(x in 1:N_AgeGroups){
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #Get Normed Betas (see Hunt, Blake )
  a <- sum(b1[1:N_AgeGroups])
  #b <- sum(b2[1:N_AgeGroups])-1
  b <- (sum(b2[1:N_AgeGroups]*b1[1:N_AgeGroups])*a)/
       inprod(b1[1:N_AgeGroups],b1[1:N_AgeGroups])

  beta[1:N_AgeGroups] <- b1[1:N_AgeGroups]/a
  betaJump[1:N_AgeGroups] <- b2[1:N_AgeGroups] - (b/a)*b1[1:N_AgeGroups]
  kappa[1:N_Year] <- a*k[1:N_Year]+b*Jtilde[1:N_Year]
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*kappa[t]+
                 betaJump[x]*(Jtilde[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

LiuLiModel_CornerOwn <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1)
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  k[1] <- 1
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  
  
  N_t[1] <- 0
  Y_t[1] <- 0
  J[1] <- 0
  
  #other time effects
  for (t in 2:N_Year){
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t] #Jump effect (deterministic)
  }
  
  p ~ dbeta(shape1 = 2, shape2 = 5)
  muY ~ T(dnorm(mean=1, sd=4), min = 0, max = 1000) #truncated normal
  sdY ~ dunif(min = 0,  max = 2)
  
  #Paramterisation for Age Effects
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
    b1[x] ~ dnorm(0, sd=5)
    b2[x] ~ dnorm(0,sd=5)
  }
  #Normalization Scheme 
  beta[1:N_AgeGroups] <- b1[1:N_AgeGroups]/sqrt(inprod(b1[1:N_AgeGroups],
                                                       b1[1:N_AgeGroups]))
  
  betaJump[1:N_AgeGroups]<- b2[1:N_AgeGroups]/sqrt(inprod(b2[1:N_AgeGroups],
                                                          b2[1:N_AgeGroups]))
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+beta[age[i]]*k[year[i]]+ #LC Model
                           betaJump[age[i]]*J[year[i]] #Jump Part
    
    y[i] ~ dnorm(mu[i], sd = sigma_eps)
    # log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) # Effect
    # lambda[i] <- eta[i]*Offset[i]
    # y[i] ~ dpois(lambda[i])
  }
})


##### OWN MODEL PROPOSITION ###################################################

OwnMod_Repara_QR <- nimbleCode({
  
  # Own Model, non differenced with QR on betas

  ##### TIME PARAMETERS ###################
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formulation
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dunif(min=-5, max=0) 
  k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  # Sum to zero constraint
  kappa[1:N_Year] <- k[1:N_Year]-mean(k[1:N_Year]) 
  
  # sigma_time ~  dunif(min = 0, max = 2) 
  # drift ~ dunif(min=-5, max=0) 
  # k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  ### Other Time Parameters ###
  N_t[1:2] <- 0 #Corner Constraint 
  R[1:N_Year,
    1:N_Year] <- RMatrix(Ntime = N_Year,a = a)
  
  for (t in 3:N_Year){ #one year longer
    N_t[t] ~ dbern(p)
  }
  
  p ~ dbeta(shape1 = 1,shape2 = 25) 

  muY ~ dchisq(df = 1)
  sdY ~ dgamma(1,1)
  #sdY ~ T(dnorm(mean=0, sd=5), min = 0, max=100) #truncated normal
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1,shape2 = 5) # low values are preferred 
  
  ###### AGE PARAMETERS########
  for(x in 1:N_AgeGroups){
    alpha[x] ~ dnorm(0, sd=5)
    b1[x]~ dnorm(0, sd=5)
    b2[x]~ dnorm(0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
   
  sigma_eps ~ dunif(min = 0, max = 2)
  
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- alpha[x]+
        QMat[x,1]*kappa[t]+
        QMat[x,2]*muY*(
          inprod(R[t,1:N_Year],N_t[1:N_Year])  #J_t
        )
      
      
      sigma[x,t] <- pow(sigma_eps,2)+
        pow(QMat[x,2]*sdY*
        inprod(R[t,1:N_Year],N_t[1:N_Year]),
        2)

      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma[x,t])
      # log(ZMat[x,t]) ~ dnorm(mu[x,t],var=sigma[x,t])
      # lambda[x,t] <- ZMat[x,t]*Offset[x,t]
      # y[x,t] ~ dpois(lambda[x,t])
    }
  }
})


## LC Jump Model (own Model) 
LCJumpOwn_SumZero_QR <- nimbleCode({
  
  # Lee Carter Model with Jumps
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  # a damping parameter  
  
  #Prior Parameterisation
  #Priors on Time Parameters
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dunif(min=-5, max=0) 
  k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  # Sum to zero constraint
  kappa[1:N_Year] <- k[1:N_Year]-mean(k[1:N_Year]) 
  
  #Other Time Parameters with corner Constraint
  J[1:2] <- 0
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  #other time effects
  for (t in 3:N_Year){
    J[t] <- J[t-1]*a+N_t[t]*Y_t[t] #Jump effect (deterministic)
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
  }
  
  p ~ dbeta(shape1 = 1,shape2 = 20) 
  muY ~ dnorm(mean=0, sd=4) #weakly informative
  sdY ~ dgamma(shape = 1,rate = 1) #weakly informative
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1,shape2 = 6) # low values are preferred 
  
  #Paramterisation for Age Effects
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+QMat[age[i],1]*kappa[year[i]]+ #LC Model
                           QMat[age[i],2]*J[year[i]] #Jump Effect
    log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) 
    lambda[i] <- eta[i]*Offset[i]
    y[i] ~ dpois(lambda[i])
  }
})

LCJumpOwn_CornerQR <- nimbleCode({
  
  # Lee Carter Model with Jumps
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  # a damping parameter  
  
  #Prior Parameterisation
  #Priors on Time Parameters
  kappa[1] <- 0 #Corner Constraint Time Index
  
  for(t in 2:(N_Year-1)){
    kappa[t] ~ dnorm(kappa[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  
  kappa[N_Year] <- 0
  sigma_time ~  dunif(0,2) 
  drift ~ dunif(min=-5, max=0)
  
  N_t[1] <- 0
  Y_t[1] <- 0
  J[1] <- 0
  for (t in 2:(N_Year)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- J[t-1]*a+N_t[t]*Y_t[t] #Jump effect (deterministic)
  }
  
  p ~ dbeta(shape1 = 1,shape2 = 3) 
  muY ~ dgamma(shape = 2,rate = 1) #weakly informative
  sdY ~ dgamma(shape = 1,rate = 1) #weakly informative
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1,shape2 = 3) # low values are preferred 
  
  #Paramterisation for Age Effects
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+QMat[age[i],1]*kappa[year[i]]+ #LC Model
                           QMat[age[i],2]*J[year[i]] #Jump Effect
    log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) 
    lambda[i] <- eta[i]*Offset[i]
    y[i] ~ dpois(lambda[i])
  }
})

#### OWN MODEL PROPOSITION WITH LOG IMPROVEMENT RATES (DIFFERENCED)#############

LCJumpOwn_Diff <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  
  #alternative for kappa
  # k[1]~dnorm(k0, sd=sigma_time)
  # for(t in 2:N_Year){
  #   k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  # }
  # sigma_time ~  dunif(0,2) 
  # drift ~ dnorm(0,sd=2) 
  # k0 ~ dnorm(0,5)
  # 
  # kappa[1:N_Year] <- SumToZero(k[1:N_Year])
  # kDiff[1:(N_Year-1)] <- kappa[2:N_Year]-kappa[1:(N_Year-1)]
  # 
  #other time effects
  N_t[1] <- 0
  Y_t[1] <- 0
  J[1] <- 0
  for (t in 2:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- J[t-1]*a+N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape1 = 2, shape2 = 5)
  muY ~ dnorm(mean = 0, sd = 10)
  sdY ~ dunif(min = 0,  max = 2)
  a ~ dunif(min = 0.01, max = 0.99) 
  
  #Dirichlet Prior
  a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
        betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

LCJumpOwn_MAQ <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  
  #alternative for kappa
  
  #other time effects
  N_t[1:Q] <- 0
  Y_t[1:Q] <- 0
  J[1:Q] <- 0
  #disadvantage of this method. Q must be known before building the model
  for (t in (Q+1):(N_Year+Q)){ #Q years longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    error[t] <- a*N_t[(t-1)]*Y_t[(t-1)]+pow(a,2)*N_t[(t-2)]*Y_t[(t-2)]+
                pow(a,3)*N_t[(t-3)]*Y_t[(t-3)]+pow(a,4)*N_t[(t-4)]*Y_t[(t-4)]+
                pow(a,Q)*N_t[(t-Q)]*Y_t[(t-Q)]
    # #Multiple assignments not allowed
    J[t] <-  N_t[t]*Y_t[t]+error[t]
  }
  
  
  
  p ~ dbeta(shape1 = 1,shape2 = 4) 
  muY ~ dnorm(mean = 0, sd = 10)
  sdY ~ dunif(min = 0,  max = 2)
  a ~ dbeta(shape1 = 1,shape2 = 3) 
  
  #Dirichlet Prior
  a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
        betaJump[x]*(J[t+Q]-J[t+(Q-1)])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

LCJumpOwn_MAQ_QR <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ T(dnorm(mean=0, sd=3), min = -1000, max = 0) #truncated normal
  
  #alternative for kappa
  
  #other time effects
  N_t[1:Q] <- 0
  Y_t[1:Q] <- 0
  J[1:Q] <- 0
  #disadvantage of this method. Q must be known before building the model
  for (t in (Q+1):(N_Year+Q)){ #Q years longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    error[t] <- a*N_t[(t-1)]*Y_t[(t-1)]+pow(a,2)*N_t[(t-2)]*Y_t[(t-2)]+
      pow(a,3)*N_t[(t-3)]*Y_t[(t-3)]+pow(a,4)*N_t[(t-4)]*Y_t[(t-4)]+
      pow(a,Q)*N_t[(t-Q)]*Y_t[(t-Q)]
    # #Multiple assignments not allowed
    J[t] <-  N_t[t]*Y_t[t]+error[t]
  }
  
  
  
  p ~ dbeta(shape1 = 1,shape2 = 3) 
  muY ~ T(dnorm(mean=1, sd=4), min = 0, max = 1000) #truncated normal
  sdY ~ dgamma(shape = 1,rate = 1) #weakly informative 
  a ~ dbeta(shape1 = 1,shape2 = 3) 
  
  #Dirichlet Prior
  for(x in 1:N_AgeGroups){
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- QMat[x,1]*k[t]+
                 QMat[x,2]*(J[t+Q]-J[t+(Q-1)])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

LCJump_QR_NonPoisson <- nimbleCode({
  
  # Own Model Proposition
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Normal Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Us
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  # Autoregressive Parameter
  
  #Prior Parameterisation
  
  #Priors on Time Parameters
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]-drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2)
  drift ~ dgamma(shape = 1, rate = 1) # positive drift turns in to negative drift
  k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  kbar <-mean(k[1:N_Year])  #mean of k vector
  # Sum to zero constraint
  for(t in 1:N_Year){
    kappa[t] <- k[t]-kbar
  }
  
  for(t in 1:(N_Year)){
      N_t[t] ~ dbern(p)
      Y_t[t] ~ dnorm(mean = muY, sd = sdY)
      J[t] <- N_t[t]*Y_t[t] #Jump effect (deterministic)
    }
  
  #Corner Constraint 
  # J[1] <- 0 #starting value AR1
  # N_t[1] <- 0 #starting values
  # Y_t[1] <- 0 #starting Values
  # 
  # for(t in 2:(N_Year)){
  #   N_t[t] ~ dbern(p)
  #   Y_t[t] ~ dnorm(mean = muY, sd = sdY)
  #   J[t] <- J[t-1]*a+N_t[t]*Y_t[t] #Jump effect (deterministic)
  # }
  
  p ~ dbeta(shape1 = 1,shape2 = 5) 
  muY ~ dunif(min = 0, max = 4) #weakly informative
  sdY ~ dgamma(shape = 2,rate = 1) #weakly informative
  
  #Works with both parameterization on real data
  #a ~ dbeta(shape1 = 1,shape2 = 5) # low values are preferred 
  
  #Paramterisation for Age Effects
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+QMat[age[i],1]*kappa[year[i]]+ #LC Model
                           QMat[age[i],2]*J[year[i]] #Jump Effect
    y[i] ~ dnorm(mu[i], sd = sigma_eps)
  }
})

####Model Tests ##############################################################
#Current Model Tests

LiuLi_Exo <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1)
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2)
  
  drift ~ dnorm(0, sd = 2)
  
  k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  kbar <-mean(k[1:N_Year])  #mean of k vector
  # Sum to zero constraint
  for(t in 1:N_Year){
    kappa[t] <- k[t]-kbar
  }
  
  #other time effects
  # for (t in 1:N_Year){
  #   N_t[t] ~ dbern(p)
  #   Y_t[t] ~ dnorm(mean = muY, sd = sdY)
  #   J[t] <- N_t[t]*Y_t[t]
  # }
  
  #Corner Constraint 
  J[1] <- 0 #starting value AR1
  N_t[1] <- 0 #starting values
  Y_t[1] <- 0 #starting Values

  for(t in 2:(N_Year)){
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- J[t-1]*a+N_t[t]*Y_t[t] #Jump effect (deterministic)
  }

  p ~ dbeta(shape1 = 1,shape2 = 6) 
  muY ~ dnorm(mean=0, sd=4) #weakly informative
  sdY ~ dgamma(shape = 1,rate = 1) #weakly informative
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1,shape2 = 4) # low values are preferred 
  
  #Paramterisation for Age Effects
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  # a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  # beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  # 
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]

  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])

  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+QMat[age[i],1]*kappa[year[i]]+ #LC Model
                           QMat[age[i],2]*J[year[i]] #Jump Effect
    y[i] ~ dnorm(mu[i], sd = sigma_eps)
    #log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) 
    #y[i] ~ dpois(lambda = lambda[i])
    
    # log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) 
    # lambda[i] <- eta[i]*Offset[i]
    # y[i] ~ dpois(lambda[i])
  }
})

LiuLi_Exo_Diff <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1)
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  drift ~ dnorm(0,sd=2) 
  
  
  
  # #other time effects
  # for (t in 2:N_Year){
  #   J[t] <- N_t[t]*Y_t[t] #Jump effect (deterministic)
  #   N_t[t] ~ dbern(p)
  #   Y_t[t] ~ dnorm(mean = muY, sd = sdY)
  # }
  
  # p ~ dbeta(shape1 = 1,shape2 = 6) 
  # #muY ~ dgamma(shape = 2,rate = 1) #weakly informative
  # muY ~ dnorm(mean=0, sd=4) #weakly informative
  # sdY ~ dgamma(shape = 1,rate = 1) #weakly informative
  
  #Paramterisation for Age Effects
  # for(x in 1:N_AgeGroups){
  #   alpha[x]~ dnorm(0, sd=5)
  #   # b1[x]~ dnorm(mean = 0, sd=5)
  #   # b2[x]~ dnorm(mean = 0, sd=5)
  # }
  
  a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])
  
  # #QR Decomposition to get orthogonal betas
  # BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  # BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  # 
  # QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  # 
  sigma_eps ~ dunif(0,2)
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+
                 betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

LiuLi_Repara_2 <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1)
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(min = 0, max = 2)
  drift ~ T(dnorm(mean=0, sd=3), min = -1000, max = 0) #truncated normal
  k0 ~ dnorm(mean = 0, sd=5) ## starting value RW 
  
  # kbar <-mean(k[1:N_Year])  #mean of k vector
  # # Sum to zero constraint
  # for(t in 1:N_Year){
  #   kappa[t] <- k[t]-kbar
  # }

  #Corner Constraint 
  N_t[1] <- 0 #starting values
  for(t in 2:(N_Year)){
    N_t[t] ~ dbern(p) # Bernoulli Jump
  }
  p ~ dbeta(shape1 = 1,shape2 = 4) 
  muY ~ dnorm(mean=0, sd=5) #weakly informative
  sdY ~ dunif(min = 0,max = 2) #weakly informative
  # 
  #Works with both parameterization on real data
  #a ~ dbeta(shape1 = 1,shape2 = 4) # low values are preferred 
  
  ###### AGE PARAMETERS########
  for(x in 1:N_AgeGroups){
    #alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(0,2)
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){

      mu[x,t] <- #alpha[x]+
                          QMat[x,1]*k[t]+
                          QMat[x,2]*N_t[t]*muY
      
      sigma[x,t] <- pow(sigma_eps,2) +
                           pow(QMat[x,2]*N_t[t]*sdY,2) #all in Variance

      #sigma[x,t] <- sigma_eps + pow(QMat[x,2]*N_t[t],2)*sdY
      
      
      ZMat[x,t] ~ dnorm(mu[x,t],var=sigma[x,t])
    }
  }
  
})

LiuLi_Repara_2_Poi <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1)
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  
  ##### TIME PARAMETERS ###################
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2)
  drift ~ dnorm(0, sd = 2)
  k0 ~ dnorm(0,sd=5) ## starting value RW 
  
  # kbar <-mean(k[1:N_Year])  #mean of k vector
  # # Sum to zero constraint
  # for(t in 1:N_Year){
  #   kappa[t] <- k[t]-kbar
  # }
  
  #Corner Constraint 
  N_t[1] <- 0 #starting values
  for(t in 2:(N_Year)){
    N_t[t] ~ dbern(p) ## Bernoulli Jump
  }
  p ~ dbeta(shape1 = 1,shape2 = 4) 
  muY ~ dnorm(mean=0, sd=4) #weakly informative
  sdY ~ dgamma(shape = 1,rate = 1) #weakly informative
  
  #Works with both parameterization on real data
  #a ~ dbeta(shape1 = 1,shape2 = 4) # low values are preferred 
  
  ###### AGE PARAMETERS########
  for(x in 1:N_AgeGroups){
    #alpha[x]~ dnorm(0, sd=5)
    b1[x]~ dnorm(mean = 0, sd=5)
    b2[x]~ dnorm(mean = 0, sd=5)
  }
  
  #QR Decomposition to get orthogonal betas
  BMat[1:N_AgeGroups,1] <- b1[1:N_AgeGroups]
  BMat[1:N_AgeGroups,2] <- b2[1:N_AgeGroups]
  
  QMat[1:N_AgeGroups,1:2] <- qr_Nimble(BMat[1:N_AgeGroups,1:2])
  
  sigma_eps ~ dunif(0,2)
  
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      
      mu[x,t] <- #alpha[x]+
        QMat[x,1]*k[t]+
        QMat[x,2]*N_t[t]*muY
      
      sigma[x,t] <- pow(sigma_eps,2) +
        pow(QMat[x,2]*N_t[t]*sdY,2) #all in Variance
      
      #sigma[x,t] <- sigma_eps + pow(QMat[x,2]*N_t[t],2)*sdY
      #sigma[x,t] <- sigma_eps+QMat[x,2]*N_t[t]*sdY #Mag er nicht
      
      log(ZMat[x,t]) ~ dnorm(mu[x,t],var=sigma[x,t])
      lambda[x,t] <- ZMat[x,t]*Offset[x,t]
      YMat[x,t] ~ dpois(lambda = lambda[x,t])
    }
  }
  
})

LCJumpOwn_S <- nimbleCode({
  
  # Own Model on Differenced Rates
  #
  # Time Effect (iid normal)+
  # Age Effect (Dirichlet via Gamma) +
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = drift,sd=sigma_time) #State space model formualtion
  }
  
  sigma_time ~  dunif(0,2) 
  
  
  #alternatively drift must be negative
  drift ~ dnorm(mean=0, sd=3) #truncated normal
  
  #other time effects (with Corner Constraint)
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0
  
  for (t in 3:(N_Year+1)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = 0, sd = 0.25)
    J[t] <- a*J[t-1]+N_t[t]*Y_t[t]
  }
  
  p ~ dbeta(shape1 = 1,shape2 = 4) #low values are preferred
  # muY ~ dnorm(mean=0, sd=4)#truncated normal
  # #muY ~ dnorm(mean = 0, sd = 5)
  # #sdY ~ dgamma(shape = 2,rate = 1) #weakly informative
  # 
  # sdY ~ dunif(min=0, max=5)
  # muY ~ dnorm(mean = 0, sd = 5)
  # sdY ~ dunif(min = 0,  max = 2)
  
  #Works with both parameterization on real data
  a ~ dbeta(shape1 = 1, shape2 = 5) # low values are preferred 
  
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
    b2[x]~ dgamma(shape = 1, rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  sigma_eps ~ dunif(min= 0, max = 2) #weakly informative 
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*k[t]+betaJump[x]*(J[t+1]-J[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=sigma_eps)
    }
  }
})

#### CURRENTLY NOT IN USE ######################################################
LiuLiModel_SumZero_TwiceDiff_Restated_NonDirichletBeta <- nimbleCode({
  
  # Liu,Li Model with Jumps (Model J1) differenced twice
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)
  
  #Prior Parameterisation
  #Priors on Time Parameters
  #twice differenced random Walk 
  for(t in 1:N_Year){
    k[t] ~ dnorm(mean = 0,sd=2*sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2) 
  
  #alternative for kappa
  # k[1]~dnorm(k0, sd=sigma_time)
  # for(t in 2:N_Year){
  #   k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  # }
  # sigma_time ~  dunif(0,2) 
  # drift ~ dnorm(0,sd=2) 
  # k0 ~ dnorm(0,5)
  # 
  # kappa[1:N_Year] <- SumToZero(k[1:N_Year])
  # kDiff[1:(N_Year-1)] <- kappa[2:N_Year]-kappa[1:(N_Year-1)]
  # 
  #other time effects
  N_t[1:2] <- 0
  Y_t[1:2] <- 0
  J[1:2] <- 0
  for (t in 3:(N_Year+2)){ #one year longer
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t]
  }
  
  JTilde[1:N_Year] <- J[3:(N_Year+2)]-2*J[2:(N_Year+1)]-J[1:N_Year]
  
  kappa[1:N_Year] <- gramschmidt_Nimble3(x = JTilde[1:N_Year],
                                         y = k[1:N_Year])
  p ~ dbeta(shape1 = 2, shape2 = 5)
  muY ~ dnorm(mean = 0, sd = 5)
  sdY ~ dunif(min = 0,  max = 2)
  
  #Dirichlet Prior
  alphaJump[1:N_AgeGroups] <- c(rep(1,3),3,3,3,rep(1,4))
  for(x in 1:N_AgeGroups){
    b1[x]~ dgamma(shape = 1, rate = 1)
    b2[x]~ dgamma(shape = alphaJump[x], rate = 1)
  }
  
  # Dirichlet Distribution is standardized Gamma Dist
  beta[1:(N_AgeGroups-1)] <- b1[1:(N_AgeGroups-1)]/sum(b1[1:N_AgeGroups])
  beta[N_AgeGroups] <- 1 - sum(beta[1:(N_AgeGroups-1)])
  
  betaJump[1:(N_AgeGroups-1)] <- b2[1:(N_AgeGroups-1)]/sum(b2[1:N_AgeGroups])
  betaJump[N_AgeGroups] <- 1 - sum(betaJump[1:(N_AgeGroups-1)])
  
  sigma_eps ~ dunif(0,2)
  
  #Putting all Parameters together
  for(x in 1:N_AgeGroups){
    for(t in 1:N_Year){
      mu[x,t] <- beta[x]*kappa[t]+
        betaJump[x]*(JTilde[t])
      ZMat[x,t] ~ dnorm(mu[x,t],sd=2*sigma_eps)
    }
  }
})

LiuLiModel_SumZero_GS <- nimbleCode({

  # Lee Carter Model with Jumps
  #
  # Time Effect (random walk with drift prior)+
  # Age Effect (Dirichtlet Prior) +
  # Age Specific Intercept (normal prior)+
  # Overdispersion Parameter (normal prior)+
  #
  # Jump Effect according to Liu, Li
  # Jump Occurence (Bernoulli Distributed)
  # Jump Size (Normally Distributed)

  #Prior Parameterisation
  #Priors on Time Parameters
  k[1] ~ dnorm(k0,sd=sigma_time)
  for(t in 2:N_Year){
    k[t] ~ dnorm(k[t-1]+drift,sd=sigma_time) #State space model formualtion
  }
  sigma_time ~  dunif(0,2)
  drift ~ dnorm(0,sd=2)
  k0 ~ dnorm(0,sd=5) ## starting value RW

  kappa[1:N_Year] <- SumToZero(k[1:N_Year])

  #other time effects
  for (t in 1:N_Year){
    N_t[t] ~ dbern(p)
    Y_t[t] ~ dnorm(mean = muY, sd = sdY)
    J[t] <- N_t[t]*Y_t[t] #Jump effect (deterministic)
  }

  #a ~ dunif(min = 0.01, max = 0.99)
  p ~ dbeta(shape1 = 2, shape2 = 5)
  muY ~ dnorm(mean = 0, sd = 10)
  sdY ~ dunif(min = 0,  max = 2)

  #Paramterisation for Age Effects
  for(x in 1:N_AgeGroups){
    alpha[x]~ dnorm(0, sd=5)
  }
  #Dirichlet Prior
  a_dirich[1:N_AgeGroups] <- rep(1,N_AgeGroups)
  beta[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])

  betaJump[1:N_AgeGroups] ~ ddirch(alpha=a_dirich[1:N_AgeGroups])

  sigma_eps ~ dunif(0,2)

  Jtilde[1:N_Year] <- gramschmidt_Nimble(x=kappa[1:N_Year],
                                          y=J[1:N_Year])

  #Putting all Parameters together
  #N is total number of cells i.e. age*year
  for(i in 1:N){
    mu[i] <- alpha[age[i]]+beta[age[i]]*kappa[year[i]]+ #LC Model
      betaJump[age[i]]*Jtilde[year[i]] #Jump Part
    log(eta[i]) ~ dnorm(mu[i],sd=sigma_eps) # Effect
    lambda[i] <- eta[i]*Offset[i]
    y[i] ~ dpois(lambda[i])
  }
})

NimbleJumpTest <- nimbleCode({
  N_t[1] <- 0
  R[1:NTime,1:NTime] <- RMatrix(N_Year = NTime,a = a)
  #J[1] <- 0
  for (t in 2:NTime) {
    N_t[t]~ dbern(p)
    
    # Funktioniert nicht ??
    # J[t] ~ dnorm(mean = N_t[t]*muY,
     #             var = pow(N_t[t]*sdY,2))
    # 
    
  }
  #J_t[1:NTime] <- R[1:NTime,1:NTime]%*%J[1:NTime]
  a ~ dbeta(shape1 = 1, shape2 = 4)
  p ~ dbeta(shape1 = 1,shape2 = 4)
  muY~ dnorm(mean = 0,sd = 4)
  sdY ~ dgamma(shape = 1,rate = 1)
  sigma_eps ~ dgamma(shape = 1,rate = 1)
  
  # Y as multivariate normal 
  for( t in 1:NTime){
    meanY[t] <- inprod(R[t,1:NTime],
                    (N_t[1:NTime]*muY))
    
    sdDat[t] <- pow(inprod(R[t,1:NTime],
                         (N_t[1:NTime]*sdY)),2)
    
    y[t]~dnorm(mean = meanY[t], 
               var = pow(sigma_eps,2) + sdDat[t])
  }
  # for( t in 1:NTime){
  #   y[t]~dnorm(mean = J_t[t], 
  #              sd =sigma_eps)
  # }
})

NimbleJumpTestDiff <- nimbleCode({
  N_t[1] <- 0
  R[1:NTime,1:NTime] <- RMatrix(Ntime = NTime,a = a)
  #J[1] <- 0
  for (t in 2:NTime) {
    N_t[t]~ dbern(p)
  }

  a ~ dbeta(shape1 = 1, shape2 = 4)
  p ~ dbeta(shape1 = 1,shape2 = 4)
  muY~ dnorm(mean = 0,sd = 4)
  sdY ~ dgamma(shape = 1,rate = 1)
  sigma_eps ~ dgamma(shape = 1,rate = 1)
  
  # Differenced Model 
  for(t in 2:NTime){
    meanY[t] <- muY*(
      inprod(R[t,1:NTime],N_t[1:NTime])-
      inprod(R[t-1,1:NTime],N_t[1:NTime])  
    )
    sdDat[t] <-
      pow(inprod(R[t,1:NTime],
                 (N_t[1:NTime]*sdY)),2)+
      pow(inprod(R[t-1,1:NTime],
                   (N_t[1:NTime]*sdY)),2)  
    
    y[t]~dnorm(mean = meanY[t], 
               var = pow(sigma_eps,2) + sdDat[t])
  }
  y[1]~dnorm(mean = 0,sd = sigma_eps)
})

NimbleJumpTestDiff_2 <- nimbleCode({
  
  #Differenced Jump Model with new Constraints
  N_t[1:2] <- 0
  R[1:(NTime+1),
    1:(NTime+1)] <- RTildeMat(Ntime = NTime,a = a)
  #J[1] <- 0
  for (t in 3:(NTime+1)) {
    N_t[t]~ dbern(p)
  }
  
  a ~ dbeta(shape1 = 1, shape2 = 4)
  p ~ dbeta(shape1 = 1,shape2 = 4)
  muY~ dnorm(mean = 0,sd = 4)
  sdY ~ dgamma(shape = 1,rate = 1)
  sigma_eps ~ dgamma(shape = 1,rate = 1)
  

  # #Warum auch immer?
  # for(t in 1:NTime){
  #   for(j in 1:t){ #calculate according to formula
  #     CovMatHelper[j] <- N_t[t+1-j]*pow_int(a,(2*j)-1)*sigma2 #scalar*scalar*scalar
  #   }
  #   Cov[t] <- sum(CovMatHelper[1:t])
  # }
  
  
  #Geht aber langsam...
  for(t in 1:NTime){
    CovVec[t] <- CovZMat(N_t = N_t[1:(NTime+1)],
                         a =  a,
                         t =  t,
                         sigma2 = pow_int(sdY,2))
  }
  
  # Differenced Model 
  for(t in 1:NTime){
    meanY[t] <- muY*(
      inprod(R[t+1,1:(NTime+1)],N_t[1:(NTime+1)])-
      inprod(R[t,1:(NTime+1)],N_t[1:(NTime+1)])  
    )
    sdDat[t] <-
      pow(inprod(R[t+1,1:(NTime+1)],
                 (N_t[1:(NTime+1)]*sdY)),2)+
      pow(inprod(R[t,1:(NTime+1)],
                 (N_t[1:(NTime+1)]*sdY)),2)-
      2*CovVec[t]
      
  
    y[t]~dnorm(mean = meanY[t], 
               var = pow(sigma_eps,2) + sdDat[t])
  }
})
