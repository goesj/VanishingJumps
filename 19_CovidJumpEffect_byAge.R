### Estimating the effect on the death rates of the Jump ######################
library(nimble);library(tidyverse); library(nimbleHMC)

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 


load(file = file.path(getwd(),"Results/SamplesSp_NoRep.RData")) # own Model
load(file = file.path(getwd(),"Results/SamplesUS_NoRep.RData")) # own Model
load(file = file.path(getwd(),"Results/SamplesIt_NoRep_2.RData")) # own Model

FutureZNoRep_Param <- function(Samples,S ,H ,OwnMod = TRUE,NAge,LC = FALSE){
  
  # @H: Forecast of H time periods ahead
  # @S: Amount of draws from the posterior
  # @NAge: Amount of age groups in the data
  # @Samples: Posterior Samples from nimble
  # @OwnMod: Binary Variable to choose either Liu,Li or Own model
  # @LC: Binary if LC model should be forecasted
  
  #Get Random indicies from posterior. All should have the same draw 
  PostDraw <- sample(1:nrow(Samples), 
                     size = S, 
                     replace = TRUE)
  
  #1) Generation of future kt's
  driftVec <- Samples[PostDraw, #Select S posterior draws
                      grep("drift",colnames(Samples))] #extract drift
  
  sigmaTime <- Samples[PostDraw,
                       grep("sigma_time",colnames(Samples))]
  
  # generation of future kt's, for each posterior draw, generate H future kt's
  #Row's are posterior draws, columns are time periods
  FutureKtMat <- sapply(1:H, function(x){
    rnorm(n = S, mean = driftVec, sd = sigmaTime)})
  
  if(LC == FALSE){ ## Lee Carter Model Forecasts
    #2) Generation of future Jt's (be very careful about indices..)
    
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
    
  } else {
    
    betaMat <- Samples[PostDraw, #Select S posterior draws
                       grep("beta", #starts with p to select "p" only 
                            colnames(Samples))][,1:NAge]
    
    sigma_epsVec <- Samples[PostDraw, #Select S posterior draws
                            grep(paste0("\\<","sigma_eps"), #starts with p to select "p" only 
                                 colnames(Samples))]
    
    FutureZArray <- array(data = 0,dim = c(NAge,H,S),
                          dimnames = list("Age"=1:NAge, "Time"=1:H, "It"=1:S))
    for(s in 1:S){
      #via Matrix multiplication (outer product)
      FutureZArray[,,s] <- 
        betaMat[s,]%*%t(FutureKtMat[s,])+ #Normal effect
        sapply(1:H, function(x){rnorm(n = NAge, 
                                      mean = 0, sd = sigma_epsVec[s])}) #error Term
    }
  }
  ReturnList <- list("Rates"=FutureZArray,
                     "Jt"=FutureJt,
                     "betaJump"=betaJumpMat)
  return(ReturnList)
}

### Calculate Future Mort Improvement Rates incl. Parameters
S <- 10000
H <- 50
set.seed(42)
FutureZSp <- FutureZNoRep_Param(H=H, OwnMod = TRUE, 
                          Samples = do.call("rbind", SamplesOwn_Sp), S=S,
                          NAge = 10, LC=FALSE)

FutureZUS <- FutureZNoRep_Param(H=H, OwnMod = TRUE, 
                                Samples = do.call("rbind", SamplesOwn_US_NoRepara), S=S,
                                NAge = 10, LC=FALSE)

FutureZIt <- FutureZNoRep_Param(H=H, OwnMod = TRUE, 
                                Samples = do.call("rbind", SamplesOwn_It), S=S,
                                NAge = 10, LC=FALSE)


FutureZList <- list(FutureZSp,
                    FutureZUS,
                    FutureZIt)

JumpEffectArray <- array(0, dim=c(10, H+1, S,3),
                       dimnames = list("Age"=1:10,
                                       "FC"=1:(H+1),
                                       "It"=1:S,
                                       "Country"=c("b) Spain","c) US","a) Italy")))
for (c in 1:3) { #over countries
  for (s in 1:S) { #over Iterations
    for(x in 1:10){ #over Age
    JumpEffectArray[x,,s,c] <- FutureZList[[c]]$betaJump[s,x]*FutureZList[[c]]$Jt[s,]  
    }
  }
}

JumpEffectQuant <- 
  reshape2::melt(JumpEffectArray, value.name = "logEffect") %>% 
  mutate(AgeGroup = factor(Age,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+"))) %>% 
  group_by(FC,AgeGroup,Country) %>% 
  ggdist::point_interval(logEffect, .width = c(0.9,0.95,0.99))

JumpEffectIt <- 
  reshape2::melt(JumpEffectArray, value.name = "logEffect") %>% 
  mutate(AgeGroup = factor(Age,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")))

Its <- sample(1:S,size = 10, replace =FALSE)
x11()
JumpEffectQuant %>% 
  filter(Country=="Spain") %>% 
  ggplot(aes(x = FC))+
  ggdist::geom_lineribbon(aes(y = logEffect, 
                              ymin = .lower, ymax = .upper,group=AgeGroup),
                          alpha = 0.5)+
  geom_line(data = filter(JumpEffectIt, It %in% Its,
                          Country=="Spain"),
            aes(x = FC, y = logEffect, group = It, col = factor(It)),
            linewidth = 1.1)+
  facet_wrap(~AgeGroup)


### small check, bene
JumpEffectQuant %>% 
  filter(AgeGroup == "30-39") %>% 
  filter(.width == 0.9) %>% 
  group_by(Country) %>% 
  summarise(mean(.upper))

#For US, the first couple years are a affected by current Jump in last two years
JumpEffectIt %>% 
  filter(AgeGroup == "30-39") %>% 
  group_by(Country,FC) %>% 
  reframe("Quant"=quantile(logEffect,0.95)) %>% 
  ggplot(aes(x = FC, y = Quant))+
  geom_line(aes(group = Country, col= Country))



#Two options for quantiles: 1: per time period never over certain value, 
#Calculation of actual,observed Covid increase
CovidIncrease <- bind_rows(mutate(LambdaVecSp, Country="b) Spain"),
                           mutate(LambdaVecUs, Country="c) US"),
                           mutate(LambdaVecIt, Country="a) Italy")) %>% 
  group_by(NewAgeInd,Country) %>% 
  filter(Year %in% c(2019,2020)) %>% 
  # filter(Year %in% c(2016:2019,2020)) %>% 
  # mutate(Period=ifelse(Year<2020,"Ref","Cov")) %>% #define reference period and Covid
  # group_by(Period, Country, NewAgeInd) %>% #group by new Period
  # reframe(Rate2 = mean(Rate)) %>%  #Caclulate the mean
  # arrange(desc(Period)) %>%  #Change order, for easier division in long format
  # group_by(NewAgeInd, Country) %>% 
  #reframe("PercentageIncrease"=exp(diff(log(Rate2)))) %>%  #calculation of mortality rate increase between mean 2016:2019 and 2020
  reframe("PercentageIncrease"=exp(diff(log(Rate)))) %>%  #calculation of mortality rate increase between 2019 and 2020
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")),
         CovidIncrease=" ")#Empty column so it can appear in legend

#Calculate the percentage for all ages
pdf(file=file.path(getwd(),"Bilder/JumpAffection_Countries.pdf"), width=16, height = 9)
JumpEffectQuant %>% 
  group_by(AgeGroup,.width,Country) %>% 
  summarise("Perc"=mean(.upper)) %>%
  ggplot(aes(x=AgeGroup))+
  geom_line(aes(y =(exp(Perc)-1)*100, group=.width, col=factor(.width),
                #linetype = factor(.width)
                ),
            linewidth = 1.2)+
  geom_line(data = CovidIncrease, aes(x =AgeGroup, 
                                      y=(PercentageIncrease-1)*100,
                                      group=Country, linewidth = "Covid Increase",
                                      linetype = "Covid Increase"))+
  scale_linetype_manual(values = 2, #set linetype of Covid Increase 
                        name ="Type")+
  scale_linewidth_manual(values = 1.2, # set line width of Covid Increase
                         guide = "none")+ #remove from legend
  scale_color_manual(values = c("#9a031e","#e36414","#023047"))+
  facet_wrap(~Country, nrow = 3, scales="free")+
  ylab("Percentage Increase")+
  xlab("Age Group")+
  labs(col="Prediction Interval")+
  theme(
    axis.text = element_text(size = 20), #change size of axis text
    axis.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.position = "bottom",
    legend.key.width = unit(1.7,"cm"),
    strip.text.x=element_text(face="bold", size = 22)
  )+
  guides(linetype = guide_legend(override.aes = list(linewidth = 1.3)))
dev.off()



### Option 2: Look at extreme effects for years with a shock
JumpEffectIt %>% 
  filter(logEffect > 0.005) %>%  #filter Iterations when there is actually a small effect 
  group_by(FC,AgeGroup,Country) %>% 
  ggdist::point_interval(logEffect, .width = c(0.9,0.95,0.99)) %>% 
  group_by(AgeGroup,.width,Country) %>% 
  summarise("Perc"=mean(.upper)) %>%
  ggplot(aes(x=AgeGroup))+
  geom_line(aes(y =(exp(Perc)-1)*100, group=.width, col=factor(.width)),linewidth = 1.1)+
  geom_line(data = CovidIncrease, aes(x =AgeGroup, y=(PercentageIncrease-1)*100,
                                      group=Country),
            linetype ="dashed", linewidth = 1.05)+
  scale_color_manual(values = c("#9a031e","#e36414","#023047"))+
  facet_wrap(~Country, nrow = 3, scales="free")+
  ylab("Percentage Increase")+
  xlab("Age Group")+
  labs(col="Prediction Interval")



Its <- sample(1:S,size = 10, replace =FALSE)
JumpEffectIt %>% 
  filter(logEffect > 0.005) %>%  #filter Iterations when there is actually a small effect 
  group_by(FC,AgeGroup,Country) %>% 
  ggdist::point_interval(logEffect, .width = c(0.9,0.95,0.99)) %>% 
  filter(Country=="Spain") %>% 
  ggplot(aes(x = FC))+
  ggdist::geom_lineribbon(aes(y = logEffect, 
                              ymin = .lower, ymax = .upper,group=AgeGroup),
                          alpha = 0.5)+
  geom_line(data = filter(JumpEffectIt, It %in% Its,
                          Country=="Spain"),
            aes(x = FC, y = logEffect, group = It, col = factor(It)),
            linewidth = 1.1)+
  facet_wrap(~AgeGroup)
