#### Script to Forecast Death Rates
library(nimble);library(tidyverse); library(nimbleHMC)

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 


load(file = file.path(getwd(),"Results/SamplesUS_NoRep.RData"))
load(file = file.path(getwd(),"Results/SamplesIt_NoRep.RData"))
load(file = file.path(getwd(),"Results/SamplesSp_NoRep.RData"))
load(file = file.path(getwd(),"Data/CovidData.RData")) # Load Data




# For Forecasting future Z_x,t we need the following
#1) Generation of future kt's 
#2) Generation of future Nt's then generation of future J'ts
#3) Generation of future error Term
# All in one Function 
FutureZ <- function(Samples,S ,H ,OwnMod = TRUE,NAge){
  
  # @H: Forecast of H time periods ahead
  # @S: Amount of draws from the posterior
  # @NAge: Amount of age groups in the data
  
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
  muJMat <- matrix(data = 0,nrow = S,ncol = H+1)
  sdJMat <- matrix(data = 0, nrow = S, ncol = H+1)
  FutureJt <- matrix(data = 0,nrow = S,ncol = H+1) #one column more
  
  
  for(s in 1:S){
    R <- RMatrix(Ntime = ncol(NtMatrix)+H,a = aVec[s]) #Creation of R Matrix
    #Written in vector notation
    muJMat[s,] <- muYVec[s]* #E(Jt) = mu*Rt*N
      as.numeric(R[ncol(NtMatrix):(ncol(NtMatrix)+H),]%*%NtTotMat[s,])
    
    #using a loop. Vector notation does not work with matrix algebra
    for(j in 1:(H+1)){ #Var(Jt)=Rt*N*sigma*N*Rt
      h <- j-1 #helper function
      sdJMat[s,j] <- sqrt(
        R[ncol(NtMatrix)+h,]%*%NtTotMat[s,]%*% #Rt*N
        pow(sdYVec[s],2)%*%t(R[ncol(NtMatrix)+h,]%*%NtTotMat[s,]) #sigma*N*Rt
        )
    }
  FutureJt[s,] <- rnorm(n = H+1, mean = muJMat[s,],sd = sdJMat[s,])
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
  return(list(FutureZArray,
               FutureNtMat))
}

FutureZNoRep <- function(Samples,S ,H ,OwnMod = TRUE,NAge,LC = FALSE){
  
  # @H: Forecast of H time periods ahead
  # @S: Amount of draws from the posterior
  # @NAge: Amount of age groups in the data
  # @Samples: Posterior Samples from nimble
  # @OwnMod: Binary Variable to choose either Liu,Li or Own model
  # @LC: Binary if LC model should be forecasted
  
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
  return(FutureZArray)
}
############## FOR US DATA #####################################################
S <- 2000
H <- 50
FutureZUS <- FutureZNoRep(H=H, OwnMod = TRUE, Samples = do.call("rbind",SamplesOwn_US_NoRepara), S=S,
                     NAge = 10)


#Future log death rates with cummulative sums
FutureLogMatArray <- array(data = 0,dim = c(10,H,S),
                           dimnames = list("NewAgeInd"=1:10, "FCHorizon"=1:H, "It"=1:S))
for(s in 1:S){
  FutureLogMatArray[,,s] <- log(LambdaMatUS)[,ncol(LambdaMatUS)]+ #Last obs rate
                                  t(apply(FutureZUS[,,s],1,cumsum)) #Cumsum
}


#Using melt, to transform array into long format for plotting
FutureLambdaVecIt_US <- 
  reshape2::melt(FutureLogMatArray, value.name = "logRate") %>% 
  mutate(Rate = exp(logRate),
         Type="FC",
         TInd = FCHorizon+41) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")))

theme_set(theme_minimal(base_size = 10))

### In Long Format ###
FutureValuesVec <- 
  FutureLambdaVecIt_US %>%  
  group_by(NewAgeInd,TInd) %>% #Calculate quantiles
  reframe("Mean"=median(Rate),
          "PiLo"=quantile(Rate, probs = c(0.25,0.1,0.05)),
          "PiUp"=quantile(Rate, probs = c(0.75,0.9,0.95))) %>% 
  mutate("Width"=rep(c("0.5","0.8","0.99"),10*H)) %>%  #Add Age Width
  rename("Rate"=Mean) %>% 
  mutate("Type"="FC")

ObservedValuesVec <- LambdaVecUs %>% 
  select(-c(Year, Deaths, Pop, ZVal)) %>% 
  mutate("PiLo"=Rate,
         "PiUp"=Rate) %>% 
  slice(rep(1:nrow(.),3)) %>% #each row 3 times
  mutate(Width = rep(c("0.5","0.8","0.99"),each=41),
         Type ="InSamp")



FutureValuesVec %>% 
  full_join(x = ObservedValuesVec, 
            y = .) %>% 
  ggplot(aes(x =TInd))+
  ggdist::geom_lineribbon(aes(y = Rate, fill = Width, 
                              ymin = PiLo, ymax = PiUp,
                              group = NewAgeInd),
                          alpha  = 0.8, linewidth = 0.9)+
  facet_wrap(~NewAgeInd, scales = "free")+
  scale_color_manual(values = c(FC="#08306b",InSamp = "black"))

  

### Option 2 with Point_interval
ObservedValuesVec2_US <- 
  LambdaVecUs %>%  #Take observed Vector and calculate PI's
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+"))) %>% 
  group_by(AgeGroup,Year) %>% 
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>% 
  mutate(TInd = match(Year, unique(Year))) %>% 
  select(-Year) %>% 
  mutate(Type ="Obs")


FutureValuesVecQuantile_US <- 
  FutureLambdaVecIt_US %>% 
  group_by(AgeGroup,TInd) %>% 
  ggdist::point_interval(Rate,.width = c(0.5,0.8,0.99)) %>% 
  mutate(Type ="FC") %>% 
  full_join(x = filter(ObservedValuesVec2_US, TInd == max(TInd)),
            y =. )

AgeFilter <- c("20-29","40-49","50-59","70-79")

ObservedValuesVec2_US %>%
  filter(AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =TInd))+
  ggdist::geom_lineribbon(aes(y = Rate,
                              ymin = .lower, ymax = .upper,
                              group = AgeGroup),
                          linewidth = 0.9)+
  ggdist::geom_lineribbon(data = filter(FutureValuesVecQuantile_US, AgeGroup %in% AgeFilter), 
                          aes(y = Rate, 
                              ymin = .lower, ymax = .upper,
                              group = AgeGroup),linewidth = 0.9,
                          col = "#400040")+
  geom_point(data = filter(FutureValuesVecQuantile_US, AgeGroup %in% AgeFilter),
             aes(y = Rate))+
  scale_fill_manual(values=c("#e0cde1","#b789b9","#98579b"))+
  facet_wrap(~AgeGroup, scales="free")


####### Including Forecasts for the Lee Carter Model ###########################
load(file = file.path(getwd(), "Results/SamplesUS_LC.RData"))
FutureZUS_LC <- FutureZNoRep(H=H, OwnMod = TRUE, 
                             Samples = do.call("rbind", Samples_LC_US), S=S,
                             NAge = 10, LC = TRUE)


FutureLogMatArrayUS_LC <- array(data = 0,dim = c(10,H,S),
                                dimnames = list("NewAgeInd"=1:10, "FCHorizon"=1:H, "It"=1:S))

for(s in 1:S){
  FutureLogMatArrayUS_LC[,,s] <- log(LambdaMatUS)[,ncol(LambdaMatUS)]+t(apply(FutureZUS_LC[,,s],1,cumsum))
}

#Get Iterations into long format
FutureLambdaVecIt_US_LC <- 
  reshape2::melt(FutureLogMatArrayUS_LC, value.name = "logRate") %>% 
  mutate(Rate = exp(logRate),
         Type="FC",
         TInd = FCHorizon+42) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")))

#Calculate Quantiles and add Model Name
FutureValuesVecQuantile_US_LC <- 
  FutureLambdaVecIt_US_LC %>% 
  group_by(AgeGroup,TInd) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>% 
  mutate(Type ="FC",
         Model="LC") 

#Too Existing Quantile frame add the one from the LC model
FutureValuesVecQuantile_US <- 
  FutureValuesVecQuantile_US %>% #use FC of Jump Model
  mutate(Model="Jump") %>% 
  bind_rows(., FutureValuesVecQuantile_US_LC) #add FC of LC model


AgeFilter <- c("40-49","50-59","60-69","70-79")

### Wörks, but why LC Model with wider PI's ? 
ObservedValuesVec2_US %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =TInd, group = AgeGroup))+
  geom_line(aes(y = Rate),linewidth = 1.2)+
  ggdist::geom_lineribbon(data = filter(FutureValuesVecQuantile_US, AgeGroup %in% AgeFilter),
                          aes(y = Rate,
                              ymin = .lower, ymax = .upper,
                              fill = interaction(.width,Model)),
                          linewidth = 1.2)+
  scale_fill_manual(values=c("#e0cde1","#b789b9","#98579b",
                             "#12492f","#f56038","#0a2f35"))+
  scale_color_manual(values =c("#806E26","blue"))+
  facet_wrap(~AgeGroup, scales = "free")

############# SPAIN ############################################################
S <- 1500
H <- 50
set.seed(42)
FutureZSp <- FutureZNoRep(H=H, OwnMod = TRUE, 
                     Samples = do.call("rbind", SamplesOwn_Sp), S=S,
                     NAge = 10, LC=TRUE)

#Future log death rates with cummulative sums
FutureLogMatArraySp <- array(data = 0,dim = c(10,H,S),
                           dimnames = list("NewAgeInd"=1:10, "FCHorizon"=1:H, "It"=1:S))
for(s in 1:S){
  FutureLogMatArraySp[,,s] <- log(LambdaMatSp)[,ncol(LambdaMatSp)]+t(apply(FutureZSp[,,s],1,cumsum))
}


#Transform FC Array into long format 
FutureLambdaVecIt_Sp <- 
  reshape2::melt(FutureLogMatArraySp, value.name = "logRate") %>% 
  mutate(Rate = exp(logRate), #Calculate Death Rate
         Type="FC", #add some idicies for plotting
         Year = FCHorizon+2022) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")))
#Get observed Values into same format. 
#This has to be done, so that there is no space between observed values and
#forecasted values when plotting
ObservedValuesVec2_Sp <- 
  LambdaVecSp %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                       levels = 1:10,
                       labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                         seq(9,89,10)),"90+"))) %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>% 
  mutate(Type ="Obs")

#Calculate Quantiles of Forecasted Death rates for line ribbon 
FutureValuesVecQuantile_Sp <- 
  FutureLambdaVecIt_Sp %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>% #Quantile Calc
  mutate(Type ="FC") %>% 
  full_join(x = filter(ObservedValuesVec2_Sp, Year == max(Year)),
            y =. )


AgeFilter <- c("40-49","50-59","60-69","70-79")

SimPaths <- c(1,398,1525)
pdf(file=file.path(getwd(),"Bilder/ForecastOwn_Spain.pdf"), width=16, height = 9)
FutureValuesVecQuantile_Sp %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =Year))+
  ggdist::geom_lineribbon(aes(y = Rate, ymin = .lower, ymax = .upper,
                              group = AgeGroup,
                              #fill = forcats::fct_rev(ordered(.width)) #standard setting
                              ), 
                          linewidth = 1.2,col = "#834333")+
  geom_line(data=filter(ObservedValuesVec2_Sp, AgeGroup %in% AgeFilter),
            aes(y=Rate, group = AgeGroup), linewidth=1.2)+
  scale_fill_manual(values=c("#f0ead5","#dbcf9c","#CCBA72"))+
  geom_vline(xintercept = 2022,linetype = "dashed", col="gray30")+
  facet_wrap(~AgeGroup, scales = "free")+
  ylab("Death Rate")+xlab("Year")+
  scale_x_continuous(breaks = c(1990,2010,2030,2050,2070))+
  theme(
    axis.text = element_text(size = 22), #change size of axis text
    axis.title = element_text(size = 25, face = "bold"),
    legend.position = "none",
    strip.text.x=element_text(face="bold", size = 22) #increase size of grouping text
)
dev.off()
  

#100,70, 30%

#Colors from "https://imagecolorpicker.com/color-code/81a88d"
#Italy Colors: Line: "#355C41", 50PI;#81a88d; 80 PI: "#a7c2af", 99 PI "#d9e5dd" (70% light)

### Option 2 single Lines

#Find out times of extreme Jumps
FutureZSp[7,,] %>% apply(., 1, which.max)

FutureZSp %>% rowSums() #find iterations with lots of jumps

LambdaVecSpHelp <- LambdaVecSp %>% 
  mutate(TInd = match(Year, unique(Year))) %>% 
  filter(NewAgeInd == 7)

#alternative, find 90% quantile
FutureZSp[7,,] %>% apply(., 1,function(x){which.min(x-quantile(x, probs=0.99))})


SimPaths2 <- c(1,398,1525)
pdf(file=file.path(getwd(),"Bilder/FutureDeathRateSpain_Path.pdf"),
    width = 16, height = 9)
FutureLambdaVecIt_Sp %>% 
  filter(NewAgeInd == 7 & It %in% SimPaths2) %>% 
  full_join(x = LambdaVecSpHelp %>% 
              filter(TInd == 42) %>% 
              slice(rep(1,3)) %>% 
              mutate(It =SimPaths2),
            y =., by = c("TInd","Rate","It","NewAgeInd")) %>% 
  ggplot(., aes(x = TInd)) +
  geom_line(aes(y = Rate, group = It, col = as.factor(It)), linewidth = 1.2)+
  scale_color_manual(labels=c("Path 1","Path 2", "Path 3"),
                              values= c("#540b0e","#132a13","#0466c8"))+ #old c("#540b0e","#e09f3e","#335c67")
  geom_line(data = LambdaVecSpHelp %>% 
              mutate(TInd = match(Year, unique(Year))), aes(x = TInd, y = Rate),
            linewidth = 1.2, linetype = "twodash")+
  geom_vline(xintercept =  42, linetype = "dashed", col="gray30")+
  scale_x_continuous(labels = c(1980,1990,2000,2010,2020,2023,2040,2050,2060,2070), 
                     breaks = c(1,11,21,31,41,51,61,71,81,91))+
  ylab("Death Rate")+
  labs(col ="")+ggtitle("Age Group 60-69")+
  theme(legend.position = "bottom",
        axis.text = element_text(size = 18), #change size of axis text
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 20))
dev.off()

###### Together with FC of the LC Model ######################################
load(file = file.path(getwd(), "Results/SamplesSp_LC.RData"))
debugonce(FutureZNoRep)
FutureZSp_LC <- FutureZNoRep(H=H, OwnMod = TRUE, 
                             Samples = do.call("rbind", Samples_LC_Sp), S=S,
                             NAge = 10, LC = TRUE)

FutureLogMatArraySp_LC <- array(data = 0,dim = c(10,H,S),
                                dimnames = list("NewAgeInd"=1:10, 
                                                "FCHorizon"=1:H, 
                                                "It"=1:S))

FutureZSp_LC[4,,] %>% apply(., 1, quantile, 0.75)

for(s in 1:S){
  FutureLogMatArraySp_LC[,,s] <- log(LambdaMatSp)[,ncol(LambdaMatSp)]+
                                  t(apply(FutureZSp_LC[,,s],1,cumsum))
}



#Get Iterations into long format
FutureLambdaVecIt_Sp_LC <- 
  reshape2::melt(FutureLogMatArraySp_LC, value.name = "logRate") %>% 
  mutate(Rate = exp(logRate),
         Type="FC",
         Year = FCHorizon+2022) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")))

#Calculate Quantiles and add Model Name
FutureValuesVecQuantile_Sp_LC <- 
  FutureLambdaVecIt_Sp_LC %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>% 
  mutate(Type ="FC",
         Model="LC") 

#Too Existing Quantile frame add the one from the LC model
FutureValuesVecQuantile_Sp <- 
  FutureValuesVecQuantile_Sp %>% #use FC of Jump Model
  mutate(Model="Jump") %>% 
  bind_rows(., FutureValuesVecQuantile_Sp_LC) #add FC of LC model
  
# AgeFilter <- c("20-29","30-39","70-79","80-89")
AgeFilter <- c("30-39")

### Wörks, but why LC Model with wider PI's ? 
ObservedValuesVec2_Sp %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  mutate(Model = "LC") %>% 
  ggplot(aes(x =Year, group = AgeGroup))+
  geom_line(aes(y = Rate),linewidth = 1.2)+
  ggdist::geom_lineribbon(data = filter(FutureValuesVecQuantile_Sp, AgeGroup %in% AgeFilter),
                          aes(y = Rate,
                              ymin = .lower, ymax = .upper,
                              fill = interaction(.width,Model),
                              color = Model),
                          linewidth = 1.2, alpha = 0.8)+
  scale_fill_manual(values=c("#f0ead5","#dbcf9c","#CCBA72",
                             "#12492f","#f56038","#0a2f35"))+
  scale_color_manual(values=c("blue","red"))+
  facet_wrap(~AgeGroup, scales = "free")


########### FOR UK WAR DATA ###################################################
load(file = file.path(getwd(), "Results/Samples_UKWar.RData"))
S <- 3000
H <- 50
set.seed(42)
FutureZUKWar <- FutureZNoRep(H=H, OwnMod = TRUE, 
                     Samples = do.call("rbind", Samples_OwnMod_OC_War), S=S,
                     NAge = 10)



#Future log death rates with cummulative sums
FutureLogMatArrayUK <- array(data = 0,dim = c(10,H,S),
                             dimnames = list("NewAgeInd"=1:10, "FCHorizon"=1:H, "It"=1:S))
for(s in 1:S){
  FutureLogMatArrayUK[,,s] <- log(LambdaMatWar)[,ncol(LambdaMatWar)]+t(apply(FutureZUKWar[,,s],1,cumsum))
}


FutureLambdaVecUK <- 
  reshape2::melt(FutureLogMatArrayUK, value.name = "logRate") %>% 
  mutate(Rate = exp(logRate),
         Type="FC",
         TInd = FCHorizon+ncol(LambdaMatWar)) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c("<1","1-4",
                                      paste0(paste0(seq(5,75,10),sep="-"),
                                             seq(14,84,10)))))

ObservedValuesVec2_UK <- 
  LambdaVecWar %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c("<1","1-4",
                                      paste0(paste0(seq(5,75,10),sep="-"),
                                             seq(14,84,10))))) %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>% 
  mutate(TInd = match(Year, unique(Year))) %>% 
  select(-Year) %>% 
  mutate(Type ="Obs")

FutureLambdaVecQuantile_UK <- 
  FutureLambdaVecUK %>% 
  group_by(AgeGroup,TInd) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>% 
  mutate(Type ="FC") %>% 
  full_join(x = filter(ObservedValuesVec2_UK, TInd == max(TInd)),
            y =. )

AgeFilter <- c("25-34","35-44","65-74")
SimPaths <- sample(1:2000,size = 3,replace = TRUE)
ObservedValuesVec2_UK %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =TInd))+
  ggdist::geom_lineribbon(aes(y = Rate,
                              ymin = .lower, ymax = .upper,
                              group = AgeGroup),
                          linewidth = 1.2)+
  ggdist::geom_lineribbon(data = filter(FutureLambdaVecQuantile_UK, 
                                        AgeGroup %in% AgeFilter), 
                          aes(y = Rate, 
                              ymin = .lower, ymax = .upper,
                              group = AgeGroup),linewidth = 1.2,
                          col = "#806E26")+
  scale_fill_manual(values=c("#f0ead5","#dbcf9c","#CCBA72"))+
  # geom_line(data = filter(FutureLambdaVecUK, It %in% SimPaths,
  #                         AgeGroup %in% AgeFilter),
  #           aes(y = Rate, group = It, col = as.factor(It)), linewidth = 0.9,
  #           alpha = .8)+
  # scale_color_manual(labels=c("Path 1","Path 2", "Path 3"),
  #                    values= c("#540b0e","#132a13","#0466c8"))+
  facet_wrap(~AgeGroup, scales = "free")



