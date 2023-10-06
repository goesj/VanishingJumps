### Goal is comparison of Jump approach with LC model ######
library(nimble);library(tidyverse); library(reshape2); library(rstan);
library(StMoMo)

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 


#load Samples
load(file = file.path(getwd(),"Results/SamplesSp.RData")) # own Model
load(file = file.path(getwd(),"Results/SamplesUS.RData")) # own Model
load(file = file.path(getwd(),"Results/SamplesIt.RData")) # own Model


#load Data
load(file = file.path(getwd(),"Data/CovidData.RData")) # Load Data

########### Comparison Own Vs StMoMo ###########################################
S <- 1500
H <- 50
set.seed(42)
#Calculate Future Mortality Improvement Rates
FutureZSp <- FutureZ(H=H, OwnMod = TRUE, 
                          Samples = do.call("rbind", SamplesOwn_Sp), S=S,
                          NAge = 10)

# Calculate Future log death rates using cummulative sums

#First create empty Array
FutureLogMatArraySp <- array(data = 0,dim = c(10,H,S),
                             dimnames = list("NewAgeInd"=1:10, 
                                             "FCHorizon"=1:H, "It"=1:S))
#Calculation of Future log death rates
for(s in 1:S){
  FutureLogMatArraySp[,,s] <- log(LambdaMatSp)[,ncol(LambdaMatSp)]+
                              t(apply(FutureZSp[,,s],1,cumsum))
}

#Transform into Long Format
FutureLambdaVecIt_Sp <- 
  reshape2::melt(FutureLogMatArraySp, value.name = "logRate") %>% 
  mutate(Rate = exp(logRate), #Calculate Death Rate
         Type="FC", #add some idicies for plotting
         Year = FCHorizon+2022) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")))



#Calculate Quantiles 
FutureValuesVecQuantile_Sp <- 
  FutureLambdaVecIt_Sp %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.9,0.95,0.99)) 


####### Compare with StMoMo ########################################
constLC <- function(ax, bx, kt, b0x, gc, wxt, ages){
  c1 <- mean(kt[1, ], na.rm = TRUE)
  c2 <- sum(bx[, 1], na.rm = TRUE)
  list(ax = ax + c1 * bx, bx = bx / c2, kt = c2 * (kt - c1))
}

LC <- StMoMo(link = "log", 
             staticAgeFun = TRUE, 
             periodAgeFun = "NP",
             constFun = constLC)


DeathMatSp <- matrix(LambdaVecSp$Deaths,
                     nrow = length(unique(LambdaVecSp$NewAgeInd)),
                     byrow = FALSE)

ExpMatSp <- matrix(LambdaVecSp$Pop,
                   nrow = length(unique(LambdaVecSp$NewAgeInd)),
                   byrow = FALSE)

### Fit Model #######
LCfit <- fit(LC, Dxt = DeathMatSp, 
             Ext = ExpMatSp)

#Calculate Forecasts
LCsim <- stats::simulate(LCfit,nsim = 500,h=50,gc.order = c(0, 1, 0))

#Put into long format
FutureLambdaVecIt_Freq <- 
  reshape2::melt(LCsim$rates, value.name = "Rate") %>%
  rename("NewAgeInd"=Var1,
         TInd = Var2, It = Var3) %>%
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")),
         Year = TInd + 1980) 

#Calculation of Quantiles
FutureValuesVecQuantile_Sp_Freq <-
  FutureLambdaVecIt_Freq %>%
  group_by(AgeGroup,Year) %>%
  mutate("logRate"=log(Rate)) %>%
  ggdist::point_interval(Rate, .width = c(0.9)) #90% PI only

#### Plot Both togethers #################
#Observed Data Frame ; for plotting needed only
ObservedValuesVec2_Sp <- 
  LambdaVecSp %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+"))) %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.9,0.95,0.99)) %>% 
  mutate(Type ="Obs")

#Add Observed Data to Quantiles of Own Model
FutureValuesVecQuantile_Sp <- 
  FutureValuesVecQuantile_Sp %>% 
  full_join(x = filter(ObservedValuesVec2_Sp, Year == max(Year)),
            y =. ) %>% #add last year, so that there is no space in plot 
  mutate(Type ="FC") %>% 
  full_join(x = filter(ObservedValuesVec2_Sp),
            y =. )

#Add Observed Values to Quantiles for Plotting
FutureValuesVecQuantile_Sp_Freq <-
  FutureValuesVecQuantile_Sp_Freq
  mutate(Type ="FC",
         Model="LC") %>% 
  full_join(x = filter(ObservedValuesVec2_Sp,
                       .width == 0.9,
                       Year == max(Year)),
            y =. ) %>% 
  mutate(.width= 90,#Use different name for Width Column, so it gets alinged a different color
         Type = "Mean FC Freq") 


AgeFilter <- c("40-49","50-59","60-69","70-79")

### Add layer after layer manually, looks nicer ....
theme_set(theme_minimal(base_size = 10))
FutureValuesVecQuantile_Sp %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =Year, group = AgeGroup))+
  ggdist::geom_lineribbon(aes(y = Rate, ymin = .lower, 
                              ymax = .upper,fill = interaction(.width),
                              color = Type,
                              linetype = Type),
                          linewidth = 1.1, alpha = 0.9)+
  geom_vline(xintercept = 2022,linetype = "dashed", col="gray30")+
  facet_wrap(~AgeGroup, scales = "free")+
  geom_ribbon(data = filter(FutureValuesVecQuantile_Sp_Freq,
                            AgeGroup %in% AgeFilter),
              aes(x = Year, ymin = .lower, ymax = .upper,
                  fill = factor(.width)), alpha = 0.7)+
  geom_line(data = filter(FutureValuesVecQuantile_Sp_Freq,
                          AgeGroup %in% AgeFilter), 
            aes(y = Rate, x = Year, col = Type, linetype = "Mean FC Freq"),
            linewidth = 1.1,
            alpha = 0.7)+
  scale_color_manual(values =c("#834333","#03045e","black"),
                     labels =c("Mean FC Own","Mean FC Freq","Observed"),
                     name = "Type")+
  scale_fill_manual(values=c("#c9b045","#CCBA72","#fdf2c5","#48cae4"),
                    labels =c("Own 90%","Own 95%","Own 99%","Freq 90%"),
                    name = "PI's")+
  scale_linetype_manual(values = c(1,3,6),
                          #c("solid","dotdash","dotted"),
                        labels = c("Mean FC Own",
                                   "Mean FC Freq",
                                   "Observed"),
                        name ="Type")+
  ylab("Death Rate")+xlab("Year")+
  scale_x_continuous(breaks = c(1990,2010,2030,2050,2070))+
  theme(
    axis.text = element_text(size = 20), #change size of axis text
    axis.title = element_text(size = 22),
    legend.position = "bottom",
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.key.width = unit(1.1,"cm"),
    plot.title = element_text(hjust = 0.5, size = 20),
    strip.text.x=element_text(face="bold", size = 22) #increase size of grouping text
  )+
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         color = guide_legend(override.aes = list(fill=NA)),
         linetype = guide_legend(override.aes = list(linewidth = 1.5))
  )


### UNITED STATES  #############################################################
S <- 2000
H <- 50
FutureZUS <- FutureZNoRep(H=H, OwnMod = TRUE, Samples = 
                            do.call("rbind",SamplesOwn_US_NoRepara), S=S,
                          NAge = 10)


#Future log death rates with cummulative sums
FutureLogMatArrayUS <- array(data = 0,dim = c(10,H,S),
                             dimnames = list("NewAgeInd"=1:10, "FCHorizon"=1:H, "It"=1:S))
for(s in 1:S){
  FutureLogMatArrayUS[,,s] <- log(LambdaMatUS)[,ncol(LambdaMatUS)]+ #Last obs rate
    t(apply(FutureZUS[,,s],1,cumsum)) #Cumsum
}


#Using melt, to transform array into long format for plotting
FutureLambdaVecIt_US <- 
  reshape2::melt(FutureLogMatArrayUS, value.name = "logRate") %>% 
  mutate(Rate = exp(logRate),
         Type="FC",
         Year = FCHorizon+2021) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")))

### Get observed Values Vector
ObservedValuesVec_US <- 
  LambdaVecUs %>%  #Take observed Vector and calculate PI's
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+"))) %>% 
  group_by(AgeGroup,Year) %>% 
  ggdist::point_interval(Rate, .width = c(0.9,0.95,0.99)) %>% 
  mutate(Type ="Obs")

### Calculate Quantiles ###
FutureValuesVecQuantile_US <- 
  FutureLambdaVecIt_US %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.9,0.95,0.99)) %>% #Quantile Calc
  full_join(x = filter(ObservedValuesVec_US, Year == max(Year)),
            y =. ) %>% #add last year, so that there is no space in plot 
  mutate(Type ="FC") %>% 
  full_join(x = filter(ObservedValuesVec_US),
            y =. )


### Freq Estimation ########
DeathMatUS <- matrix(LambdaVecUs$Deaths,
                     nrow = length(unique(LambdaVecUs$NewAgeInd)),
                     byrow = FALSE)

ExpMatUS <- matrix(LambdaVecUs$Pop,
                   nrow = length(unique(LambdaVecUs$NewAgeInd)),
                   byrow = FALSE)

LCfitUS <- fit(LC, Dxt = DeathMatUS, 
               Ext = ExpMatUS)

LCsimUS <- stats::simulate(LCfitUS,nsim = 500,h=50,gc.order = c(0, 1, 0))

FutureLambdaVecIt_Freq_US <- 
  reshape2::melt(LCsimUS$rates, value.name = "Rate") %>%
  rename("NewAgeInd"=Var1,
         TInd = Var2, It = Var3) %>%
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")),
         Year = TInd + 1980) 

FutureValuesVecQuantile_US_Freq <-
  FutureLambdaVecIt_Freq_US %>%
  group_by(AgeGroup,Year) %>%
  mutate("logRate"=log(Rate)) %>%
  ggdist::point_interval(Rate, .width = c(0.9)) %>%
  mutate(Type ="FC",
         Model="LC")%>% 
  full_join(x = filter(ObservedValuesVec_US,
                       .width == 0.9,
                       Year == max(Year)),
            y =. )

### Add layer after layer manually, looks nicer ....
AgeFilter <- c("20-29","30-39","40-49","50-59","60-69","70-79")
FutureValuesVecQuantile_US %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =Year, group = AgeGroup))+
  ggdist::geom_lineribbon(aes(y = Rate, ymin = .lower, 
                              ymax = .upper,fill = interaction(.width),
                              color = Type),
                          linewidth = 1.1, alpha = 0.8)+
  scale_fill_manual(values=c("#400040","#98579B","#b789b9"),
                    name =  "Own Model PI", labels = c("90%","95%", "99%"))+
  geom_vline(xintercept = 2022,linetype = "dashed", col="gray30")+
  facet_wrap(~AgeGroup, scales = "free")+
  geom_ribbon(data = filter(FutureValuesVecQuantile_US_Freq,
                            AgeGroup %in% AgeFilter),
              aes(x = Year, ymin = .lower, ymax = .upper),
              fill = "#48cae4", alpha = 0.5)+
  geom_line(data = filter(FutureValuesVecQuantile_US_Freq,
                          AgeGroup %in% AgeFilter), 
            aes(y = Rate, x = Year), col = "#03045e",
            linetype ="dashed", linewidth = 0.75,
            alpha = 0.7)+
  scale_color_manual(values =c("#2b012d","black","#03045e"),
                     labels =c("FC","Obs","LC"))+
  ylab("Death Rate")+xlab("Year")+
  scale_x_continuous(breaks = c(1990,2010,2030,2050,2070))+
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         colors = guide_legend(override.aes = list(fill=NA)))




######## ITALY  ################################################################ 
load(file = file.path(getwd(),"Results/SamplesIt_NoRep_2.RData")) # own Model

S <- 2000
H <- 50
FutureZIT <- FutureZNoRep(H=H, OwnMod = TRUE, Samples = 
                            do.call("rbind",SamplesOwn_It), S=S,
                          NAge = 10)


#Future log death rates with cummulative sums
FutureLogMatArrayIT <- array(data = 0,dim = c(10,H,S),
                             dimnames = list("NewAgeInd"=1:10, "FCHorizon"=1:H, "It"=1:S))
for(s in 1:S){
  FutureLogMatArrayIT[,,s] <- log(LambdaMatIt)[,ncol(LambdaMatIt)]+ #Last obs rate
    t(apply(FutureZIT[,,s],1,cumsum)) #Cumsum
}


#Using melt, to transform array into long format for plotting
FutureLambdaVecIt_IT <- 
  reshape2::melt(FutureLogMatArrayIT, value.name = "logRate") %>% 
  mutate(Rate = exp(logRate),
         Type="FC",
         Year = FCHorizon+2021) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")))

### Get observed Values Vector
ObservedValuesVec_IT <- 
  LambdaVecIt %>%  #Take observed Vector and calculate PI's
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+"))) %>% 
  group_by(AgeGroup,Year) %>% 
  ggdist::point_interval(Rate, .width = c(0.9,0.95,0.99)) %>% 
  mutate(Type ="Obs")

### Calculate Quantiles ###
FutureValuesVecQuantile_IT <- 
  FutureLambdaVecIt_IT %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.9,0.95,0.99)) %>% #Quantile Calc
  full_join(x = filter(ObservedValuesVec_IT, Year == max(Year)),
            y =. ) %>% #add last year, so that there is no space in plot 
  mutate(Type ="FC") %>% 
  full_join(x = filter(ObservedValuesVec_IT),
            y =. )


### Freq Estimation ########
DeathMatIT <- matrix(LambdaVecIt$Deaths,
                     nrow = length(unique(LambdaVecIt$NewAgeInd)),
                     byrow = FALSE)

ExpMatIT <- matrix(LambdaVecIt$Pop,
                   nrow = length(unique(LambdaVecIt$NewAgeInd)),
                   byrow = FALSE)

LCfitIT <- fit(LC, Dxt = DeathMatIT, 
               Ext = ExpMatIT)

LCsimIT <- stats::simulate(LCfitIT,nsim = 500,h=49,gc.order = c(0, 1, 0))

FutureLambdaVecIt_Freq_IT <- 
  reshape2::melt(LCsimIT$rates, value.name = "Rate") %>%
  rename("NewAgeInd"=Var1,
         TInd = Var2, It = Var3) %>%
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")),
         Year = TInd + 1980) 

FutureValuesVecQuantile_IT_Freq <-
  FutureLambdaVecIt_Freq_IT %>%
  group_by(AgeGroup,Year) %>%
  mutate("logRate"=log(Rate)) %>%
  ggdist::point_interval(Rate, .width = c(0.9)) %>%
  mutate(Type ="FC",
         Model="LC")%>% 
  full_join(x = filter(ObservedValuesVec_IT,
                       .width == 0.9,
                       Year == max(Year)),
            y =. )

### Add layer after layer manually, looks nicer ....
FutureValuesVecQuantile_IT %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =Year, group = AgeGroup))+
  ggdist::geom_lineribbon(aes(y = Rate, ymin = .lower, 
                              ymax = .upper,fill = interaction(.width),
                              color = Type),
                          linewidth = 1.1, alpha = 0.8)+
  scale_fill_manual(values=c("#355C41","#81A88D","#a7c2af"),
                    name =  "Own Model PI", labels = c("90%","95%", "99%"))+
  geom_vline(xintercept = 2022,linetype = "dashed", col="gray30")+
  facet_wrap(~AgeGroup, scales = "free")+
  geom_ribbon(data = filter(FutureValuesVecQuantile_IT_Freq,
                            AgeGroup %in% AgeFilter),
              aes(x = Year, ymin = .lower, ymax = .upper),
              fill = "#ef233c", alpha = 0.5)+
  geom_line(data = filter(FutureValuesVecQuantile_IT_Freq,
                          AgeGroup %in% AgeFilter), 
            aes(y = Rate, x = Year), col = "#9d0208",
            linetype ="dashed", linewidth = 0.9,
            alpha = 0.7)+
  scale_color_manual(values =c("#182b1e","black","#9d0208"),
                     labels =c("FC","Obs","LC"))+
  #scale_color_manual(values = "#03045e", name ="MeanFC")+
  ylab("Death Rate")+xlab("Year")+
  scale_x_continuous(breaks = c(1990,2010,2030,2050,2070))+
  guides(fill = guide_legend(override.aes = list(linetype = 0)),
         colors = guide_legend(override.aes = list(fill=NA)))





#### SPAIN OLD STUFF ###########################################################
## Comparison Jump, with Jump Component Removed #############################
#Same model but w/o the jump component included in the forecast
FutureZSp_noJump <- FutureZNoRep(H=H, OwnMod = TRUE, 
                                 Samples = do.call("rbind", SamplesOwn_Sp), S=S,
                                 NAge = 10, LC=TRUE)



#Future log death rates with cummulative sums
FutureLogMatArraySp_noJump <- array(data = 0,dim = c(10,H,S),
                             dimnames = list("NewAgeInd"=1:10, "FCHorizon"=1:H, "It"=1:S))
for(s in 1:S){
  FutureLogMatArraySp_noJump[,,s] <- log(LambdaMatSp)[,ncol(LambdaMatSp)]+
                                     t(apply(FutureZSp_noJump[,,s],1,cumsum))
}

#Transform into Long Format
FutureLambdaVecIt_Sp_noJump <- 
  reshape2::melt(FutureLogMatArraySp_noJump, value.name = "logRate") %>% 
  mutate(Rate = exp(logRate), #Calculate Death Rate
         Type="FC", #add some idicies for plotting
         Year = FCHorizon+2022) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")))

#Calculate Quantiles
FutureValuesVecQuantile_Sp_noJump <- 
  FutureLambdaVecIt_Sp_noJump %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(logRate, .width = c(0.5,0.8,0.99)) %>% #Quantile Calc
  mutate(Type ="FC") %>% 
  mutate(Model="LC")

#add Both together
FutureValuesQuantile_Both <- 
  FutureValuesVecQuantile_Sp %>% 
  mutate(Model = "Jump") %>% 
  bind_rows(., FutureValuesVecQuantile_Sp_noJump)


FutureValuesQuantile_Both %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =Year, group = AgeGroup))+
  ggdist::geom_lineribbon(aes(y = logRate,
                              ymin = .lower, ymax = .upper,
                              fill = interaction(.width,Model),
                              color = Model),
                          linewidth = 1.2, alpha = 0.8)+
  scale_fill_manual(values=c("#f0ead5","#dbcf9c","#CCBA72",
                             "#f48c06","#d00000","#370617"))+
  geom_vline(xintercept = 2022,linetype = "dashed", col="gray30")+
  facet_wrap(~AgeGroup, scales = "free")+
  ylab("Death Rate")+xlab("Year")+
  scale_x_continuous(breaks = c(1990,2010,2030,2050,2070))
  

### Compare the quantiles of both mortality improvement rates
FutureZSp_noJump[5:8,,] %>% apply(., c(1,2), quantile, 0.95) %>% 
  rowMeans()

FutureZSp[5:8,,] %>% apply(., c(1,2), quantile, 0.95) %>% 
  rowMeans()

## the 95% quantile is slightly higher, though 99% quantile substantially higher
FutureZSp_noJump[5:8,,] %>% apply(., c(1,2), quantile, 0.99) %>% 
  apply(., 1, summary)

FutureZSp[5:8,,] %>% apply(., c(1,2), quantile, 0.99) %>% 
  apply(., 1, summary)

#however, what about the minimum i.e. the 1% quantile is substantially lower,
#since whenever there is a jump is has to go down afterwards
FutureZSp_noJump[5:8,,] %>% apply(., c(1,2), quantile, 0.01) %>% 
  apply(., 1, summary)

FutureZSp[5:8,,] %>% apply(., c(1,2), quantile, 0.01) %>% 
  apply(., 1, summary)

##They might cancel out 

### Comparison with STAN MOMO ##############################
library(StanMoMo)
options(mc.cores = parallel::detectCores())
LC_Stan <- lc_stan(death = matrix(LambdaVecSp$Deaths,
                                  nrow = 10, byrow = FALSE),
                   exposure =  matrix(LambdaVecSp$Pop, 
                                      nrow = 10, byrow = FALSE),
                   forecast = 50, family = "poisson")

ForecastsLC_Stan <- rstan::extract(LC_Stan)$mufor

FutureValuesIt_StanLC <- 
  ForecastsLC_Stan %>% t() %>% 
  data.frame() %>% 
  mutate("NewAgeInd"=rep(1:10, 50),
         "FCHorizon"=rep(1:50,each = 10),
         .before=1) %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "It",values_to = "Rate") %>% 
  mutate("It"=rep(1:4000, 500)) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")),
         "Year"=FCHorizon + 2022)

FutureValuesQuantile_Sp_Stan <- 
  FutureValuesIt_StanLC %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(logRate, .width = c(0.5,0.8,0.99)) %>% 
  mutate(Type ="FC",
         Model="Stan") 

#add Both together
FutureValuesQuantile_Both2 <- 
  FutureValuesVecQuantile_Sp %>% 
  mutate(Model = "Jump") %>% 
  bind_rows(., FutureValuesQuantile_Sp_Stan)


FutureValuesQuantile_Both2 %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =Year, group = AgeGroup))+
  ggdist::geom_lineribbon(aes(y = logRate,
                              ymin = .lower, ymax = .upper,
                              fill = interaction(.width,Model),
                              color = Model),
                          linewidth = 1.2, alpha = 0.8)+
  scale_fill_manual(values=c("#f0ead5","#dbcf9c","#CCBA72",
                             "#f48c06","#d00000","#370617"))+
  geom_vline(xintercept = 2022,linetype = "dashed", col="gray30")+
  facet_wrap(~AgeGroup, scales = "free")+
  ylab("Death Rate")+xlab("Year")+
  scale_x_continuous(breaks = c(1990,2010,2030,2050,2070))

### Same shit different story....


#### Compare with SVD Estimation #######################################
#Calcualate Parameters via SVD 
DeathMatSp <- matrix(LambdaVecSp$Deaths,
                     nrow = length(unique(LambdaVecSp$NewAgeInd)),
                     byrow = FALSE)

ExpMatSp <- matrix(LambdaVecSp$Pop,
                   nrow = length(unique(LambdaVecSp$NewAgeInd)),
                   byrow = FALSE)

InSampRange <- 1:42

#### SVD Estimation ##########
MortMatSp <- DeathMatSp[,InSampRange]/ExpMatSp[,InSampRange]

h <- ncol(MortMatSp)
ax <- apply(MortMatSp, 1, function(x) log(prod(x^(1/h))))
ZMat_SVD_Sp <- log(MortMatSp)-ax

#Do SVD 
SVD <- svd(ZMat_SVD_Sp)

# Extract first principal component
sumu <- sum(SVD$u[,1])

#Parameter Estimates
bx <- SVD$u[,1]/sumu
kt <- SVD$d[1] * SVD$v[,1] * sumu #sigma*singularvalue*sumbx

RWDModel <- forecast::Arima(kt, order = c(0,1,0), include.drift = TRUE)

#RWD Forecast
KTMat <- matrix(data=0,nrow = 50, ncol = 500)
for( s in 1:500){
  Error <- rnorm(n = 50,mean = 0, sd = sqrt(RWDModel$sigma2))
  KTMat[,s] <- RWDModel$x[42]+cumsum(RWDModel$coef+Error)
}

LogMatArrayFC <- array(data = 0,dim = c(10,50,500),
                       dimnames = list("NewAgeInd"=1:10,
                                       "FCHorizon"=1:50,
                                       "It"=1:500)) 
for(h in 1:50){
  for(a in 1:10){
    LogMatArrayFC[a,h,] <- ax[a]+bx[a]*KTMat[h,]
  }
}

FutureLambdaVecIt_SVD <- 
  reshape2::melt(LogMatArrayFC, value.name = "logRate") %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+"))) %>% 
  mutate(Type="FC",
         Year = FCHorizon+2022,
         Rate = exp(logRate))

FutureValuesVecQuantile_Sp_SVD <- 
  FutureLambdaVecIt_SVD %>%
  group_by(AgeGroup,Year) %>%
  mutate("logRate"=log(Rate)) %>%
  ggdist::point_interval(logRate, .width = c(0.5,0.8,0.99)) %>%
  mutate(Type ="FC",
         Model="SVD")


FutureValuesQuantile_Both3 <- 
  FutureValuesVecQuantile_Sp %>% 
  mutate(Model = "Jump") %>% 
  bind_rows(., FutureValuesVecQuantile_Sp_SVD)


FutureValuesQuantile_Both3 %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =Year, group = AgeGroup))+
  ggdist::geom_lineribbon(aes(y = logRate,
                              ymin = .lower, ymax = .upper,
                              fill = interaction(.width,Model),
                              color = Model),
                          linewidth = 1.2, alpha = 0.8)+
  scale_fill_manual(values=c("#806E26","#CCBA72","#f0ead5",
                             "#574f7d","#95adbe","#e0f0ea"))+
  geom_vline(xintercept = 2022,linetype = "dashed", col="gray30")+
  facet_wrap(~AgeGroup, scales = "free")+
  ylab("Death Rate")+xlab("Year")+
  scale_x_continuous(breaks = c(1990,2010,2030,2050,2070))


