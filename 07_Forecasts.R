### Script to Forecast future Death Rates ######################################
library(nimble);library(tidyverse); library(reshape2); library(rstan);
library(StMoMo);library(ggdist)

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 


#load Samples
load(file = file.path(getwd(),"Results/SamplesSp.RData")) # own Model
load(file = file.path(getwd(),"Results/SamplesUS.RData")) # own Model
load(file = file.path(getwd(),"Results/SamplesIt.RData")) # own Model


#load Data
load(file = file.path(getwd(),"Data/CovidData.RData")) # Load Data

#  Exemplary for Country of Spain
## For selection of other country, change Samples and Observed Data 

########### 1. Own Model only ##################################################

########## 1.1 Forecast Mortaltiy Improvemnt ###################################
S <- 1500  #select number of samples from Posterior predictive
H <- 50    #Foreacst ahead period

#Calculate Future Mortality Improvement Rates
FutureZSp <- FutureZ(H=H, OwnMod = TRUE,  #binary for own model vs Liu,Li
                     Samples = do.call("rbind", SamplesOwn_Sp), 
                     S=S,
                     NAge = 10) #Number of age groups

###### 1.2. Future log Death rates #############################################

#Calculation of future log death rates via cummulative sum 

#First create empty Array
FutureLogMatArraySp <- array(data = 0,dim = c(10,H,S),
                             dimnames = list("NewAgeInd"=1:10, 
                                             "FCHorizon"=1:H, "It"=1:S))
#Calculation of Future log death rates
for(s in 1:S){
  FutureLogMatArraySp[,,s] <- log(LambdaMatSp)[,ncol(LambdaMatSp)]+
    t(apply(FutureZSp[,,s],1,cumsum))
}

#Transform Forecasts into Long Format
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
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) 

###### 1.3 Plot Results ########################################################
#Get observed Values into same format. 
#This has to be done, so that there is no space between observed values and
#forecasted values when plotting
ObservedValuesVec_Sp <- 
  LambdaVecSp %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+"))) %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>% 
  mutate(Type ="Obs")

#Add Observed values to Quantiles
FutureValuesVecQuantile_Sp <- 
  FutureValuesVecQuantile_Sp %>% 
  mutate(Type ="FC") %>% 
  full_join(x = filter(ObservedValuesVec_Sp, Year == max(Year)),
            y =. )

#Select suitable age range
AgeFilter <- c("40-49","50-59","60-69","70-79")

#Plot results
FutureValuesVecQuantile_Sp %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =Year))+
  ggdist::geom_lineribbon(aes(y = Rate, ymin = .lower, ymax = .upper,
                              group = AgeGroup,
                              #fill = forcats::fct_rev(ordered(.width)) #standard setting
  ), 
  linewidth = 1.2,col = "#834333")+
  geom_line(data=filter(ObservedValuesVec_Sp, AgeGroup %in% AgeFilter),
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



# ######### 2. Comparison Own Vs StMoMo ########################################

#### 2.1 Calculate Estimates of StMoMo #########################################
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

# Fit Model
LCfit <- fit(LC, Dxt = DeathMatSp, 
             Ext = ExpMatSp)

######## 2.2. Calcualte Forecasts of StMoMo ####################################
LCsim <- stats::simulate(LCfit,nsim = 1500,
                         h=50,gc.order = c(0, 1, 0))

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

########### 2.3 Plot Both togethers ############################################
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

### Add layer after layer, 
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


