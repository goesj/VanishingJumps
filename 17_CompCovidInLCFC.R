## Test: is Covid Jump in 95% Interval of Lee Carter Model? ##
source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 

load(file = file.path(getwd(),"Data/CovidData.RData")) # Load Data

##### SPAIN ###################################################################
# first test on frequentist Model 
DeathMatSp <- matrix(LambdaVecSp$Deaths,
                     nrow = length(unique(LambdaVecSp$NewAgeInd)),
                     byrow = FALSE)

ExpMatSp <- matrix(LambdaVecSp$Pop,
                   nrow = length(unique(LambdaVecSp$NewAgeInd)),
                   byrow = FALSE)

InSampRange <- 1:42
OOSRange <- tail(InSampRange,1):42

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

### Own Simulation of Future Rates
forecast::forecast(forecast::Arima(kt, order = c(0,1,0), 
                                   include.drift = TRUE),
                   bootstrap = TRUE, h = 50)

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

FutureLambdaVecIt_Freq <- 
  reshape2::melt(LogMatArrayFC, value.name = "logRate") %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+"))) %>% 
  mutate(Type="FC",
        Year = FCHorizon+2022,
        Rate = exp(logRate))

FutureValuesVecQuantile_Sp_Freq <- 
  FutureLambdaVecIt_Freq %>%
    group_by(AgeGroup,Year) %>%
    mutate("logRate"=log(Rate)) %>%
    ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>%
    mutate(Type ="FC",
           Model="Freq")

# ## Frequentist Estimation
# library(StMoMo)
# constLC <- function(ax, bx, kt, b0x, gc, wxt, ages){
#   c1 <- mean(kt[1, ], na.rm = TRUE)
#   c2 <- sum(bx[, 1], na.rm = TRUE)
#   list(ax = ax + c1 * bx, bx = bx / c2, kt = c2 * (kt - c1))
# }
# LC <- StMoMo(link = "log", staticAgeFun = TRUE, 
#              periodAgeFun = "NP",
#             constFun = constLC)
# 
# LCfit <- fit(LC, Dxt = DeathMatSp[,InSampRange], 
#              Ext = ExpMatSp[,InSampRange])
# 
# LCsim <- stats::simulate(LCfit,nsim = 500,h=10,gc.order = c(0, 1, 0))
# 
# FutureLambdaVecIt_Freq <- 
#   reshape2::melt(LCsim$rates, value.name = "Rate") %>%
#   rename("NewAgeInd"=Var1,
#          TInd = Var2, It = Var3) %>%
#   mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
#                            levels = 1:10,
#                            labels = c(paste0(paste0(seq(0,80,10),sep="-"),
#                                              seq(9,89,10)),"90+")),
#          Year = TInd + 1980) 
# 
# FutureValuesVecQuantile_Sp_Freq <-
#   FutureLambdaVecIt_Freq %>%
#   group_by(AgeGroup,Year) %>%
#   mutate("logRate"=log(Rate)) %>%
#   ggdist::point_interval(logRate, .width = c(0.5,0.8,0.99)) %>%
#   mutate(Type ="FC",
#          Model="Freq")

# 
# qxt <- DeathMatSp / ExpMatSp
# 
# plot(InSampRange, qxt[7,InSampRange], xlim = range(LCfit$years, LCsim$years),
#      ylim = range(qxt[7, ], LCsim$rates[7, , 1:50]), type = "l")
# lines(OOSRange, qxt[7,OOSRange],col="red")
# matlines(LCsim$years, LCsim$rates[7, , 1:100], type = "l", lty = 1)
# 
# library("fanplot")
# probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
# matplot(LCfit$years, t(qxt[c(5:7), InSampRange]),
#         xlim = range(LCfit$years, LCsim$years), ylim = range(0.001, 0.02), 
#         pch = 20, col = "black",
#         log = "y", xlab = "year", ylab = "mortality rate")
# fan(t(LCsim$rates[5, , ]), start = OOSRange[2], probs = probs, n.fan = 4,
#        fan.col = colorRampPalette(c("black", "white")), ln = NULL)
# fan(t(LCsim$rates[6, , ]), start = OOSRange[2], probs = probs, n.fan = 4,
#        fan.col = colorRampPalette(c("red", "white")), ln = NULL)
# fan(t(LCsim$rates[7, , ]), start = OOSRange[2], probs = probs, n.fan = 4,
#        fan.col = colorRampPalette(c("blue", "white")), ln = NULL)
# lines(OOSRange,qxt[5, OOSRange],col="darkgreen", lwd =2)
# lines(OOSRange,qxt[6, OOSRange],col="darkgreen", lwd =2)
# lines(OOSRange,qxt[7, OOSRange],col="darkgreen", lwd =2)
# 
# LCfit <- fit(LC, Dxt = DeathMatSp[,InSampRange], Ext = ExpMatSp[,InSampRange])
# LCsim <- stats::simulate(LCfit,nsim = 1000,h=50, gc.order = c(0, 1, 0))
# 
# FutureLambdaVecIt_Freq <- 
#   reshape2::melt(LCsim$rates, value.name = "Rate") %>% 
#   rename("NewAgeInd"=Var1,
#          TInd = Var2, It = Var3) %>% 
#   mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
#                            levels = 1:10,
#                            labels = c(paste0(paste0(seq(0,80,10),sep="-"),
#                                              seq(9,89,10)),"90+")))
# FutureValuesVecQuantile_Sp_Freq <- 
# FutureLambdaVecIt_Freq %>% 
#   group_by(AgeGroup,TInd) %>% 
#   mutate("logRate"=log(Rate)) %>% 
#   ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>% 
#   mutate(Type ="FC",
#          Model="Freq") 

#### Compare FC of Diff LC model with those of Freq LC Model ###################
FutureValuesVecQuantile_Sp <- 
  FutureValuesVecQuantile_Sp %>% 
  mutate(Model="Jump") %>% 
  bind_rows(., FutureValuesVecQuantile_Sp_Freq)


AgeFilter <- c("20-29","40-49","50-59","60-69")

ObservedValuesVec2_Sp %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  mutate(Model = "LC") %>% 
  ggplot(aes(x =TInd, group = AgeGroup))+
  geom_line(aes(y = Rate),linewidth = 1.2)+
  ggdist::geom_lineribbon(data = filter(FutureValuesVecQuantile_Sp, AgeGroup %in% AgeFilter),
                          aes(y = Rate,
                              ymin = .lower, ymax = .upper,
                              fill = interaction(.width,Model)),
                          linewidth = 1.2)+
  scale_fill_manual(values=c("#f0ead5","#dbcf9c","#CCBA72",
                             "#12492f","#f56038","#0a2f35"))+
  facet_wrap(~AgeGroup, scales = "free")
  

#### Comparison with own LC Model on Diff Data #################################
load(file = file.path(getwd(), "Results/SamplesSp_LC.RData"))

FutureZSp_LC <- FutureZNoRep(H=H, OwnMod = TRUE, 
                             Samples = do.call("rbind", Samples_LC_Sp), S=S,
                             NAge = 10, LC = TRUE)


FutureLogMatArraySp_LC <- array(data = 0,dim = c(10,H,S),
                                dimnames = list("NewAgeInd"=1:10, "FCHorizon"=1:H, "It"=1:S))

for(s in 1:S){
  FutureLogMatArraySp_LC[,,s] <- log(LambdaMatSp)[,42]+t(apply(FutureZSp_LC[,,s],1,cumsum))
}


InSampRange <- 1:42

FutureLambdaVecIt_LC <- 
  reshape2::melt(FutureLogMatArraySp_LC, value.name = "logRate") %>%
  mutate(Rate = exp(logRate),
         Type="FC",
         TInd = FCHorizon+42) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")))


FutureValuesVecQuantile_Sp_LC <-
  FutureLambdaVecIt_LC %>%
  group_by(AgeGroup,TInd) %>%
  mutate("logRate"=log(Rate)) %>%
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>%
  mutate(Type ="FC",
         Model="Freq")

FutureValuesVecQuantile_Sp <- 
  FutureValuesVecQuantile_Sp %>% 
  filter(Model =="Jump") %>% 
  bind_rows(., FutureValuesVecQuantile_Sp_LC)


AgeFilter <- c("20-29","40-49","50-59","60-69")

ObservedValuesVec2_Sp %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  mutate(Model = "LC") %>% 
  ggplot(aes(x =TInd, group = AgeGroup))+
  geom_line(aes(y = Rate),linewidth = 1.2)+
  ggdist::geom_lineribbon(data = filter(FutureValuesVecQuantile_Sp, AgeGroup %in% AgeFilter),
                          aes(y = Rate,
                              ymin = .lower, ymax = .upper,
                              fill = interaction(.width,Model)),
                          linewidth = 1.2)+
  scale_fill_manual(values=c("#f0ead5","#dbcf9c","#CCBA72",
                             "#12492f","#f56038","#0a2f35"))+
  scale_color_manual(values =c("#806E26","blue"))+
  facet_wrap(~AgeGroup, scales = "free")



#### Comparison with Stan ################
#### in Stan, with FC's #####
StanDatSp <- list("T"=length(unique(LambdaVecSp$Year)),
                  "A"=max(LambdaVecSp$NewAgeInd),
                  "y"=LambdaVecSp$Deaths,
                  "E"=LambdaVecSp$Pop,
                  "TFor"=50)


options(mc.cores = parallel::detectCores())
StanLC_Sp <- rstan::stan(file.path(getwd(),"Stan Code/LeeCarter_SumZeroTime.stan"),
                         data=StanDatSp,
                         chains = 4, iter=2500,warmup=1500, 
                         save_warmup=FALSE, thin=1,
                         control = list(adapt_delta = 0.81))

summary(StanLC_Sp)$summary %>% round(., 4)

LogDeathRatesStanFC <- rstan::extract(StanLC_Sp)$mufor

FutureValuesIt_StanLCSp <- 
  LogDeathRatesStanFC %>% t() %>% 
  data.frame() %>% 
  mutate("NewAgeInd"=rep(1:10, 50),
         "FCHorizon"=rep(1:50,each = 10),
         .before=1) %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "It",values_to = "logRate") %>% 
  mutate("It"=rep(1:2000, 500)) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")),
         "TInd"=FCHorizon + 42)

FutureValuesVecQuantile_Sp_Stan <- 
  FutureValuesIt_StanLCSp %>% 
  group_by(AgeGroup,TInd) %>% 
  mutate("Rate"=exp(logRate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>% 
  mutate(Type ="FC",
         Model="Stan") 

FutureValuesVecQuantile_Sp <- 
  FutureValuesVecQuantile_Sp %>% 
  mutate(Model="Jump") %>% 
  bind_rows(., FutureValuesVecQuantile_Sp_Stan)


AgeFilter <- c("20-29","40-49","50-59","60-69")

ObservedValuesVec2_Sp %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  mutate(Model = "LC") %>% 
  ggplot(aes(x =TInd, group = AgeGroup))+
  geom_line(aes(y = Rate),linewidth = 1.2)+
  ggdist::geom_lineribbon(data = filter(FutureValuesVecQuantile_Sp, AgeGroup %in% AgeFilter),
                          aes(y = Rate,
                              ymin = .lower, ymax = .upper,
                              fill = interaction(.width,Model)),
                          linewidth = 1.2)+
  scale_fill_manual(values=c("#f0ead5","#dbcf9c","#CCBA72",
                             "#12492f","#f56038","#0a2f35"))+
  scale_color_manual(values =c("#806E26","blue"))+
  facet_wrap(~AgeGroup, scales = "free")


### Karims Code #############################################################
library(StanMoMo)
options(mc.cores = parallel::detectCores())
LC_Stan <- lc_stan(death = matrix(LambdaVecSp$Deaths,
                                  nrow = 10, byrow = FALSE),
                   exposure =  matrix(LambdaVecSp$Pop, 
                                      nrow = 10, byrow = FALSE),
                   forecast = 50, family = "poisson")

summary(LC_Stan)$summary %>% round(., 4)

ForecastsLC_Stan <- rstan::extract(LC_Stan)$mufor

FutureValuesIt_StanLC <- 
  ForecastsLC_Stan %>% t() %>% 
  data.frame() %>% 
  mutate("NewAgeInd"=rep(1:10, 50),
         "FCHorizon"=rep(1:50,each = 10),
         .before=1) %>% 
  pivot_longer(., cols = 3:ncol(.), names_to = "It",values_to = "logRate") %>% 
  mutate("It"=rep(1:2000, 500)) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                             seq(9,89,10)),"90+")),
         "TInd"=FCHorizon + 42)

FutureValuesVecQuantile_Sp_Stan <- 
  FutureValuesIt_StanLCSp %>% 
  group_by(AgeGroup,TInd) %>% 
  mutate("Rate"=exp(logRate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.5,0.8,0.99)) %>% 
  mutate(Type ="FC",
         Model="Stan") 

FutureValuesVecQuantile_Sp <- 
  FutureValuesVecQuantile_Sp %>% 
  mutate(Model="Jump") %>% 
  bind_rows(., FutureValuesVecQuantile_Sp_Stan)


AgeFilter <- c("20-29","40-49","50-59","60-69")

ObservedValuesVec2_Sp %>% 
  filter(AgeGroup %in% AgeFilter) %>% 
  mutate(Model = "LC") %>% 
  ggplot(aes(x =TInd, group = AgeGroup))+
  geom_line(aes(y = Rate),linewidth = 1.2)+
  ggdist::geom_lineribbon(data = filter(FutureValuesVecQuantile_Sp, AgeGroup %in% AgeFilter),
                          aes(y = Rate,
                              ymin = .lower, ymax = .upper,
                              fill = interaction(.width,Model)),
                          linewidth = 1.2)+
  scale_fill_manual(values=c("#f0ead5","#dbcf9c","#CCBA72",
                             "#12492f","#f56038","#0a2f35"))+
  scale_color_manual(values =c("#806E26","blue"))+
  facet_wrap(~AgeGroup, scales = "free")

