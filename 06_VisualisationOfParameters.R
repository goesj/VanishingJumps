### Visualisation of Results ##############
library(nimble);library(tidyverse);library(ggdist); library(gghighlight);
library(cowplot)


source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 


load(file = file.path(getwd(),"Results/SamplesUS.RData"))
load(file = file.path(getwd(),"Results/SamplesIt.RData"))
load(file = file.path(getwd(),"Results/SamplesSp.RData"))


########### SAMPLING RESULTS####################################################
theme_set(theme_minimal(base_size = 10)) #set global theme
#### 1.1.1 Data Preparation ####################################################
SamplesItSingle <- do.call(rbind, 
                           SamplesOwn_It) %>% 
  data.frame() %>% 
  mutate(Country="Italy")

SamplesUSSingle <- do.call(rbind, 
                           SamplesOwn_US) %>% 
  data.frame() %>% 
  mutate(Country="United States")

SamplesSpSingle <- do.call(rbind, 
                           SamplesOwn_Sp) %>% 
  data.frame() %>% 
  mutate(Country = "Spain")


SamplesTotal <- bind_rows(SamplesItSingle,
                          SamplesUSSingle,
                          SamplesSpSingle)

######## 1.2. Parameter Density Plot ###########################################
SampleDensPlot <- function(Samples, Parameter = p){
  if(!Parameter %in% c("p","a","muY","sdY")){
    print("please select one of p, a, muY or sdY")
  } else  {
    Plot <- Samples %>% 
      select(a,p,muY,sdY,Country) %>% 
      #Attention: Factor level is turned around, so that last factor (Italy) is on top
      mutate("Country" = factor(Country, 
                                levels = c("United States","Spain","Italy"))) %>% 
      ggplot(aes(x = a, y = Country))+    #change here parameter to sample
      ggdist::stat_slabinterval(aes(fill_ramp=after_stat(level),
                                    fill=Country,color=Country),
                                geom="slabinterval",
                                point_interval =  "mean_qi",
                                slab_type = "pdf", 
                                .width = c(.50, .85, 0.99),
                                alpha=0.8,
                                slab_linewidth=0.6,
                                point_size = 2, 
                                linewidth=2,
                                shape=21)+
      scale_fill_manual(values =c("#98579B","#c9b045","#81A88D"))+ 
      ggdist::scale_fill_ramp_discrete(range = c(0.2, 0.95))+
      scale_color_manual(values=c("#400040","#806E26","#355C41"))+
      ggdist::scale_slab_color_discrete(c("#400040","#806E26","#355C41"),
                                        aesthetics = "slab_color",
                                        na.translate = FALSE)+
      theme(axis.text = element_text(size = 24), #change size of axis text
            axis.title = element_text(size = 25, face = "bold"),
            legend.position = "none",
      )+
      ylab("")
    return(Plot)
  }
}


#Create plot select one of a, p, muY, and sdY for Parameter
SampleDensPlot(Samples = SamplesTotal, Parameter =  "p")



########## 1.3. Jump Visualisation #############################################
#Plot for Visualisation of Jump Occurence
SamplesTotal %>% 
  select(contains("N_t"),Country) %>%
  mutate(Country = factor(Country, levels = c("Italy","Spain","United States"),
                          labels = c("a) Italy","b) Spain","c) United States"))) %>% 
  group_by(Country) %>% 
  summarise(across(.cols = where(is.numeric),.fns = mean)) %>% 
  pivot_longer(., cols = 2:ncol(.), names_to="Time", values_to = "PostMean") %>% 
  mutate(TimeFac = factor(Time, levels = unique(Time), ordered = TRUE)) %>% 
  ggplot(aes(x=TimeFac, group = Country, col = Country))+
  geom_segment(aes(x = TimeFac, xend = TimeFac, y = 0, yend = PostMean), linewidth = 1.8)+
  geom_point(aes(y = PostMean), size = 3, pch = 18)+
  facet_wrap(~Country, ncol = 1)+
  scale_color_manual(values=c("#81A88D","#c9b045","#98579B"))+
  scale_x_discrete(labels = 1981:2022)+
  ylab("Posterior Mean")+xlab("")+
  theme(legend.position = "none",
        axis.text = element_text(size =20), #change size of axis text
        axis.title = element_text(size = 22),
        strip.text.x=element_text(face="bold", size = 22),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
        )


############ 1.4. Age Pattern of Mortality Jump ################################

#Plot to create the age pattern of Mortality Jump for all three Countries
#https://r-graph-gallery.com/web-line-chart-small-multiple-all-group-greyed-out.html
font <- "serif"

SamplesTotal %>% 
  select_if(grepl(paste0(c("Country","betaJump"), collapse="|"), names(.))) %>% 
  mutate(Country = factor(Country, levels = c("Italy","Spain","United States"),
                          labels = c("a) Italy","b) Spain","c) United States"))) %>% 
  pivot_longer(., cols = 1:10, names_to = "BetaVal",values_to = "Val") %>%  #transform into long 
  mutate(BetaVal = fct_inorder(as.factor(BetaVal),ordered = NA)) %>% 
  group_by(Country,BetaVal) %>%  #group by
  reframe(mean = mean(Val),
          PIL = quantile(Val, probs=0.25),
          PIU = quantile(Val, probs=0.75)) %>% 
  ggplot(aes(x = BetaVal, group = Country)) +
  geom_line(aes(y =mean, col=Country), linewidth = 0.75)+
  geom_ribbon(aes(ymin = PIL, ymax = PIU, fill = Country), alpha = 0.3)+
  gghighlight::gghighlight(use_direct_label = FALSE,
              unhighlighted_params = list(colour = alpha("grey85", 1)))+
  facet_wrap(~  Country, nrow = 3)+
  scale_fill_manual(values=c("#81A88D","#c9b045","#98579B"))+
  scale_color_manual(values=c("#182b1e","#834333","#2b012d"))+
  ylab("Posterior Value")+
  xlab("Age Group")+
  scale_x_discrete(labels = c(paste0(paste0(seq(0,80,10),sep="-"),
                                     seq(9,89,10)),"90+"))+  
  theme(
        axis.text = element_text(size = 20), #change size of axis text
        axis.title = element_text(size = 22),
        strip.text.x=element_text(face="bold", size = 22), #Text of Countries
        legend.position = "none",
        axis.text.x.bottom = element_text(vjust=2)
  )


############ 1.5 All Parameters in one #########################################

## Creation of plot that visualizes all Parameters into one plot

##### 1.5.1 United States ######################################################
#Starting with the Age Effects beta_x, beta_x^J

# 1.5.1.1 Beta Parameters 
AgeVecCovid <-  c(paste0(paste0(seq(0,80,10),sep="-"),seq(9,89,10)),"90+")
BetaParamUS <- SamplesUSSingle %>% 
  select(grep("beta",colnames(.))) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate("AgeGroup"=rep(AgeVecCovid,2*nrow(SamplesUSSingle)),#Create Age Group Variable
         "Type"=rep(c(rep("Norm",10),rep("Jump",10)),nrow(SamplesUSSingle))
  ) %>% mutate("Type"=factor(Type, 
                        levels = c("Norm","Jump"),
                        labels = c(expression(beta[x]),expression(beta[x]^{(J)}))),
          "AgeGroup"=factor(AgeGroup, levels = AgeVecCovid)) %>% #Change to factor 
  ggplot(aes(x = AgeGroup, y = Val,group = Type))+
  ggdist::stat_slabinterval(
    aes(group = after_stat(level), fill = after_stat(level), size = NULL),
    geom = "lineribbon",
    .width = c(0.5, 0.8, .95), show_slab = FALSE,
    show.legend = NA, color ="#2b012d"   
  )+
  facet_wrap(~Type,nrow = 2,labeller = label_parsed)+
  scale_fill_manual(values = c("#f7ebfd","#b789b9","#98579b")) +
  ylab("Posterior Mean")+xlab("Age Group")+
  theme(legend.position = "none",
        axis.text = element_text(size = 20), #change size of axis text
        axis.title = element_text(size = 22),
        axis.title.y = element_text(vjust=1.3),
        axis.text.x = element_text(vjust=2),
        strip.text.x=element_text(face="bold", size = 22),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


#  1.5.1.2 Other Parameters 

#select Other Parameters to plot
OtherParams <- 
  c("muY","sdY","a","p","sigma_time","drift","sigma_eps")

#Get Indicies of correct Parameters
Ind <- sapply(OtherParams, 
              function(y){grep(paste0("\\<",y), #match beginning of word
                               x = colnames(SamplesUSSingle))}) %>% unlist() %>% 
  unique() 

OtherParamPlotUS <- 
  SamplesUSSingle %>% 
  select(all_of(Ind)) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate("Type" = factor(rep(
    c(rep("Jump Intensity",2),rep("Shock Params",2),rep("Error RW",1),
      rep("Drift",1),"Error Term"),
    nrow(SamplesUSSingle)),
    levels=c("Jump Intensity","Shock Params","Error RW","Drift",
             "Error Term"))) %>%
  mutate(Param = factor(Param,
                        levels = unique(Param),
         labels = c(expression(mu[Y]),expression(sigma[Y]), expression(" a "),
                    expression(p), expression(sigma[xi]),
                    expression(d),expression(sigma[r])))) %>% 
  # # add y min a y max by group
  mutate("ymin"=rep(c(rep(0.25,2),rep(0,2),rep(0.1,1), 
                      rep(-0.2,1),rep(0.02,1)),nrow(SamplesUSSingle)),
        "ymax"=rep(c(rep(2.25,2),rep(0.5,2),rep(0.3,1), 
                     rep(0,1),rep(0.025,1)),nrow(SamplesUSSingle))) %>%
  ggplot(aes(y = Val, x = Param, group = Type))+
  ggdist::stat_slabinterval(
    aes(y = Val,fill =Param, colour = after_stat(level), size = NULL),
    geom = "interval",
    show_point = FALSE, .width = c(0.5, 0.8, 0.95), show_slab = FALSE,
    show.legend = NA
  )+
  stat_summary( #Add mean to plot
    geom = "point",fun = "mean",col = "#200020",size = 3 #200020
  )+
  facet_wrap(~Type, scales = "free",shrink = TRUE, nrow = 1)+
  scale_x_discrete(labels = function(l) parse(text=l))+ #transforms text into expression
  scale_color_manual(values = c("#eaddeb","#b789b9","#98579b"))+
  ggdist::scale_slab_alpha_discrete(range = c(0.1, 1))+
  ylab("Posterior Mean")+ xlab("Other Parameters")+
  geom_blank(aes(y = ymin))+geom_blank(aes(y = ymax))+
  theme(legend.position = "none",
        axis.text = element_text(size = 25), #change size of axis text
        axis.title = element_text(size = 22),
        axis.title.y = element_text(vjust=1.3),
        axis.text.x = element_text(vjust=1),
        strip.text.x = element_blank(), #remove headers of facet wrap
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
        )


#1.5.1.3 Add both plots together 
## Add both plots together via Cowplot
cowplot::plot_grid(BetaParamUS, OtherParamPlotUS, byrow = FALSE, nrow = 2,
                   rel_heights = c(1.4,1)) #making first plot 1.4 times as high


######### 1.5.2 Spain ##########################################################
## 1.5.2.1 Beta Parameters 
AgeVecCovid <-  c(paste0(paste0(seq(0,80,10),sep="-"),seq(9,89,10)),"90+")
BetaParamSp <- SamplesSpSingle %>% 
  select(grep("beta",colnames(.))) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate("AgeGroup"=rep(AgeVecCovid,2*nrow(SamplesSpSingle)),#Create Age Group Variable
         "Type"=rep(c(rep("Norm",10),rep("Jump",10)),nrow(SamplesSpSingle))
  ) %>% mutate("Type"=factor(Type, 
                             levels = c("Norm","Jump"),
                             labels = c(expression(beta[x]),expression(beta[x]^{(J)}))),
               "AgeGroup"=factor(AgeGroup, levels = AgeVecCovid)) %>% #Change to factor 
  ggplot(aes(x = AgeGroup, y = Val,group = Type))+
  ggdist::stat_slabinterval(
    aes(group = after_stat(level), fill = after_stat(level), size = NULL),
    geom = "lineribbon",
    .width = c(0.5, 0.8, .95), show_slab = FALSE,
    show.legend = NA, color ="#834333"
  )+
  facet_wrap(~Type,nrow = 2,labeller = label_parsed)+
  scale_fill_manual(values = c("#fdf2c5","#ecdc94","#c9b045")) +
  ylab("Posterior Mean")+xlab("Age Group")+
  theme(legend.position = "none",
        axis.text = element_text(size = 20), #change size of axis text
        axis.title = element_text(size = 22),
        axis.title.y = element_text(vjust=1.3),
        axis.text.x = element_text(vjust=2),
        strip.text.x=element_text(face="bold", size = 22),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


## 1.5.2.2 Other Parameters 
OtherParams <- 
  c("muY","sdY","a","p","sigma_time","drift","sigma_eps")
#Get Indicies of correct Parameters
Ind <- sapply(OtherParams, 
              function(y){grep(paste0("\\<",y), #match beginning of word
                               x = colnames(SamplesSpSingle))}) %>% unlist() %>% 
  unique() 

OtherParamPlotSp <- 
  SamplesSpSingle %>% 
  select(all_of(Ind)) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate("Type" = factor(rep(
    c(rep("Jump Intensity",2),rep("Shock Params",2),rep("Error RW",1),
      rep("Drift",1),"Error Term"),
    nrow(SamplesUSSingle)),
    levels=c("Jump Intensity","Shock Params","Error RW","Drift","Error Term"))) %>%
  mutate(Param = factor(Param,
                        levels = unique(Param),
                        labels = c(expression(mu[Y]),expression(sigma[Y]), 
                                   expression(" a "),
                                   expression(p), expression(sigma[xi]),
                                   expression(" d "),expression(sigma[r])))) %>% 
  # # add y min a y max by group 
  mutate("ymin"=rep(c(rep(0.25,2),rep(0,2),rep(0.15,1), 
                      rep(-0.4,1),rep(0.025,1)),nrow(SamplesUSSingle)),
         "ymax"=rep(c(rep(2.25,2),rep(0.5,2),rep(0.45,1), 
                      rep(-0.1,1),rep(0.04,1)),nrow(SamplesUSSingle))) %>%
  ggplot(aes(y = Val, x = Param, group = Type))+
  ggdist::stat_slabinterval(
    aes(y = Val,fill =Param, colour = after_stat(level), size = NULL),
    geom = "interval",
    show_point = FALSE, .width = c(0.5, 0.8, 0.95), show_slab = FALSE,
    show.legend = NA
  )+
  stat_summary( #Add mean to plot
    geom = "point",fun = "mean",col = "#834333", size = 3
  )+
  facet_wrap(~Type, scales = "free",shrink = TRUE, nrow = 1)+
  scale_x_discrete(labels = function(l) parse(text=l))+ #transforms text into expression
  scale_color_manual(values = c("#fdf2c5","#ecdc94","#c9b045"))+
  ggdist::scale_slab_alpha_discrete(range = c(0.1, 1))+
  ylab("Posterior Mean")+ xlab("Other Parameters")+
  geom_blank(aes(y = ymin))+geom_blank(aes(y = ymax))+
  theme(legend.position = "none",
        axis.text = element_text(size = 25), #change size of axis text
        axis.title = element_text(size = 22),
        axis.title.y = element_text(vjust=1.3),
        axis.text.x = element_text(vjust=1),
        strip.text.x = element_blank(), #remove headers of facet wrap
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())



## 1.5.2.3 Add plots together using Cowplot 
cowplot::plot_grid(BetaParamSp, OtherParamPlotSp, byrow = FALSE, nrow = 2,
                   rel_heights = c(1.4,1)) #making first plot 1.4 times as high



#####  1.5.3 Italy  #########################################################
## 1.5.3.1 Beta Parameters 

AgeVecCovid <-  c(paste0(paste0(seq(0,80,10),sep="-"),seq(9,89,10)),"90+")

BetaParamIt <- SamplesItSingle %>% 
  select(grep("beta",colnames(.))) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate("AgeGroup"=rep(AgeVecCovid,2*nrow(SamplesItSingle)),#Create Age Group Variable
         "Type"=rep(c(rep("Norm",10),rep("Jump",10)),nrow(SamplesSpSingle))
  ) %>% mutate("Type"=factor(Type, 
                             levels = c("Norm","Jump"),
                             labels = c(expression(beta[x]),expression(beta[x]^{(J)}))),
               "AgeGroup"=factor(AgeGroup, levels = AgeVecCovid)) %>% #Change to factor 
  ggplot(aes(x = AgeGroup, y = Val,group = Type))+
  ggdist::stat_slabinterval(
    aes(group = after_stat(level), fill = after_stat(level), size = NULL),
    geom = "lineribbon",
    .width = c(0.5, 0.8, .95), show_slab = FALSE,
    show.legend = NA, color ="#182b1e"
  )+
  facet_wrap(~Type,nrow = 2,labeller = label_parsed)+
  scale_fill_manual(values = c("#d9e5dd","#a7c2af","#81a88d")) +
  #scale_interval_color_discrete(values = "#EBF90E")+ #EBF90E ,#400040
  ylab("Posterior Mean")+xlab("Age Group")+
  theme(legend.position = "none",
        axis.text = element_text(size = 20), #change size of axis text
        axis.title = element_text(size = 22),
        axis.title.y = element_text(vjust=1.3),
        axis.text.x = element_text(vjust=2),
        strip.text.x=element_text(face="bold", size = 22),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

## 1.5.3.2 Other Parameters 
OtherParams <- 
  c("muY","sdY","a","p","sigma_time","drift","sigma_eps")
#Get Indicies of correct Parameters
Ind <- sapply(OtherParams, 
              function(y){grep(paste0("\\<",y), #match beginning of word
                               x = colnames(SamplesItSingle))}) %>% unlist() %>% 
  unique() 

OtherParamPlotIt <- 
  SamplesItSingle %>% 
  select(all_of(Ind)) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate("Type" = factor(rep(
    c(rep("Jump Intensity",2),rep("Shock Params",2),rep("Error RW",1),
      rep("Drift",1),"Error Term"),
    nrow(SamplesUSSingle)),
    levels=c("Jump Intensity","Shock Params","Error RW","Drift","Error Term"))) %>%
  mutate(Param = factor(Param,
                        levels = unique(Param),
                        labels = c(expression(mu[Y]),expression(sigma[Y]), 
                                   expression(" a "),
                                   expression(p), expression(sigma[xi]),
                                   expression(" d "),expression(sigma[r])))) %>% 
  mutate("ymin"=rep(c(rep(0,2),rep(0,2),rep(0.15,1), 
                      rep(-0.2,1),rep(0.0275,1)),nrow(SamplesUSSingle)),
         "ymax"=rep(c(rep(0.75,2),rep(0.35,2),rep(0.35,1), 
                      rep(-0.1,1),rep(0.035,1)),nrow(SamplesUSSingle))) %>%
  ggplot(aes(y = Val, x = Param, group = Type))+
  ggdist::stat_slabinterval(
    aes(y = Val,fill =Param, colour = after_stat(level), size = NULL),
    geom = "interval",
    show_point = FALSE, .width = c(0.5, 0.8, 0.95), show_slab = FALSE,
    show.legend = NA
  )+
  stat_summary( #Add mean to plot
    geom = "point",fun = "mean",col = "#182b1e", size = 3
  )+
  facet_wrap(~Type, scales = "free",shrink = TRUE, nrow = 1)+
  scale_x_discrete(labels = function(l) parse(text=l))+ #transforms text into expression
  scale_color_manual(values = c("#d9e5dd","#a7c2af","#81a88d"))+
  ggdist::scale_slab_alpha_discrete(range = c(0.1, 1))+
  ylab("Posterior Mean")+ xlab("Other Parameters")+
  geom_blank(aes(y = ymin))+geom_blank(aes(y = ymax))+
  theme(legend.position = "none",
        axis.text = element_text(size = 25), #change size of axis text
        axis.title = element_text(size = 22),
        axis.title.y = element_text(vjust=1.3),
        axis.text.x = element_text(vjust=1),
        strip.text.x = element_blank(), #remove headers of facet wrap
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

## 1.5.3.3 Other Parameters 
cowplot::plot_grid(BetaParamIt, OtherParamPlotIt, byrow = FALSE, nrow = 2,
                   rel_heights = c(1.4,1)) #making first plot 1.4 times as high

##### 2. Visualization UK WAR DATA ###############################################
load(file = file.path(getwd(),"Results/Samples_UKWar.RData"))

#Total Samples
SamplesLiuLi_Tot <- do.call(rbind, 
                            SamplesLiuLi_OC_War) 

######### 2.1 Comparison Liu-Li Freq vs. Liu-Li Bayesian estimates ############
#Parameters to visualizse
Params <- c("beta","betaJump","muY","sdY","p",
            "drift","sigma_time","sigma_eps")

#Get Indicies of correct Parameters
Ind <- sapply(Params, 
              function(y){grep(paste0("\\<",y), #match beginning of word
                               x = colnames(SamplesLiuLi_Tot))}) %>% unlist() %>% 
  unique() 

#Create Age Group Vector
AgeGroupVec <- c("<1","1-4",paste0(paste0(seq(5,75,10),sep="-"),seq(14,84,10)))
#Create data set with mean estimates of Liu,Li(2015) (frequentist estimates)
FreqEstLiu <- data.frame("Param"=colnames(SamplesLiuLi_Tot)[Ind],
                        "MeanEst"=c(0.1125,0.2566,0.1522,0.0665,0.072, #beta 1-5
                                    0.0864,0.0783,0.0609,0.0579,0.0563, #beta 6-10
                                    0,0.0722,0.0806,0.3896,0.3081, #betaJ 1-5
                                    0.1108,0.0285,0.0106,0.0064,0, #betaJ 6-10
                                    1.9984,1.3205,0.1338,-0.2038,0.459,0.0437),
                        "Type"=c(rep("Norm",10),rep("Jump",10), #Needed for Plots
                                 rep("JumpParam",2),rep("Other",3),
                                 "Error Term"),
                        "AgeGroup"=c(rep(AgeGroupVec,2),        #Needed for Plots
                                     rep("Other",6))) %>% 
  mutate(Type = factor(Type,   #Change type Variable to Factor
                       levels = c("Norm","Jump","JumpParam","Other","Error Term"),
                       labels = c(expression(beta[x]),
                                  expression(beta[x]^{(J)}),
                                  "Jump Parameters","Other", "Error Term")),
         "AgeGroup"=factor(AgeGroup, levels = c(AgeGroupVec,"Other")),
         Param = factor(Param,
                  levels = Param,
                  labels = c(paste0("beta",1:10),paste0("betaJump",1:10),
                             expression(mu[Y]),expression(sigma[Y]),
                           expression(p), expression(d), expression(sigma[xi]),
                           expression(sigma[r])))) 

#Create Long Data Frame for Plotting
SamplesLiuLi_Tot_Long <- SamplesLiuLi_Tot %>% 
  as.data.frame() %>% #transform into Datafram
  select(all_of(Ind)) %>% #Select all relevant Parameters
  pivot_longer(cols = 1:ncol(.), names_to = "Param", values_to = "Val")
               
  
### Beta Parametsr with geom Ribbon ####
BetaParamPlotRibbon <- SamplesLiuLi_Tot_Long %>% 
  filter(grepl("beta",Param)) %>% 
  mutate("AgeGroup"=rep(AgeGroupVec,2*nrow(SamplesLiuLi_Tot)),#Create Age Group Variable
         "Type"=rep(c(rep("Norm",10),rep("Jump",10)),nrow(SamplesLiuLi_Tot))
  ) %>%  #transform Type into Factor Variable to allow greek facet labels
  mutate("Type"=factor(Type, 
                       levels = c("Norm","Jump"),
                       labels = c(expression(beta[x]),expression(beta[x]^{(J)}))),
         "AgeGroup"=factor(AgeGroup, levels = AgeGroupVec)) %>% #Change to factor 
  ggplot(aes(x = AgeGroup, y = Val,group = Type))+
  ggdist::stat_slabinterval(
    aes(group = after_stat(level), fill = after_stat(level), size = NULL),
    geom = "lineribbon",
    .width = c(0.5, 0.8, .95), show_slab = FALSE,
    show.legend = NA,
  )+
  facet_wrap(~Type,nrow = 2,labeller = label_parsed)+
  scale_fill_manual(values = c("#c6dbef","#6baed6","#08519c")) +
  geom_point(data = FreqEstLiu[1:20,], #Add Point Estimates
             aes(y = MeanEst, x = AgeGroup,group = Type), col="red")+ 
  ylab("Posterior Mean")+xlab("Age Group")+
  theme(legend.position = "none",
        axis.text = element_text(size = 20), #change size of axis text
        axis.title = element_text(size = 22),
        axis.title.y = element_text(vjust=1.3),
        axis.text.x = element_text(vjust=2),
        strip.text.x=element_text(face="bold", size = 22),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

  
#Other Parameters
OtherParamPlot <- 
SamplesLiuLi_Tot_Long %>% 
  filter(!grepl("beta",Param)) %>% 
  mutate("Type" = factor(rep(c(rep("Jump Parameters",2),rep("Other",3),"Error Term"),
                             nrow(SamplesLiuLi_Tot)),
                         levels=c("Jump Parameters","Other","Error Term"))) %>%
  mutate(Param = factor(Param,
                        levels = unique(Param),
                        labels = c(expression(mu[Y]),expression(sigma[Y]),
                                   expression(p), expression(d), expression(sigma[xi]),
                                   expression(sigma[r])))) %>% 
  # add y min a y max by group 
  mutate("ymin"=rep(c(1,1,-0.25,-0.25,-0.25,0.04),nrow(SamplesLiuLi_Tot)),
         "ymax"=rep(c(3.5,3.5,0.55,0.55,0.55,0.05),nrow(SamplesLiuLi_Tot))) %>%  
  ggplot(aes(y = Val, x = Param, group = Type))+
  ggdist::stat_interval(
    aes(fill = Param, color_ramp = after_stat(level)),
    .width = c(.50, .80, .95),
  )+
  geom_point(data = FreqEstLiu[21:26,], aes(y = MeanEst, x = Param, group = Type), col="red")+
  facet_wrap(~Type, scales = "free",shrink = TRUE)+
  scale_color_manual(values = c("#4292c6","#2171b5","#08519c"))+ 
  ylab("Posterior Mean")+ xlab("Other Parameters")+
  geom_blank(aes(y = ymin))+geom_blank(aes(y = ymax))+
  scale_x_discrete(labels = function(l) parse(text=l))+ #transforms text into expression
  theme(legend.position = "none",
        axis.text = element_text(size = 25), #change size of axis text
        axis.title = element_text(size = 22),
        axis.title.y = element_text(vjust=1.3),
        axis.text.x = element_text(vjust=1),
        strip.text.x = element_blank(), #remove headers of facet wrap
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())




cowplot::plot_grid(BetaParamPlotRibbon, OtherParamPlot, byrow = FALSE, nrow = 2,
                   rel_heights = c(1.4,1)) #making first plot 1.4 times as high


######### 1.6.2. Jump Occurence WAR ############################################
#Circular Plot for Jump occurence ##########
#Code taken from R graph gallery
DataCircle <- SamplesLiuLi_Tot %>% 
  data.frame() %>% 
  select(contains("N_t")) %>%
  summarise(across(.cols = where(is.numeric),.fns = mean)) %>% 
  pivot_longer(., cols = 1:ncol(.), names_to="Time", values_to = "PostMean") %>%
  mutate(Year = as.factor(1901:2010))


#Set Number of empty labels
emptyLabels <- 5
to_add <- matrix(NA, emptyLabels, ncol(DataCircle))
colnames(to_add) <- colnames(DataCircle)
DataCircle <- rbind(DataCircle, to_add)
DataCircle$id <- 1:nrow(DataCircle)

number_of_bar <- nrow(DataCircle)
# substract 0.5 because the letter must have the angle of the center of the bars. 
angle <-  90 - 360 * (DataCircle$id-0.5) /number_of_bar     

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
DataCircle$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
DataCircle$angle<-ifelse(angle < -90, angle+180, angle)

## Add a small constant
Constant <- 0.5
DataCircle <- DataCircle %>% 
  mutate("PostMean"=PostMean+Constant) %>% 
  mutate("Jump"=ifelse(PostMean>0.5+Constant,"yes","no"))

#Remove every other year
DataCircle$Year[!DataCircle$Year %in% c(seq(1901,1913,2),1914:1919,seq(1921,1939,2),1940:1945,seq(1947,2010,2))] <- NA

  ggplot(DataCircle, aes(x=id, y=PostMean)) +  
  geom_bar(aes(fill ="UK"),stat="identity",position = "dodge", width = 0.5)+
  scale_fill_manual(name="", values="#08519c")+
  geom_rect(aes(
  xmin = 0, xmax = 111, ymin = 0, ymax = Constant-0.1),
  fill = "grey97", color = "grey97"
  )+
  ylim(-1,1.6) +
  geom_hline(aes(yintercept = Constant), color = "grey60") +
  geom_hline(aes(yintercept = Constant+1), color = "grey60") +
  ggtitle("a.) Jump Occurence War Data")+
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.5,0.5),
    legend.text = element_text(size = 18),
    plot.title = element_text(hjust = 0.5, size = 20,
                              vjust = -10), #adjust position of title
    plot.margin=unit(c(-1.5,-1.5,-1.5,-1.5), "cm") #remove Margins from plot
  ) +
  coord_polar(start = 0) + # This makes the coordinate polar instead of cartesian.
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=DataCircle, aes(x=id, y=PostMean+0.1, 
                                 label=Year, 
                                 hjust=hjust), 
            color="black", fontface="bold",alpha=0.7, size=6, 
            angle= DataCircle$angle, inherit.aes = FALSE)
dev.off()
