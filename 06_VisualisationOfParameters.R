### Visualisation  ##############
pacman::p_load("nimble","tidyverse","ggdist","MCMCvis","ggpubr",
               "cowplot","gghighlight","reshape2")

source("01_NimbleModels.R")
source("02_Functions.R") 


load(file = file.path(getwd(),"Results/SamplesUS.RData"))
load(file = file.path(getwd(),"Results/SamplesPl.RData"))
load(file = file.path(getwd(),"Results/SamplesSp.RData"))


########### SAMPLING RESULTS####################################################
theme_set(theme_minimal(base_size = 10)) #set global theme
#### 1.1.1 Data Preparation ####################################################
SamplesPlSingle <- 
  MCMCvis::MCMCchains(SamplesOwn_Pl_MA1) %>% 
  data.frame() %>% 
  mutate(Country="Poland")

SamplesUSSingle <- 
  MCMCvis::MCMCchains(SamplesOwn_US_MA1) %>% 
  data.frame() %>% 
  mutate(Country="United States")

SamplesSpSingle <- 
  MCMCvis::MCMCchains(SamplesOwn_Sp_MA1) %>% 
  data.frame() %>% 
  mutate(Country = "Spain")


SamplesTotal <- bind_rows(SamplesUSSingle,
                          SamplesSpSingle,
                          SamplesPlSingle)

######## 1.2. Parameter Density Plot ###########################################
SampleDensPlot <- function(Samples, Parameter = p){
  if(!Parameter %in% c("p","a","muY","sdY","b")){
    print("please select one of p, a, b,muY or sdY")
  } else  {
    Plot <- Samples %>% 
      #select(a,p,b,muY,sdY,Country) %>% 
      #Attention: Factor level is turned around, so that last factor (US) is on top
      mutate("Country" = factor(Country, 
                                levels = c("Poland","Spain","United States"))) %>% 
      ggplot(aes(x = .data[[Parameter]],   #select parameter for x 
                 y = Country))+    #change here parameter to sample
      ggdist::stat_slabinterval(aes(fill_ramp=after_stat(level),
                                    fill=Country),
                                geom="slab",
                                show_point = FALSE,
                                slab_type = "pdf", 
                                .width = c(.50, .80, 0.95))+
      ggdist::stat_pointinterval(aes(color = Country), #point interval seperate
                                 slab_linewidth=0.6,
                                 point_size = 2, 
                                 linewidth=2)+
      scale_fill_manual(values =c("#81A88D","#c9b045","#98579B"),
                        guide ="none")+ #Colors for Countries
      scale_fill_ramp_discrete(na.translate = FALSE, 
                               range = c(0.2,1))+ #drops the unnecessary NA
      scale_color_manual(values=c("#355C41","#806E26","#400040"), #color of mean estimate
                         guide ="none")+ # no legend
      labs(fill_ramp = "Credible Interval")+
      guides(fill_ramp =guide_legend(override.aes = 
                                       list(fill = "#222222"), #darken the input of the fill aes in legend
                                     reverse = TRUE))+
      theme(axis.text = element_text(size = 24), #change size of axis text
            axis.title = element_text(size = 25, face = "bold"),
            axis.title.x = element_text(margin = margin(t = 10)),
            legend.position = "bottom",
            legend.key.height = unit(0.3,"cm"), #adjust size of legend keys
            legend.key.width = unit(0.8,"cm"),
            legend.title = element_text(size = 20),
            legend.text = element_text(size = 20,
                                       margin = margin(r = 10, unit = "pt",
                                                       l = 10)), #increase spacing of text
            #plot.title = element_text(size = 25, face ="bold", hjust = 0.5)
      )+
      #ggtitle(label = Parameter)+xlab("")
      ylab("")
    return(Plot)
  }
}


#Create plot select one of a, p, muY, and sdY for Parameter
SampleDensPlot(Samples = SamplesTotal, Parameter =  "b")


########## 1.3. Jump Visualisation #############################################
#Plot for Visualisation of Jump Occurence
SamplesTotal %>% 
  select(contains("N_t"),Country) %>%
  mutate(Country = factor(Country, levels = c("United States","Spain","Poland"),
                          labels = c("a) United States" ,"b) Spain","c) Poland"))) %>% 
  group_by(Country) %>% 
  summarise(across(.cols = where(is.numeric),.fns = mean)) %>% 
  pivot_longer(., cols = 2:ncol(.), names_to="Time", values_to = "PostMean") %>% 
  mutate(TimeFac = factor(Time, levels = unique(Time), ordered = TRUE)) %>% 
  ggplot(aes(x=TimeFac, group = Country, col = Country))+
  geom_segment(aes(x = TimeFac, xend = TimeFac, y = 0, yend = PostMean), linewidth = 1.8)+
  geom_point(aes(y = PostMean), size = 3, pch = 18)+
  facet_wrap(~Country, ncol = 1)+
  scale_color_manual(values=c("#98579B","#c9b045","#81A88D"))+
  scale_x_discrete(labels = 1991:2023)+
  ylab("Posterior Mean")+xlab("")+
  theme(legend.position = "none",
        axis.text = element_text(size =20), #change size of axis text
        axis.title = element_text(size = 22),
        strip.text.x=element_text(face="bold", size = 22),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major.y = element_line(colour = "gray87"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
  )


############ 1.4. Age Pattern of Mortality Jump ################################

#Plot to create the age pattern of Mortality Jump for all three Countries
#https://r-graph-gallery.com/web-line-chart-small-multiple-all-group-greyed-out.html

SamplesTotal %>% 
  select_if(grepl(paste0(c("Country","betaJump"), collapse="|"), names(.))) %>% 
  mutate(Country = factor(Country, levels = c("United States","Spain","Poland"),
                          labels = c("a) United States" ,"b) Spain","c) Poland"))) %>% 
  pivot_longer(., cols = 1:10, names_to = "BetaVal",values_to = "Val") %>%  #transform into long 
  mutate(BetaVal = fct_inorder(as.factor(BetaVal),ordered = NA)) %>% 
  group_by(Country,BetaVal) %>%  #group by
  reframe(mean = mean(Val),
          PIL = quantile(Val, probs=0.1),
          PIU = quantile(Val, probs=0.9)) %>% 
  ggplot(aes(x = BetaVal, group = Country)) +
  geom_line(aes(y =mean, col=Country), linewidth = 0.75)+
  geom_ribbon(aes(ymin = PIL, ymax = PIU, fill = Country), alpha = 0.3)+
  gghighlight::gghighlight(use_direct_label = FALSE,
                           unhighlighted_params = list(colour = alpha("grey85", 1)))+
  facet_wrap(~  Country, nrow = 3)+
  scale_fill_manual(values=c("#98579B","#c9b045","#81A88D"))+
  scale_color_manual(values=c("#2b012d","#834333","#182b1e"))+
  ylab("Posterior Value")+
  xlab("Age Group")+
  scale_x_discrete(labels = c("0-4",paste0(paste0(seq(5,75,10),sep="-"),seq(14,84,10)),"85+"))+  
  theme(#panel.background = element_rect(fill='gray98', colour='white'),
    axis.text = element_text(size = 20), #change size of axis text
    axis.title = element_text(size = 22),
    strip.text.x=element_text(face="bold", size = 22), #Text of Countries
    legend.position = "none",
    axis.text.x.bottom = element_text(vjust=2)
  )

########### 1.5. Jump Parameter Visualization ##################################
# 1.5.1 Starting with hyper paramters of Y_t
SingleParamPl <- MCMCvis::MCMCchains(SamplesOwn_Pl_MA1, params = c("muY","sdY")) %>% 
  as.data.frame() %>% 
  pivot_longer(., cols = 1:ncol(.), names_to = "Param", values_to = "Val") %>% 
  mutate(Country = "Poland")

SingleParamUS <- MCMCvis::MCMCchains(SamplesOwn_US_MA1, params = c("muY","sdY")) %>% 
  as.data.frame() %>% 
  pivot_longer(., cols = 1:ncol(.), names_to = "Param", values_to = "Val") %>% 
  mutate(Country = "US")

SingleParamSp <- MCMCvis::MCMCchains(SamplesOwn_Sp_MA1, params = c("muY","sdY")) %>% 
  as.data.frame() %>% 
  pivot_longer(., cols = 1:ncol(.), names_to = "Param", values_to = "Val") %>% 
  mutate(Country = "Spain")

SingleParamTotal <- rbind(SingleParamUS,
                          SingleParamSp,
                          SingleParamPl) %>% 
  mutate(Param = factor(Param, 
                        levels = c("muY","sdY"), 
                        labels = c(expression(mu[Y]), expression(sigma[Y])))) %>%   
  mutate("Country" = factor(Country, 
                            levels = c("US","Spain","Poland")))

BarPlotJumpParameters <- SingleParamTotal %>% 
  ggplot(aes(x = Val, y = Param, group = Country))+
  ggdist::stat_slabinterval(
    aes(x = Val, color = Country, color_ramp = after_stat(level),
        size = NULL),
    #side = "both",
    geom = "interval",
    show_point = FALSE, .width = c(0.5, 0.8, 0.95), show_slab = FALSE,
    show.legend = TRUE, position = position_dodge(width = 0.3, preserve = "total"), 
  )+
  stat_summary( #Add mean to plot
    color =c("#2b012d","#2b012d",
             "#834333","#834333",
             "#182b1e","#182b1e"),
    geom = "point",fun = "mean", size = 2.5,
    pos = position_dodge(width = 0.3, preserve = "total"),
  )+
  scale_color_manual(values = c("#98579b","#c9b045","#81a88d"))+
  labs(color_ramp = "Credible Interval")+ #change name of "level" title in legend
  guides(color_ramp = guide_legend(reverse = TRUE))+ #change order of apperance in legend
  xlab("")+ylab("")+
  ggtitle("b)")+
  scale_y_discrete(labels = function(l) parse(text=l))+ #transforms text into expression
  theme(axis.text = element_text(size = 25), #change size of axis text
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        legend.position = "bottom")+coord_flip()
# 

### 1.5.2 Next is Jump Intensity Y_t###
SingleParamPl_Jump <- MCMCvis::MCMCchains(SamplesOwn_Pl_MA1, params = c("Y_t")) %>% 
  as.data.frame() %>% select(30,31) %>% 
  pivot_longer(., cols = 1:ncol(.), names_to = "Param", values_to = "Val") %>% 
  mutate(Country = "Poland")

SingleParamUS_Jump <- MCMCvis::MCMCchains(SamplesOwn_US_MA1, params = c("Y_t")) %>% 
  as.data.frame() %>% select(30,31) %>% 
  pivot_longer(., cols = 1:ncol(.), names_to = "Param", values_to = "Val") %>% 
  mutate(Country = "US")

SingleParamSp_Jump <- MCMCvis::MCMCchains(SamplesOwn_Sp_MA1, params = c("Y_t")) %>% 
  as.data.frame() %>% select(30) %>% 
  pivot_longer(., cols = 1:ncol(.), names_to = "Param", values_to = "Val") %>% 
  mutate(Country = "Spain")

SingleParamTotal_Jump <- rbind(SingleParamUS_Jump,
                               SingleParamSp_Jump,
                               SingleParamPl_Jump) %>% 
  mutate(Param = factor(Param, 
                        levels = c("Y_t[30]","Y_t[31]") , 
                        labels = c(expression(Y[2020]), expression(Y[2021])))) %>% 
  mutate("Country" = factor(Country, 
                            levels = c("US","Spain","Poland")))

DensityPlotJumpIntensity <- 
  SingleParamTotal_Jump %>% 
  ggplot(aes(x = Val, y = Param, group = Country))+
  ggdist::stat_slabinterval(
    aes(x = Val,
        fill = Country, # color by Country
        fill_ramp = after_stat(level), #Fill by Country as well
        color = Country), # color of line underneath dist
    side = "right",.width = c(0.5, 0.8, 0.95),
    show_point = TRUE, show_slab = TRUE, point_interval = "mean_qi",
    position = position_dodge(width = 0.66, preserve = "total"), 
  )+
  scale_fill_manual(values = c("#98579b","#c9b045","#81a88d"))+
  scale_color_manual(values = c("#4e2b50","#655823","#3f4e44"),
                     guide ="none")+
  xlim(c(0,2))+
  xlab("Posterior Value")+ylab("")+
  ggtitle("a)")+
  guides(level ="none")+
  scale_y_discrete(labels = function(l) parse(text=l))+ #transforms text into expression
  theme(axis.text = element_text(size = 25), #change size of axis text
        axis.title = element_text(size = 22),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        plot.title = element_text(hjust = 0.5, size = 25, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")+coord_flip()

#Both Togehter
ComLegend <- ggpubr::get_legend(BarPlotJumpParameters)

# ggarrange for common legend
ggpubr::ggarrange(DensityPlotJumpIntensity,
                  BarPlotJumpParameters,
                  common.legend = TRUE, 
                  legend.grob = ComLegend,
                  legend = "bottom", align = "h")


############ 1.5 All Parameters in one #########################################

## Creation of plot that visualizes all Parameters into one plot

##### 1.5.1 United States ######################################################
#Starting with the Age Effects beta_x, beta_x^J

# 1.5.1.1 Beta Parameters 
AgeVecCovid <-  c("0-4",paste0(paste0(seq(5,75,10),sep="-"),seq(14,84,10)),"85+")
BetaParamUS <- SamplesUSSingle %>% 
  select(grep("beta",colnames(.))) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate("AgeGroup"=rep(AgeVecCovid,2*nrow(SamplesUSSingle)),#Create Age Group Variable
         "Type"=rep(c(rep("Norm",10),rep("Jump",10)),nrow(SamplesUSSingle))
  ) %>% mutate("Type"=factor(Type, #Plot Variables by Type
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
  ylab("Posterior Value")+xlab("Age Group")+
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
  mutate("Type" = factor(rep( #Plot Variables by Type (for facet_wrap)
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
  mutate("ymin"=rep(c(rep(0,2),rep(0,2),rep(0.1,1), rep(-0.2,1),rep(0.02,1)),
                    nrow(SamplesUSSingle)),
         "ymax"=rep(c(rep(3,2),rep(0.6,2),rep(0.3,1), rep(0,1),rep(0.025,1)),
                    nrow(SamplesUSSingle))) %>%
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
  ylab("Posterior Value")+ xlab("Other Parameters")+
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
AgeVecCovid <-  c("0-4",paste0(paste0(seq(5,75,10),sep="-"),seq(14,84,10)),"85+")
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
  ylab("Posterior Value")+xlab("Age Group")+
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
  mutate("ymin"=rep(c(rep(0,2),rep(0,2),rep(0.2,1), rep(-0.4,1),rep(0.025,1)),
                    nrow(SamplesUSSingle)),
         "ymax"=rep(c(rep(3,2),rep(0.5,2),rep(0.45,1), rep(-0.1,1),rep(0.04,1)),
                    nrow(SamplesUSSingle))) %>% 
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
  ylab("Posterior Value")+ xlab("Other Parameters")+
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


#### 1.5.3. Poland ############################################################
AgeVecCovid <-  c("0-4",paste0(paste0(seq(5,75,10),sep="-"),seq(14,84,10)),"85+")
BetaParamPoland <- SamplesPlSingle %>% 
  select(grep("beta",colnames(.))) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate("AgeGroup"=rep(AgeVecCovid,2*nrow(SamplesPlSingle)),#Create Age Group Variable
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
  ylab("Posterior Value")+xlab("Age Group")+
  theme(legend.position = "none",
        axis.text = element_text(size = 20), #change size of axis text
        axis.title = element_text(size = 22),
        axis.title.y = element_text(vjust=1.3),
        axis.text.x = element_text(vjust=2),
        strip.text.x=element_text(face="bold", size = 22),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


OtherParams <- 
  c("muY","sdY","a","p","sigma_time","drift","sigma_eps")
#Get Indicies of correct Parameters
Ind <- sapply(OtherParams, 
              function(y){grep(paste0("\\<",y), #match beginning of word
                               x = colnames(SamplesPlSingle))}) %>% unlist() %>% 
  unique() 

OtherParamPlotPl <- 
  SamplesPlSingle %>% 
  select(all_of(Ind)) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate("Type" = factor(rep(
    c(rep("Jump Intensity",2),rep("Shock Params",1),rep("Error RW",1),rep("Drift",1),"Error Term"),
    nrow(SamplesPlSingle)),
    levels=c("Jump Intensity","Shock Params","Error RW","Drift","Error Term"))) %>%
  ggplot(aes(y = Val, x = Param, group = Type))+
  ggdist::stat_slabinterval(
    aes(y = Val,fill = Param, colour = after_stat(level), size = NULL),
    geom = "interval",
    show_point = FALSE, .width = c(0.5, 0.8, 0.95), show_slab = FALSE,
    show.legend = NA
  )+
  stat_summary( #Add mean to plot
    geom = "point",fun = "mean",col = "#182b1e", size = 3
  )+
  facet_wrap(~Type, scales = "free",shrink = TRUE, nrow = 1)+
  scale_x_discrete(labels = function(l) parse(text=l))+ #transforms text into expression
  scale_color_manual(values = c("#d9e5dd","#a7c2af","#81a88d"))+#Color Palette Blues of R Color Brewer
  ggdist::scale_slab_alpha_discrete(range = c(0.1, 1))+
  ylab("Posterior Value")+ xlab("Other Parameters")+
  geom_blank(aes(y = ymin))+geom_blank(aes(y = ymax))+
  theme(legend.position = "none",
        axis.text = element_text(size = 25), #change size of axis text
        axis.title = element_text(size = 22),
        axis.title.y = element_text(vjust=1.3),
        axis.text.x = element_text(vjust=1),
        strip.text.x = element_blank(), #remove headers of facet wrap
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

### 2. Visualization UK WAR Forecasts  #########################################
load(file = file.path(getwd(),"Data/UKWARData.RData")) # Load Data
load(file.path(getwd(),"Big_results/Samples_UKWar_80.RData"))

#Total Samples
pacman::p_load("MCMCvis","scoringRules","truncnorm")
S <- 20000
H <- 30
LastYearObs <- 1980

#Create Forecasts 
FutureZWAR_AR <- FutureZ(H=H, Mod = "AR", Samples = SamplesAR_War_80, S=S)
FutureZWar_LC <- FutureZ(H=H, Mod = "Liu", LC = TRUE,
                         Samples = SamplesLC_War_80, S=S)

LastYearObsInd <- which(1901:2010 == LastYearObs)

FutureLogMatArrayWar_AR <- FutureLogMatArrayWar_Liu <-  FutureLogMatArrayWar_MA <- 
  FutureLogMatArrayWar_LC <- 
  array(data = 0,dim = c(10,H,S),
        dimnames = list("NewAgeInd"=1:10, "FCHorizon"=1:H, "It"=1:S))

for(s in 1:S){
  FutureLogMatArrayWar_AR[,,s] <- log(LambdaMatWar)[,LastYearObsInd]+ #Last obs rate
    t(apply(FutureZWAR_AR$Rates[,,s],1,cumsum)) #Cumsum
  FutureLogMatArrayWar_LC[,,s] <- log(LambdaMatWar)[,LastYearObsInd]+ #Last obs rate
    t(apply(FutureZWar_LC$Rates[,,s],1,cumsum)) #Cumsum
}

# Create Helper Data Files for Plotting 
FutureLambdaVecIt_UK_LC <- 
  reshape2::melt(FutureLogMatArrayWar_LC, value.name = "logRate") %>% 
  mutate(Rate = exp(logRate),
         Type="FC",
         Year = FCHorizon+LastYearObs) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c("0-4",paste0(paste0(seq(5,75,10),sep="-"),
                                                   seq(14,84,10)),"85+")))


#variance of Liu,Li model higher resulting in wider PI's
FutureLambdaVecIt_UK_AR <- 
  reshape2::melt(FutureLogMatArrayWar_AR, value.name = "logRate") %>% 
  mutate(Rate = exp(logRate),
         Type="FC",
         Year = FCHorizon+LastYearObs) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c("0-4",paste0(paste0(seq(5,75,10),sep="-"),
                                                   seq(14,84,10)),"85+")))

### Get observed Values Vector
ObservedValuesVec_UK <- 
  LambdaVecWar %>%  #Take observed Vector and calculate PI's
  filter(Year <= (LastYearObs+ H)) %>% 
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels =c("0-4",paste0(paste0(seq(5,75,10),sep="-"),
                                                  seq(14,84,10)),"85+"))) %>% 
  group_by(AgeGroup,Year) %>% 
  ggdist::point_interval(Rate, .width = c(0.95)) %>% 
  mutate(Type ="Obs")

### Calculate Quantiles of each forecast and add observed values to vector #####
FutureValuesVecQuantile_UK_AR <- 
  FutureLambdaVecIt_UK_AR %>% 
  group_by(AgeGroup,Year) %>% 
  mutate("logRate"=log(Rate)) %>% 
  ggdist::point_interval(Rate, .width = c(0.95)) %>% 
  #ggdist::point_interval(Rate, .width = c(0.8,0.95)) %>% #Quantile Calc
  full_join(x = filter(ObservedValuesVec_UK, Year == LastYearObs),
            y =. ) %>% #add last year, so that there is no space in plot 
  mutate(Type ="FC",
         Model ="AR") %>% 
  full_join(x = filter(ObservedValuesVec_UK),
            y =. )

FutureValuesVecQuantile_UK_LC <-
  FutureLambdaVecIt_UK_LC %>%
  group_by(AgeGroup,Year) %>%
  ggdist::point_interval(Rate, .width = c(0.95)) %>%
  mutate(Type ="FC",
         Model="LC") %>% 
  full_join(x = filter(ObservedValuesVec_UK,
                       .width == 0.95,
                       Year == LastYearObs),
            y =. ) %>% 
  mutate(.width= 95,##Use different name for Width Column, so it gets alinged a different color
         Type = "Mean LC") 


##### FINAL PLOT ###############################################################
AgeFilter <- c("45-54","55-64","65-74","75-84")
YearFilter <- c(1960:2010) 
FutureValuesVecQuantile_UK_LC %>% 
  filter(Year %in% YearFilter,
         AgeGroup %in% AgeFilter) %>% 
  ggplot(aes(x =Year, group = AgeGroup))+
  ggdist::geom_lineribbon(aes(y = Rate, ymin = .lower,
                              ymax = .upper,
                              fill = interaction(.width),
                              color = Type,
                              linetype = Type),
                          linewidth = 1.1, alpha = 0.65)+
  geom_vline(xintercept = LastYearObs,linetype = "dashed", col="gray30")+
  ggdist::geom_lineribbon(data = filter(FutureValuesVecQuantile_UK_AR,
                                        Year %in% YearFilter,
                                        AgeGroup %in% AgeFilter),
                          aes(y = Rate, ymin = .lower, 
                              ymax = .upper, 
                              fill = interaction(.width),
                              color = Type,
                              linetype = Type),linewidth = 1.1, alpha = 0.65)+
  geom_line(aes(y = Rate, x = Year, col = Type, linetype = "Mean LC"),
            linewidth = 1.1,
            alpha = 1)+
  facet_wrap(~AgeGroup, scales = "free")+
  scale_color_manual(values =c("#67000D","#011f4b","black"),
                     labels =c("Mean AR","Mean LC","Observed"),
                     name = "Type")+
  scale_fill_manual(values = c("#00496F","#ED8B00"),
                    labels = c("LC 95%","AR 95%"),
                    name = "PI's")+
  scale_linetype_manual(
    values = c(1,1,6),
    labels = c("Mean AR",
               "Mean LC",
               "Observed"),
    name ="Type")+
  ylab("Death Rate")+xlab("Year")+
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



# 3. Visualization of Percentage Increase ######################################

#First load forecasts of own Model   
load(file = file.path(getwd(),"Results/SamplesSp.RData")) # own Model
load(file = file.path(getwd(),"Results/SamplesUS.RData")) # own Model
load(file = file.path(getwd(),"Results/SamplesPl.RData")) # own Model

#Load observed Data
load(file = file.path(getwd(),"Data/CovidData.RData")) # Load Data
 
  
### Calculate Future Mort Improvement Rates incl. Parameters
S <- 10000
H <- 50
set.seed(42)

FutureZUS <-  FutureZ(H=H, Mod = "AR", Samples = SamplesOwn_US_AR, S=S)

FutureZSp <- FutureZ(H=H, Mod = "AR", Samples = SamplesOwn_Sp_AR, S=S)

FutureZPl <- FutureZ(H=H, Mod = "AR", Samples = SamplesOwn_Pl_AR, S=S)


FutureZList <- list( FutureZUS,
                     FutureZSp,
                     FutureZPl)


JumpEffectArray <- array(0, dim=c(10, H+1, S,3),
                         dimnames = list("Age"=1:10,
                                         "FC"=1:(H+1),
                                         "It"=1:S,
                                         "Country"=c("a) US", "b) Spain","c) Poland")))
for (c in 1:3) { #over countries
  for (s in 1:S) { #over Iterations
    for(x in 1:10){ #over Age
      JumpEffectArray[x,,s,c] <- FutureZList[[c]]$betaJump[s,x]*FutureZList[[c]]$Jt[s,]  
    }
  }
}

#Transform big Array into long format
JumpEffectQuant <- 
  reshape2::melt(JumpEffectArray, value.name = "logEffect") %>% 
  mutate(AgeGroup = factor(Age,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c("0-4",paste0(paste0(seq(5,75,10),sep="-"),
                                                   seq(14,84,10)),"85+"))) %>% 
  group_by(FC,AgeGroup,Country) %>% 
  ggdist::point_interval(logEffect, .width = c(0.9,0.95,0.99))

#Jump Effect File by iteration
JumpEffectIt <- 
  reshape2::melt(JumpEffectArray, value.name = "logEffect") %>% 
  mutate(AgeGroup = factor(Age,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c("0-4",paste0(paste0(seq(5,75,10),sep="-"),
                                                   seq(14,84,10)),"85+")))


#Calculate the actual, observed Covid Increase
CovidIncrease <- bind_rows(mutate(LambdaVecUs, Country="a) US"),
                           mutate(LambdaVecSp, Country="b) Spain"),
                           mutate(LambdaVecPl, Country="c) Poland")) %>% 
  group_by(NewAgeInd,Country) %>% 
  filter(Year %in% c(2016:2019,2020)) %>%
  mutate(Period=ifelse(Year<2020,"Ref","Cov")) %>% #define reference period and Covid
  group_by(Period, Country, NewAgeInd) %>% #group by new Period
  reframe(Rate2 = mean(Rate)) %>%  #Caclulate the mean
  arrange(desc(Period)) %>%  #Change order, for easier division in long format
  group_by(NewAgeInd, Country) %>%
  reframe("PercentageIncrease"=exp(diff(log(Rate2)))) %>%  #calculation of mortality rate increase between mean 2016:2019 and 202
  mutate(AgeGroup = factor(NewAgeInd,   #Change type Variable to Factor
                           levels = 1:10,
                           labels = c("0-4",paste0(paste0(seq(5,75,10),sep="-"),
                                                   seq(14,84,10)),"85+")),
         CovidIncrease=" ")#Empty column so it can appear in legend

#Plot the Results (however in this plot legend is not correct)
AreaPlot <- JumpEffectQuant %>% 
  group_by(AgeGroup,.width,Country) %>% 
  summarise("Perc"=mean(.upper)) %>%
  mutate("PercIncrease"=(exp(Perc)-1)*100) %>% 
  mutate(".width"=factor(.width, levels = c(0.99,0.95,0.9))) %>% #change order for geom area
  ggplot(aes(x=AgeGroup))+
  geom_area(aes(y = PercIncrease, group = .width, 
                fill = interaction(factor(.width),Country),
                col = interaction(factor(.width),Country)),
            position = "identity")+
  geom_line(data = CovidIncrease, aes(x =AgeGroup, 
                                      y=(PercentageIncrease-1)*100,
                                      group=Country, 
                                      linewidth = "Covid Increase",
                                      linetype = "Covid Increase"))+
  scale_linetype_manual(values = 2, #set linetype of Covid Increase 
                        name ="Type")+
  scale_linewidth_manual(values = 1.2, # set line width of Covid Increase
                         guide = "none")+ #remove from legend
  scale_fill_manual(values = c("#EBDCEB","#C299C2","#98579b",
                               "#F7EFD9","#E3CF90","#c9b045",
                               "#E5EDE7","#B2CAB9","#81a88d"))+
  scale_color_manual(values = c("#EBDCEB","#C299C2","#98579b",
                                "#F7EFD9","#E3CF90","#c9b045",
                                "#E5EDE7","#B2CAB9","#81a88d"),
                     guide ="none")+ #remove from legend
  facet_wrap(~Country, nrow = 3, scales="free")+
  labs(fill = "Credible Interval")+
  ylab("Percentage Increase")+
  xlab("Age Group")+
  labs(col="Prediction Interval")+
  theme(
    axis.text = element_text(size = 20), #change size of axis text
    axis.title = element_text(size = 22),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 20),
    legend.position = "none",
    legend.key.width = unit(1.7,"cm"),
    strip.text.x=element_text(face="bold", size = 22)
  )+
  guides(linetype = guide_legend(override.aes = list(linewidth = 1.3)))

## Create an interval helper plot to extract the legend
IntervalPlot <- JumpEffectIt %>% 
  group_by(AgeGroup,Country) %>% 
  mutate("PercIncrease"=(exp(logEffect)-1)*100) %>% 
  ggplot(aes(x = AgeGroup, y = PercIncrease, group = Country))+
  ggdist::stat_slabinterval(
    aes(color = Country, color_ramp = after_stat(level),
        size = NULL),
    #side = "both",
    geom = "interval",
    show_point = FALSE, .width = c(0.9, 0.95, 0.99), show_slab = FALSE
    #position = position_dodge(width = 0.3, preserve = "total"), 
  )+
  geom_point(data = CovidIncrease, aes(x = AgeGroup, 
                                       y = (PercentageIncrease-1)*100,
                                       group = Country))+
  geom_line(data = CovidIncrease, aes(x =AgeGroup,
                                      y=(PercentageIncrease-1)*100,
                                      group=Country,
                                      linewidth = "Covid Increase",
                                      linetype = "Covid Increase"))+
  scale_linetype_manual(values = 2, #set linetype of Covid Increase
                        name ="Type")+
  scale_linewidth_manual(values = 1.2, # set line width of Covid Increase
                         guide = "none")+ #remove from legend
  scale_color_manual(values = c("#98579b","#c9b045","#81a88d"),
                     guide = "none")+
  facet_wrap(~Country, nrow = 3, scales="free")+
  labs(fill = "Credible Interval")+
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
  labs(color_ramp = "Credible Interval")+
  guides(linetype = guide_legend(override.aes = list(linewidth = 1.3))) 

#extract the legend
Legend <- ggpubr::get_legend(IntervalPlot)

## Add both together
ggpubr::ggarrange(AreaPlot,
                  common.legend = TRUE, 
                  legend.grob = Legend,
                  legend = "bottom")


