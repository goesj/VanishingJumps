### Visualization of Results ##############
library(nimble);library(tidyverse);library(ggdist)

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 


load(file = file.path(getwd(),"Results/SamplesUS_NoRep.RData"))
load(file = file.path(getwd(),"Results/SamplesIt_NoRep.RData"))
load(file = file.path(getwd(),"Results/SamplesSp_NoRep.RData"))

###### 1.1 Data Description ###################################################
load(file = file.path(getwd(),"Data/CovidData.RData")) # Load Data
#### 1.1.1 US Data ###########################################################
pdf(file=file.path(getwd(),
                   "Bilder/LogDeathRateUS.pdf"),
    width=16, height = 9)
LambdaVecUs %>% 
  ungroup() %>% 
  mutate("AgeGroups"=rep(
    c(paste0(
      paste0(seq(0,80,10),sep="-"),
      seq(10,90,10)),"90+"),
    42)) %>%
  ggplot(data=., aes(x=Year,group=AgeGroups))+
  geom_line(aes(y = log(Rate), col=AgeGroups))+
  ylab("Log Death Rate")+
  facet_wrap(~AgeGroups,scales = "free")+
  theme(panel.background = element_rect(fill='gray95', colour='white'),
        axis.text = element_text(size = 17), #change size of axis text
        axis.title = element_text(size = 20),
        legend.text=element_text(size=17),
        legend.title = element_text(size=20)
  ) #change size of axis title)+
dev.off()


####  1.1.2 Spain Data #######################################################
pdf(file=file.path(getwd(),"Bilder/LogDeathRatesSpain.pdf"),width=16, height = 9)
LambdaVecSp %>% 
  ungroup() %>% 
  mutate("AgeGroups"=rep(
    c(paste0(
      paste0(seq(0,80,10),sep="-"),
      seq(9,89,10)),"90+"),
    42)) %>%
  filter(NewAgeInd %in% c(5:10)) %>% 
  ggplot(data=., aes(x=Year,group=AgeGroups))+
  geom_line(aes(y = log(Rate), col=AgeGroups), linewidth = 1.2)+
  geom_vline(xintercept = c(1999,2003,2015), linetype ="dashed", col = "gray50")+
  ylab("Log Death Rate")+
  facet_wrap(~AgeGroups,scales = "free")+
  theme(#panel.background = element_rect(fill='gray95', colour='white'),
        legend.text=element_text(size=17),
        legend.title = element_text(size=20),
        axis.text = element_text(size = 18), #change size of axis text
        axis.title = element_text(size = 18),
        strip.text.x=element_text(face="bold", size = 22),
        legend.position = "none"
  ) #change size of axis title)+
dev.off()

LambdaVecSp %>% 
  ungroup() %>% 
  mutate("AgeGroups"=rep(
    c(paste0(
      paste0(seq(0,80,10),sep="-"),
      seq(9,89,10)),"90+"),
    42)) %>%
  filter(NewAgeInd %in% c(5:10)) %>% 
  ggplot(data=., aes(x=Year,group=AgeGroups))+
  geom_line(aes(y = ZVal, col=AgeGroups), linewidth = 1.2)+
  geom_vline(xintercept = c(1999,2003,2015), linetype ="dashed", col = "gray50")+
  ylab("Log Death Rate")


##### 1.1.3 Italian Data #####################################################
pdf(file=file.path(getwd(),"Bilder/LogDeathRatesItaly.pdf"),width=16, height = 9)
LambdaVecIt %>% 
  ungroup() %>% 
  mutate("AgeGroups"=rep(
    c(paste0(
      paste0(seq(0,80,10),sep="-"),
      seq(9,89,10)),"90+"),
    42)) %>%
  filter(NewAgeInd %in% c(5:8)) %>% 
  ggplot(data=., aes(x=Year,group=AgeGroups))+
  geom_line(aes(y = log(Rate), col=AgeGroups), linewidth = 1.2)+
  geom_vline(xintercept = c(1983,2003,2015,2020), linetype ="dashed", col = "gray50")+
  ylab("Log Death Rate")+
  facet_wrap(~AgeGroups,scales = "free")+
  theme(legend.text=element_text(size=17),
        legend.title = element_text(size=20),
        axis.text = element_text(size = 18), #change size of axis text
        axis.title = element_text(size = 18),
        strip.text.x=element_text(face="bold", size = 22),
        legend.position = "none"
  ) #change size of axis title)+
dev.off()
# Shiny App with all results (Ask Karim for a Website)
pdf(file=file.path(getwd(),"Bilder/MortImprovementsItaly.pdf"),width=16, height = 9)
LambdaVecIt %>% 
  ungroup() %>% 
  mutate("AgeGroups"=rep(
    c(paste0(
      paste0(seq(0,80,10),sep="-"),
      seq(9,89,10)),"90+"),
    42)) %>%
  filter(NewAgeInd %in% c(5:10)) %>% 
  ggplot(data=., aes(x=Year,group=AgeGroups))+
  geom_line(aes(y = ZVal, col=AgeGroups), linewidth = 1.2)+
  geom_vline(xintercept = c(1983,2003,2015,2020), 
             linetype ="dashed", col = "gray50")+
  geom_hline(yintercept = 0, col = "black")+
  ylab("Mortality Improvements")+
  theme(legend.text=element_text(size=17),
        legend.title = element_text(size=20),
        axis.text = element_text(size = 18), #change size of axis text
        axis.title = element_text(size = 18),
        strip.text.x=element_text(face="bold", size = 22)
  )
dev.off()


########### SAMPLING RESULTS####################################################
library(ggdist)

#1.1.1 Add all Samples into one big Data frame 
SamplesItSingle <- do.call(rbind, 
                           SamplesOwn_It) %>% 
  data.frame() %>% 
  mutate(Country="Italy")

SamplesUSSingle <- do.call(rbind, 
                           SamplesOwn_US_NoRepara) %>% 
  data.frame() %>% 
  mutate(Country="United States")

# SamplesEWSingle <- do.call(rbind, 
#                            SamplesOwn_EW_Rep) %>% 
#   data.frame() %>% 
#   mutate(Country = "UK")


SamplesSpSingle <- do.call(rbind, 
                           SamplesOwn_Sp) %>% 
  data.frame() %>% 
  mutate(Country = "Spain")


SamplesTotal <- bind_rows(SamplesItSingle,
                          SamplesUSSingle,
                          SamplesSpSingle)

######## 1.2. Parameter Density Plot ###########################################
theme_set(theme_minimal(base_size = 10))

x11()
pdf(file=file.path(getwd(),"Bilder/VanishingParamComp_NoRep.pdf"), width=16, height = 9)
SamplesTotal %>% 
  select(a,p,muY,sdY,Country) %>% 
  #Attention: Factor level is turned around, so that last factor (Italy) is on top
  mutate("Country" = factor(Country, levels = c("United States","Spain","Italy"))) %>% 
  ggplot(aes(x = a, y = Country))+
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
  scale_fill_manual(values =c("#98579B","#c9b045","#81A88D"))+ #Colors for Countries
  ggdist::scale_fill_ramp_discrete(range = c(0.2, 0.95))+
  scale_color_manual(values=c("#400040","#806E26","#355C41"))+
  ggdist::scale_slab_color_discrete(c("#400040","#806E26","#355C41"),
                                    aesthetics = "slab_color",
                                    na.translate = FALSE)+
  theme(#panel.background = element_rect(fill='grey98', colour='white'), #change color?
        axis.text = element_text(size = 24), #change size of axis text
        axis.title = element_text(size = 25, face = "bold"),
        legend.position = "none",
  )+
  ylab("")
# xlim(c(0,2.5))+
# xlab(expression(sigma[Y]))
dev.off()


#Combining multiple legends: https://stackoverflow.com/questions/12410908/combine-legends-for-color-and-shape-into-a-single-legend
# scale_fill_manual(name = "Country PI",
#                   labels = as.vector(outer(c("United States", "Spain","Italy"),
#                                            c("0.99","0.85","0.5","NA"),
#                                            paste, sep=" ")),
#                   values =rep(c("#98579B","#c9b045","#81A88D"),4))+




#Coefficient of Variation
SamplesTotal %>% 
  select(a,p,muY,sdY,Country) %>% 
  #Attention: Factor level is turned around, so that last factor (Italy) is on top
  mutate("Country" = factor(Country, levels = c("United States","Spain","Italy"))) %>% 
  group_by(Country) %>% 
  summarise("CV"=sd(muY)/mean(muY))
  



########## 1.3. Jump Visualisation #############################################
## Time Series Plot##
SamplesTotal %>% 
  select(contains("N_t"),Country) %>%
  group_by(Country) %>% 
  select(-Country) %>% 
  summarise(across(.cols = where(is.numeric),.fns = mean)) %>% 
  pivot_longer(., cols = 2:ncol(.), names_to="Time", values_to = "PostMean") %>% 
  mutate(TimeFac = factor(Time, levels = unique(Time), ordered = TRUE)) %>% 
  ggplot(aes(x=TimeFac, group = Country))+
  geom_line(aes(y  = PostMean, col = Country))+
  facet_wrap(~Country, ncol = 1)+
  scale_color_manual(values=c("#81A88D","#c9b045","#98579B"))+
  scale_x_discrete(labels = 1981:2022)



#with Geom segment and point
pdf(file=file.path(getwd(),"Bilder/PosteriorMeanJump_Countries_noRep.pdf"), width=16, height = 9)
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
dev.off()



## Same plot but in a Circle

#Set Number of empty labels
emptyLabels <- 3
Constant <- 0.5

DataCircleTot <- 
  SamplesTotal %>% 
  select(contains("N_t"),Country)  %>% 
  group_by(Country) %>% 
  summarise(across(.cols = where(is.numeric),.fns = mean)) %>% 
  pivot_longer(., cols = 2:ncol(.), names_to="Time", values_to = "PostMean") %>%
  group_by(Country) %>% 
  group_modify(~ add_row(.data=., Time=rep(NA,emptyLabels))) %>% #Add Empty Rows
  mutate("Year"=1980+(1:n())) %>% #get Year for Label
  mutate("id"=1:n()) %>%  #Create an ID Variable
  mutate("PostMean"=PostMean+Constant) 


#Create Helper file for geom_text
HelpDataText <- 
  DataCircleTot %>%
  filter(Year < 2023) %>% 
  group_by(Year) %>% 
  reframe("MP"=max(PostMean,na.rm = TRUE)) %>% #get maximum posterior mean for placement of Year 
  add_row(.data=., Year=rep(NA,emptyLabels)) %>% #Add 5 Empty Rows
  mutate("id"=1:n()) %>% 
  mutate("angle"= 90 - 360 *(id-0.5)/n()) %>% # I substract 0.5 because the letter must have the angle of the center of the bars.
  mutate("hjust"=ifelse(angle< -90,1,0)) %>% # If I am on the left part of the plot, my labels have currently an angle < -90
  mutate("angle"=ifelse(angle< -90,angle+180,angle))


#With geom_bar instead of geom_segment for grouped bar plot
pdf(file=file.path(getwd(),"Bilder/PosteriorMeanJump_Countries_Circle.pdf"), width=10, height = 10,
    paper = "special")
BarEurope <- ggplot(DataCircleTot, aes(x=Year,
                          y=PostMean, group = Country, 
                          fill= Country)) +  
  geom_bar(stat ="identity", position = "dodge", width = 0.7) +
  geom_rect(aes(
    xmin = 1980, xmax = 2023, ymin = 0, ymax = Constant-0.1),
    fill = "grey97", color = "grey97"
  )+
  ylim(-1,1.6) +
  geom_hline(aes(yintercept = Constant), color = "grey60") +
  #geom_hline(aes(yintercept = Constant+0.5), color = "grey60") +
  geom_hline(aes(yintercept = Constant+1), color = "grey60") +
  ggtitle("b.) Jump Occurence COVID Data")+
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.5, 0.5), # move legend to the center
    legend.text = element_text(size = 18),
    #legend.title = element_text(size = 20),
    legend.title =element_blank(),
    plot.title = element_text(hjust = 0.5, size = 20,
                              vjust = -10), #adjust position of title
    plot.margin=unit(c(-1.5,-1.5,-1.5,-1.5), "cm") #remove Margins from plot
  ) +
  coord_polar(start = 0) + # This makes the coordinate polar instead of cartesian.
  scale_fill_manual(values=c("#81A88D","#c9b045","#98579B"))+
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=HelpDataText, 
                        aes(x=Year, y=MP+0.1, 
                                 label=Year, 
                                 hjust=hjust), 
            color="black", fontface="bold",alpha=0.7, size=6, 
            angle= HelpDataText$angle, 
            inherit.aes = FALSE)
dev.off()




############ 1.4. Age Pattern of Mortality Jump ###############
# #Option 1
# SamplesTotal %>% 
#   select_if(grepl(paste0(c("Country","betaJump"), collapse="|"), names(.))) %>% 
#   pivot_longer(., cols = 1:10, names_to = "BetaVal",values_to = "Val") %>%  #transform into long 
#   mutate(BetaVal = fct_inorder(as.factor(BetaVal),ordered = NA)) %>% 
#   group_by(Country,BetaVal) %>%  #group by
#   reframe(mean = mean(Val),
#           PIL = quantile(Val, probs=0.25),
#           PIU = quantile(Val, probs=0.75)) %>% #calcualte mean
#   ggplot(aes(x = BetaVal, group = Country))+
#   geom_line(aes(y =mean, col=Country), linewidth = 1.2)+
#   geom_ribbon(aes(ymin = PIL, ymax = PIU, fill = Country), alpha = 0.3)+
#   ylab("Posterior Value")+
#   xlab("Beta Jump")+
#   scale_fill_manual(values=c("#81A88D","#CCBA72","#98579B"))+
#   scale_color_manual(values=c("#355C41","#806E26","#400040"))+
#   scale_x_discrete(labels =as.expression(paste0("beta[",1:10,"]")))
# 
# ### Option 2
# SamplesTotal %>% 
#   select_if(grepl(paste0(c("Country","betaJump"), collapse="|"), names(.))) %>% 
#   pivot_longer(., cols = 1:10, names_to = "BetaVal",values_to = "Val") %>%  #transform into long 
#   mutate(BetaVal = fct_inorder(as.factor(BetaVal),ordered = NA)) %>% 
#   group_by(Country,BetaVal) %>%  #group by
#   reframe(mean = mean(Val),
#           PIL = quantile(Val, probs=0.25),
#           PIU = quantile(Val, probs=0.75)) %>% #calcualte mean
#   ggplot(aes(x = BetaVal, group = Country))+
#   geom_line(aes(y =mean, col=Country), linewidth = 1.2)+
#   geom_ribbon(aes(ymin = PIL, ymax = PIU, fill = Country), alpha = 0.3)+
#   facet_wrap(~  Country, nrow = 3)+
#   scale_fill_manual(values=c("#81A88D","#CCBA72","#98579B"))+
#   scale_color_manual(values=c("#355C41","#806E26","#400040"))+
#   ylab("Posterior Value")+
#   xlab("Mortality Effect")+
#   scale_x_discrete(labels = c(expression(beta[1]^{(J)},beta[2]^{(J)},beta[3]^{(J)},
#                                          beta[4]^{(J)} ,beta[5]^{(J)},beta[6]^{(J)},
#                                          beta[7]^{(J)},beta[8]^{(J)},beta[9]^{(J)},
#                                          beta[10]^{(J)})))+
#   scale_y_continuous(position = 'left', sec.axis = dup_axis())+ #duplicated y axis
#   theme(#panel.background = element_rect(fill='gray98', colour='white'),
#         axis.text = element_text(size = 22), #change size of axis text
#         axis.title = element_text(size = 22),
#         strip.text.x=element_text(face="bold", size = 22), #Text of Countries
#         legend.position = "none",
#         axis.title.x = element_blank(),
#         axis.title.y.left  = element_blank(), #remoce y axis title on left
#         axis.text.y.right = element_blank(),  #remove y axis text on right
#         axis.title.y.right = element_text(vjust = -5), #move left axis title to the left
#         #axis.text.y.left = element_text(hjust = 10)
#   )

#Option 3 (With other estimtates in Behind)
pacman::p_load("gghighlight")
pacman::p_load("ggtext")

#https://r-graph-gallery.com/web-line-chart-small-multiple-all-group-greyed-out.html
font <- "serif"


pdf(file=file.path(getwd(),"Bilder/AgeMortalityPattern.pdf"), width=16, height = 9)
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
  #scale_color_manual(values=c("#355C41","#806E26","#400040"))+
  scale_color_manual(values=c("#182b1e","#834333","#2b012d"))+
  ylab("Posterior Value")+
  xlab("Age Group")+
  scale_x_discrete(labels = c(paste0(paste0(seq(0,80,10),sep="-"),seq(9,89,10)),"90+"))+  
  #scale_y_continuous(position = 'left', sec.axis = dup_axis())+ #duplicated y axis
  theme(#panel.background = element_rect(fill='gray98', colour='white'),
        axis.text = element_text(size = 20), #change size of axis text
        axis.title = element_text(size = 22),
        strip.text.x=element_text(face="bold", size = 22), #Text of Countries
        legend.position = "none",
        #axis.title.y.left  = element_blank(), #remoce y axis title on left
        #axis.text.y.right = element_blank(),  #remove y axis text on right
        #axis.title.y.right = element_text(vjust = -5), #move left axis title to the left
        axis.text.x.bottom = element_text(vjust=2)
  )
dev.off()

### Same plot idea but with beta x behind in gray, or in light colors
# !!!!! DIFFICULT FOR INTERPRETATION !!!!! 
## However, difficult to interpret. Since kappa t is usually negative and jump effect positve
# Hence, a high value of beta is good, while a high value of betaJump bad (in terms of death rates)
BetaJumpFrame <- SamplesTotal %>% 
  select_if(grepl(paste0(c("Country","betaJump"), collapse="|"), names(.))) %>% 
  pivot_longer(., cols = 1:10, names_to = "Beta",values_to = "Val") %>%  #transform into long 
  mutate(Beta = fct_inorder(as.factor(Beta),ordered = NA)) %>% 
  group_by(Country,Beta) %>%  #group by
  reframe(mean = mean(Val),
          PIL = quantile(Val, probs=0.25),
          PIU = quantile(Val, probs=0.75))%>% 
  mutate("AgeID" = rep(1:10,3),
         "Type"="Jump") 

BetaFrame <- SamplesTotal %>% 
  select_if(grepl(paste0(c("Country","beta"), collapse="|"), names(.))) %>% 
  select(c(1:10,21)) %>% #only select beta, not beta Jump aswell
  pivot_longer(., cols = 1:10, names_to = "Beta",values_to = "Val") %>%  #transform into long 
  mutate(Beta = fct_inorder(as.factor(Beta),ordered = NA)) %>% 
  group_by(Country,Beta) %>%  #group by
  reframe(mean = mean(Val),
          PIL = quantile(Val, probs=0.25),
          PIU = quantile(Val, probs=0.75)) %>% 
  mutate("AgeID" = rep(1:10,3),
         "Type"="Normal") 

TotFrame <- rbind(BetaFrame,
                  BetaJumpFrame)

ggplot(data=TotFrame, aes(group=Country, x = AgeID))+
  geom_line(data = . %>% filter(Type =="Normal"), 
            aes(y = mean),col = "gray85", alpha= 0.8 ,linewidth = 0.75)+
  geom_ribbon(data = . %>% filter(Type =="Normal"),
              aes(ymin = PIL, ymax = PIU),col = "gray85", alpha = 0.3)+
  geom_line(data = . %>% filter(Type =="Jump"), 
            aes(y = mean, col = Country), linewidth = 0.75)+
  geom_ribbon(data = . %>% filter(Type =="Jump"),
              aes(ymin = PIL, ymax = PIU, fill = Country), alpha = 0.3)+
  facet_wrap(~Country, nrow = 3)+
  scale_fill_manual(values=c("#81A88D","#CCBA72","#98579B"))+
  scale_color_manual(values=c("#355C41","#806E26","#400040"))+
  ylab("Posterior Value")+
  xlab("Mortality Effect")+
  scale_x_discrete(labels = c(expression(beta[1]^{(J)},beta[2]^{(J)},beta[3]^{(J)},
                                         beta[4]^{(J)} ,beta[5]^{(J)},beta[6]^{(J)},
                                         beta[7]^{(J)},beta[8]^{(J)},beta[9]^{(J)},
                                         beta[10]^{(J)})))


##### 1.5 All Parameters in one ################################################
##### 1.5.1 United States #######################################################
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
    show.legend = NA, color ="#2b012d"   #400040
  )+
  #theme()+
  facet_wrap(~Type,nrow = 2,labeller = label_parsed)+
  scale_fill_manual(values = c("#f7ebfd","#b789b9","#98579b")) +
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
    c(rep("Jump Intensity",2),rep("Shock Params",2),rep("Error RW",1),rep("Drift",1),"Error Term"),
    nrow(SamplesUSSingle)),
    levels=c("Jump Intensity","Shock Params","Error RW","Drift","Error Term"))) %>%
  mutate(Param = factor(Param,
                        levels = unique(Param),
         labels = c(expression(mu[Y]),expression(sigma[Y]), expression(" a "),
                    expression(p), expression(sigma[xi]),
                    expression(d),expression(sigma[r])))) %>% 
  # # add y min a y max by group
  #https://stackoverflow.com/questions/40137060/how-to-set-limits-on-rounded-facet-wrap-y-axis/42590252#42590252
  mutate("ymin"=rep(c(rep(0.25,2),rep(0,2),rep(0.1,1), rep(-0.2,1),rep(0.02,1)),nrow(SamplesUSSingle)),
        "ymax"=rep(c(rep(2.25,2),rep(0.5,2),rep(0.3,1), rep(0,1),rep(0.025,1)),nrow(SamplesUSSingle))) %>%
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
  scale_color_manual(values = c("#eaddeb","#b789b9","#98579b"))+#Color Palette Blues of R Color Brewer
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



pdf(file=file.path(getwd(),"Bilder/AllParamsOwn_US_NoRep.pdf"), width=16, height = 9)
cowplot::plot_grid(BetaParamUS, OtherParamPlotUS, byrow = FALSE, nrow = 2,
                   rel_heights = c(1.4,1)) #making first plot 1.4 times as high
dev.off()


#Looks awful, and is difficult to interpret. Thus not added to the plot
TimeParametersPlot <- 
  SamplesUSSingle %>% 
  select(grep("k",colnames(.))) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate(Type="Delta Kappa",
         Year=rep(1982:2021, nrow(SamplesUSSingle))) %>% 
  ggplot(aes(x = Year, y = Val))+
  stat_slabinterval(
    aes(group = after_stat(level), fill = after_stat(level), size = NULL),
    geom = "lineribbon",
    .width = c(0.5, 0.8, .95), show_slab = FALSE,
    show.legend = NA, color ="#2b012d"   #400040
  )


######### 1.5.2 Spain #########################################################
AgeVecCovid <-  c(paste0(paste0(seq(0,80,10),sep="-"),seq(9,89,10)),"90+")
BetaParamSpain <- SamplesSpSingle %>% 
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
    c(rep("Jump Intensity",2),rep("Shock Params",2),rep("Error RW",1),rep("Drift",1),"Error Term"),
    nrow(SamplesUSSingle)),
    levels=c("Jump Intensity","Shock Params","Error RW","Drift","Error Term"))) %>%
  mutate(Param = factor(Param,
                        levels = unique(Param),
                        labels = c(expression(mu[Y]),expression(sigma[Y]), expression(" a "),
                                   expression(p), expression(sigma[xi]),
                                   expression(" d "),expression(sigma[r])))) %>% 
  # # add y min a y max by group 
  # #https://stackoverflow.com/questions/40137060/how-to-set-limits-on-rounded-facet-wrap-y-axis/42590252#42590252
  mutate("ymin"=rep(c(rep(0.25,2),rep(0,2),rep(0.15,1), rep(-0.4,1),rep(0.025,1)),nrow(SamplesUSSingle)),
         "ymax"=rep(c(rep(2.25,2),rep(0.5,2),rep(0.45,1), rep(-0.1,1),rep(0.04,1)),nrow(SamplesUSSingle))) %>%
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
  scale_color_manual(values = c("#fdf2c5","#ecdc94","#c9b045"))+#Color Palette Blues of R Color Brewer
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



pdf(file=file.path(getwd(),"Bilder/AllParamsOwn_Sp_NoRep.pdf"), width=16, height = 9)
cowplot::plot_grid(BetaParamSpain, OtherParamPlotSp, byrow = FALSE, nrow = 2,
                   rel_heights = c(1.4,1)) #making first plot 1.4 times as high
dev.off()


#####  1.5.3 Italy  #########################################################
#Italy colors: Line #355C41", 50PI;#81a88d; 80 PI: "#a7c2af", 99 PI "#d9e5dd" (70% light)
AgeVecCovid <-  c(paste0(paste0(seq(0,80,10),sep="-"),seq(9,89,10)),"90+")
BetaParamItaly <- SamplesItSingle %>% 
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
    c(rep("Jump Intensity",2),rep("Shock Params",2),rep("Error RW",1),rep("Drift",1),"Error Term"),
    nrow(SamplesUSSingle)),
    levels=c("Jump Intensity","Shock Params","Error RW","Drift","Error Term"))) %>%
  mutate(Param = factor(Param,
                        levels = unique(Param),
                        labels = c(expression(mu[Y]),expression(sigma[Y]), expression(" a "),
                                   expression(p), expression(sigma[xi]),
                                   expression(" d "),expression(sigma[r])))) %>% 
  mutate("ymin"=rep(c(rep(0,2),rep(0,2),rep(0.15,1), rep(-0.2,1),rep(0.0275,1)),nrow(SamplesUSSingle)),
         "ymax"=rep(c(rep(0.75,2),rep(0.35,2),rep(0.35,1), rep(-0.1,1),rep(0.035,1)),nrow(SamplesUSSingle))) %>%
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
  scale_color_manual(values = c("#d9e5dd","#a7c2af","#81a88d"))+#Color Palette Blues of R Color Brewer
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

pdf(file=file.path(getwd(),"Bilder/AllParamsOwn_It_NoRep.pdf"), width=16, height = 9)
cowplot::plot_grid(BetaParamItaly, OtherParamPlotIt, byrow = FALSE, nrow = 2,
                   rel_heights = c(1.4,1)) #making first plot 1.4 times as high
dev.off()
##### Visualization UK WAR DATA ###############################################
load(file = file.path(getwd(),"Results/Samples_UKWar.RData"))

#Total Samples
SamplesLiuLi_Tot <- do.call(rbind, 
                            SamplesLiuLi_OC_War) 

######### 1.6 Comparison Liu-Li Freq vs. Liu-Li Bayesian estimates ############
#Parameters to visualizse
library(ggdist)
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
#pdf(file=file.path(getwd(),"Bilder/AgePatternUKWar.pdf"), width=16, height = 9)
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
#dev.off()
  
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
  #https://stackoverflow.com/questions/40137060/how-to-set-limits-on-rounded-facet-wrap-y-axis/42590252#42590252
  mutate("ymin"=rep(c(1,1,-0.25,-0.25,-0.25,0.04),nrow(SamplesLiuLi_Tot)),
         "ymax"=rep(c(3.5,3.5,0.55,0.55,0.55,0.05),nrow(SamplesLiuLi_Tot))) %>%  
  ggplot(aes(y = Val, x = Param, group = Type))+
  ggdist::stat_interval(
    aes(fill = Param, color_ramp = after_stat(level)),
    .width = c(.50, .80, .95),
  )+
  geom_point(data = FreqEstLiu[21:26,], aes(y = MeanEst, x = Param, group = Type), col="red")+
  facet_wrap(~Type, scales = "free",shrink = TRUE)+
  scale_color_manual(values = c("#4292c6","#2171b5","#08519c"))+ #Color Palette Blues of R Color Brewer
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


#Option 2 with Stat eye (remove sigma_eps sice Variance is too small)
#Beta Parameter
BetaParamPlotInt <- SamplesLiuLi_Tot_Long %>% 
  filter(grepl("beta",Param)) %>% 
  mutate("AgeGroup"=rep(AgeGroupVec,2*nrow(SamplesLiuLi_Tot)),#Create Age Group Variable
         "Type"=rep(c(rep("Norm",10),rep("Jump",10)),4000)
  ) %>%  #transform Type into Factor Variable to allow greek facet labels
  mutate("Type"=factor(Type, 
                       levels = c("Norm","Jump"),
                       labels = c(expression(beta[x]),expression(beta[x]^{(J)}))),
         "AgeGroup"=factor(AgeGroup, levels = AgeGroupVec)) %>% 
  ggplot(aes(x = AgeGroup, y = Val,group = Type))+
  ggdist::stat_interval(
    aes(fill = Param, color_ramp = after_stat(level)),
    .width = c(.50, .80, .95),
  )+
  geom_point(data = FreqEstLiu[1:20,], aes(y = MeanEst, x = AgeGroup,group = Type), col="red")+
  scale_color_manual(values = c("#4292c6","#2171b5","#08519c"))+
  facet_wrap(~Type,nrow = 2,labeller = label_parsed)+
  ylab("Posterior Mean")+xlab("Age Group")+
  theme(legend.position = "none",
        axis.text = element_text(size = 20), #change size of axis text
        axis.title = element_text(size = 22),
        axis.title.y = element_text(vjust=1.3),
        axis.text.x = element_text(vjust=2),
        strip.text.x=element_text(face="bold", size = 22),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

OtherParamPlot2 <- SamplesLiuLi_Tot_Long %>% 
  filter(!grepl("beta",Param)) %>% 
  mutate("Type" = factor(rep(c(rep("Jump Parameters",2),rep("Other",4)),nrow(SamplesLiuLi_Tot)))) %>% 
  filter(Param!="sigma_eps") %>% 
  ggplot(aes(y = Val, x = Param,group = Type))+
  stat_eye(aes(fill=after_stat(level)),side="left",.width = c(0.5,0.8,0.95), 
           position = "dodge")+
  geom_point(data = FreqEstLiu[21:25,], aes(y = MeanEst, x = Param, group = Type), col="red")+
  facet_wrap(~Type, scales = "free")+
  scale_fill_brewer()+
  ylab("Posterior Mean")+ xlab("Parameters")+
  theme(legend.position = "none")



pdf(file=file.path(getwd(),"Bilder/AllParamLiuLi_UKWar.pdf"), width=16, height = 9)
cowplot::plot_grid(BetaParamPlotRibbon, OtherParamPlot, byrow = FALSE, nrow = 2,
                   rel_heights = c(1.4,1)) #making first plot 1.4 times as high
dev.off()


pdf(file=file.path(getwd(),"Bilder/AllParamLiuLi_Option2_UKWar.pdf"), width=16, height = 9)
cowplot::plot_grid(BetaParamPlotInt, OtherParamPlot2, byrow = FALSE, nrow = 2,
                   rel_heights = c(1.4,1)) #making first plot 1.4 times as high
dev.off()

######### 1.6.2. Jump Occurence WAR ############################################
pdf(file=file.path(getwd(),"Bilder/JumpOccurenceUK.pdf"), width=16, height = 9)
SamplesLiuLi_Tot %>% 
  data.frame() %>% 
  select(contains("N_t")) %>%
  summarise(across(.cols = where(is.numeric),.fns = mean)) %>% 
  pivot_longer(., cols = 1:ncol(.), names_to="Time", values_to = "PostMean") %>%
  mutate(Year = 1901:2010) %>% 
  ggplot(aes(x=Year))+
  geom_segment(aes(x = Year, xend = Year, y = 0, yend = PostMean), linewidth = 1.8,
               color = "#08519c")+
  geom_point(aes(y = PostMean), size = 3, pch = 18,col = "#08519c")+
  scale_x_continuous(breaks =  c(1901,1914,1919,1940,1945,seq(1955,2005,10)))+
  ylab("Posterior Mean")+xlab("")+
  theme(legend.position = "none",
        axis.text = element_text(size =20), #change size of axis text
        axis.title = element_text(size = 22),
        strip.text.x=element_text(face="bold", size = 22),
        axis.text.x = element_text(angle = 45, hjust = 1),
        #panel.grid.major.x = element_blank(),
        #panel.grid.minor.x = element_blank()
  )
dev.off()


### Circular Plot ##########
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
# I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
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
# pdf(file=file.path(getwd(),"Bilder/JumpOccurenceUK_Circle.pdf"), width=10, height = 10,
#     paper = "special")
BarUK <- 
  ggplot(DataCircle, aes(x=id, y=PostMean)) +  
  # geom_point(aes(y = PostMean), size = 3, pch = 18,col = "#08519c")+
  # geom_segment(
  #   aes(x = id, xend = id, y = 0, yend = PostMean),
  #   size = 1.2,
  #   #color = ifelse(DataCircle$Jump=="yes","#08519c","gray80")
  #   color = "#08519c"
  # ) +
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

pdf(file=file.path(getwd(),"Bilder/JumpOccurenceAll.pdf"), width=16, height = 9)
cowplot::plot_grid(BarUK, BarEurope, byrow = FALSE, nrow = 1)
dev.off()

####### 2. Comparison of Reparametization Type #################################
load(file = file.path(getwd(),"Results/SamplesUS_NoRepara.RData"))
load(file = file.path(getwd(),"Results/SamplesUS_Rep2.RData"))

load(file = file.path(getwd(),"Results/SamplesSp_Rep2.RData"))


SamplesUSSingle <- do.call(rbind, 
                           SamplesOwn_Sp) %>% 
  data.frame() %>% 
  mutate(Country="United States")

SamplesUSSingleRep <- 
  do.call(rbind, SamplesOwn_Sp_Rep) %>% 
  data.frame()  


OtherParams <- 
  c("beta","muY","sdY","a","p","sigma_time","drift","sigma_eps")

#Get Indicies of correct Parameters
IndUS <- sapply(OtherParams, 
              function(y){grep(paste0("\\<",y), #match beginning of word
                               x = colnames(SamplesUSSingle))}) %>% unlist() %>% 
  unique() 

IndUSRep <-  sapply(OtherParams, 
                    function(y){grep(paste0("\\<",y), #match beginning of word
                                     x = colnames(SamplesUSSingleRep))}) %>% unlist() %>% 
  unique()

AgeVecCovid <-  c(paste0(paste0(seq(0,80,10),sep="-"),seq(9,89,10)),"90+")

SamplesUSRepBeta <- SamplesUSSingleRep %>% 
  select(all_of(IndUSRep)) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate(Model = "Rep") %>% 
  filter(grepl("beta",Param)) %>% 
  mutate("AgeGroup"=rep(AgeVecCovid,2*nrow(SamplesUSSingleRep)),#Create Age Group Variable
         "Type"=rep(c(rep("Norm",10),rep("Jump",10)),nrow(SamplesUSSingleRep)))%>% 
  mutate("Type"=factor(Type, levels = c("Norm","Jump"),
                       labels = c(expression(beta[x]),expression(beta[x]^{(J)}))),
         "AgeGroup"=factor(AgeGroup, levels = AgeVecCovid))

SamplesUSOther <- SamplesUSSingleRep %>% 
  select(all_of(IndUSRep)) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate(Model = "Rep") %>% 
  filter(!grepl("beta",Param))


SamplesUSSingle %>% 
  select(all_of(IndUS)) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate(Model = "NoRep") %>% 
  filter(grepl("beta",Param)) %>% 
  mutate("AgeGroup"=rep(AgeVecCovid,2*nrow(SamplesUSSingle)),#Create Age Group Variable
         "Type"=rep(c(rep("Norm",10),rep("Jump",10)),nrow(SamplesUSSingle)))%>% 
  mutate("Type"=factor(Type, levels = c("Norm","Jump"),
                       labels = c(expression(beta[x]),expression(beta[x]^{(J)}))),
         "AgeGroup"=factor(AgeGroup, levels = AgeVecCovid)) %>% #Change to factor 
  bind_rows(x = ., y = SamplesUSRepBeta) %>% 
  ggplot(aes(x = AgeGroup, y = Val,fill = Model))+
  geom_boxplot()+
  facet_wrap(~Type, nrow = 2,labeller = label_parsed)


SamplesUSSingle %>% 
  select(all_of(IndUS)) %>% 
  pivot_longer(cols=1:ncol(.), names_to="Param", values_to = "Val") %>% 
  mutate(Model = "NoRep") %>% 
  filter(!grepl("beta",Param)) %>% 
  bind_rows(x = ., y = SamplesUSOther) %>% 
  ggplot(aes(x = Param, y = Val,fill = Model))+
  geom_boxplot()






