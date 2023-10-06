library(nimble);library(tidyverse); library(nimbleHMC)

source("01_NimbleModels.R")
source("02_Functions.R") #loading stan 
load(file = file.path(getwd(),"Results/SamplesUS_NoRep.RData"))
load(file = file.path(getwd(),"Results/SamplesIt_NoRep.RData"))
load(file = file.path(getwd(),"Results/SamplesSp_NoRep.RData"))

###### US DATA #################################################################
SumOutUSOwn <- SummaryOutput(SamplesOwn_US_NoRepara,
                             params = c("beta","betaJump",
                                        "drift","sigma_time","sigma_eps",
                                        "p","a",
                                        "muY","sdY")) %>% 
  as.data.frame()

#https://tex.stackexchange.com/questions/358462/including-greek-letters-when-saving-matrix-from-r-to-tex-table
ColNamesLatex <- c("Parameter","Mean","Sd","10\\%-PI","90\\%-PI",
                   paste0("split-","$\\hat{R}$"),"Tail\\_ESS","Bulk\\_ESS")

RowNamesLatex <-  c(paste0("$\\beta_{",1:10,"}$"),
                    paste0("$\\beta^{(J)}_{",1:10,"}$"),
                    "$d$","$\\sigma_\\xi$","$\\sigma_{r}$",
                    "$p$","$a$","$\\mu_{Y}$",
                    "$\\sigma_{Y}$")


colnames(SumOutUSOwn) <- ColNamesLatex

#Change Names of Parameter for Latex 
rownames(SumOutUSOwn) <- RowNamesLatex

hline <- c(-1,0,nrow(SumOutUSOwn))
htype <- c("\\toprule ", "\\midrule ", "\\botrule ")

print(xtable::xtable(SumOutUSOwn[,-1], 
                     caption = "Posterior Estimates of the United States", 
                     digits = c(0,4,4,4,4,2,0,0),
                     align=c("c",rep("l",7))),
      add.to.row = list(pos = as.list(hline), #add top rule etc.
                        command = htype),
      hline.after = NULL,
      include.rownames=TRUE, 
      sanitize.colnames.function = identity, #needed so that backslahes are correctly shown
      sanitize.rownames.function = identity, #needed for backslashes
      caption.placement = "top") 

###### Spanish Data ############################################################
SumOutSpOwn <- SummaryOutput(SamplesOwn_Sp,
                             params = c("beta","betaJump",
                                        "drift","sigma_time","sigma_eps",
                                        "p","a",
                                        "muY","sdY")) %>% 
  as.data.frame()

#Change Names of Parameter for Latex 
colnames(SumOutSpOwn) <- ColNamesLatex
  
rownames(SumOutSpOwn) <- RowNamesLatex

print(xtable::xtable(SumOutSpOwn[,-1], 
                     caption = "Posterior Estimates of the Spain", 
                     digits = c(0,4,4,4,4,2,0,0),
                     align=c("c",rep("l",7))),
      add.to.row = list(pos = as.list(hline), #add top rule etc.
                        command = htype),
      hline.after = NULL,
      include.rownames=TRUE, 
      sanitize.colnames.function = identity, #needed so that backslahes are correctly shown
      sanitize.rownames.function = identity, #needed for backslashes
      caption.placement = "top") 

###### Italian DATA ############################################################
SumOutItOwn <- SummaryOutput(SamplesOwn_It,
                             params = c("beta","betaJump",
                                        "drift","sigma_time","sigma_eps",
                                        "p","a",
                                        "muY","sdY")) %>% 
  as.data.frame()

#Change Names of Parameter for Latex 
colnames(SumOutItOwn) <- ColNamesLatex
rownames(SumOutItOwn) <- RowNamesLatex

print(xtable::xtable(SumOutItOwn[,-1], 
                     caption = "Posterior Estimates of Italy", 
                     digits = c(0,4,4,4,4,2,0,0),
                     align=c("c",rep("l",7))),
      add.to.row = list(pos = as.list(hline), #add top rule etc.
                        command = htype),
      hline.after = NULL,
      include.rownames=TRUE, 
      sanitize.colnames.function = identity, #needed so that backslahes are correctly shown
      sanitize.rownames.function = identity, #needed for backslashes
      caption.placement = "top") 


##### UK WAR DATA ##############################################################
load(file = file.path(getwd(),"Results/Samples_UKWar.RData"))

SumOutOwnUK <- SummaryOutput(Samples_OwnMod_OC_War, 
                           params = c("beta","betaJump", 
                                      "drift","sigma_time","sigma_eps",
                                      "p","a","muY","sdY")) %>% 
  as.data.frame()

colnames(SumOutOwnUK) <- c("Parameter","Mean","Sd","10\\%-PI","90\\-PI",
                         paste0("split-","$\\hat{R}$"),"Tail\\_ESS","Bulk\\_ESS")


#Change Names of Parameter for Latex 
rownames(SumOutOwnUK) <- c(paste0("$\\beta_{",1:10,"}$"),
                         paste0("$\\beta^{(J)}_{",1:10,"}$"),
                         "$d$","$\\sigma_\\xi$","$\\sigma_{r}$",
                         "$p$","$a$","$\\mu_{Y}$","$\\sigma_{Y}$"
                         )

hlineUK <- c(-1,0,nrow(SumOutOwnUK))
print(xtable::xtable(SumOutOwnUK[,-1], 
                     caption = "Posterior Estimates of England and Wales", 
                     digits = c(0,4,4,4,4,2,0,0),
                     align=c("c",rep("l",7))),
      add.to.row = list(pos = as.list(hlineUK), #add top rule etc.
                        command = htype),
      hline.after = NULL,
      include.rownames=TRUE, 
      sanitize.colnames.function = identity, #needed so that backslahes are correctly shown
      sanitize.rownames.function = identity, #needed for backslashes
      caption.placement = "top")


################ WAIC AND LOO-CV TABLE #########################################
library(xtable)
load(file.path(getwd(),"Results/WAIC_US_NoRep.RData"))
load(file.path(getwd(),"Results/WAIC_It_NoRep.RData"))
load(file.path(getwd(),"Results/WAIC_Sp_NoRep.RData"))
load(file.path(getwd(),"Results/WAIC_War.RData"))




print(xtable::xtable(cbind(CompDataFrameIt[,c(1,3)],
                   CompDataFrameSp[,c(3)],
                   CompDataFrameUS[,c(3)])),
      include.rownames = FALSE)

print(xtable::xtable(
  cbind(CompDataFrameIt[,c(1,2)],
      CompDataFrameSp[,2],
      CompDataFrameUS[,2])),
      include.rownames = FALSE)


###### WAIC for UK War Data ####################################################
load(file.path(getwd(),"Results/WAIC_War"))


print(xtable::xtable(CompDataFrameWar),
      include.rownames = FALSE)


