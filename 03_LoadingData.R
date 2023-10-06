#### Script to load data for empirical analysis #########
library(tidyverse);library(openxlsx);library(rstan)

source("02_Functions.R") #needs stan 

############## COVID DATA ######################################################
############# 1.1 UNITED STATES ################################################
USDeaths <- read.table(file = file.path(getwd(),"Data/Deaths_5x1_USA.txt"),
                       header = TRUE)
USExp <- read.table(file = file.path(getwd(),"Data/Exposures_5x1_USA.txt"),
                    header = TRUE)
YearMin <- 1980

#Put Data into correct Format
#Create Data Frame
#HMD Data is Age groups of 5, including 0. For Sake of comparison, create age groups of 10

#Helper to Create age Groups of 0-10,10-20,20-30,...,90+
AgeWidthType <- 1 #Age with type 1 for COVID DATA
Helper <- AgeLabFun_HMD(Type = AgeWidthType)

TotalDataUS <- data.frame("Y"=USDeaths$Total,
                          "Offset"=USExp$Total,
                          "Year"=USDeaths$Year,
                          "AgeGroup"=USDeaths$Age) %>% 
  filter(Year > YearMin) %>% #only 1980 onward
  mutate("AInd"=match(AgeGroup, unique(AgeGroup))) %>% #Get Age ID
  mutate("NewAgeInd"=Helper$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise(across(c(Y,Offset),sum)) %>% #sum over Age ID
  mutate("TInd"=match(Year, unique(Year))) %>% 
  mutate("Y"=round(Y,0), #why are there decimal points in data?
         "Offset"=round(Offset,0)) %>% 
  rename("Deaths"=Y) %>% 
  rename("Pop"=Offset) %>% 
  arrange(Year)


#Create Vector of Death Rates and Mort Improvement Rates
LambdaVecUs <- TotalDataUS %>% mutate("Rate"=Deaths/Pop) %>% 
  mutate("ZVal"=c(NA,diff(log(Rate))))%>% 
  arrange(Year)

LambdaMatUS <- matrix(LambdaVecUs$Rate,
                      nrow = length(unique(LambdaVecUs$NewAgeInd)), 
                      byrow = FALSE)

#Differenced Log Death Rates
ZMatUS <- apply(log(LambdaMatUS), 1, diff) %>% t()


######### 1.2. Italy ############################################################
## Problem: HMD DATA only available until 2019, though Eurostat has Data for 2022.
# Solution: Merge both DataSets

#Eurostat Death Data starts from 1985
#Eurostat Population Data starts from 1992

#HMD Data from 1981 - 1992 respectively 1981 - 1985 for Pop/Death
library(openxlsx)
Deaths_It_Eu <- read.xlsx(xlsxFile = 
                            file.path(getwd(),"Data/DeathsItaly_Eurostat.xlsx"),
                          startRow = 8, sheet = 3)

DeathsIt_Eu22_We <- read.xlsx(xlsxFile = 
                                file.path(getwd(),"Data/WeeklyDeaths2022_It_Fr_Sp.xlsx"),
                              startRow = 25, sheet = "ItalyTotal")

PopIt_EU <- read.xlsx(xlsxFile = 
                        file.path(getwd(),"Data/PopulationItaly_Eurostat.xlsx"),
                      startRow = 8, sheet = 3)

#Helper to Create age Groups of <5,5-14,15-24,...,85+
HelperAgeLab <- AgeLabFun_EU(AgeLabels = Deaths_It_Eu$`AGE.(Labels)`)


DeathsIt_EU21 <- Deaths_It_Eu %>% 
  filter(`AGE.(Labels)`!= "Total") %>% #remove total age group 
  pivot_longer(.,                      #transform into long format 
               cols = 3:ncol(Deaths_It_Eu), 
               names_to = "Year", 
               values_to = "Deaths") %>% 
  mutate("AInd"=match(`AGE.(Labels)`, unique(`AGE.(Labels)`))) %>% #Get Age ID
  mutate("NewAgeInd"=HelperAgeLab$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise(across(Deaths,sum)) %>% 
  filter(Year >1984) %>% # data is available from 1985 onward
  mutate(Year = as.numeric(Year)) # transform year into numeric


## Add Deaths of Year 2022 of provisional Yearly Data from EuroStat
DeathsIt_Eu22 <- DeathsIt_Eu22_We %>% 
  mutate("NewAgeInd"= AgeLabFun_EU22()) %>% 
  group_by(NewAgeInd) %>% 
  summarise("Deaths"=sum(Total)) %>% 
  mutate("Year" = 2022)


DeathsGroupedIt_EU <- dplyr::full_join(x=DeathsIt_EU21,
                                       y = DeathsIt_Eu22, 
                                       by=c("NewAgeInd","Year","Deaths"),
                                       na_matches="never")

## Get HMD Data of deaths and combine 
DeathsIt_HMD <- read.table(file = file.path(getwd(),"Data/Deaths_5x1_Italy.txt"),
                           header = TRUE)
YearMin <- 1980

AgeWidthType <- 1 #Age with type 1 for COVID DATA
HelperHMD <- AgeLabFun_HMD(Type = AgeWidthType)

DeathsGrouped_It_HMD <- DeathsIt_HMD %>% filter(Year > YearMin) %>% 
  select(Year,Age,Total) %>% 
  mutate("AInd"=match(Age, unique(Age))) %>% #Get Age ID
  mutate("NewAgeInd"=HelperHMD$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise(across(Total,sum)) %>% 
  filter(Year < 1985) %>% 
  mutate(Deaths = Total) %>% 
  select(NewAgeInd, Year, Deaths)

DeathsIt_Tot <- full_join(x= DeathsGrouped_It_HMD, 
                          y= DeathsGroupedIt_EU, 
                          by=c("NewAgeInd","Year","Deaths"),
                          na_matches="never")


## Get Exposure. Approximate by Population Estimate
PopIt_EU <- read.xlsx(xlsxFile = 
                        file.path(getwd(),"Data/PopulationItaly_Eurostat.xlsx"),
                      startRow = 8, sheet = 3)
PopIt_HMD <- read.table(file = file.path(getwd(),"Data/Exposures_5x1_Italy.txt"),
                        header = TRUE)

#Select Year from 1993 onward
PopGrouped_EU_It <- PopIt_EU %>% 
  select(c(1,2,as.character(1993:2022))) %>% #from 1993 onward there is data available for 90+ years
  filter(`AGE.(Labels)`!= "Total") %>% #remove total age group
  pivot_longer(.,                      #transform into long format 
               cols = 3:ncol(.), 
               names_to = "Year", 
               values_to = "Pop") %>% 
  mutate("AInd"=match(`AGE.(Labels)`, unique(`AGE.(Labels)`))) %>% #Get Age ID
  mutate("NewAgeInd"=HelperAgeLab$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise(across(Pop,sum)) %>% 
  mutate("Year"=as.numeric(Year)) # transform year into integer

PopGrouped_HMD_It <- 
  PopIt_HMD %>% filter(Year > YearMin & Year < 1993) %>% 
  select(Year,Age,Total) %>% 
  mutate("AInd"=match(Age, unique(Age))) %>% #Get Age ID
  mutate("NewAgeInd"=HelperHMD$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise(across(Total,sum)) %>% 
  mutate(Pop = Total) %>% 
  select(-3)

PopIt_Tot <- rbind(PopGrouped_HMD_It, 
                   PopGrouped_EU_It)


#### Total Data (add both datasets together)
TotDataIt <- dplyr::left_join(x = DeathsIt_Tot,
                              y = PopIt_Tot, 
                              by =c("NewAgeInd","Year"),
                              na_matches="never")

# Only take the oldest Age Groups 
LambdaVecIt <- TotDataIt %>% 
  mutate("Rate"=Deaths/Pop) %>% 
  mutate("ZVal"=c(NA,diff(log(Rate)))) %>% 
  arrange(Year)

LambdaMatIt <- matrix(LambdaVecIt$Rate,
                      nrow = length(unique(LambdaVecIt$NewAgeInd)),
                      byrow = FALSE)

#Differenced Log Death Rates (Creation of Z Matrix )
ZMatIt <- apply(log(LambdaMatIt), 1, diff) %>% t()



############ 1.3. Spain ########################################################
#For Spain all the needed Data is available at Eurostat. No HMD needed
Deaths_Sp_EU <- openxlsx::read.xlsx(xlsxFile = 
                                      file.path(getwd(),"Data/DeathsSpain_Eurostat.xlsx"),
                                    startRow = 9, sheet = 3, na.strings = ":")

DeathsSp_Eu22_We <- openxlsx::read.xlsx(xlsxFile = 
                                          file.path(getwd(),"Data/WeeklyDeaths2022_It_Fr_Sp.xlsx"),
                                        startRow = 25, sheet = "SpainTotal")

PopIt_EU <- openxlsx::read.xlsx(xlsxFile = 
                                  file.path(getwd(),"Data/PopulationSpain_Eurostat.xlsx"),
                                startRow = 9, sheet = 3, na.strings = ":")

HelperAgeLab <- AgeLabFun_EU(AgeLabels = Deaths_Sp_EU$`AGE.(Labels)`)


DeathsSP_EU21 <- Deaths_Sp_EU %>% 
  filter(`AGE.(Labels)`!= "Total") %>% #remove total age group 
  pivot_longer(.,                      #transform into long format 
               cols = 3:ncol(Deaths_Sp_EU), 
               names_to = "Year", 
               values_to = "Deaths") %>% 
  mutate("AInd"=match(`AGE.(Labels)`, unique(`AGE.(Labels)`))) %>% #Get Age ID
  mutate("NewAgeInd"=HelperAgeLab$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise(across(Deaths,sum)) %>% 
  filter(Year >1980) %>% 
  mutate(Year = as.numeric(Year)) # transform year into numeric


## Add Deaths of Year 2022 of provisional Yearly Data from EuroStat
DeathsSp_Eu22 <- DeathsSp_Eu22_We %>% 
  mutate("NewAgeInd"= AgeLabFun_EU22()) %>% 
  group_by(NewAgeInd) %>% 
  summarise("Deaths"=sum(Total)) %>% 
  mutate("Year" = 2022)


DeathsGroupedSp_EU <- dplyr::full_join(x=DeathsSP_EU21,
                                       y = DeathsSp_Eu22, 
                                       by=c("NewAgeInd","Year","Deaths"),
                                       na_matches="never")

## Get Exposure. Value is based on annual (January 1st) population estimates
PopSp_EU <- openxlsx::read.xlsx(xlsxFile = 
                                  file.path(getwd(),"Data/PopulationSpain_Eurostat.xlsx"),
                                startRow = 9, sheet = 3, na.strings = ":")

#Select Year from 1981 onward
PopGrouped_EU_Sp <- PopSp_EU %>% 
  select(c(1,2,as.character(1981:2022))) %>% 
  filter(`AGE.(Labels)`!= "Total") %>% #remove total age group
  pivot_longer(.,                      #transform into long format 
               cols = 3:ncol(.), 
               names_to = "Year", 
               values_to = "Pop") %>% 
  mutate("AInd"=match(`AGE.(Labels)`, unique(`AGE.(Labels)`))) %>% #Get Age ID
  mutate("NewAgeInd"=HelperAgeLab$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise(across(Pop,sum)) %>% 
  mutate("Year"=as.numeric(Year)) # transform year into integer


#### Total Data (add both datasets together)
TotDataSp <- dplyr::left_join(x = DeathsGroupedSp_EU,
                              y = PopGrouped_EU_Sp, 
                              by =c("NewAgeInd","Year"),
                              na_matches="never")

LambdaVecSp <- TotDataSp %>% 
  mutate("Rate"=Deaths/Pop) %>% 
  mutate("ZVal"=c(NA,diff(log(Rate)))) %>% 
  arrange(Year)


## Creation of Z Matrix 
LambdaMatSp <- matrix(LambdaVecSp$Rate,
                      nrow = length(unique(LambdaVecSp$NewAgeInd)),
                      byrow = FALSE)

#Differenced Log Death Rates
ZMatSp <- apply(log(LambdaMatSp), 1, diff) %>% t()


### 1.4 Save all data ##########################################################
save(LambdaMatUS, ZMatUS, LambdaVecUs, TotalDataUS,
     LambdaMatIt,  ZMatIt, LambdaVecIt, TotDataIt,
     LambdaMatSp, ZMatSp, LambdaVecSp, TotDataSp,
     file = file.path(getwd(),"Data/CovidData.RData"))



######################### UK WAR DATA ##########################################
#Liu, Li Data
EnglandDeaths <- read.table(file = 
                              file.path(getwd(),"Data/Deaths_5x1_EnglandWales.txt"),
                            header = TRUE)
EnglandExp <- read.table(file = 
                           file.path(getwd(),"Data/Exposures_5x1_EnglandWales.txt"),
                         header = TRUE)


YearMax <- 2011 
YearMin <-  1900

#Helper Data.frame for combination of Age groups as done by Liu,Li
# Ages are given in groups of 10 instead of 5 (see Liu, Li p.137 bottom)
#<1, 1-4, 5-15, 15-25,...
AgeWidthType <- 2 #Age with type 1 for COVID DATA
Helper <- AgeLabFun_HMD(Type = AgeWidthType)


TotalDataWar <- data.frame("Deaths"=EnglandDeaths$Total,
                        "Pop"=EnglandExp$Total,
                        "Year"=EnglandDeaths$Year,
                        "AgeGroup"=EnglandDeaths$Age) %>% 
  mutate("AInd"=match(AgeGroup, unique(AgeGroup))) %>% #Get Age ID
  mutate("NewAgeInd"=Helper$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise(across(c(Deaths,Pop),sum)) %>% #sum over Age ID
  filter(Year > YearMin & Year < YearMax) %>% #only 1901 onward
  filter(NewAgeInd %in% c(1:10)) %>% #only first 10 Age Groups (till 75-85)
  mutate("TInd"=match(Year, unique(Year))) %>% 
  mutate("Deaths"=round(Deaths,0), #why are there decimal points in data?
         "Pop"=round(Pop,0)) %>% 
  arrange(Year)

#Log Death Rates
LambdaVecWar <- TotalDataWar %>% 
  mutate("Rate"=Deaths/Pop) %>% 
  mutate("ZVal"=c(NA,diff(log(Rate))))%>% 
  arrange(Year)


LambdaMatWar <- matrix(TotalDataWar$Deaths/TotalDataWar$Pop,
                       nrow = length(unique(TotalDataWar$NewAgeInd)), 
                       byrow = FALSE)

#Differenced Log Death Rates
ZMatWar <- apply(log(LambdaMatWar), 1, diff) %>% t()


save(LambdaMatWar, ZMatWar, LambdaVecWar, TotalDataWar,
     file = file.path(getwd(),"Data/UKWARData.RData"))
