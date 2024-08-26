#### Script to load data for empirical analysis #########
pacman::p_load("tidyverse","openxlsx")

source("02_Functions.R") #needs stan 

############## COVID DATA ######################################################
############# 1.1 UNITED STATES ################################################
USDeaths <- read.table(file = file.path(getwd(),"Data/Deaths_5x1_USA.txt"),
                       header = TRUE)
USExp <- read.table(file = file.path(getwd(),"Data/Exposures_5x1_USA.txt"),
                    header = TRUE)
YearMin <- 1990

#Put Data into correct Format
#Create Data Frame
#HMD Data is Age groups of 5, including 0. For Sake of comparison, create age groups of 10

#Helper to Create age Groups of <5,5-14,15-24,...,90+
Helper <- AgeLabFun_HMD()

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


# Provisional Data United States #####
## New Data on Wonder http://wonder.cdc.gov/mcd-icd10-provisional.html
USProvisionalData <- read.table(file = file.path(getwd(), "Data/Provisional_Data_US_22_23_New.txt"),
                                header = TRUE,sep = "", dec = ".", fill = TRUE, na.strings = "")

USProvisionalData <- USProvisionalData %>% 
  select(YearCode, Ten.Year_AgeGroups, Deaths, Population) %>% 
  mutate("NewAgeInd"=rep(c(rep(1,2),2:10),2))

## Summarise by age 
US_22_23_Dat <- 
  USProvisionalData %>% 
  group_by(NewAgeInd,YearCode) %>% 
  reframe(Deaths=sum(Deaths),
          Pop = sum(Population)) %>% #sum over age groups 
  rename("Year"=YearCode) %>% 
  arrange(Year) %>% 
  mutate("TInd"=rep(c(41,42),each = 10))


TotalDataUS <- rbind(TotalDataUS,
                     US_22_23_Dat) #add 2022 Data

LambdaVecUs <- TotalDataUS %>% mutate("Rate"=Deaths/Pop) %>% 
  mutate("ZVal"=c(NA,diff(log(Rate))))%>% 
  arrange(Year)

LambdaMatUS <- matrix(LambdaVecUs$Rate,
                      nrow = length(unique(LambdaVecUs$NewAgeInd)), 
                      byrow = FALSE)

#Differenced Log Death Rates
ZMatUS <- apply(log(LambdaMatUS), 1, diff) %>% t()



############ 1.2. Spain ########################################################
#For Spain all the needed Data is available at Eurostat. No HMD needed
Deaths_Sp_EU <- openxlsx::read.xlsx(xlsxFile = 
                                      file.path(getwd(),"Data/DeathsSpain_Eurostat.xlsx"),
                                    startRow = 9, sheet = 3, na.strings = ":") %>% 
  slice(-nrow(.)) #remove "unknown" age group

DeathsSp_Eu23_We <- openxlsx::read.xlsx(xlsxFile = 
                                          file.path(getwd(),"Data/WeeklyDeaths2023_It_Sp.xlsx"),
                                        startRow = 47, sheet = "SpainTotal") %>% 
  slice(-nrow(.))

PopSp_EU <- openxlsx::read.xlsx(xlsxFile = 
                                  file.path(getwd(),"Data/PopulationSpain_Eurostat.xlsx"),
                                startRow = 9, sheet = 3, na.strings = ":")

HelperAgeLab <- AgeLabFun_EU(AgeLabels = Deaths_Sp_EU$`AGE.(Labels)`)

DeathsSP_EU22 <- Deaths_Sp_EU %>% 
  filter(`AGE.(Labels)`!= "Total") %>% #remove total age group 
  pivot_longer(.,                      #transform into long format 
               cols = 3:ncol(Deaths_Sp_EU), 
               names_to = "Year", 
               values_to = "Deaths") %>% 
  mutate("AInd"=match(`AGE.(Labels)`, unique(`AGE.(Labels)`))) %>% #Get Age ID
  mutate("NewAgeInd"=HelperAgeLab$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise(across(Deaths,sum)) %>% 
  filter(Year >YearMin) %>% 
  mutate(Year = as.numeric(Year)) # transform year into numeric


## Add Deaths of Year 2022 of provisional Yearly Data from EuroStat
DeathsSp_Eu23 <- DeathsSp_Eu23_We %>% 
  mutate("NewAgeInd"= AgeLabFun_EU23()) %>% 
  group_by(NewAgeInd) %>% 
  summarise("Deaths"=sum(Total)) %>% 
  mutate("Year" = 2023)


DeathsGroupedSp_EU <- dplyr::full_join(x=DeathsSP_EU22,
                                       y = DeathsSp_Eu23, 
                                       by=c("NewAgeInd","Year","Deaths"),
                                       na_matches="never")


#Select Year from 1981 onward
PopGrouped_EU_Sp <- PopSp_EU %>% 
  select(c(1,2,as.character(1991:2023))) %>% 
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

######### 1.3. Poland ##########################################################
#Provisional Population data from : 
#https://stat.gov.pl/en/topics/population/population/population-size-and-structure-and-vital-statistics-in-poland-by-territorial-division-as-of-31-december,3,33.html
Deaths_Pl_Eu <- openxlsx::read.xlsx(xlsxFile = 
                                      file.path(getwd(),"Data/DeathsPoland_Eurostat.xlsx"),
                                    startRow = 8, sheet = 3, na.strings = ":") %>% slice(-nrow(.))

DeathsPl_Eu23_We <- openxlsx::read.xlsx(xlsxFile = 
                                          file.path(getwd(),"Data/WeeklyDeaths_2023_Pl.xlsx"),
                                        startRow = 1,
                                        sheet = "PolandTotal") %>% 
  slice(-nrow(.))


#Helper to Create age Groups of <5,5-14,15-24,...,85+
HelperAgeLab <- AgeLabFun_EU(AgeLabels = Deaths_Pl_Eu$`AGE.(Labels)`)


DeathsPl_EU22 <- Deaths_Pl_Eu %>% 
  filter(`AGE.(Labels)`!= "Total") %>% #remove total age group 
  pivot_longer(.,                      #transform into long format 
               cols = 3:ncol(Deaths_Pl_Eu), 
               names_to = "Year", 
               values_to = "Deaths") %>% 
  mutate("AInd"=match(`AGE.(Labels)`, unique(`AGE.(Labels)`))) %>% #Get Age ID
  mutate("NewAgeInd"=HelperAgeLab$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise("Deaths"=sum(Deaths, na.rm =TRUE)) %>% 
  filter(Year >YearMin) %>% # data is available from 1998 onward
  mutate(Year = as.numeric(Year)) # transform year into numeric


## Add Deaths of Year 2022 of provisional Yearly Data from EuroStat
DeathsPl_Eu23 <- DeathsPl_Eu23_We %>% 
  mutate("NewAgeInd"= AgeLabFun_EU23()) %>% 
  group_by(NewAgeInd) %>% 
  summarise("Deaths"=sum(Total)) %>% 
  mutate("Year" = 2023)


DeathsPl_Tot <- dplyr::full_join(x=DeathsPl_EU22,
                                       y = DeathsPl_Eu23, 
                                       by=c("NewAgeInd","Year","Deaths"),
                                       na_matches="never")

## Get Exposure. Approximate by Population Estimate
PopPl_EU <- openxlsx::read.xlsx(xlsxFile = 
                                  file.path(getwd(),"Data/PopulationPoland_Eurostat.xlsx"),
                                startRow = 8, sheet = 3, na.strings = ":")

# PopPl_HMD <- read.table(file = file.path(getwd(),"Data/Exposures_5x1_Poland.txt"),
#                         header = TRUE)

#Select Year from 1993 onward
PopGrouped_EU_Pl <- PopPl_EU %>% 
  select(c(1,2,as.character(1991:2020))) %>% #from 1990 onward there is data available for 90+ years
  filter(`AGE.(Labels)`!= "Total") %>% #remove total age group
  pivot_longer(.,                      #transform into long format 
               cols = 3:ncol(.), 
               names_to = "Year", 
               values_to = "Pop") %>% 
  mutate("AInd"=match(`AGE.(Labels)`, unique(`AGE.(Labels)`))) %>% #Get Age ID
  mutate("NewAgeInd"=HelperAgeLab$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise("Pop"=sum(Pop, na.rm =TRUE)) %>% 
  mutate("Year"=as.numeric(Year)) # transform year into integer


#### Get Population for 2021, 2022, and 2023 from Stat.Pl 
PopPl_StatPl <- openxlsx::read.xlsx(xlsxFile = 
                                      file.path(getwd(),"Data/PopulationPoland_StatPl.xlsx"),
                                    startRow = 2, sheet = "Pop21-23", na.strings = ":")

PopGrouped_StatPl <- PopPl_StatPl %>% 
  filter("Age"!= "Total") %>% #remove total age group
  pivot_longer(.,                      #transform into long format 
               cols = 2:ncol(.), 
               names_to = "Year", 
               values_to = "Pop") %>% 
  mutate("AInd"=match(Age, unique(Age))) %>% #Get Age ID
  mutate("NewAgeInd"=HelperAgeLab$AgeNew[AInd]) %>% #add new Age ID
  group_by(NewAgeInd,Year) %>% #group by new Age ID
  summarise("Pop"=sum(Pop, na.rm =TRUE)) %>% 
  mutate("Year"=as.numeric(Year))

PopPl_Tot <- rbind(PopGrouped_EU_Pl,
                   PopGrouped_StatPl)

#### Total Data (add both datasets together)
TotDataPl <- dplyr::left_join(x = DeathsPl_Tot,
                              y = PopPl_Tot, 
                              by =c("NewAgeInd","Year"),
                              na_matches="never")

# Only take the oldest Age Groups 
LambdaVecPl <- TotDataPl %>% 
  mutate("Rate"=Deaths/Pop) %>% 
  mutate("ZVal"=c(NA,diff(log(Rate)))) %>% 
  arrange(Year)

LambdaMatPl <- matrix(LambdaVecPl$Rate,
                      nrow = length(unique(LambdaVecPl$NewAgeInd)),
                      byrow = FALSE)

#Differenced Log Death Rates (Creation of Z Matrix )
ZMatPl <- apply(log(LambdaMatPl), 1, diff) %>% t()



### 1.4 Save all data ##########################################################
save(LambdaMatUS, ZMatUS, LambdaVecUs, TotalDataUS,
     LambdaMatSp, ZMatSp, LambdaVecSp, TotDataSp,
     LambdaMatPl, ZMatPl, LambdaVecPl, TotDataPl,
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
