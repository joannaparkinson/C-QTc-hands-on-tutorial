## Created by: pharmacometrician
## Date: dd/mm/yyyy
## Deriving analysis data for exploratory plots and model

library(dplyr)
library(readr)

rm(list=ls())

## Import dataset
input.data <- read.csv("https://physionet.org/physiobank/database/ecgrdvq/SCR-002.Clinical.Data.csv")

# prepare analysis ready dataset
qtpk <- input.data %>% 
  
  ## calculate QTcF, QTcB, and HR (for each triplicate measurement)
       mutate(QTcF_3 = QT/((RR/1000)^(1/3)),       
              QTcB_3 = QT/((RR/1000)^0.5),         
              HR_3   = 60/(RR/1000)) %>%   
  #
  # Calculate mean of triplicate 
  group_by(RANDID, EXTRT, TPT) %>%   
  mutate(QTcF = mean(QTcF_3, na.rm=T),  
         QTcB = mean(QTcB_3, na.rm=T),
         QTm  = mean(QT, na.rm=T),
         RRm  = mean(RR, na.rm=T), 
         HR   = mean(HR_3, na.rm=T)) %>%
  ## Only keep first values by ID, TREATMENT, and TIME 
  distinct(RANDID, EXTRT, TPT, .keep_all=T) %>%  
  
  ## Extract individual baselines, calculate change from baseline for QTcF and HR 
  group_by(RANDID, EXTRT) %>%
  mutate(QTcF.B        = QTcF[BASELINE=="Y"],   
         QTcF.CFB      = QTcF-QTcF.B,     
         HR.B          = HR[BASELINE=="Y"],      
         HR.CFB        = HR-HR.B) %>%  
  
  ## Calculate ddQTcF and ddHR (parallel design) 
  group_by(RANDID, TPT) %>%
  mutate(Placebo.dQTcF      = QTcF.CFB[EXTRT=="Placebo"],
         Placebo.dHR        = HR.CFB[EXTRT=="Placebo"]) %>%
  group_by(TPT) %>%
  mutate(Placebo.dQTcF_mean = mean(Placebo.dQTcF, na.rm=T), 
         ddQTcF             = QTcF.CFB-Placebo.dQTcF_mean,
         Placebo.dHR_mean   = mean(Placebo.dHR, na.rm=T), 
         ddHR               = HR.CFB-Placebo.dHR_mean
         ) %>% 
  ungroup() %>%
  
  ## Impute 0 conc for placebo observations and predose observations that are NA   
  mutate(PCSTRESN = ifelse(EXTRT=="Placebo", 0, PCSTRESN),      
         PCSTRESN = ifelse(TPT<0 & is.na(PCSTRESN), 0, PCSTRESN),
         ## Add indicator for active treatment (1==ACTIVE)
         ACTIVE   = as.factor(ifelse(EXTRT=="Placebo", 0, 1))) %>%
  
  ## Rename variables to more convenient names 
  rename(TREAT     = EXTRT,          
         TIME      = TPT, 
         CONC.DRUG = PCSTRESN) %>%
  
  ## include only relevant arms used in this example
  filter(TREAT=="Dofetilide" | TREAT=="Placebo" | TREAT=="Verapamil HCL") 

## Assign new ID to each subject in each arm so make the study parallel design
qtpk<-mutate(qtpk, USUBJID=rep(seq(from=1, to=length(unique(qtpk$RANDID))*length(unique(qtpk$TREAT))),each=length(unique(qtpk$TIME)))) %>%
  select(USUBJID, everything())

## create dofetilide dataset

qtpk_dofetilide <- qtpk %>%
  filter(TREAT %in% c("Dofetilide","Placebo")) %>%
  
  ## calculate population mean baseline QTcF
  mutate(QTcF.mB   = mean(QTcF[BASELINE=="Y"]),
         CONC      = CONC.DRUG/1000 # changing concentration units to ng/ml
         ) %>%
  #select relevant columns
  select(USUBJID,TREAT,ACTIVE,TIME,QTcF,QTcB,QTcF.CFB,HR.CFB,ddHR,ddQTcF,CONC,QTm,RRm,QTcF.mB,QTcF.B)

## create verapamil dataset

qtpk_verapamil <- qtpk %>%
  filter(TREAT %in% c("Verapamil HCL","Placebo")) %>%
  
  ## calculate population mean baseline QTcF
  mutate(QTcF.mB   = mean(QTcF[BASELINE=="Y"]),
         CONC      = CONC.DRUG) %>%
  #select relevant columns
  select(USUBJID,TREAT,ACTIVE,TIME,QTcF,QTcB,QTcF.CFB,HR.CFB,ddHR,ddQTcF,CONC,QTm,RRm,QTcF.mB,QTcF.B)

## save data - here, the data is saved in the "data/derived" subfolder
write.csv(qtpk_dofetilide, "data/derived/qtpk_dofetilide.csv")
write.csv(qtpk_verapamil, "data/derived/qtpk_verapamil.csv")

