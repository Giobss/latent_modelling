rm(list=ls())

## INSTALL PACKAGES - if needed
#install.packages("pracma")
#install.packages("lmerTest")
#install.packages("emmeans")

library(lme4)
library(mediation)
library(readxl)
library(pracma)
library(lmerTest)
library(ggplot2)
library(emmeans)

## set working directory

targetfolder<-"/project/3022035.01/analysis_AMvol"
setwd(targetfolder)

## IMPORT DATA 

#(1) PRS

prs_dataset <- read.table('/project/3022035.01/PRS/v2-noid/PRSice.best',header=T) 
prs_specs <- read.table('/project/3022035.01/PRS/v2/PRSice.prsice',header=T) 
prs_summary<-read.table('/project/3022035.01/PRS/v2/PRSice.summary',header=T) 

## NB PRSice files imported as factors
#prs_dataset$FID<-as.numeric(as.character(prs_dataset$FID))
#prs_dataset$IID<-as.numeric(as.character(prs_dataset$IID))
#prs_dataset$In_Regression<-as.character(prs_dataset$In_Regression)

#(2) clinical data

clinical_data<-read_excel('/project/3022035.01/clinical/LEAP_t1_Core clinical variables_19-07-19-withvalues.xlsx')

#(3) subcortical data

subcortical_data<-read.csv("/project/3022035.01/anatomy/fs_all_output/SubVolumes.csv", header=T)

## REGRESSION OF CLINICAL OUTCOME ON GENETICS

# Match IDs across datasets

prsinclude<-match(clinical_data$subjects,prs_dataset$IID)

clininclude<-which(!is.na(prsinclude))
prsinclude<-prsinclude[clininclude]

clin<-clinical_data[clininclude,]
prs<-prs_dataset[prsinclude,]

# outcome variable of interest: VABS social score
availsub_vabs<-which(clin$t1_vabsdscoress_dss!=999)

# outcome variable of interest: SRS
availsub_srs<-which(clin$t1_srs_rawscore_combined!=999)

## CHOSEN OUTCOME SRS
availsub<-availsub_vabs

clin<-clin[availsub,]
prs<-prs[availsub,]

# covariates for the model: full scale IQ (fsiq)
noiq<-which(clin$t1_fsiq==999)

# final dataset
dataset<-data.frame(cbind(prs,clin))
names(dataset)<-c(names(prs),names(clin))

## CHECK NORMALITY OF INPUT VARIABLES AND OUTCOME

shapiro.test(dataset$t1_vabsdscoress_dss)
shapiro.test(dataset$t1_srs_rawscore_combined)
shapiro.test(dataset$PRS)
shapiro.test(dataset$t1_ageyrs)
shapiro.test(dataset$t1_fsiq)
shapiro.test(dataset$t1_site)

# scale PRS

dataset$PRS<-dataset$PRS*1000

# group factor - NEEDED????
dataset$group<-as.factor(dataset$group)

# mixed model
# test the interaction age:sex, age*sex
# add site
lin_model<-lmer(t1_srs_rawscore_combined~t1_ageyrs+t1_sex+PRS+t1_fsiq+t1_site+(1|group),dataset)

# to compare models use anova(,)
# to check the coefficients anova()

# check residuals
qqnorm(residuals(lin_model))
qqline(residuals(lin_model))

# check model quality
plot(lin_model)

# plot

p<-ggplot(data,aes(x=PRS,y=vabsdscoress_dss))+geom_point(aes(y=predict(lin_model),group=IID,colour=group))+
  scale_colour_manual(name="regression",values=c("green","red","blue","black"),labels=c("TD","ASD","ID","ID-ASD"))+geom_smooth(method=lm, se=FALSE, fullrange=TRUE)+labs(x="PRS*1000", y="VABS Social Score")+theme_bw()+
  scale_y_continuous(breaks = c(85,100,115), labels = c(85,100,115))+geom_smooth(method="lm")
  
print(p)

## MEDIATION ANALYSIS: B~G+A, where G=genetics; A=amygdala; B=behaviour.

results <- mediate(model.A, model.B, treat='PRS', mediator='amygdala',
                   boot=TRUE, sims=500)
summary(results)
  




