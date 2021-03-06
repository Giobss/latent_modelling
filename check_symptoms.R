load("~/VABS_lcga/R/latentVABS_classspecific_diag.RData")
data$class<-m3q$pprob$class
nnz(which(data$class==1))
length(which(data$class==1))
length(which(data$class==2))
length(which(data$class==3))
adosdata<-read.xlsx('C:/Users/giorgia/Documents/VABS_lcga/R/ados_SSARRB.xlsx')
load('xlsx')
library('xlsx')
adosdata<-read.xlsx('C:/Users/giorgia/Documents/VABS_lcga/R/ados_SSARRB.xlsx')
adosdata<-read.xlsx('C:/Users/giorgia/Documents/VABS_lcga/R/ados_SSARRB.xlsx',1)
match(adosdata$id,data$id)
check<-match(adosdata$id,data$id)
check<-match(data$id,adosdata$id)
check
ados_included<-adosdata[check,]
View(ados_included)
ados_included$class<-data$class
ados_included_clean<-ados_included
ados_included_clean[which(ados_included$ssa>776),]<-[]
ados_included_clean<-ados_included_clean[-which(ados_included$ssa>776),]
ados_included_clean<-ados_included_clean[-which(ados_included$rrb>776),]
View(ados_included_clean)
aov()
aov(ssa~class,ados_included_clean)
ssa_aov<-aov(ssa~class,ados_included_clean)
summary(ssa_aov)
TukeyHSD(ssa_aov)
library(multcomp)
summary(glht(ssa_aov, linfct = mcp(group = "Tukey")))
summary(glht(ssa_aov, linfct = mcp(class = "Tukey")))
ados_included_clean$class<-factor(ados_included_clean$class)
ssa_aov<-aov(ssa~class,ados_included_clean)
summary(ssa_aov)
aov(ssa~class,ados_included_clean)
rrb_aov<-aov(rrb~class,ados_included_clean)
aov(rrb~class,ados_included_clean)
summary(rrb_aov)
summary(glht(ssa_aov, linfct = mcp(class = "Tukey")))
summary(glht(rrb_aov, linfct = mcp(class = "Tukey")))
library(car)
leveneTest(weight~class,ados_included_clean)
leveneTest(class,ados_included_clean)
leveneTest(weight ~ class,ados_included_clean)
leveneTest(ssa ~ class,ados_included_clean)
leveneTest(rrb ~ class,ados_included_clean)
shapiro.test(residuals(ssa_aov))
shapiro.test(residuals(rrb_aov))
kruskal.test(rrb~class,ados_included_clean)
kruskal.test(ssa~class,ados_included_clean)
