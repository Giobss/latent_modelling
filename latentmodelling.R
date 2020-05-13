library(lcmm)
library(ggplot2)
library(gridExtra)

data<-read.table('C:/Users/giorgia/Documents/VABS_lcga/dataset_mselvabs_complete.dat')
names(data)<-c('id','phase','group','gender','elc1','elc2','elc3','elc4','age1','GM1','VR1','FM1','RL1','EL1','age2','GM2','VR2','FM2','RL2','EL2','age3','GM3','VR3','FM3','RL3','EL3','age4','VR4','FM4','RL4','EL4','ageV1','comm1','DL1','soc1','mot1','ageV2','comm2','DL2','soc2','mot2','ageV3','comm3','DL3','soc3','mot3','ageV4','comm4','DL4','soc4','mot4','outcome')

data$id<-as.numeric(as.character(data$id))
data$id[1]<-1

nooutcome<-which(data$outcome>776)
to_exclude<-c(nooutcome[-2],167,199)
data<-data[-to_exclude,]

comm<-c(data$comm1,data$comm2,data$comm3,data$comm4)
dl<-c(data$DL1,data$DL2,data$DL3,data$DL4)
mot<-c(data$mot1,data$mot2,data$mot3,data$mot4)
soc<-c(data$soc1,data$soc2,data$soc3,data$soc4)
age<-c(data$ageV1,data$ageV2,data$ageV3,data$ageV4)
gender<-c(data$gender,data$gender,data$gender,data$gender)
id<-c(data$id,data$id,data$id,data$id)

dataset<-as.data.frame(cbind(id,gender,age,comm,dl,soc,mot))

elc<-c(data$elc1,data$elc2,data$elc3,data$elc4)
elcnaindx<-which(elc>776)
elc[elcnaindx]<-NA

agenaindx<-which(dataset$age>776)
commnaindx<-which(dataset$comm>776)
dlnaindx<-which(dataset$dl>776)
motnaindx<-which(dataset$mot>776)
socnaindx<-which(dataset$soc>776)

dataset$age[agenaindx]<-NA
dataset$comm[commnaindx]<-NA
dataset$dl[dlnaindx]<-NA
dataset$mot[motnaindx]<-NA
dataset$soc[socnaindx]<-NA

#naindx<-which(is.na(dataset$age))
#dataset<-dataset[-naindx,]

# initial model with ng=1 for the random initial values
m1l <- multlcmm(comm+dl+mot+soc~age+gender,random=~age+gender,subject='id',ng=1,nwg=FALSE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1)

# gridsearch with 10 iterations from 50 random departures
m2l <- gridsearch(rep = 50, maxiter = 100, minit = m1l,multlcmm(comm+dl+mot+soc~age+gender,random=~age+gender,subject='id',mixture=~age+gender,ng=2,classmb=~age+gender,nwg=TRUE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1))
m3l <- gridsearch(rep = 50, maxiter = 100, minit = m1l,multlcmm(comm+dl+mot+soc~age+gender,random=~age+gender,subject='id',mixture=~age+gender,ng=3,classmb=~age+gender,nwg=TRUE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1))
m4l <- gridsearch(rep = 50, maxiter = 100, minit = m1l,multlcmm(comm+dl+mot+soc~age+gender,random=~age+gender,subject='id',mixture=~age+gender,ng=4,classmb=~age+gender,nwg=TRUE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1))
m5l <- gridsearch(rep = 50, maxiter = 100, minit = m1l,multlcmm(comm+dl+mot+soc~age+gender,random=~age+gender,subject='id',mixture=~age+gender,ng=5,classmb=~age+gender,nwg=TRUE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1))
m6l <- gridsearch(rep = 50, maxiter = 100, minit = m1l,multlcmm(comm+dl+mot+soc~age+gender,random=~age+gender,subject='id',mixture=~age+gender,ng=6,classmb=~age+gender,nwg=TRUE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1))

# initial model with ng=1 for the random initial values
m1q <- multlcmm(comm+dl+mot+soc~age+I(age^2)+gender,random=~age+I(age^2)+gender,subject='id',ng=1,nwg=FALSE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1)

# gridsearch with 10 iterations from 50 random departures
m2q <- gridsearch(rep = 50, maxiter = 100, minit = m1q,multlcmm(comm+dl+mot+soc~age+I(age^2)+gender,random=~age+I(age^2)+gender,subject='id',mixture=~age+I(age^2)+gender,ng=2,classmb=~age+I(age^2)+gender,nwg=TRUE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1))
m3q <- gridsearch(rep = 50, maxiter = 100, minit = m1q,multlcmm(comm+dl+mot+soc~age+I(age^2)+gender,random=~age+I(age^2)+gender,subject='id',mixture=~age+I(age^2)+gender,ng=3,classmb=~age+I(age^2)+gender,nwg=TRUE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1))
m4q <- gridsearch(rep = 50, maxiter = 100, minit = m1q,multlcmm(comm+dl+mot+soc~age+I(age^2)+gender,random=~age+I(age^2)+gender,subject='id',mixture=~age+I(age^2)+gender,ng=4,classmb=~age+I(age^2)+gender,nwg=TRUE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1))
m5q <- gridsearch(rep = 50, maxiter = 100, minit = m1q,multlcmm(comm+dl+mot+soc~age+I(age^2)+gender,random=~age+I(age^2)+gender,subject='id',mixture=~age+I(age^2)+gender,ng=5,classmb=~age+I(age^2)+gender,nwg=TRUE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1))
m6q <- gridsearch(rep = 50, maxiter = 100, minit = m1q,multlcmm(comm+dl+mot+soc~age+I(age^2)+gender,random=~age+I(age^2)+gender,subject='id',mixture=~age+I(age^2)+gender,ng=6,classmb=~age+I(age^2)+gender,nwg=TRUE,idiag=TRUE,data=dataset,link="linear",nsim=500,na.action=1))

dataset$id<- as.character(dataset$id)
people2 <- as.data.frame(m3q$pprob[,1:2])
dataset$group2 <- factor(people2$class[sapply(data$id, function(x) which(people2$id==x))])
dataset$elc <- elc

p1 <- ggplot(dataset, aes(age, soc, group=id, colour=group2)) + geom_point() + geom_smooth(aes(group=group2), method="loess", size=2, se=F)  +  labs(x="age (months)",y="Soc (standard score)",colour="Latent Class") + labs(title="Social Skills - Raw") + theme_minimal()
p2 <- ggplot(dataset, aes(age, soc, group=id, colour=group2)) + geom_point()+ geom_smooth(aes(group=id, colour=group2),size=0.3, se=F) + geom_smooth(aes(group=group2), method="loess", size=1.5, se=T)  + labs(x="age (months)",y="Soc (standard score)",colour="Group") + labs(title="Socialization Skills", legend.position="none") + theme_minimal()

p3 <- ggplot(dataset, aes(age, comm, group=id, colour=group2)) + geom_point() + geom_smooth(aes(group=group2), method="loess", size=2, se=F)  +  labs(x="age (months)",y="Comm (standard score)",colour="Latent Class") + labs(title="Communication Skills - Raw") + theme_minimal()
p4 <- ggplot(dataset, aes(age, comm, group=id, colour=group2)) + geom_point()+ geom_smooth(aes(group=id, colour=group2),size=0.3, se=F) + geom_smooth(aes(group=group2), method="loess", size=1.5, se=T)  + labs(x="age (months)",y="Comm (standard score)",colour="Group") + labs(title="Communication Skills", legend.position="none") + theme_minimal()

p5 <- ggplot(dataset, aes(age, dl, group=id, colour=group2)) + geom_point() + geom_smooth(aes(group=group2), method="loess", size=2, se=F)  +  labs(x="age (months)",y="DL (standard score)",colour="Latent Class") + labs(title="Daily Living Skills - Raw") + theme_minimal()
p6 <- ggplot(dataset, aes(age, dl, group=id, colour=group2)) + geom_point()+ geom_smooth(aes(group=id, colour=group2),size=0.3, se=F) + geom_smooth(aes(group=group2), method="loess", size=1.5, se=T)  + labs(x="age (months)",y="DL (standard score)",colour="Group") + labs(title="Daily Living Skills", legend.position="none") + theme_minimal()

p7 <- ggplot(dataset, aes(age, mot, group=id, colour=group2)) + geom_point() + geom_smooth(aes(group=group2), method="loess", size=2, se=F)  +  labs(x="age (months)",y="Mot (standard score)",colour="Latent Class") + labs(title="Motor Skills - Raw") + theme_minimal()
p8 <- ggplot(dataset, aes(age, mot, group=id, colour=group2)) + geom_point()+ geom_smooth(aes(group=id, colour=group2),size=0.3, se=F) + geom_smooth(aes(group=group2), method="loess", size=1.5, se=T)  + labs(x="age (months)",y="Mot (standard score)",colour="Group") + labs(title="Motor Skills", legend.position="none") + theme_minimal()

grid.arrange(p2,p4,p6,p8)

pelc<-ggplot(dataset, aes(age, elc, group=id, colour=group2)) + geom_point()+ geom_smooth(aes(group=id, colour=group2),size=0.3, se=F) + geom_smooth(aes(group=group2), method="loess", size=1.5, se=T)  + labs(x="age (months)",y="ELC (standard-score)",colour="Latent Class") + labs(title="Early Learning Composite Score - Smoothed", legend.position="none") + theme_minimal()

#+ scale_y_continuous(limits = c(13,37)) 