rm(list=ls()) #Clear workspace

file.choose()
dat<-data.frame(read.csv("/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/photogramm_measure_2015.csv"))


#####Functions#####
func <- function(x, y, par) {
  B = par[1]
  D1 = par[2]
  D2 = par[3]
  D3 = par[4]
  C1 = par[5]
  C2 = par[6]
  C3 = par[7]
  yhat<-B*(1 - 1/(1 +((x/D1)^C1)+((x/D2)^C2)+((x/D3)^C3))) #JPPS equation
  res<-y-yhat #residuals
  sum((res)^2) # optimize by minimizing residual sums of squares
}

reg.func<-function(par,a,b){
  b0<-par[1]
  b1<-par[2]
  bhat<-b0 + b1 * a # linear model
  res<-b-bhat #residuals
  sum(abs(res)^2) #residual sums of squares
}

weight.func<-function(x, y, par, w){
  B = par[1]
  D1 = par[2]
  D2 = par[3]
  D3 = par[4]
  C1 = par[5]
  C2 = par[6]
  C3 = par[7]
  yhat<-B*(1 - 1/(1 +((x/D1)^C1)+((x/D2)^C2)+((x/D3)^C3)))
  res<-y-yhat
  sum(w*(res)^2) # minimizing sums of squares with weighting matrix w
}


####### Size Bootstrap #######

#To standardize where we guess the parameter B (the assymptote), use the mean of >20 years#adults range >15 
adult.mean<-mean(dat$Avg.Size[dat$Age > 20]) 
adult.mean.female<-mean(dat$Avg.Size[dat$Age > 20 & dat$Sex=='0'])
adult.mean.male<-mean(dat$Avg.Size[dat$Age > 20 & dat$Sex=='1'])

####### Weighted Size Bootstrap #######
#First set of "original" estimates generated using simulated annealing (SANN)
mod.new<- optim(par=c(B=adult.mean, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age, y=dat$Avg.Size,method="SANN",control=list(maxit=10000,trace=T)); mod.new 
plot(Avg.Size~Age,data=dat,main='Average Body Size')
curve(mod.new$par[1]*(1-1/(1 +(x/mod.new$par[2])^mod.new$par[5]+(x/mod.new$par[3])^mod.new$par[6]+(x/mod.new$par[4])^mod.new$par[7])), add=TRUE,lwd=2,lty=2,col='black')
mod.new<-mod.old

w<-rep(1,length(dat$Age)) #initial weighting covariance matrix is all 1s, so it will have no effect (like multiplying any number by 1)

tol=1e-9 #tolerance
maxit=500 #max iterations for one pass of the iteratively weighted optimization
iter=1000 #max number of total samples in the bootstrap

mod<-matrix(nrow=iter,ncol=7) #parameter matrix 
res.tot<-matrix(nrow=iter,ncol=length(dat$Age)) #residual matrix

#Double for loop to generate bootstrap estimates
for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new<- optim(par=c(B=adult.mean, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age, y=dat$Avg.Size,method="SANN",control=list(maxit=10000,trace=F))
  for (i in 1:maxit){ #here's the iterative weighting for loop
    oldcof<-as.vector(mod.new$par)
    B<-oldcof[1]
    D1<-oldcof[2]
    D2<-oldcof[3]
    D3<-oldcof[4]
    C1<-oldcof[5]
    C2<-oldcof[6]
    C3<-oldcof[7]
    yhat<-B*(1 - 1/(1 +((dat$Age/D1)^C1)+((dat$Age/D2)^C2)+((dat$Age/D3)^C3))) #JPPS eq.
    res<-dat$Avg.Size-yhat #residuals
    mod.res<-optim(par=c(.5,.5),reg.func,a=dat$Avg.Size,b=res,control=list(abstol=tol),method="Nelder") #this models any linear trend to the residuals; essentially the same as lm() function
    gam0<-mod.res$par[1] 
    gam1<-mod.res$par[2]
    w.new<-abs(1/(gam0 + gam1 * dat$Avg.Size)) # estimate new weights  
    mod.new<-optim(par=as.vector(mod.new$par),weight.func,x=dat$Age,y=dat$Avg.Size,w=w.new,control=list(abstol=tol),method="SANN") #run optimization again with new weights
    newcof<-as.vector(mod.new$par) 
    dif<-sum(sqrt((oldcof-newcof)^2))/length(newcof) #compare the difference; check tolerance with last pass of iterative weighting
    dif<-sqrt((oldcof[1]-newcof[1])^2) 
    if(dif<=tol){break}  
  }
  mod[j,]<-mod.new$par #save jth model parameters
  res.tot[j,]<-res #save jth residuals
}

mod<-matrix(mod[complete.cases(mod)],ncol=7);mod
mod.mean<-apply( mod, 2, median, na.rm=T);mod.mean
mod.sd<-apply( mod, 2, sd, na.rm=T);mod.sd

plot(dat$Avg.Size~dat$Age, main='Male and Female Iteratively Weighted',pch=19,col=c('pink','blue')[dat$Sex],xlab='Age',ylab='Avg. Size'); legend(x=30,y=900, legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2))
curve(mod.mean[1]*(1-1/(1 +(x/mod.mean[2])^mod.mean[5]+(x/mod.mean[3])^mod.mean[6]+(x/mod.mean[4])^mod.mean[7])), add=TRUE,lwd=1,lty=1,col='black')

curve(mod.old$par[1]*(1-1/(1 +(x/mod.old$par[2])^mod.old$par[5]+(x/mod.old$par[3])^mod.old$par[6]+(x/mod.old$par[4])^mod.old$par[7])), add=TRUE,lwd=1,lty=2,col='black')

curve(mod.mean[1] * (((x/mod.mean[2])^(mod.mean[5] - 1) * (mod.mean[5] * (1/mod.mean[2])) + (x/mod.mean[3])^(mod.mean[6] - 1) * (mod.mean[6] * (1/mod.mean[3])) + (x/mod.mean[4])^(mod.mean[7] - 1) * (mod.mean[7] * (1/mod.mean[4])))/(1 + ((x/mod.mean[2])^mod.mean[5]) + ((x/mod.mean[3])^mod.mean[6]) + ((x/mod.mean[4])^mod.mean[7]))^2),lty=1, xlim=c(.5,25),main='Male and Female Iteratively Weighted',xlab='Age',ylab='Growth Velocity')

curve(mod.old$par[1] * (((x/mod.old$par[2])^(mod.old$par[5] - 1) * (mod.old$par[5] * (1/mod.old$par[2])) + (x/mod.old$par[3])^(mod.old$par[6] - 1) * (mod.old$par[6] * (1/mod.old$par[3])) + (x/mod.old$par[4])^(mod.old$par[7] - 1) * (mod.old$par[7] * (1/mod.old$par[4])))/(1 + ((x/mod.old$par[2])^mod.old$par[5]) + ((x/mod.old$par[3])^mod.old$par[6]) + ((x/mod.old$par[4])^mod.old$par[7]))^2), add=T, lty=2);legend(x=15,y=120,legend=c('Bootstrap','Original'),lty=c(1,2))

#Export
exp.weighted.data<-as.data.frame(mod)  #create data frame to export
colnames(exp.weighted.data)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.weighted.data)# check first three rows of data frame

write.table(exp.weighted.data, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/male&female_combined_weightedbootstrap.csv", sep = ",", row.names=F) #Export

exp.weighted.resid<-as.data.frame(res.tot)  #create data frame to export
head(exp.weighted.resid)# check first three rows of data frame

write.table(exp.weighted.resid, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/weightedresid.csv", sep = ",", row.names=F) #Export

###### Non-weighted Size Bootstrap #######
tol=1e-9
maxit=500
iter=1000

mod.2<-matrix(nrow=iter,ncol=7)
res.tot.2<-matrix(nrow=iter,ncol=length(dat$Age))

for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new<- optim(par=c(B=adult.mean, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age, y=dat$Avg.Size,method="SANN",control=list(maxit=10000,trace=F))
    oldcof<-as.vector(mod.new$par)
    B<-oldcof[1]
    D1<-oldcof[2]
    D2<-oldcof[3]
    D3<-oldcof[4]
    C1<-oldcof[5]
    C2<-oldcof[6]
    C3<-oldcof[7]
    yhat<-B*(1 - 1/(1 +((dat$Age/D1)^C1)+((dat$Age/D2)^C2)+((dat$Age/D3)^C3)))
    res<-dat$Avg.Size-yhat
  mod.2[j,]<-mod.new$par
  res.tot.2[j,]<-res
}

mod.2<-matrix(mod.2[complete.cases(mod.2)],ncol=7);mod.2
mod.2.mean<-apply( mod.2, 2, median, na.rm=T);mod.2.mean
mod.2.sd<-apply( mod.2, 2, sd, na.rm=T);mod.2.sd

plot(dat$Avg.Size~dat$Age, main='Male and Female Non-weighted',pch=19,col=c('pink','blue')[dat$Sex],xlab='Age',ylab='Avg. Size'); legend(x=30,y=900, legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2),lwd=c(NA,NA,1,2))
curve(mod.2.mean[1]*(1-1/(1 +(x/mod.2.mean[2])^mod.2.mean[5]+(x/mod.2.mean[3])^mod.2.mean[6]+(x/mod.2.mean[4])^mod.2.mean[7])), add=TRUE,lwd=1,lty=1)

curve(mod.old$par[1]*(1-1/(1 +(x/mod.old$par[2])^mod.old$par[5]+(x/mod.old$par[3])^mod.old$par[6]+(x/mod.old$par[4])^mod.old$par[7])), add=TRUE,lwd=2,lty=2,col='black')

curve(mod.2.mean[1] * (((x/mod.2.mean[2])^(mod.2.mean[5] - 1) * (mod.2.mean[5] * (1/mod.2.mean[2])) + (x/mod.2.mean[3])^(mod.2.mean[6] - 1) * (mod.2.mean[6] * (1/mod.2.mean[3])) + (x/mod.2.mean[4])^(mod.2.mean[7] - 1) * (mod.2.mean[7] * (1/mod.2.mean[4])))/(1 + ((x/mod.2.mean[2])^mod.2.mean[5]) + ((x/mod.2.mean[3])^mod.2.mean[6]) + ((x/mod.2.mean[4])^mod.2.mean[7]))^2), xlim=c(.5,25),main='Male and Female Non-Weighted',xlab='Age',ylab='Growth Velocity')

curve(mod.old$par[1] * (((x/mod.old$par[2])^(mod.old$par[5] - 1) * (mod.old$par[5] * (1/mod.old$par[2])) + (x/mod.old$par[3])^(mod.old$par[6] - 1) * (mod.old$par[6] * (1/mod.old$par[3])) + (x/mod.old$par[4])^(mod.old$par[7] - 1) * (mod.old$par[7] * (1/mod.old$par[4])))/(1 + ((x/mod.old$par[2])^mod.old$par[5]) + ((x/mod.old$par[3])^mod.old$par[6]) + ((x/mod.old$par[4])^mod.old$par[7]))^2), add=T,lwd=1,lty=2);legend(x=15,y=120,legend=c('Bootstrap','Original'),lty=c(1,2))

#Export
exp.data<-as.data.frame(mod.2)  #create data frame to export
colnames(exp.data)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.data)# check first three rows of data frame

write.table(exp.data, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/combined_male&female_non-weighted_bootstrap_Dec1.csv", sep = ",", row.names=F) #Export

exp.resid<-as.data.frame(res.tot.2)  #create data frame to export
head(exp.resid)# check first three rows of data frame

write.table(exp.resid, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/combined_male&female_nonweighted_resid_Dec1.csv", sep = ",", row.names=F) #Export


####### Female Iteratively Weighted Size Bootstrap #######

mod.new.female<- optim(par=c(B=adult.mean.female, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age[which(dat$Sex=='0')], y=dat$Avg.Size[which(dat$Sex=='0')], method="SANN",control=list(maxit=10000,trace=T));mod.new.female 
plot(Avg.Size~Age,data=dat[which(dat$Sex=='0'),],main='Average Body Size',col='red')
curve(mod.new.female$par[1]*(1-1/(1 +(x/mod.new.female$par[2])^mod.new.female$par[5]+(x/mod.new.female$par[3])^mod.new.female$par[6]+(x/mod.new.female$par[4])^mod.new.female$par[7])), add=TRUE,lwd=2,lty=2)
mod.old.female<-mod.new.female

w<-rep(1,length(dat$Age[which(dat$Sex=='0')]))

tol=1e-9
maxit=500
iter=1000

mod.female<-matrix(nrow=iter,ncol=7)
res.tot.female<-matrix(nrow=iter,ncol=length(dat$Age[which(dat$Sex=='0')]))

for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new.female<- optim(par=c(B=adult.mean.female, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age[which(dat$Sex=='0')], y=dat$Avg.Size[which(dat$Sex=='0')], method="SANN",control=list(maxit=10000,trace=F));mod.new.female 
  for (i in 1:maxit){
    oldcof<-as.vector(mod.new.female$par)
    B<-oldcof[1]
    D1<-oldcof[2]
    D2<-oldcof[3]
    D3<-oldcof[4]
    C1<-oldcof[5]
    C2<-oldcof[6]
    C3<-oldcof[7]
    fem.age<-dat$Age[which(dat$Sex=='0')]
    fem.size<-dat$Avg.Size[which(dat$Sex=='0')]
    yhat<-B*(1 - 1/(1 +((fem.age/D1)^C1)+((fem.age/D2)^C2)+((fem.age/D3)^C3)))
    res<- fem.size-yhat
    mod.res<-optim(par=c(.5,.5),reg.func,a=fem.size,b=res,control=list(abstol=tol),method="Nelder")
    gam0<-mod.res$par[1]
    gam1<-mod.res$par[2]
    w.new<-abs(1/(gam0 + gam1 * fem.size)) # estimate new weights  
    mod.new<-optim(par=as.vector(mod.new.female$par),weight.func,x=fem.age,y=fem.size,w=w.new,control=list(abstol=tol),method="SANN")
    newcof<-as.vector(mod.new.female$par)
    dif<-sum(sqrt((oldcof-newcof)^2))/length(newcof)
    dif<-sqrt((oldcof[1]-newcof[1])^2) # only a single parameter
    if(dif<=tol){break}  
  }
  mod.female[j,]<-mod.new.female$par
  res.tot.female[j,]<-res
}

mod.female<-matrix(mod.female[complete.cases(mod.female)],ncol=7);mod.female
mod.female.mean<-apply(mod.female, 2, median, na.rm=T);mod.female.mean
mod.female.sd<-apply( mod.female, 2, sd, na.rm=T);mod.female.sd

plot(dat$Avg.Size~dat$Age, main='Iteratively Weighted',pch=19,col=c('pink','blue')[dat$Sex],xlab='Age',ylab='Avg. Size'); legend(x=30,y=900, legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2),lwd=c(NA,NA,1,2))

curve(mod.old.female$par[1]*(1-1/(1 +(x/mod.old.female$par[2])^mod.old.female$par[5]+(x/mod.old.female$par[3])^mod.old.female$par[6]+(x/mod.old.female$par[4])^mod.old.female$par[7])), add=TRUE,lwd=2,lty=2,col='pink')

curve(mod.female.mean[1]*(1-1/(1 +(x/mod.female.mean[2])^mod.female.mean[5]+(x/mod.female.mean[3])^mod.female.mean[6]+(x/mod.female.mean[4])^mod.female.mean[7])), add=TRUE,lwd=1,lty=1,col='pink')

curve(mod.female.mean[1] * (((x/mod.female.mean[2])^(mod.female.mean[5] - 1) * (mod.female.mean[5] * (1/mod.female.mean[2])) + (x/mod.female.mean[3])^(mod.female.mean[6] - 1) * (mod.female.mean[6] * (1/mod.female.mean[3])) + (x/mod.female.mean[4])^(mod.female.mean[7] - 1) * (mod.female.mean[7] * (1/mod.female.mean[4])))/(1 + ((x/mod.female.mean[2])^mod.female.mean[5]) + ((x/mod.female.mean[3])^mod.female.mean[6]) + ((x/mod.female.mean[4])^mod.female.mean[7]))^2), xlim=c(0,25),,col='pink',ylab='Growth Velocity',xlab='Age',main='Iteratively Weighted Growth Curve');legend(x=15,y=120,legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2),lwd=c(NA,NA,1,2))

curve(mod.old.female$par[1] * (((x/mod.old.female$par[2])^(mod.old.female$par[5] - 1) * (mod.old.female$par[5] * (1/mod.old.female$par[2])) + (x/mod.old.female$par[3])^(mod.old.female$par[6] - 1) * (mod.old.female$par[6] * (1/mod.old.female$par[3])) + (x/mod.old.female$par[4])^(mod.old.female$par[7] - 1) * (mod.old.female$par[7] * (1/mod.old.female$par[4])))/(1 + ((x/mod.old.female$par[2])^mod.old.female$par[5]) + ((x/mod.old.female$par[3])^mod.old.female$par[6]) + ((x/mod.old.female$par[4])^mod.old.female$par[7]))^2), lty=2,lwd=2,col='pink',add=T)

#Export
exp.female.data<-as.data.frame(mod.female)  #create data frame to export
colnames(exp.female.data)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.female.data)# check first three rows of data frame

write.table(exp.female.data, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/female_weighted_bootstrap_Aug3.csv", sep = ",", row.names=F) #Export

exp.female.resid<-as.data.frame(res.tot.female)  #create data frame to export
head(exp.female.resid)# check first three rows of data frame

write.table(exp.female.resid, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/female_weighted_resid_Aug3.csv", sep = ",", row.names=F) #Export


###### Female Non-weighted Size Bootstrap #######
tol=1e-9
maxit=500
iter=1000

mod.female.2<-matrix(nrow=iter,ncol=7)
res.tot.2.female<-matrix(nrow=iter,ncol=length(dat$Age[which(dat$Sex=='0')]))

for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new.female<- optim(par=c(B=adult.mean.female, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age[which(dat$Sex=='0')], y=dat$Avg.Size[which(dat$Sex=='0')], method="SANN",control=list(maxit=10000,trace=F));mod.new.female 
  oldcof<-as.vector(mod.new.female$par)
  B<-oldcof[1]
  D1<-oldcof[2]
  D2<-oldcof[3]
  D3<-oldcof[4]
  C1<-oldcof[5]
  C2<-oldcof[6]
  C3<-oldcof[7]
  fem.age<-dat$Age[which(dat$Sex=='0')]
  fem.size<-dat$Avg.Size[which(dat$Sex=='0')]
  yhat<-B*(1 - 1/(1 +((fem.age/D1)^C1)+((fem.age/D2)^C2)+((fem.age/D3)^C3)))
  res<-fem.size-yhat
  mod.female.2[j,]<-mod.new.female$par
  res.tot.2.female[j,]<-res
}

mod.female.2<-matrix(mod.female.2[complete.cases(mod.female.2)],ncol=7);mod.female.2
mod.female.2.mean<-apply( mod.female.2, 2, median, na.rm=T);mod.female.2.mean
mod.female.2.sd<-apply( mod.female.2, 2, sd, na.rm=T);mod.female.2.sd

plot(dat$Avg.Size~dat$Age, main='Non-Weighted female length',pch=19,col=c('pink','blue')[dat$Sex],xlab='Age',ylab='Avg. Size'); legend(x=30,y=900, legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2),lwd=c(NA,NA,1,2))
curve(mod.female.2.mean[1]*(1-1/(1 +(x/mod.female.2.mean[2])^mod.female.2.mean[5]+(x/mod.female.2.mean[3])^mod.female.2.mean[6]+(x/mod.female.2.mean[4])^mod.female.2.mean[7])), add=TRUE,lwd=1,lty=1,col='pink')

curve(mod.old.female$par[1]*(1-1/(1 +(x/mod.old.female$par[2])^mod.old.female$par[5]+(x/mod.old.female$par[3])^mod.old.female$par[6]+(x/mod.old.female$par[4])^mod.old.female$par[7])), add=TRUE,lwd=2,lty=2,col='pink')

curve(mod.female.2.mean[1] * (((x/mod.female.2.mean[2])^(mod.female.2.mean[5] - 1) * (mod.female.2.mean[5] * (1/mod.female.2.mean[2])) + (x/mod.female.2.mean[3])^(mod.female.2.mean[6] - 1) * (mod.female.2.mean[6] * (1/mod.female.2.mean[3])) + (x/mod.female.2.mean[4])^(mod.female.2.mean[7] - 1) * (mod.female.2.mean[7] * (1/mod.female.2.mean[4])))/(1 + ((x/mod.female.2.mean[2])^mod.female.2.mean[5]) + ((x/mod.female.2.mean[3])^mod.female.2.mean[6]) + ((x/mod.female.2.mean[4])^mod.female.2.mean[7]))^2),col='pink', xlim=c(.5,25),ylim=c(0,150),ylab='Growth Velocity',xlab='Age',main='Non-weighted Growth Curve');legend(x=15,y=120,legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2),lwd=c(NA,NA,1,2))

curve(mod.old.female$par[1] * (((x/mod.old.female$par[2])^(mod.old.female$par[5] - 1) * (mod.old.female$par[5] * (1/mod.old.female$par[2])) + (x/mod.old.female$par[3])^(mod.old.female$par[6] - 1) * (mod.old.female$par[6] * (1/mod.old.female$par[3])) + (x/mod.old.female$par[4])^(mod.old.female$par[7] - 1) * (mod.old.female$par[7] * (1/mod.old.female$par[4])))/(1 + ((x/mod.old.female$par[2])^mod.old.female$par[5]) + ((x/mod.old.female$par[3])^mod.old.female$par[6]) + ((x/mod.old.female$par[4])^mod.old.female$par[7]))^2), col='pink',lwd=2,lty=2,add=T)

#Export
exp.female.data2<-as.data.frame(mod.female.2)  #create data frame to export
colnames(exp.female.data2)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.female.data2)# check first three rows of data frame

write.table(exp.female.data2, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/non-weighted_female_bootstrap_Aug3.csv", sep = ",", row.names=F) #Export

exp.female.resid2<-as.data.frame(res.tot.2.female)  #create data frame to export
head(exp.female.resid2)# check first three rows of data frame

write.table(exp.female.resid2, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/nonweighted_female_resid_Aug3.csv", sep = ",", row.names=F) #Export


####### Male Iteratively Weighted Size Bootstrap #######

mod.new.male <-optim(par=c(B=adult.mean.male, D1= 9.19 , D2=9.97, D3=12.12, C1= 4.79, C2= 19.03, C3=1.22), fn=func, x=dat$Age[which(dat$Sex=='1')], y=dat$Avg.Size[which(dat$Sex=='1')], method="SANN",control=list(maxit=10000,trace=T));mod.new.male
plot(Avg.Size~Age,data=dat[which(dat$Sex=='1'),],main='Average Body Size',col='blue')
curve(mod.new.male$par[1]*(1-1/(1 +(x/mod.new.male$par[2])^mod.new.male$par[5]+(x/mod.new.male$par[3])^mod.new.male$par[6]+(x/mod.new.male$par[4])^mod.new.male$par[7])), add=TRUE,lwd=2,lty=2)
mod.old.male<-mod.new.male

w<-rep(1,length(dat$Age[which(dat$Sex=='1')]))

tol=1e-9
maxit=500
iter=1000

mod.male<-matrix(nrow=iter,ncol=7)
res.tot.male<-matrix(nrow=iter,ncol=length(dat$Age[which(dat$Sex=='1')]))

for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new.male<- optim(par=c(B=adult.mean.male, D1= 9.19 , D2=9.97, D3=12.12, C1= 4.79, C2= 19.03, C3=1.22), fn=func, x=dat$Age[which(dat$Sex=='1')], y=dat$Avg.Size[which(dat$Sex=='1')], method="SANN",control=list(maxit=10000,trace=F));mod.new.male 
  for (i in 1:maxit){
    oldcof<-as.vector(mod.new.male$par)
    B<-oldcof[1]
    D1<-oldcof[2]
    D2<-oldcof[3]
    D3<-oldcof[4]
    C1<-oldcof[5]
    C2<-oldcof[6]
    C3<-oldcof[7]
    male.age<-dat$Age[which(dat$Sex=='1')]
    male.size<-dat$Avg.Size[which(dat$Sex=='1')]
    yhat<-B*(1 - 1/(1 +((male.age/D1)^C1)+((male.age/D2)^C2)+((male.age/D3)^C3)))
    res<- male.size-yhat
    mod.res<-optim(par=c(.5,.5),reg.func,a=male.size,b=res,control=list(abstol=tol),method="Nelder")
    gam0<-mod.res$par[1]
    gam1<-mod.res$par[2]
    w.new<-abs(1/(gam0 + gam1 * male.size)) # estimate new weights  
    mod.new<-optim(par=as.vector(mod.new.male$par),weight.func,x=male.age,y=male.size,w=w.new,control=list(abstol=tol),method="SANN")
    newcof<-as.vector(mod.new.male$par)
    dif<-sum(sqrt((oldcof-newcof)^2))/length(newcof)
    dif<-sqrt((oldcof[1]-newcof[1])^2) # only a single parameter
    if(dif<=tol){break}  
  }
  mod.male[j,]<-mod.new.male$par
  res.tot.male[j,]<-res
}

mod.male<-matrix(mod.male[complete.cases(mod.male)],ncol=7);mod.male
mod.male.mean<-apply(mod.male, 2, median, na.rm=T);mod.male.mean
mod.male.sd<-apply( mod.male, 2, sd, na.rm=T);mod.male.sd

plot(Avg.Size~Age,data=dat[which(dat$Sex=='1'),],main='Average Body Size',col='blue')

curve(mod.old.male$par[1]*(1-1/(1 +(x/mod.old.male$par[2])^mod.old.male$par[5]+(x/mod.old.male$par[3])^mod.old.male$par[6]+(x/mod.old.male$par[4])^mod.old.male$par[7])), add=TRUE,lwd=2,lty=2,col='blue')

curve(mod.male.mean[1]*(1-1/(1 +(x/mod.male.mean[2])^mod.male.mean[5]+(x/mod.male.mean[3])^mod.male.mean[6]+(x/mod.male.mean[4])^mod.male.mean[7])), add=TRUE,lwd=1,lty=1,col='blue')

curve(mod.male.mean[1] * (((x/mod.male.mean[2])^(mod.male.mean[5] - 1) * (mod.male.mean[5] * (1/mod.male.mean[2])) + (x/mod.male.mean[3])^(mod.male.mean[6] - 1) * (mod.male.mean[6] * (1/mod.male.mean[3])) + (x/mod.male.mean[4])^(mod.male.mean[7] - 1) * (mod.male.mean[7] * (1/mod.male.mean[4])))/(1 + ((x/mod.male.mean[2])^mod.male.mean[5]) + ((x/mod.male.mean[3])^mod.male.mean[6]) + ((x/mod.male.mean[4])^mod.male.mean[7]))^2),col='blue',add=T)

curve(mod.old.male$par[1] * (((x/mod.old.male$par[2])^(mod.old.male$par[5] - 1) * (mod.old.male$par[5] * (1/mod.old.male$par[2])) + (x/mod.old.male$par[3])^(mod.old.male$par[6] - 1) * (mod.old.male$par[6] * (1/mod.old.male$par[3])) + (x/mod.old.male$par[4])^(mod.old.male$par[7] - 1) * (mod.old.male$par[7] * (1/mod.old.male$par[4])))/(1 + ((x/mod.old.male$par[2])^mod.old.male$par[5]) + ((x/mod.old.male$par[3])^mod.old.male$par[6]) + ((x/mod.old.male$par[4])^mod.old.male$par[7]))^2), lwd=2,lty=2,col='blue',add=T)

#Export
exp.male.data<-as.data.frame(mod.male)  #create data frame to export
colnames(exp.male.data)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.male.data)# check first three rows of data frame

write.table(exp.male.data, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/weighted_male_bootstrap_Aug3.csv", sep = ",", row.names=F) #Export

exp.male.resid<-as.data.frame(res.tot.male)  #create data frame to export
head(exp.male.resid)# check first three rows of data frame

write.table(exp.male.resid, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/weighted_male_resid_Aug3.csv", sep = ",", row.names=F) #Export

###### Male Non-weighted Size Bootstrap #######
tol=1e-9
maxit=500
iter=1000

mod.male.2<-matrix(nrow=iter,ncol=7)
res.tot.2.male<-matrix(nrow=iter,ncol=length(dat$Age[which(dat$Sex=='1')]))

for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new.male<- optim(par=c(B=adult.mean.male, D1= 9.19 , D2=9.97, D3=12.12, C1= 4.79, C2= 19.03, C3=1.22), fn=func, x=dat$Age[which(dat$Sex=='1')], y=dat$Avg.Size[which(dat$Sex=='1')], method="SANN",control=list(maxit=10000,trace=F));mod.new.male 
  oldcof<-as.vector(mod.new.male$par)
  B<-oldcof[1]
  D1<-oldcof[2]
  D2<-oldcof[3]
  D3<-oldcof[4]
  C1<-oldcof[5]
  C2<-oldcof[6]
  C3<-oldcof[7]
  male.age<-dat$Age[which(dat$Sex=='1')]
  male.size<-dat$Avg.Size[which(dat$Sex=='1')]
  yhat<-B*(1 - 1/(1 +((male.age/D1)^C1)+((male.age/D2)^C2)+((male.age/D3)^C3)))
  res<-male.size-yhat
  mod.male.2[j,]<-mod.new.male$par
  res.tot.2.male[j,]<-res
}

mod.male.2<-matrix(mod.male.2[complete.cases(mod.male.2)],ncol=7);mod.male.2
mod.male.2.mean<-apply( mod.male.2, 2, median, na.rm=T);mod.male.2.mean
mod.male.2.sd<-apply( mod.male.2, 2, sd, na.rm=T);mod.male.2.sd

plot(Avg.Size~Age,data=dat[which(dat$Sex=='1'),],main='Average Body Size',col='red')
curve(mod.male.2.mean[1]*(1-1/(1 +(x/mod.male.2.mean[2])^mod.male.2.mean[5]+(x/mod.male.2.mean[3])^mod.male.2.mean[6]+(x/mod.male.2.mean[4])^mod.male.2.mean[7])), add=TRUE,lwd=1,lty=1,col='blue')

curve(mod.old.male$par[1]*(1-1/(1 +(x/mod.old.male$par[2])^mod.old.male$par[5]+(x/mod.old.male$par[3])^mod.old.male$par[6]+(x/mod.old.male$par[4])^mod.old.male$par[7])), add=TRUE,lwd=2,lty=2,col='blue')

curve(mod.male.2.mean[1] * (((x/mod.male.2.mean[2])^(mod.male.2.mean[5] - 1) * (mod.male.2.mean[5] * (1/mod.male.2.mean[2])) + (x/mod.male.2.mean[3])^(mod.male.2.mean[6] - 1) * (mod.male.2.mean[6] * (1/mod.male.2.mean[3])) + (x/mod.male.2.mean[4])^(mod.male.2.mean[7] - 1) * (mod.male.2.mean[7] * (1/mod.male.2.mean[4])))/(1 + ((x/mod.male.2.mean[2])^mod.male.2.mean[5]) + ((x/mod.male.2.mean[3])^mod.male.2.mean[6]) + ((x/mod.male.2.mean[4])^mod.male.2.mean[7]))^2),col='blue',add=T)

curve(mod.old.male$par[1] * (((x/mod.old.male$par[2])^(mod.old.male$par[5] - 1) * (mod.old.male$par[5] * (1/mod.old.male$par[2])) + (x/mod.old.male$par[3])^(mod.old.male$par[6] - 1) * (mod.old.male$par[6] * (1/mod.old.male$par[3])) + (x/mod.old.male$par[4])^(mod.old.male$par[7] - 1) * (mod.old.male$par[7] * (1/mod.old.male$par[4])))/(1 + ((x/mod.old.male$par[2])^mod.old.male$par[5]) + ((x/mod.old.male$par[3])^mod.old.male$par[6]) + ((x/mod.old.male$par[4])^mod.old.male$par[7]))^2), lwd=2,lty=2,col='blue',add=T)

#Export
exp.male.data2<-as.data.frame(mod.male.2)  #create data frame to export
colnames(exp.male.data2)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.male.data2)# check first three rows of data frame

write.table(exp.male.data2, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/nonweighted_male_bootstrap_Aug3.csv", sep = ",", row.names=F) #Export

exp.male.resid2<-as.data.frame(res.tot.2.male)  #create data frame to export
head(exp.male.resid2)# check first three rows of data frame

write.table(exp.male.resid2, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/nonweighted_male_resid_Aug3.csv", sep = ",", row.names=F) #Export


#######  Length Bootstrap #######
dat$Avg.Size<-dat$Avg.Length
#dat$Avg.Length<-read.table()
adult.mean<-mean(dat$Avg.Length[dat$Age > 20])
adult.mean.female<-mean(dat$Avg.Length[dat$Age > 20 & dat$Sex=='0'])
adult.mean.male<-mean(dat$Avg.Length[dat$Age > 20 & dat$Sex=='1'])

####### Weighted Length Bootstrap #######

mod.new<- optim(par=c(B=adult.mean, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age, y=dat$Avg.Length,method="SANN",control=list(maxit=10000,trace=T));mod.new
mod.old<-mod.new
plot(Avg.Length~Age,data=dat,main='Male and Female Weighted Average Body Length')
curve(mod.new$par[1]*(1-1/(1 +(x/mod.new$par[2])^mod.new$par[5]+(x/mod.new$par[3])^mod.new$par[6]+(x/mod.new$par[4])^mod.new$par[7])), add=TRUE,lwd=2,lty=2,col='black', ylab='Body Length')

w<-rep(1,length(dat$Age))

tol=1e-9
maxit=500
iter=1000

mod.l<-matrix(nrow=iter,ncol=7)
res.tot.l<-matrix(nrow=iter,ncol=length(dat$Age))

for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new<- optim(par=c(B=adult.mean, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age, y=dat$Avg.Length,method="SANN",control=list(maxit=10000,trace=F))
  for (i in 1:maxit){
    oldcof<-as.vector(mod.new$par)
    B<-oldcof[1]
    D1<-oldcof[2]
    D2<-oldcof[3]
    D3<-oldcof[4]
    C1<-oldcof[5]
    C2<-oldcof[6]
    C3<-oldcof[7]
    yhat<-B*(1 - 1/(1 +((dat$Age/D1)^C1)+((dat$Age/D2)^C2)+((dat$Age/D3)^C3)))
    res<-dat$Avg.Length-yhat
    mod.res<-optim(par=c(.5,.5),reg.func,a=dat$Avg.Length,b=res,control=list(abstol=tol),method="Nelder")
    gam0<-mod.res$par[1]
    gam1<-mod.res$par[2]
    w.new<-abs(1/(gam0 + gam1 * dat$Avg.Length)) # estimate new Length  
    mod.new<-optim(par=as.vector(mod.new$par),weight.func,x=dat$Age,y=dat$Avg.Length,w=w.new,control=list(abstol=tol),method="SANN")
    newcof<-as.vector(mod.new$par)
    dif<-sum(sqrt((oldcof-newcof)^2))/length(newcof)
    dif<-sqrt((oldcof[1]-newcof[1])^2) # only a single parameter
    if(dif<=tol){break}  
  }
  mod.l[j,]<-mod.new$par
  res.tot.l[j,]<-res
}

mod.l<-matrix(mod.l[complete.cases(mod.l)],ncol=7);mod.l
mod.mean.l<-apply(mod.l, 2, median, na.rm=T);mod.mean.l
mod.sd.l<-apply(mod.l, 2, sd, na.rm=T);mod.sd.l

plot(dat$Avg.Length~dat$Age, main='Male and Female Iteratively Weighted for Length',pch=19,col=c('pink','blue')[dat$Sex],xlab='Age',ylab='Avg.Length'); legend(x=60,y=50, legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2))
curve(mod.mean.l[1]*(1-1/(1 +(x/mod.mean.l[2])^mod.mean.l[5]+(x/mod.mean.l[3])^mod.mean.l[6]+(x/mod.mean.l[4])^mod.mean.l[7])), add=TRUE,lwd=1,lty=1,col='black')
#NFB testing to put 95% CI on graph
a<-predict(mod.l, interval="confidence")
lines(dat$Avg.Length, a[,2], lty=2) #Upper limit
lines(dat$Avg.Length, a[,3], lty=2) #Lower limit

curve(mod.old$par[1]*(1-1/(1 +(x/mod.old$par[2])^mod.old$par[5]+(x/mod.old$par[3])^mod.old$par[6]+(x/mod.old$par[4])^mod.old$par[7])), add=TRUE,lwd=1,lty=2,col='black')

curve(mod.mean.l[1] * (((x/mod.mean.l[2])^(mod.mean.l[5] - 1) * (mod.mean.l[5] * (1/mod.mean.l[2])) + (x/mod.mean.l[3])^(mod.mean.l[6] - 1) * (mod.mean.l[6] * (1/mod.mean.l[3])) + (x/mod.mean.l[4])^(mod.mean.l[7] - 1) * (mod.mean.l[7] * (1/mod.mean.l[4])))/(1 + ((x/mod.mean.l[2])^mod.mean.l[5]) + ((x/mod.mean.l[3])^mod.mean.l[6]) + ((x/mod.mean.l[4])^mod.mean.l[7]))^2),lty=1, xlim=c(.5,25),main='Male and Female Iteratively Weighted',xlab='Age',ylab='Growth Velocity')

curve(mod.old$par[1] * (((x/mod.old$par[2])^(mod.old$par[5] - 1) * (mod.old$par[5] * (1/mod.old$par[2])) + (x/mod.old$par[3])^(mod.old$par[6] - 1) * (mod.old$par[6] * (1/mod.old$par[3])) + (x/mod.old$par[4])^(mod.old$par[7] - 1) * (mod.old$par[7] * (1/mod.old$par[4])))/(1 + ((x/mod.old$par[2])^mod.old$par[5]) + ((x/mod.old$par[3])^mod.old$par[6]) + ((x/mod.old$par[4])^mod.old$par[7]))^2), add=T, lty=2);legend(x=30,y=10,legend=c('Bootstrap','Original'),lty=c(1,2))



#Export
exp.weighted.data<-as.data.frame(mod.l)  #create data frame to export
colnames(exp.weighted.data)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.weighted.data)# check first three rows of data frame

write.table(exp.weighted.data, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/weightedlengthbootstrap_combined_LENGTH_july6.csv", sep = ",", row.names=F) #Export

exp.weighted.resid<-as.data.frame(res.tot.l)  #create data frame to export
head(exp.weighted.resid)# check first three rows of data frame

write.table(exp.weighted.resid, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/weightedlengthresid_combined_LENGTH_july6.csv", sep = ",", row.names=F) #Export

###### Non-weighted Length Bootstrap #######
tol=1e-9
maxit=500
iter=1000

mod.2.l<-matrix(nrow=iter,ncol=7)
res.tot.2.l<-matrix(nrow=iter,ncol=length(dat$Age))

for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new<- optim(par=c(B=adult.mean, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age, y=dat$Avg.Length,method="SANN",control=list(maxit=10000,trace=F))
  oldcof<-as.vector(mod.new$par)
  B<-oldcof[1]
  D1<-oldcof[2]
  D2<-oldcof[3]
  D3<-oldcof[4]
  C1<-oldcof[5]
  C2<-oldcof[6]
  C3<-oldcof[7]
  yhat<-B*(1 - 1/(1 +((dat$Age/D1)^C1)+((dat$Age/D2)^C2)+((dat$Age/D3)^C3)))
  res<-dat$Avg.Length-yhat
  mod.2.l[j,]<-mod.new$par
  res.tot.2.l[j,]<-res
}

mod.2.l<-matrix(mod.2.l[complete.cases(mod.2.l)],ncol=7);mod.2.l
mod.2.mean.l<-apply( mod.2.l, 2, median, na.rm=T);mod.2.mean.l
mod.2.sd.l<-apply( mod.2.l, 2, sd, na.rm=T);mod.2.sd.l

plot(dat$Avg.Length~dat$Age, main='Male and Female Non-weighted Length',pch=19,col=c('pink','blue')[dat$Sex],xlab='Age',ylab='Avg.Length'); legend(x=60, y=50, legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2),lwd=c(NA,NA,1,2))
curve(mod.2.mean.l[1]*(1-1/(1 +(x/mod.2.mean.l[2])^mod.2.mean.l[5]+(x/mod.2.mean.l[3])^mod.2.mean.l[6]+(x/mod.2.mean.l[4])^mod.2.mean.l[7])), add=TRUE,lwd=1,lty=1)

curve(mod.old$par[1]*(1-1/(1 +(x/mod.old$par[2])^mod.old$par[5]+(x/mod.old$par[3])^mod.old$par[6]+(x/mod.old$par[4])^mod.old$par[7])), add=TRUE,lwd=2,lty=2,col='black')

curve(mod.2.mean.l[1] * (((x/mod.2.mean.l[2])^(mod.2.mean.l[5] - 1) * (mod.2.mean.l[5] * (1/mod.2.mean.l[2])) + (x/mod.2.mean.l[3])^(mod.2.mean.l[6] - 1) * (mod.2.mean.l[6] * (1/mod.2.mean.l[3])) + (x/mod.2.mean.l[4])^(mod.2.mean.l[7] - 1) * (mod.2.mean.l[7] * (1/mod.2.mean.l[4])))/(1 + ((x/mod.2.mean.l[2])^mod.2.mean.l[5]) + ((x/mod.2.mean.l[3])^mod.2.mean.l[6]) + ((x/mod.2.mean.l[4])^mod.2.mean.l[7]))^2), xlim=c(.5,25),main='Male and Female Non-Weighted Length Curve',xlab='Age',ylab='Growth Velocity')

curve(mod.old$par[1] * (((x/mod.old$par[2])^(mod.old$par[5] - 1) * (mod.old$par[5] * (1/mod.old$par[2])) + (x/mod.old$par[3])^(mod.old$par[6] - 1) * (mod.old$par[6] * (1/mod.old$par[3])) + (x/mod.old$par[4])^(mod.old$par[7] - 1) * (mod.old$par[7] * (1/mod.old$par[4])))/(1 + ((x/mod.old$par[2])^mod.old$par[5]) + ((x/mod.old$par[3])^mod.old$par[6]) + ((x/mod.old$par[4])^mod.old$par[7]))^2), add=T,lwd=1,lty=2);legend(x=15,y=120,legend=c('Bootstrap','Original'),lty=c(1,2))

#Export
exp.data<-as.data.frame(mod.2.l)  #create data frame to export
colnames(exp.data)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.data)# check first three rows of data frame

write.table(exp.data, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/non-weighted_bootstrap_length_aug10.csv", sep = ",", row.names=F) #Export

exp.resid<-as.data.frame(res.tot.2.l)  #create data frame to export
head(exp.resid)# check first three rows of data frame

write.table(exp.resid, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve//nonweighted_resid_length_july6.csv", sep = ",", row.names=F) #Export


####### Female Iteratively Weighted Length Bootstrap #######

mod.new.female<- optim(par=c(B=adult.mean.female, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age[which(dat$Sex=='0')], y=dat$Avg.Length[which(dat$Sex=='0')], method="SANN",control=list(maxit=10000,trace=T));mod.new.female 
plot(Avg.Length~Age,data=dat[which(dat$Sex=='0'),],main='Average Length for Female',col='red')
curve(mod.new.female$par[1]*(1-1/(1 +(x/mod.new.female$par[2])^mod.new.female$par[5]+(x/mod.new.female$par[3])^mod.new.female$par[6]+(x/mod.new.female$par[4])^mod.new.female$par[7])), add=TRUE,lwd=2,lty=2)
mod.old.female<-mod.new.female

w<-rep(1,length(dat$Age[which(dat$Sex=='0')]))

tol=1e-9
maxit=500
iter=1000

mod.female.l<-matrix(nrow=iter,ncol=7)
res.tot.female.l<-matrix(nrow=iter,ncol=length(dat$Age[which(dat$Sex=='0')]))

for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new.female<- optim(par=c(B=adult.mean.female, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age[which(dat$Sex=='0')], y=dat$Avg.Length[which(dat$Sex=='0')], method="SANN",control=list(maxit=10000,trace=F));mod.new.female 
  for (i in 1:maxit){
    oldcof<-as.vector(mod.new.female$par)
    B<-oldcof[1]
    D1<-oldcof[2]
    D2<-oldcof[3]
    D3<-oldcof[4]
    C1<-oldcof[5]
    C2<-oldcof[6]
    C3<-oldcof[7]
    fem.age<-dat$Age[which(dat$Sex=='0')]
    fem.size<-dat$Avg.Length[which(dat$Sex=='0')]
    yhat<-B*(1 - 1/(1 +((fem.age/D1)^C1)+((fem.age/D2)^C2)+((fem.age/D3)^C3)))
    res<- fem.size-yhat
    mod.res<-optim(par=c(.5,.5),reg.func,a=fem.size,b=res,control=list(abstol=tol),method="Nelder")
    gam0<-mod.res$par[1]
    gam1<-mod.res$par[2]
    w.new<-abs(1/(gam0 + gam1 * fem.size)) # estimate new weights  
    mod.new<-optim(par=as.vector(mod.new.female$par),weight.func,x=fem.age,y=fem.size,w=w.new,control=list(abstol=tol),method="SANN")
    newcof<-as.vector(mod.new.female$par)
    dif<-sum(sqrt((oldcof-newcof)^2))/length(newcof)
    dif<-sqrt((oldcof[1]-newcof[1])^2) # only a single parameter
    if(dif<=tol){break}  
  }
  mod.female.l[j,]<-mod.new.female$par
  res.tot.female.l[j,]<-res
}

mod.female.l<-matrix(mod.female.l[complete.cases(mod.female.l)],ncol=7);mod.female.l
mod.female.mean.l<-apply(mod.female.l, 2, median, na.rm=T);mod.female.mean.l
mod.female.sd.l<-apply( mod.female.l, 2, sd, na.rm=T);mod.female.sd.l

plot(dat$Avg.Length~dat$Age, main='Iteratively Weighted for length for Females',pch=19,col=c('pink','blue')[dat$Sex],xlab='Age',ylab='Avg.Length'); legend(x=60,y=50, legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2),lwd=c(NA,NA,1,2))
 
curve(mod.old.female$par[1]*(1-1/(1 +(x/mod.old.female$par[2])^mod.old.female$par[5]+(x/mod.old.female$par[3])^mod.old.female$par[6]+(x/mod.old.female$par[4])^mod.old.female$par[7])), add=TRUE,lwd=2,lty=2,col='pink')

curve(mod.female.mean.l[1]*(1-1/(1 +(x/mod.female.mean.l[2])^mod.female.mean.l[5]+(x/mod.female.mean.l[3])^mod.female.mean.l[6]+(x/mod.female.mean.l[4])^mod.female.mean.l[7])), add=TRUE,lwd=1,lty=1,col='pink')

curve(mod.female.mean.l[1] * (((x/mod.female.mean.l[2])^(mod.female.mean.l[5] - 1) * (mod.female.mean.l[5] * (1/mod.female.mean.l[2])) + (x/mod.female.mean.l[3])^(mod.female.mean.l[6] - 1) * (mod.female.mean.l[6] * (1/mod.female.mean.l[3])) + (x/mod.female.mean.l[4])^(mod.female.mean.l[7] - 1) * (mod.female.mean.l[7] * (1/mod.female.mean.l[4])))/(1 + ((x/mod.female.mean.l[2])^mod.female.mean.l[5]) + ((x/mod.female.mean.l[3])^mod.female.mean.l[6]) + ((x/mod.female.mean.l[4])^mod.female.mean.l[7]))^2), xlim=c(0,40),,col='pink',ylab='Growth Velocity',xlab='Age',main='Iteratively Weighted Growth Curve for females');legend(x=50,y=20,legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2),lwd=c(NA,NA,1,2))

curve(mod.old.female$par[1] * (((x/mod.old.female$par[2])^(mod.old.female$par[5] - 1) * (mod.old.female$par[5] * (1/mod.old.female$par[2])) + (x/mod.old.female$par[3])^(mod.old.female$par[6] - 1) * (mod.old.female$par[6] * (1/mod.old.female$par[3])) + (x/mod.old.female$par[4])^(mod.old.female$par[7] - 1) * (mod.old.female$par[7] * (1/mod.old.female$par[4])))/(1 + ((x/mod.old.female$par[2])^mod.old.female$par[5]) + ((x/mod.old.female$par[3])^mod.old.female$par[6]) + ((x/mod.old.female$par[4])^mod.old.female$par[7]))^2), lty=2,lwd=2,col='pink',add=T)

#Export
exp.female.data<-as.data.frame(mod.female.l)  #create data frame to export
colnames(exp.female.data)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.female.data)# check first three rows of data frame

write.table(exp.female.data, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve//female_weighted_length_bootstrap_aug10.csv", sep = ",", row.names=F) #Export

exp.female.resid<-as.data.frame(res.tot.female.l)  #create data frame to export
head(exp.female.resid)# check first three rows of data frame

write.table(exp.female.resid, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/female_weighted_resid_length_aug10.csv", sep = ",", row.names=F) #Export

###### Female Non-weighted Length Bootstrap #######
tol=1e-9
maxit=500
iter=1000

mod.female.2.l<-matrix(nrow=iter,ncol=7)
res.tot.2.female.l<-matrix(nrow=iter,ncol=length(dat$Age[which(dat$Sex=='0')]))

for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new.female<- optim(par=c(B=adult.mean.female, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age[which(dat$Sex=='0')], y=dat$Avg.Length[which(dat$Sex=='0')], method="SANN",control=list(maxit=10000,trace=F));mod.new.female 
  oldcof<-as.vector(mod.new.female$par)
  B<-oldcof[1]
  D1<-oldcof[2]
  D2<-oldcof[3]
  D3<-oldcof[4]
  C1<-oldcof[5]
  C2<-oldcof[6]
  C3<-oldcof[7]
  fem.age<-dat$Age[which(dat$Sex=='0')]
  fem.size<-dat$Avg.Length[which(dat$Sex=='0')]
  yhat<-B*(1 - 1/(1 +((fem.age/D1)^C1)+((fem.age/D2)^C2)+((fem.age/D3)^C3)))
  res<-fem.size-yhat
  mod.female.2.l[j,]<-mod.new.female$par
  res.tot.2.female.l[j,]<-res
}

mod.female.2.l<-matrix(mod.female.2.l[complete.cases(mod.female.2.l)],ncol=7);mod.female.2.l
mod.female.2.mean.l<-apply( mod.female.2.l, 2, median, na.rm=T);mod.female.2.mean.l
mod.female.2.sd.l<-apply( mod.female.2.l, 2, sd, na.rm=T);mod.female.2.sd.l

plot(dat$Avg.Length~dat$Age, main='Non-Weighted Female Length',pch=19,col=c('pink','blue')[dat$Sex],xlab='Age',ylab='Avg.Length'); legend(x=40,y=50, legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2),lwd=c(NA,NA,1,2))
curve(mod.female.2.mean.l[1]*(1-1/(1 +(x/mod.female.2.mean.l[2])^mod.female.2.mean.l[5]+(x/mod.female.2.mean.l[3])^mod.female.2.mean.l[6]+(x/mod.female.2.mean.l[4])^mod.female.2.mean.l[7])), add=TRUE,lwd=1,lty=1,col='pink')

curve(mod.old.female$par[1]*(1-1/(1 +(x/mod.old.female$par[2])^mod.old.female$par[5]+(x/mod.old.female$par[3])^mod.old.female$par[6]+(x/mod.old.female$par[4])^mod.old.female$par[7])), add=TRUE,lwd=2,lty=2,col='pink')

curve(mod.female.2.mean.l[1] * (((x/mod.female.2.mean.l[2])^(mod.female.2.mean.l[5] - 1) * (mod.female.2.mean.l[5] * (1/mod.female.2.mean.l[2])) + (x/mod.female.2.mean.l[3])^(mod.female.2.mean.l[6] - 1) * (mod.female.2.mean.l[6] * (1/mod.female.2.mean.l[3])) + (x/mod.female.2.mean.l[4])^(mod.female.2.mean.l[7] - 1) * (mod.female.2.mean.l[7] * (1/mod.female.2.mean.l[4])))/(1 + ((x/mod.female.2.mean.l[2])^mod.female.2.mean.l[5]) + ((x/mod.female.2.mean.l[3])^mod.female.2.mean.l[6]) + ((x/mod.female.2.mean.l[4])^mod.female.2.mean.l[7]))^2),col='pink', xlim=c(.5,25),ylim=c(0,150),ylab='Growth Velocity',xlab='Age',main='Non-weighted Growth Curve');legend(x=15,y=120,legend=c('Female','Male','Bootstrap','Original'),col=c('pink','blue','black','black'),pch=c(19,19,NA,NA),lty=c(NA,NA,1,2),lwd=c(NA,NA,1,2))

curve(mod.old.female$par[1] * (((x/mod.old.female$par[2])^(mod.old.female$par[5] - 1) * (mod.old.female$par[5] * (1/mod.old.female$par[2])) + (x/mod.old.female$par[3])^(mod.old.female$par[6] - 1) * (mod.old.female$par[6] * (1/mod.old.female$par[3])) + (x/mod.old.female$par[4])^(mod.old.female$par[7] - 1) * (mod.old.female$par[7] * (1/mod.old.female$par[4])))/(1 + ((x/mod.old.female$par[2])^mod.old.female$par[5]) + ((x/mod.old.female$par[3])^mod.old.female$par[6]) + ((x/mod.old.female$par[4])^mod.old.female$par[7]))^2), col='pink',lwd=2,lty=2,add=T)

#Export
exp.female.data2<-as.data.frame(mod.female.2.l)  #create data frame to export
colnames(exp.female.data2)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.female.data2)# check first three rows of data frame

write.table(exp.female.data2, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/non-weighted_female_length_bootstrap_aug10.csv", sep = ",", row.names=F) #Export

exp.female.resid2<-as.data.frame(res.tot.2.female.l)  #create data frame to export
head(exp.female.resid2)# check first three rows of data frame

write.table(exp.female.resid2, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/nonweighted_female_resid_length_aug10.csv", sep = ",", row.names=F) #Export


####### Male Iteratively Weighted Length Bootstrap #######

mod.new.male <-optim(par=c(B=adult.mean.male, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age[which(dat$Sex=='1')], y=dat$Avg.Length[which(dat$Sex=='1')], method="SANN",control=list(maxit=10000,trace=T));mod.new.male
plot(Avg.Length~Age,data=dat[which(dat$Sex=='1'),],main='Male Original Body Length',col='blue')
curve(mod.new.male$par[1]*(1-1/(1 +(x/mod.new.male$par[2])^mod.new.male$par[5]+(x/mod.new.male$par[3])^mod.new.male$par[6]+(x/mod.new.male$par[4])^mod.new.male$par[7])), add=TRUE,lwd=2,lty=2)
mod.old.male<-mod.new.male

w<-rep(1,length(dat$Age[which(dat$Sex=='1')]))

tol=1e-9
maxit=500
iter=1000

mod.male.l<-matrix(nrow=iter,ncol=7)
res.tot.male.l<-matrix(nrow=iter,ncol=length(dat$Age[which(dat$Sex=='1')]))

for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new.male<- optim(par=c(B=adult.mean.male, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age[which(dat$Sex=='1')], y=dat$Avg.Length[which(dat$Sex=='1')], method="SANN",control=list(maxit=10000,trace=F));mod.new.male 
  for (i in 1:maxit){
    oldcof<-as.vector(mod.new.male$par)
    B<-oldcof[1]
    D1<-oldcof[2]
    D2<-oldcof[3]
    D3<-oldcof[4]
    C1<-oldcof[5]
    C2<-oldcof[6]
    C3<-oldcof[7]
    male.age<-dat$Age[which(dat$Sex=='1')]
    male.size<-dat$Avg.Length[which(dat$Sex=='1')]
    yhat<-B*(1 - 1/(1 +((male.age/D1)^C1)+((male.age/D2)^C2)+((male.age/D3)^C3)))
    res<- male.size-yhat
    mod.res<-optim(par=c(.5,.5),reg.func,a=male.size,b=res,control=list(abstol=tol),method="Nelder")
    gam0<-mod.res$par[1]
    gam1<-mod.res$par[2]
    w.new<-abs(1/(gam0 + gam1 * male.size)) # estimate new lengths  
    mod.new<-optim(par=as.vector(mod.new.male$par),weight.func,x=male.age,y=male.size,w=w.new,control=list(abstol=tol),method="SANN")
    newcof<-as.vector(mod.new.male$par)
    dif<-sum(sqrt((oldcof-newcof)^2))/length(newcof)
    dif<-sqrt((oldcof[1]-newcof[1])^2) # only a single parameter
    if(dif<=tol){break}  
  }
  mod.male.l[j,]<-mod.new.male$par
  res.tot.male.l[j,]<-res
}

mod.male.l<-matrix(mod.male.l[complete.cases(mod.male.l)],ncol=7);mod.male.l
mod.male.mean.l<-apply(mod.male.l, 2, median, na.rm=T);mod.male.mean.l
mod.male.sd.l<-apply( mod.male.l, 2, sd, na.rm=T);mod.male.sd.l
apply(mod.male.l, MARGIN = 2, FUN = quantile, probs = .975)
apply(mod.male.l, MARGIN = 2, FUN = quantile, probs = .025)


plot(Avg.Length~Age,data=dat[which(dat$Sex=='1'),],main='Weighted Bootstrap Length for Male',col='blue')

curve(mod.old.male$par[1]*(1-1/(1 +(x/mod.old.male$par[2])^mod.old.male$par[5]+(x/mod.old.male$par[3])^mod.old.male$par[6]+(x/mod.old.male$par[4])^mod.old.male$par[7])), add=TRUE,lwd=2,lty=2,col='blue')

curve(mod.male.mean.l[1]*(1-1/(1 +(x/mod.male.mean.l[2])^mod.male.mean.l[5]+(x/mod.male.mean.l[3])^mod.male.mean.l[6]+(x/mod.male.mean.l[4])^mod.male.mean.l[7])), add=TRUE,lwd=1,lty=1,col='blue')

curve(mod.male.mean.l[1] * (((x/mod.male.mean.l[2])^(mod.male.mean.l[5] - 1) * (mod.male.mean.l[5] * (1/mod.male.mean.l[2])) + (x/mod.male.mean.l[3])^(mod.male.mean.l[6] - 1) * (mod.male.mean.l[6] * (1/mod.male.mean.l[3])) + (x/mod.male.mean.l[4])^(mod.male.mean.l[7] - 1) * (mod.male.mean.l[7] * (1/mod.male.mean.l[4])))/(1 + ((x/mod.male.mean.l[2])^mod.male.mean.l[5]) + ((x/mod.male.mean.l[3])^mod.male.mean.l[6]) + ((x/mod.male.mean.l[4])^mod.male.mean.l[7]))^2),
      col='blue',add=F, xlim=c(.5,25),ylim=c(0,20),ylab='Growth Velocity',xlab='Age',main='Weighted Male Length Growth Velocity Curve')

curve(mod.old.male$par[1] * (((x/mod.old.male$par[2])^(mod.old.male$par[5] - 1) * (mod.old.male$par[5] * (1/mod.old.male$par[2])) + (x/mod.old.male$par[3])^(mod.old.male$par[6] - 1) * (mod.old.male$par[6] * (1/mod.old.male$par[3])) + (x/mod.old.male$par[4])^(mod.old.male$par[7] - 1) * (mod.old.male$par[7] * (1/mod.old.male$par[4])))/(1 + ((x/mod.old.male$par[2])^mod.old.male$par[5]) + ((x/mod.old.male$par[3])^mod.old.male$par[6]) + ((x/mod.old.male$par[4])^mod.old.male$par[7]))^2), lwd=2,lty=2,col='blue',add=T)

#Export
exp.male.data<-as.data.frame(mod.male.l)  #create data frame to export
colnames(exp.male.data)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.male.data)# check first three rows of data frame

write.table(exp.male.data, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/weighted_male_bootstrap_LENGTH_aug10.csv", sep = ",", row.names=F) #Export

exp.male.resid<-as.data.frame(res.tot.male.l)  #create data frame to export
head(exp.male.resid)# check first three rows of data frame

write.table(exp.male.resid, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/weighted_male_resid_LENGTH_aug10.csv", sep = ",", row.names=F) #Export

###### Male Non-weighted Length Bootstrap #######
tol=1e-9
maxit=500
iter=1000

mod.male.2.l<-matrix(nrow=iter,ncol=7)
res.tot.2.male.l<-matrix(nrow=iter,ncol=length(dat$Age[which(dat$Sex=='1')]))

for (j in 1:iter){
  cat(j, "iterations","\n")
  mod.new.male<- optim(par=c(B=adult.mean.male, D1= 8.48 , D2=8.93, D3=44.81, C1= 1.92, C2= 6.95, C3=0.79), fn=func, x=dat$Age[which(dat$Sex=='1')], y=dat$Avg.Length[which(dat$Sex=='1')], method="SANN",control=list(maxit=10000,trace=F));mod.new.male 
  oldcof<-as.vector(mod.new.male$par)
  B<-oldcof[1]
  D1<-oldcof[2]
  D2<-oldcof[3]
  D3<-oldcof[4]
  C1<-oldcof[5]
  C2<-oldcof[6]
  C3<-oldcof[7]
  male.age<-dat$Age[which(dat$Sex=='1')]
  male.size<-dat$Avg.Length[which(dat$Sex=='1')]
  yhat<-B*(1 - 1/(1 +((male.age/D1)^C1)+((male.age/D2)^C2)+((male.age/D3)^C3)))
  res<-male.size-yhat
  mod.male.2.l[j,]<-mod.new.male$par
  res.tot.2.male.l[j,]<-res
}

mod.male.2.l<-matrix(mod.male.2.l[complete.cases(mod.male.2.l)],ncol=7);mod.male.2.l
mod.male.2.mean.l<-apply( mod.male.2.l, 2, median, na.rm=T);mod.male.2.mean.l
mod.male.2.sd.l<-apply( mod.male.2.l, 2, sd, na.rm=T);mod.male.2.sd.l

plot(Avg.Length~Age,data=dat[which(dat$Sex=='1'),],main='Average Body Length for Males',col='blue')
curve(mod.male.2.mean.l[1]*(1-1/(1 +(x/mod.male.2.mean.l[2])^mod.male.2.mean.l[5]+(x/mod.male.2.mean.l[3])^mod.male.2.mean.l[6]+(x/mod.male.2.mean.l[4])^mod.male.2.mean.l[7])), add=TRUE,lwd=1,lty=1,col='blue')

curve(mod.old.male$par[1]*(1-1/(1 +(x/mod.old.male$par[2])^mod.old.male$par[5]+(x/mod.old.male$par[3])^mod.old.male$par[6]+(x/mod.old.male$par[4])^mod.old.male$par[7])), add=TRUE,lwd=2,lty=2,col='blue')

curve((mod.male.2.mean.l[1] * (((x/mod.male.2.mean.l[2])^(mod.male.2.mean.l[5] - 1) * (mod.male.2.mean.l[5] * (1/mod.male.2.mean.l[2])) + (x/mod.male.2.mean.l[3])^(mod.male.2.mean.l[6] - 1) * (mod.male.2.mean.l[6] * (1/mod.male.2.mean.l[3])) + (x/mod.male.2.mean.l[4])^(mod.male.2.mean.l[7] - 1) * (mod.male.2.mean.l[7] * (1/mod.male.2.mean.l[4])))/(1 + ((x/mod.male.2.mean.l[2])^mod.male.2.mean.l[5]) + ((x/mod.male.2.mean.l[3])^mod.male.2.mean.l[6]) + ((x/mod.male.2.mean.l[4])^mod.male.2.mean.l[7]))^2)),
      col='blue',add=F, from=0, to=30, xlim=c(.5,30),ylim=c(0,20),ylab='Growth Velocity',xlab='Age',main='Non-Weighted Male Length Growth Velocity Curve')

curve(mod.old.male$par[1] * (((x/mod.old.male$par[2])^(mod.old.male$par[5] - 1) * (mod.old.male$par[5] * (1/mod.old.male$par[2])) + (x/mod.old.male$par[3])^(mod.old.male$par[6] - 1) * (mod.old.male$par[6] * (1/mod.old.male$par[3])) + (x/mod.old.male$par[4])^(mod.old.male$par[7] - 1) * (mod.old.male$par[7] * (1/mod.old.male$par[4])))/(1 + ((x/mod.old.male$par[2])^mod.old.male$par[5]) + ((x/mod.old.male$par[3])^mod.old.male$par[6]) + ((x/mod.old.male$par[4])^mod.old.male$par[7]))^2), lwd=2,lty=2,col='green',add=T)

#Export
exp.male.data2<-as.data.frame(mod.male.2.l)  #create data frame to export
colnames(exp.male.data2)<-c('B','D1','D2','D3','C1','C2','C3') #name columns of data frame
head(exp.male.data2)# check first three rows of data frame

write.table(exp.male.data2, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/nonweighted_male_bootstrap_length_aug10.csv", sep = ",", row.names=F) #Export

exp.male.resid2<-as.data.frame(res.tot.2.male.l)  #create data frame to export
head(exp.male.resid2)# check first three rows of data frame

write.table(exp.male.resid2, file = "/Users/nbrazeau/Dropbox/ZPM&NFB Folder, July 2015/JPPS_Curve/nonweighted_male_resid_length_agu10.csv", sep = ",", row.names=F) #Export