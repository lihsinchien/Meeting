
##################################
### fit model ####################
##################################
rm(list=ls())

## library, function ###############
#install.packages("survival")
library(survival) #age-matched case control, clogit()

f.mu2<-function(x,assir.age,xb.b1){ #x=mu30, assir.age=P30, xb.b1=beta2*x2+beta3*x3
  diff<-sum(exp(xb.b1+x)/(1+exp(xb.b1+x)))-assir.age*length(xb.b1) 
  #	exp(xb.b1+x)/(1+exp(xb.b1+x) 估計每個年齡的母體發病率
  return(diff)
}

## load data
load("data_mcc.RData")

##fit age-matched case control data
# stage 1

r <- clogit(case~x1+x2+x3+x4+x5+strata(group), data=data.s)
summary(r)

b1 <- r$coefficients
#b1.cov <- r$var

# stage 2

x.s <- as.matrix(data.s[,7:11]) 

tmp.index1 <- which(data.s$case==0) #抓control的位置
xb.b1.tmp1 <- x.s[tmp.index1,]%*%b1 #control對應到的xbeta的值
assir <- pop_data[,3]/pop_data[,2] #各年齡層的發生率

age <- (30:76)/10
mu.1 <- mu.2 <- NULL #mu.1: all control, mu.2: age-specific control
i=1
for(i in 1:length(age)){ #估每個年齡層的mu
  mu.tmp1 <- uniroot(f.mu2, assir.age=assir[i], xb.b1=xb.b1.tmp1, interval=c(-10,10))$root
  mu.1 <- c(mu.1, mu.tmp1)
  
  tmp.index2 <- intersect(which(data.s$case==0), which(data.s$age==age[i]))
  if(length(tmp.index2) >= 10){ #總共多少人(樣本數夠大則符合條件)
    xb.b1.tmp2 <- x.s[tmp.index2,]%*%b1 #算出xbeta(=beta2*x2+beta3*x3)
    mu.tmp2 <- uniroot(f.mu2, assir.age=assir[i], xb.b1=xb.b1.tmp2, interval=c(-10,10))$root
    mu.2 <- c(mu.2, mu.tmp2)
  }else{
    mu.2 <- c(mu.2,NA)
  }
}

#cbind(age, mu.1, mu.2) 

age.2 <- age^2
age.3 <- age^3
coe.age1 <- lm(mu.1~age+age.2+age.3)$coefficients
coe.age2 <- lm(mu.2~age+age.2+age.3)$coefficients

coe.out<-c(coe.age1,coe.age2,b1)

## bootstrap #################

n<-1000 # group size
n.boot<-20 # number of bootstrap sample

coe.boot<-NULL

i.boot=1
tic.boot<-Sys.time() #start time

for(i.boot in 1:n.boot){
  g.boot<-sample(1:n,n,replace=T)
  data.boot<-rbind(data.s[6*(g.boot-1)+1,],data.s[6*(g.boot-1)+2,],data.s[6*(g.boot-1)+3,],data.s[6*(g.boot-1)+4,],data.s[6*(g.boot-1)+5,],data.s[6*(g.boot),])
  
  r.boot<-clogit(case~x1+x2+x3+x4+x5+strata(group),data=data.boot)
  b1.boot<-r.boot$coefficients
  
  ##### stage 2
  x.boot<-as.matrix(data.boot[,7:11])
  
  tmp.index1<-which(data.boot$case==0)
  xb.b1.tmp1.boot<-x.boot[tmp.index1,]%*%b1.boot
  
  age<-30:76/10
  mu.1.boot<-mu.2.boot<-NULL
  i=1
  for(i in 1:length(age)){
    mu.tmp1<-uniroot(f.mu2,assir.age=assir[i],xb.b1=xb.b1.tmp1.boot,interval=c(-10,10))$root
    mu.1.boot<-c(mu.1.boot,mu.tmp1)
    
    tmp.index2<-intersect(which(data.boot$case==0),which(data.boot$age==age[i]))
    if(length(tmp.index2)>=10){
      xb.b1.tmp2.boot<-x.boot[tmp.index2,]%*%b1.boot
      mu.tmp2<-uniroot(f.mu2,assir.age=assir[i],xb.b1=xb.b1.tmp2.boot,interval=c(-10,10))$root
      mu.2.boot<-c(mu.2.boot,mu.tmp2)
    }else{
      mu.2.boot<-c(mu.2.boot,NA)
    }
  }
  tmp<-age
  tmp2<-tmp^2
  tmp3<-tmp^3
  
  coe.age1.boot<-lm(mu.1.boot~tmp+tmp2+tmp3)$coefficients
  coe.age2.boot<-lm(mu.2.boot~tmp+tmp2+tmp3)$coefficients
  coe.out.boot<-c(coe.age1.boot,coe.age2.boot,b1.boot)
  
  coe.boot<-rbind(coe.boot,coe.out.boot)
#  if(i.boot%%1000==0){
#    print(i.boot)
#    print(Sys.time()-tic.boot)
#  }
}

coe.out.all<-rbind(coe.out,coe.boot)
colnames(coe.out.all)[1:8]<-c("const_all","age_all","age2_all","age3_all","const_sub","age_sub","age2_sub","age3_sub")



write.table(coe.out.all,file="boot_test.txt",row.names=F)

  
  
  
  
  
  



