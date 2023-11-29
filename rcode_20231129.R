
setwd("D:/cycu/meeting/CJH/20231129")

save(data.s,pop_data,coe.true, file="data_mcc.RData")

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
ls()
dim(data.s) #1000 matched group, 1:5 match
head(data.s)


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

table.out<-rbind(c(coe.age1,b1),c(coe.age2,b1))
round(table.out,digits=2)

