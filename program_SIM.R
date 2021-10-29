# setwd("C:/Users/Cati/Desktop/Ricardo/Cati_Montse/datos")



install.packages("RGCCA")
library(RGCCA)

install.packages("e1071")
library(e1071)

install.packages("MultiSkew")
library(MultiSkew)

install.packages("robustbase")
library(robustbase)
# Reading data
data<-read.table("data_cost.csv",header=TRUE,sep=",")

# data<-data[(data$total_cos<100000),] #To eliminate extreme

cost<-data$total_cost/data$Nclaims
hist(cost)
plot(density(cost))
# Explanatory variables
mix1<-data$Nclaims
mix2<-data$age 
mix3<-data$gender
mix4<-data$agelic    
mix5<-data$agecar 
mix6<-data$parking 
mix7<-data$tkm/1000
mix8<-data$nightkm
mix8<-(mix8/data$tkm)*100
mix9<-data$urbankm 
mix10<-data$speedkm


# Dependent variable
cost<-cost/1000
Log_cost<-log(cost)
hist(Log_cost)
data<-as.matrix(Log_cost)
n<-nrow(data)
summary(data)
boots<-matrix(0,n,1)
set.seed(123456)
for (i in 1:500){xb<-sample(Log_cost,n,replace = TRUE)
 boots[i]<-skewness(xb)}
z<-((sum(boots)/500)-skewness(Log_cost))/sd(boots)
sum((boots/500))
skewness(Log_cost)
quantile(boots,0.95)
SkewBoot(data, 500, 400,"Partial")
data2<-data.frame(as.matrix(cbind(cost,Log_cost,mix2,mix4,mix5,mix6,mix7,mix8,mix9,mix10)))

# Descriptive statistics
Means<-colMeans(data2)
Variances<-diag(var(data2))
Desviacion<-sqrt(Variances)
Q<-apply(data2,2,quantile)
# Los concatenamos por filas
Estadisticos<-rbind(Means,Desviacion,Q)
Estadisticos<-t(Estadisticos)
colnames(Estadisticos)<-c('Means','STD','Min','Q25','Median','Q75','Max')

knitr::kable(Estadisticos,digits=8, caption="Estad??sticos Descriptivos")
hist(data2$mix10)
sum(data2$mix10>15)
quantile(data2$mix10,c(0.8,0.9,0.91,0.95,0.99))
# Scatterplots
layout(matrix(c(1,2,3,4,5,6,7,0), 2, 4, byrow = TRUE))
plot(mix2,Log_cost,pch=19, cex=.3, xlab="Age", ylab="log(cost)")
plot(mix4,Log_cost,pch=19, cex=.3, xlab="Age of driving licence", ylab="log(cost)")
plot(mix5,Log_cost,pch=19, cex=.3, xlab="Age of car", ylab="log(cost)")
# plot(mix6,lcoste_medio, xlab="=1 if car is parked in the garage at night or $=0$ on the contrary ", ylab="Average cost per claim")
plot(mix7,Log_cost,pch=19, cex=.3, xlab="Annual kilometers", ylab="log(cost)")
plot(mix8,Log_cost,pch=19, cex=.3, xlab="Percentage of night kilometers", ylab="log(cost)")
plot(mix9,Log_cost,pch=19, cex=.3, xlab="Percentage of kilometers on urban roads", ylab="log(cost)")
plot(mix10,Log_cost,pch=19, cex=.3, xlab="Percentage of kilometers with speeding ", ylab="log(cost)")

# Covariates in the trhree estimated models: TO SELECT ONE VECTOR
# dataX<-as.matrix(cbind(mix10,mix2,mix4,mix5,mix6,mix7,mix8,mix9)) #All variables
# dataX<-as.matrix(cbind(mix10,mix7,mix8,mix9)) #Only telematics
dataX<-as.matrix(cbind(mix2,mix4,mix5,mix6)) #Only no telematics

dimX<-ncol(dataX)
nn<-nrow(dataX)

Y<-Log_cost
summary(Y)

# Initial parameters
ini<-lm(Y~dataX)
summary(ini)

attributes(ini)
ini<-as.matrix(ini$coefficients)
ini<-ini[-1]
ini<-as.matrix(ini/ini[1])
ini<-t(ini)

# kernel function, its derivatives and its integral; Gaussian

Ke<-function(x){return(dnorm(x))}

dKe<-function(x){
  r<--dnorm(x)*x
  return(r)}

ddKe<-function(x){
  r<-dnorm(x)*(x-1)*(x+1)
  return(r)}

pKe<-function(x){return(pnorm(x))}

# R rounds all numbers smaller than c1 to 0. 
# That's why I introduce c1 so that the likelihood is well defined.
c1=5e-324

# nn=sample size

singleindex<-function(nn,t2=ini[-1], c1=5e-324){
  dim1<-length(t2)
  one<-1
  # save results for estimated parameters
  # when two selected bandwidths  
    tediffh<-array(0, dim=c(dim1,1))
  
 
    Xv<-dataX
    
    #define theta, while theta1=1 (always! Model restriction)
    
    thetav<-c(one,t2)
    
    print(thetav)
    
    # generate data
    
    score1<-colSums(thetav*t(Xv))
    Y1<-Y
    
    # likelihood function. It returns -loglikelihood and hence we have minimization problem
    #Theta=argumn
    #h1 is associated with Theta'*x and h2 with y.
    likelih<-function(argumn,h1,h2){
      
      Likelihood<-rep(1,times=nn)
      
      for(k in 1:nn){
        
        # to define conditional density and survival function we need following help functions
        
        help1<-c(Ke(c(sum(c(one,argumn)*Xv[k,])-colSums(c(one,argumn)*t(Xv)))/h1)*Ke((Y1[k]-Y1)/h2))
        
        help2<-c(Ke(c(sum(c(one,argumn)*Xv[k,])-colSums(c(one,argumn)*t(Xv)))/h1))
        
        
        term1<-sum(help1[-k])/(h2*sum(help2[-k]))
        
        if(sum(help1[-k])==0){term1<-0}
        
        logDensity<-log(max(term1,c1))
        
        Likelihood[k]<-logDensity
      }
      return(-sum(Likelihood))
    }
    
    
    # likelihhod function needed to optimize in (h1,h2)    
    
    likeliarg<-function(hh,ar){
      return(likelih(argumn=ar,h1=hh[1],h2=hh[2]))
    }
    
    # likelihood function needed to optimize in (h sd(X*theta),h*sd(Z))    
    likelihnew<-function(ar,hh){
      return(likelih(argumn=ar,h1=hh*sd(colSums(c(one,ar)*t(Xv))),h2=hh*sdkm(Y1=Y1)))
    }
    
    # standard deviation of Z
    
    sdkm<-function(Y1=Y1){
      return(sqrt((sum(Y1^2)/nn)-(sum(Y1)/nn)^2))
    }
    
    ########################################################
    
    # choosing start parameters for optimization
    
    h1start<-nn^(-2/13)
    h2start<-nn^(-4/13)
    
    h1A<-1
    h2A<-1
    teA<-rep(0,times=dim1)
    te<-rep(1,times=dim1)
    
    count<-0
      while((count<50)&((sum(abs(teA-te)/dim1)>0.000001)|(abs(h1A-h1start)>0.001)|
                      (abs(h2A-h2start)>0.1))){
      h1A<-h1start
      h2A<-h2start
      te<-teA
      
      teA<-optim(par=c(te), fn=likelih,h1=h1start,h2=h2start)$par
      
      # this next line is "just in case". It checks if the new found 
      # parameter is smaller than the previous one (recall that we seach for minimum)
      
      if(likelih(argumn=teA,h1=h1A,h2=h2A)>likelih(argumn=te,h1=h1A,h2=h2A)){
        teA<-te
      }
      
      # minimize -loglikelihood in (h1,h2)
      
      
      sdX1<-sd(colSums(ini%*%t(Xv)))
      sdX2<-IQR(colSums(ini%*%t(Xv)))/1.349
      sdX<-min(sdX1,sdX2)
      lowX<-0.1*sdX*nn^(-2/13)
      upX<-4*sdX*nn^(-2/13)
      sdY1<-sdkm(Y1=Y1)
      sdY2<-IQR(Y1)/1.349
      sdY<-min(sdY1,sdY2)
      lowY<-0.1*sdY*nn^(-4/13)
      upY<-4*sdY*nn^(-4/13)
      
      ha<-optim(par=c(h1start,h2start), fn=likeliarg,lower=c(lowX,upX),upper=c(lowY,upY),method="L-BFGS-B",ar=teA)
      
      
      h1start<-ha$par[1]
      h2start<-ha$par[2]
      
    
      print("count="); print(count)
      count<-count+1
      
    }
    
    tediffh<-c(teA)
    

  ##################### printing results ########################

  print("estimated parameter different h:")
  result<-c(h1start,h2start,tediffh)  
  return(result)
}

tiempo.ini<-Sys.time()
par<-singleindex(nn=nn)
tiempo.fin<-Sys.time()
tiempo<-tiempo.fin-tiempo.ini
tiempo

# Printing the estimated parameter
par

# Likelihood function
likeli<-function(Xv,Y1,argumn,h1,h2){
  argumn<-argumn[-1]
  Likelihood<-rep(1,times=nn)
  
  for(k in 1:nn){
    
    # to define conditional density and survival function we need following help functions
    
    help1<-c(Ke(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1)*Ke((Y1[k]-Y1)/h2))
    
    help2<-c(Ke(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1))
    
    
    term1<-sum(help1[-k])/(h2*sum(help2[-k]))
    
    if(sum(help1[-k])==0){term1<-0}
    
    logDensity<-log(max(term1,c1))
    
    Likelihood[k]<-logDensity
  }
  return(sum(Likelihood))
}


# Optimal parameters
theta<-c(1,par[2:dimX+1])
theta<-as.matrix(theta)

h1<-par[1]
h2<- par[2]
h1
h2
# Maximum of Log-likelihood
logLik<-likeli(dataX,Y,theta,h1,h2)
logLik


# Ploting index
index<-dataX%*%theta
par(mfrow=c(1,1))
plot(density(index))
Xv<-dataX
Y1<-Y

#Conditionated density
condens<-function(argumn,h1,h2){
  argumn<-argumn[-1]
  cdens<-rep(0,times=nn)
  for(k in 1:nn){
    
    # to define conditional density and survival function we need following help functions
    
    help1<-c(Ke(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1)*Ke((Y1[k]-Y1)/h2))
    
    help2<-c(Ke(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1))
    
    dens<-sum(help1[-k])/(h2*(sum(help2[-k]+c1)))
    cdens[k]<-(max(c1,dens))
  }
  return(cdens)
}

# First derivative of likelihood
deriv1<-function(argumn,h1,h2,j){
  argumn<-argumn[-1]
  ss<-rep(0,times=nn)
  for(k in 1:nn){
    # to define conditional density and survival function we need following help functions
    help1<-c(Ke((Y1[k]-Y1)/h2)*dKe(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1)*((Xv[k,j]-Xv[,j])/h1))
    help2<-c(Ke(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1))
    help3<-c(Ke(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1)*Ke((Y1[k]-Y1)/h2))
    help4<-c(dKe(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1)*((Xv[k,j]-Xv[,j])/h1))
    num<-(1/h2)*((sum(help1[-k])*sum(help2[-k]))-(sum(help3[-k])*sum(help4[-k])))
    den<-(sum(help2[-k]))**2
    ss[k]<-num/max(c1,den)
  }
  return(ss)
}

# Gradient
fderiv<-matrix(0,nn,dimX-1)
sderiv<-array(0,dim=c(dimX-1,dimX-1))
for(j in 1:dimX-1){
  fderiv[,j]<-deriv1(theta,h1,h2,j+1)/condens(theta,h1,h2)
}

# Sigma1
Sigma1<-(1/nn)*t(fderiv)%*%fderiv
det(Sigma1)

# Second derivative of likelihood
deriv2<-function(argumn,h1,h2,j1,j2){
  argumn<-argumn[-1]
  ss<-rep(0,times=nn)
  for(k in 1:nn){
    # to define conditional density and survival function we need following help functions
    help1<-c(Ke((Y1[k]-Y1)/h2)*ddKe(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1)*((Xv[k,j1]-Xv[,j1])/h1)*((Xv[k,j2]-Xv[,j2])/h1))
    help2<-c(Ke(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1))
    help3<-c(Ke((Y1[k]-Y1)/h2)*dKe(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1)*((Xv[k,j1]-Xv[,j1])/h1))
    help4<-c(dKe(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1)*((Xv[k,j2]-Xv[,j2])/h1))
    help5<-c(Ke((Y1[k]-Y1)/h2)*Ke(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1))
    help6<-c(Ke((Y1[k]-Y1)/h2)*dKe(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1)*((Xv[k,j2]-Xv[,j2])/h1))
    help7<-c(dKe(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1)*((Xv[k,j1]-Xv[,j1])/h1))
    help8<-c(ddKe(c(sum(c(1,argumn)*Xv[k,])-colSums(c(1,argumn)*t(Xv)))/h1)*((Xv[k,j1]-Xv[,j1])/h1)*((Xv[k,j2]-Xv[,j2])/h1))
    num<-(1/h2)*((((sum(help1[-k]))*(sum(help2[-k]))+(sum(help3[-k]))*(sum(help4[-k])))*(sum(help2[-k])**2)-(2*sum(help3[-k])*(sum(help2[-k])**2)*sum(help4[-k]))-
                    (((sum(help6[-k]))*(sum(help7[-k]))+(sum(help5[-k]))*(sum(help8[-k])))*(sum(help2[-k])**2)-(2*sum(help5[-k])*sum(help7[-k])*
                                                                                                                  sum(help2[-k])*sum(help4[-k])))))
    
    
    den<-(sum(help2[-k]))**4
    ss[k]<-num/max(c1,den)
  }
  return(ss)
}

# Hessian

id2<-(condens(theta,h1,h2))**2    
w1<-1/max(c1,id2)
id1<-condens(theta,h1,h2)
w2<-1/max(c1,id1)
for(j in 1:dimX-1){
  for(k in 1:dimX-1){
    dd<--w1*deriv1(theta,h1,h2,j+1)*deriv1(theta,h1,h2,k+1)+w2*deriv2(theta,h1,h2,j+1,k+1)
    sderiv[j,k]<-sum(dd)
  }
}

det(sderiv)
# Sigma2
Sigma2<-solve(sderiv)

# Variance-covariance matrix of parameters 
Sigma<-Sigma2%*%Sigma1%*%Sigma2

# Standard errors
STE<-sqrt(diag(Sigma))

# Standard Normal ratio
z<-theta[-1]/STE
STE
z
# P-values
pval<-1-pnorm(abs(z))

# RESULTS
print("Estimation Results")
Results<-cbind(theta,c(0,STE),c(0,z),c(0,pval))

colnames(Results)<-c('Coefficients','STE','Z','p-value')

knitr::kable(Results,digits=6, caption="Estimation Results")

# Predictive analysis
# Conditionated cummulative distribution function of y given x
condcdf<-function(y,x,ydata,xdata,theta,nr,h1,h2){
  help1<-c(Ke(c(sum(theta*x)-colSums(theta*t(xdata)))/h1)*pKe((y-ydata)/h2))
  help2<-c(Ke(c(sum(theta*x)-colSums(theta*t(xdata)))/h1))
  if(sum(help2)>=c1){
    rr<-sum(help1)/((sum(help2)))}
  if(sum(help2)<c1){
    rr<-sum(help1)/((sum(help2+c1)))}
  # if(rr>1){rr<-1}
  return(rr)}

# Conditionated probability distribution function of y given the vector x
condpdf<-function(y,x,ydata,xdata,theta,nr,h1,h2){
  help1<-c(Ke(c(sum(theta*x)-colSums(theta*t(xdata)))/h1)*Ke((y-ydata)/h2))
  help2<-c(Ke(c(sum(theta*x)-colSums(theta*t(xdata)))/h1))
  if(sum(help2)>=c1){
    dens<-sum(help1)/(h2*(sum(help2)))}
  if(sum(help2)<c1){
    dens<-sum(help1)/(h2*(sum(help2+c1)))}
  rr<-(max(c1,dens))
  return(rr)
}


# Conditionated quantile given vector x
FUN<- function(y,x,ydata,xdata,theta,nr,h1,h2,alpha){
  cdf<-condcdf(y,x,ydata,xdata,theta,nr,h1,h2)
  rr=cdf-alpha
  return(rr)
}

newton<-function(y,x,ydata,xdata,theta,nr,h1,h2,alpha,converge){
  #z is a scalar, v is a vector
  ff<-FUN(y,x,ydata,xdata,theta,nr,h1,h2,alpha)
  while(abs(ff)>converge){
    dd<-condpdf(y,x,ydata,xdata,theta,nr,h1,h2)
    delta=-ff/dd
    y=y+delta
    ff=FUN(y,x,ydata,xdata,theta,nr,h1,h2,alpha)}
  return(y)
}


# Firts derivative
deriv_ind<-function(y,x,ydata,xdata,theta,nr,h1,h2){
  help1<-c(Ke(c(sum(theta*x)-colSums(theta*t(xdata)))/h1)*pKe((y-ydata)/h2))
  help2<-c(Ke(c(sum(theta*x)-colSums(theta*t(xdata)))/h1))
  cdf<-sum(help1)/sum(help2)
  help3<-c(dKe(c(sum(theta*x)-colSums(theta*t(xdata)))/h1)*pKe((y-ydata)/h2))
  help4<-c(dKe(c(sum(theta*x)-colSums(theta*t(xdata)))/h1))
  rr<-(1/h1)*((sum(help3)/sum(help2))-cdf*(sum(help4)/sum(help2)))
  return(rr)  
}

# Conditionated cummulative distribution function of y given the value of index
condcdf_index<-function(y,ind,ydata,index,nr,h1,h2){
  #z is a scalar, v is a vector, x is a vector of variable values.
  help1<-c(Ke(c(ind-index)/h1)*pKe((y-ydata)/h2))
  help2<-c(Ke(c(ind-index)/h1))
  if(sum(help2)>=c1){
    rr<-sum(help1)/sum(help2)}
  if(sum(help2)<c1){
    rr<-sum(help1)/sum(help2+c1)}
  return(rr)}

# Conditionated cummulative probability function of y given the value of index
condpdf_index<-function(y,ind,ydata,index,nr,h1,h2){
  # to define conditional density and survival function we need following help functions
  help1<-c(Ke(c(ind-index)/h1)*Ke((y-ydata)/h2))
  help2<-c(Ke(c(ind-index)/h1))
  if(sum(help2)>=c1){
    dens<-sum(help1)/(h2*(sum(help2)))}
  if(sum(help2)<c1){
    dens<-sum(help1)/(h2*(sum(help2+c1)))}
  rr<-(max(c1,dens))
  return(rr)
}

# Conditionated quantile given value of index
FUN_index<- function(y,ind,ydata,index,nr,h1,h2,alpha){
  cdf<-condcdf_index(y,ind,ydata,index,nr,h1,h2)
  rr=cdf-alpha
  return(rr)
}
newton_index<-function(y,ind,ydata,index,nr,h1,h2,alpha,converge){
  #z is a scalar, v is a vector
  ff<-FUN_index(y,ind,ydata,index,nr,h1,h2,alpha)
  # print(ff)
  while(abs(ff)>converge){
    dd<-condpdf_index(y,ind,ydata,index,nr,h1,h2)
    # print(dd)
    if(dd==0){dd=c1}
    delta=-ff/dd
    y=y+delta
    # print(delta)
    ff=FUN_index(y,ind,ydata,index,nr,h1,h2,alpha)}
  return(y)
}

FUN_index_y<- function(y){
  #z is a scalar, v is a vector
  cdf<-condcdf_index(y,ind,ydata,index,nr,h1,h2)
  rr=cdf-alpha
  return(rr)
}

index<-dataX%*%theta
summary(index)
theta<-t(theta)
thetac<-c(theta)
icost<-cbind(index,Y)
sdX1<-sd(index)
sdX2<-IQR(index)/1.349
sdX<-min(sdX1,sdX2)
h1<- as.numeric(1*sdX*nn**(-1/5))*1.5
sdY1<-sd(Y)
sdY2<-IQR(Y)/1.349
sdY<-min(sdY1,sdY2)
h2<- as.numeric(3.572*sdY1*nn**(-1/3))

plot(density(index))
pd<-density(index)
# h1<-pd$bw  

mi<-min(index)
ma<-max(index)
grid<-as.matrix((mi*100):(ma*100)/100)

alpha=0.50
y<-0
ind<-grid[1]
ydata<-Y
index<-index
nr<-nn
FUN_index_y(100)
FUN_index_y(-10)
uniroot(FUN_index_y,c(-1000,10000))$root
CM_Index50<-matrix(0,nrow(grid),1)
for(i in 1:nrow(grid)){
  ind<-grid[i]
  ydata<-Y
  index<-index
  nr<-nn
  CM_Index50[i]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
} 

alpha=0.90
y<-0

CM_Index90<-matrix(0,nrow(grid),1)
for(i in 1:nrow(grid)){
  ind<-grid[i]
  ydata<-Y
  index<-index
  nr<-nn
  CM_Index90[i]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
} 

alpha=0.95
y<-0

CM_Index95<-matrix(0,nrow(grid),1)
for(i in 1:nrow(grid)){
  ind<-grid[i]
  ydata<-Y
  index<-index
  nr<-nn
  CM_Index95[i]<-as.numeric(uniroot(FUN_index_y,c(-100,500))$root)
} 


alpha=0.99
CM_Index99<-matrix(0,nrow(grid),1)  
for(i in 1:nrow(grid)){
  ind<-grid[i]
  ydata<-Y
  index<-index
  nr<-nn
  CM_Index99[i]<-as.numeric(uniroot(FUN_index_y,c(-10,500))$root)
}

CM_Index50_or<-exp(CM_Index50)
CM_Index90_or<-exp(CM_Index90)
CM_Index95_or<-exp(CM_Index95)
CM_Index99_or<-exp(CM_Index99)

# No paramétrico
h1b<-h1*1

CE_Index<-matrix(0,nrow(grid),1)
for(i in 1:nrow(grid)){
  argum<-(grid[i]-index)/h1b
  num<-sum(Ke(argum)*Y)
  dem<-sum(Ke(argum))
  if(dem==0){dem<-c1}
  CE_Index[i]<-num/dem} 

# Plotting conditional mean and quantile
plot(index,exp(Y), pch=1, cex=.7,xlim = c(mi,34),ylim=c(0,25),xlab='Index',ylab='Fitted Cost')
lines(grid,CM_Index90_or,lty=2,lwd=1)
lines(grid,CM_Index95_or,lty=3,lwd=1)
lines(grid,CM_Index99_or,lty=4,lwd=1)
lines(grid,exp(CE_Index),lty=1,lwd=1)
lines(grid,CM_Index50_or,lty=5,lwd=1)
legend('topright',c('Mean','Median','0.90 quantile','0.95 quantile','0.99 quantile' ),lty =c(1,5,2,3,4),lwd = c(1,1,1,1,1),cex = 0.60,bty = "n" )


# Marginal Effect
summary(dataX)
nrow(dataX)
index<-dataX%*%theta
summary(index)
theta<-t(theta)
thetac<-c(theta)
icost<-cbind(index,Y)
sdX1<-sd(index)
sdX2<-IQR(index)/1.349
sdX<-min(sdX1,sdX2)
h1<- as.numeric(1*sdX*nn**(-1/5))*1.5
sdY1<-sd(Y)
sdY2<-IQR(Y)/1.349
sdY<-min(sdY1,sdY2)
h2<- as.numeric(3.572*sdY1*nn**(-1/3))


c1=5e-324

# Marginal effect of the skm
gridskm<-as.matrix((0:150)/10)
thetac<-c(theta)
CQ_Index_skm<-matrix(0,nrow(gridskm),2)
X<-as.matrix(dataX)
# At the q50
meskm50<-matrix(0,nrow(gridskm),2)
alpha=0.5
# x<-as.matrix(colMedians(dataX))
x<-as.matrix(apply(dataX,2,min))
# x<-as.matrix(colMeans(dataX))
x[2]<-20.59
x[5]<-1
x
for(i in 1:nrow(gridskm)){
  x[1]<-gridskm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_skm[i,1]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  meskm50[i,1]<--(1/condpdf(CQ_Index_skm[i,1],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_skm[i,1],x,Y,X,thetac,nn,h1,h2)*(thetac[1]/CQ_Index_skm[i,1])
} 

x[2]<-34.07
x[5]<-1
x
for(i in 1:nrow(gridskm)){
  x[1]<-gridskm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_skm[i,2]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  meskm50[i,2]<--(1/condpdf(CQ_Index_skm[i,2],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_skm[i,2],x,Y,X,thetac,nn,h1,h2)*(thetac[1]/CQ_Index_skm[i,2])
} 
# At the q95

meskm95<-matrix(0,nrow(gridskm),2)
alpha=0.95
# x<-as.matrix(colMedians(dataX))
x<-as.matrix(apply(dataX,2,min))
# x<-as.matrix(colMeans(dataX))
x[2]<-20.59
x[5]<-1
x
for(i in 1:nrow(gridskm)){
  x[1]<-gridskm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_skm[i,1]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  meskm95[i,1]<--(1/condpdf(CQ_Index_skm[i,1],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_skm[i,1],x,Y,X,thetac,nn,h1,h2)*(thetac[1]/CQ_Index_skm[i,1])
} 

x[2]<-34.07
x[5]<-1
x
for(i in 1:nrow(gridskm)){
  x[1]<-gridskm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_skm[i,2]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  meskm95[i,2]<--(1/condpdf(CQ_Index_skm[i,2],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_skm[i,2],x,Y,X,thetac,nn,h1,h2)*(thetac[1]/CQ_Index_skm[i,2])
} 

summary(meskm50)
summary(meskm95)
par(mfrow=c(1,2))
plot(gridskm,meskm50[,1],xlim=c(0,15),ylim=c(-0.40,0.14), cex.main=0.8,type="l",xlab="Speed %Km",ylab="Marginal Effects",main="Median")
lines(gridskm,meskm50[,2],col=1,lty=2)

plot(gridskm,meskm95[,1],xlim=c(0,15),ylim=c(-0.014,0.05), cex.main=0.8,type="l",xlab="Speed %Km",ylab="Marginal Effects",main="0.95 quantile")
lines(gridskm,meskm95[,2],col=1,lty=2)



c1=5e-324

# Marginal effect of the tkm
gridtkm<-as.matrix((0:350)/10)
thetac<-c(theta)
CQ_Index_tkm<-matrix(0,nrow(gridtkm),2)
X<-as.matrix(dataX)
# At the q50
metkm50<-matrix(0,nrow(gridtkm),2)
alpha=0.5
# x<-as.matrix(colMedians(dataX))
x<-as.matrix(apply(dataX,2,min))
# x<-as.matrix(colMeans(dataX))
x[2]<-20.59
x[5]<-1
x
for(i in 1:nrow(gridtkm)){
  x[6]<-gridtkm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_tkm[i,1]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  metkm50[i,1]<--(1/condpdf(CQ_Index_tkm[i,1],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_tkm[i,1],x,Y,X,thetac,nn,h1,h2)*(thetac[6]/CQ_Index_tkm[i,1])
} 
summary(metkm50)
x[2]<-34.07
x[5]<-1
x
for(i in 1:nrow(gridtkm)){
  x[6]<-gridtkm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_tkm[i,2]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  metkm50[i,2]<--(1/condpdf(CQ_Index_tkm[i,2],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_tkm[i,2],x,Y,X,thetac,nn,h1,h2)*(thetac[6]/CQ_Index_tkm[i,2])
} 
# At the q95

metkm95<-matrix(0,nrow(gridtkm),2)
alpha=0.95
# x<-as.matrix(colMedians(dataX))
x<-as.matrix(apply(dataX,2,min))
# x<-as.matrix(colMeans(dataX))
x[2]<-20.59
x[5]<-1
x
for(i in 1:nrow(gridtkm)){
  x[6]<-gridtkm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_tkm[i,1]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  metkm95[i,1]<--(1/condpdf(CQ_Index_tkm[i,1],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_tkm[i,1],x,Y,X,thetac,nn,h1,h2)*(thetac[6]/CQ_Index_tkm[i,1])
} 

x[2]<-34.07
x[5]<-1
x
for(i in 1:nrow(gridtkm)){
  x[6]<-gridtkm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_tkm[i,2]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  metkm95[i,2]<--(1/condpdf(CQ_Index_tkm[i,2],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_tkm[i,2],x,Y,X,thetac,nn,h1,h2)*(thetac[6]/CQ_Index_tkm[i,2])
} 
summary(metkm95)
par(mfrow=c(1,2))
plot(gridtkm,metkm50[,1],xlim=c(0,35),ylim=c(-0.004,0.0014), cex.main=0.8,type="l",xlab="Total Km",ylab="Marginal Effects",main="Median")
lines(gridtkm,metkm50[,2],col=1,lty=2)

plot(gridtkm,metkm95[,1],xlim=c(0,35),ylim=c(-0.0010,-0.0012), cex.main=0.8,type="l",xlab="Total Km",ylab="Marginal Effects",main="0.95 quantile")
lines(gridtkm,metkm95[,2],col=1,lty=2)

c1=5e-324

# Marginal effect of the nkm
gridnkm<-as.matrix((0:400)/10)
thetac<-c(theta)
CQ_Index_nkm<-matrix(0,nrow(gridnkm),2)
X<-as.matrix(dataX)
# At the q50
menkm50<-matrix(0,nrow(gridnkm),2)
alpha=0.5
# x<-as.matrix(colMedians(dataX))
x<-as.matrix(apply(dataX,2,min))
# x<-as.matrix(colMeans(dataX))
x[2]<-20.59
x[5]<-1
x
for(i in 1:nrow(gridnkm)){
  x[7]<-gridnkm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_nkm[i,1]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  menkm50[i,1]<--(1/condpdf(CQ_Index_nkm[i,1],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_nkm[i,1],x,Y,X,thetac,nn,h1,h2)*(thetac[7]/CQ_Index_nkm[i,1])
} 
summary(menkm50)
x[2]<-34.07
x[5]<-1
x
for(i in 1:nrow(gridnkm)){
  x[7]<-gridnkm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_nkm[i,2]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  menkm50[i,2]<--(1/condpdf(CQ_Index_nkm[i,2],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_nkm[i,2],x,Y,X,thetac,nn,h1,h2)*(thetac[7]/CQ_Index_nkm[i,2])
} 
# At the q95

menkm95<-matrix(0,nrow(gridnkm),2)
alpha=0.95
# x<-as.matrix(colMedians(dataX))
x<-as.matrix(apply(dataX,2,min))
# x<-as.matrix(colMeans(dataX))
x[2]<-20.59
x[5]<-1
x
for(i in 1:nrow(gridnkm)){
  x[7]<-gridnkm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_nkm[i,1]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  menkm95[i,1]<--(1/condpdf(CQ_Index_nkm[i,1],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_nkm[i,1],x,Y,X,thetac,nn,h1,h2)*(thetac[7]/CQ_Index_nkm[i,1])
} 

x[2]<-34.07
x[5]<-1
x
for(i in 1:nrow(gridnkm)){
  x[7]<-gridnkm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_nkm[i,2]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  menkm95[i,2]<--(1/condpdf(CQ_Index_nkm[i,2],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_nkm[i,2],x,Y,X,thetac,nn,h1,h2)*(thetac[7]/CQ_Index_nkm[i,2])
} 
summary(menkm50)
summary(menkm95)
par(mfrow=c(1,2))
plot(gridnkm,menkm50[,1],xlim=c(0,35),ylim=c(0.015,0.002), cex.main=0.8,type="l",xlab="Night %Km",ylab="Marginal Effects",main="Median")
lines(gridnkm,menkm50[,2],col=1,lty=2)

plot(gridnkm,menkm95[,1],xlim=c(0,35),ylim=c(0.0032,0.0016), cex.main=0.8,type="l",xlab="Night %Km",ylab="Marginal Effects",main="0.95 quantile")
lines(gridnkm,menkm95[,2],col=1,lty=2)


# Marginal effect of the ukm
gridukm<-as.matrix((50:700)/10)
thetac<-c(theta)
CQ_Index_ukm<-matrix(0,nrow(gridukm),2)
X<-as.matrix(dataX)
# At the q50
meukm50<-matrix(0,nrow(gridukm),2)
alpha=0.5
# x<-as.matrix(colMedians(dataX))
x<-as.matrix(apply(dataX,2,min))
# x<-as.matrix(colMeans(dataX))
x[2]<-20.59
x[5]<-1
x
for(i in 1:nrow(gridukm)){
  x[8]<-gridukm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_ukm[i,1]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  meukm50[i,1]<--(1/condpdf(CQ_Index_ukm[i,1],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_ukm[i,1],x,Y,X,thetac,nn,h1,h2)*(thetac[8]/CQ_Index_ukm[i,1])
} 
summary(meukm50)
x[2]<-34.07
x[5]<-1
x
for(i in 1:nrow(gridukm)){
  x[8]<-gridukm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_ukm[i,2]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  meukm50[i,2]<--(1/condpdf(CQ_Index_ukm[i,2],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_ukm[i,2],x,Y,X,thetac,nn,h1,h2)*(thetac[8]/CQ_Index_ukm[i,2])
} 
# At the q95

meukm95<-matrix(0,nrow(gridukm),2)
alpha=0.95
# x<-as.matrix(colMedians(dataX))
x<-as.matrix(apply(dataX,2,min))
# x<-as.matrix(colMeans(dataX))
x[2]<-20.59
x[5]<-1
x
for(i in 1:nrow(gridukm)){
  x[8]<-gridukm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_ukm[i,1]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  meukm95[i,1]<--(1/condpdf(CQ_Index_ukm[i,1],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_ukm[i,1],x,Y,X,thetac,nn,h1,h2)*(thetac[8]/CQ_Index_ukm[i,1])
} 

x[2]<-34.07
x[5]<-1
x
for(i in 1:nrow(gridukm)){
  x[8]<-gridukm[i]
  ind<-sum(x*thetac)
  xc<-c(x)
  ydata<-Y
  index<-index
  nr<-nn
  y<-0
  CQ_Index_ukm[i,2]<-as.numeric(uniroot(FUN_index_y,c(-10,10))$root)
  meukm95[i,2]<--(1/condpdf(CQ_Index_ukm[i,2],xc,Y,X,thetac,nn,h1,h2))*deriv_ind(CQ_Index_ukm[i,2],x,Y,X,thetac,nn,h1,h2)*(thetac[8]/CQ_Index_ukm[i,2])
} 
summary(meukm50)
summary(meukm95)
summary(menkm50)
summary(menkm95)
par(mfrow=c(1,2))
plot(gridskm,meskm50[,1],xlim=c(0,15),ylim=c(-0.40,0.14), cex.main=0.8,type="l",xlab="Speed %Km",ylab="Marginal Effects",main="Median")
lines(gridskm,meskm50[,2],col=1,lty=2)

plot(gridskm,meskm95[,1],xlim=c(0,15),ylim=c(-0.014,0.05), cex.main=0.8,type="l",xlab="Speed %Km",ylab="Marginal Effects",main="0.95 quantile")
lines(gridskm,meskm95[,2],col=1,lty=2)

par(mfrow=c(1,2))
plot(gridtkm,metkm50[,1],xlim=c(0,35),ylim=c(-0.004,0.0014), cex.main=0.8,type="l",xlab="Total Km",ylab="Marginal Effects",main="Median")
lines(gridtkm,metkm50[,2],col=1,lty=2)

plot(gridtkm,metkm95[,1],xlim=c(0,35),ylim=c(-0.0010,-0.0012), cex.main=0.8,type="l",xlab="Total Km",ylab="Marginal Effects",main="0.95 quantile")
lines(gridtkm,metkm95[,2],col=1,lty=2)


par(mfrow=c(1,2))
plot(gridnkm,menkm50[,1],xlim=c(0,35),ylim=c(0.002,0.015), cex.main=0.8,type="l",xlab="Night %Km",ylab="Marginal Effects",main="Median")
lines(gridnkm,menkm50[,2],col=1,lty=2)

plot(gridnkm,menkm95[,1],xlim=c(0,35),ylim=c(0.0016,0.0032), cex.main=0.8,type="l",xlab="Night %Km",ylab="Marginal Effects",main="0.95 quantile")
lines(gridnkm,menkm95[,2],col=1,lty=2)

par(mfrow=c(1,2))
plot(gridukm,meukm50[,1],xlim=c(5,70),ylim=c(-0.011,0.018), cex.main=0.8,type="l",xlab="Urban %Km",ylab="Marginal Effects",main="Median")
lines(gridukm,meukm50[,2],col=1,lty=2)

plot(gridukm,meukm95[,1],xlim=c(5,70),ylim=c(-0.002,0.0040), cex.main=0.8,type="l",xlab="Urban %Km",ylab="Marginal Effects",main="0.95 quantile")
lines(gridukm,meukm95[,2],col=1,lty=2)