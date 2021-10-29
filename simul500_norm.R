# setwd("D:\\RiskCenter\\Ricardo\\Cati_Montse\\Simulaciones")
setwd("D:\\RiskCenter\\Ricardo\\Cati_Montse\\Simulaciones\\MetRicardo\\n500")

# install.packages("RGCCA")
library(RGCCA)

# NUMBER OF REPLICAS
nrep<-500

# Sample size
n<-500
n
# Explanatory variables (trivariate)
# install.packages("mvtnorm")
library(mvtnorm)
# Number of explanatory variables
k=3
x3<-matrix(0,n,(nrep*k))
sig0<-diag(k)

set.seed(123)
for (i in 1:nrep){x3[,(k*(i-1)+1):(k*i)]<-rmvnorm(n, mean = rep(0, k), sigma = sig0)}

# Index
index<-matrix(0,n,nrep)
dep_norm0<-matrix(0,n,nrep)
b0<-as.matrix(c(1,1.3,0.5))
for (j in 1:nrep){
for (i in 1:n){
  index[i,j]<-x3[i,(k*(j-1)+1):(k*j)]%*%b0
  dep_norm0[i,j]<-rnorm(1, mean = index[i,j], sd = abs(index[i,j]))
  }
}
# Nonparametric index estimation

b_ini<-matrix(1,nrep,k)
for (j in 1:nrep){
Y<-dep_norm0[,j]
datosX<-x3[,(k*(j-1)+1):(k*j)]
Yr<-Y-datosX[,1]
datosXr<-datosX[,2:k]
ini<-lm(Yr~datosXr)
ini<-as.matrix(ini$coefficients)
ini<-ini[-1]
# ini<-as.matrix(ini/ini[1])
ini<-as.matrix(ini)
b_ini[j,2:k]<-t(ini)
}

b_ini2<-b_ini[,2:k]

bias_lm<-colSums(b_ini2-b0[2:k])/nrep-b0[2:k]
bias_lm
S_lm<-var(b_ini2)
MSE_lm<-diag(S_lm)+bias_lm^2
MSE_lm
nrep
n
singleindex<-function(nn,t2=b_ini[,2:k], type, c1=5e-324, MM=nrep){
  
  dim1<-ncol(t2)
  
  # save results for estimated parameters
  # when two selected bandwidths  
  
  tediffh<-array(0, dim=c(MM,dim1))
  dd<-array(0, dim=c(MM,3))
  mh<-array(0, dim=c(MM,2))
  Ke<-function(x){return(dnorm(x))}
  
  dKe<-function(x){
    r<--dnorm(x)*x
    return(r)}
  
  ddKe<-function(x){
    r<-dnorm(x)*(x-1)*(x+1)
    return(r)}
  
  pKe<-function(x){return(pnorm(x))}
  # Ke<-function(x){
  #   r<-3/4*(1-x^2)*(abs(x)<=1)
  #   return(r)}
  # dKe<-function(x){
  #   r<--3/2*x*(abs(x)<=1)
  #   return(r)}
  # ddKe<-function(x){
  #   r<--3/2*(abs(x)<=1)
  #   return(r)}
  # 
  # pKe<-function(x){
  #   r<-(3/4*x-(x^3)/4+1/2)*(abs(x)<=1)+(x>1)
  #   return(r)}
  
  
  # likelihood function. It returns -loglikelihood and hence we have minimization problem
  #Theta=argumn
  #h1 is associated with Theta'*x and h2 with y.
  likelih<-function(argumn,h1,h2){
    
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
    return(-sum(Likelihood))
  }
  
  
  # likelihhod function needed to optimize in (h1,h2)    
  
  likeliarg<-function(hh,ar){
    return(likelih(argumn=ar,h1=hh[1],h2=hh[2]))
  }
  
  # likelihood function needed to optimize in (h sd(X*theta),h*sd(Z))    
  likelihnew<-function(ar,hh){
    return(likelih(argumn=ar,h1=hh*sd(colSums(c(1,ar)*t(Xv))),h2=hh*sdkm(Y1=Y1)))
  }
  
  # standard deviation of Z
  
  sdkm<-function(Y1=Y1){
    return(sqrt((sum(Y1^2)/nn)-(sum(Y1)/nn)^2))
  }
  l<-1
  
  # can be used for(l in 1:MM); I cannot recall why used "while" but it does not make a difference
  
  while(l<=MM){
    datosX<-x3[,(k*(l-1)+1):(k*l)]
    
    Xv<-datosX
    
    #define theta, while theta1=1 (always! Model restriction)
    
    # thetav<-c(1,t2[l,])
    
    # print(thetav)
    
    # generate data
    
    # score1<-colSums(thetav*t(Xv))
    Y1<-dep_norm0[,l]
    
    # kernel function, its derivatives and its integral; Gaussian
    

    
    ########################################################
    
    # choosing start parameters for optimization
    
    h1start<-nn^(-2/13)
    h2start<-nn^(-4/13)
    
    h1A<-1
    h2A<-1
    teA<-t2[l,]
    # teA<-t2
    # teA<-rep(0,times=dim1)
    te<-rep(1,times=dim1)
    
    count<-0
    iter<-count
    # we optimize first in the parameter theta then in bandwidths
    # the optimization goes until the distances between consecutive steps <0.01 or 50 steps were reached
    # in our censored case average count was 3.    
    
    #     while((count<50)&((abs(teA-te)[1]>0.000001)|(abs(teA-te)[2]>0.000001)|(abs(teA-te)[3]>0.000001)|(abs(teA-te)[4]>0.000001)|
    #                       (abs(teA-te)[5]>0.000001)|(abs(teA-te)[6]>0.000001)|(abs(teA-te)[7]>0.000001)|(abs(teA-te)[8]>0.000001)|
    #                       (abs(teA-te)[9]>0.000001)|(abs(teA-te)[10]>0.000001)|(abs(h1A-h1start)>0.001)|
    #                       (abs(h2A-h2start)>0.001))){
    
    #while((count<50)&((abs(teA-te)[1]>0.000001)|(abs(teA-te)[2]>0.000001)|(abs(teA-te)[3]>0.000001)|(abs(teA-te)[4]>0.000001)|
    #                  (abs(teA-te)[5]>0.000001)|(abs(teA-te)[6]>0.000001)|(abs(teA-te)[7]>0.000001)|(abs(teA-te)[8]>0.000001)|
    #                  (abs(teA-te)[9]>0.000001)|(abs(h1A-h1start)>0.0001)|
    #                  (abs(h2A-h2start)>0.0001))){
    
    while((count<100)&((sum(abs(teA-te)/dim1)>0.001)|(abs(h1A-h1start)>0.001)|
                      (abs(h2A-h2start)>0.001))){
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
      b<-c(1,teA)
      sdX<-sd(colSums(b%*%t(Xv)))
      lowX<-0.01*sdX*nn^(-2/13)
      upX<-10*sdX*nn^(-2/13)
      sdY<-sdkm(Y1=Y1)
      #sdY<-IQR(Y1)/1.349
      lowY<-0.01*sdY*nn^(-4/13)
      upY<-10*sdY*nn^(-4/13)
      
      ha<-optim(par=c(h1start,h2start), fn=likeliarg,lower=c(lowX,lowY),upper=c(upX,upY),method="L-BFGS-B",ar=teA)
      
      
      h1start<-ha$par[1]
      h2start<-ha$par[2]
      
      
      #     print(c(h1start,h2start,teA))
      
      print("count="); print(count)
      count<-count+1
      d1<-sum(abs(teA-te)/dim1)
      d2<-abs(h1A-h1start)
      d3<-abs(h1A-h1start)
    }
    
    # print(c(h1start,h2start,teA))
    
    tediffh[l,]<-c(teA)
    mh[l,]<-c(h1start,h2start)
    print("l="); print(l)
    dd[l,]<-as.matrix(c(d1,d2,d3))
    l<-l+1
  }
  
  
  ##################### printing results ########################
  print(iter)
  mtheta<-as.matrix(tediffh)
  mh<-as.matrix(mh)
  result<-cbind(mtheta,mh,dd)
  write.csv(result,file = "result_norm0_500.txt")
  # print("mean(estimated parameter) different h:")
  
  # print(rowSums(tediffh)/MM)
  # 
   # print("bias:")
  # 
   bias1<-colSums(tediffh)/MM-b0[2:k]
  # 
   print(bias1)
  # 
  # 
 sthe<-var(tediffh)
  # te2covariance<-array(vect1, dim=c(3,3))
  # 
  print("covariance of theta:")
  print(sthe)
  # 
  print("MSE")
  # 
  mse1<-diag(sthe)+bias1^2
  # 
  print(mse1)
  
}

tiempo.ini<-Sys.time()
singleindex(nn=n,type=1)
tiempo.fin<-Sys.time()
tiempo<-tiempo.fin-tiempo.ini
tiempo