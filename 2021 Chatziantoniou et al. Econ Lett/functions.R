if (!require("MASS")) install.packages("MASS")
if (!require("urca")) install.packages("urca")
if (!require("moments")) install.packages("moments")
if (!require("openxlsx")) install.packages("openxlsx")
if (!require("parallel")) install.packages("parallel")
if (!require("quantreg")) install.packages("quantreg")
if (!require("WeightedPortTest")) install.packages("WeightedPortTest")

QVAR = function(y, p, tau=0.5) {
  y = as.matrix(y)
  k = ncol(y)

  resid = B = NULL
  for (i in 1:k) {
    x = embed(y, p+1)
    rq_est = quantreg::rq(y[-c(1:p), i] ~ x[, -c(1:k)], tau=tau)
    B = rbind(B, coef(rq_est)[-1]) 
    resid = cbind(resid, rq_est$residuals)
  }  
  Q = (t(resid)%*%resid)/nrow(resid)
  results = list(B=B, Q=Q)
}
GFEVD = function(Phi, Sigma, n.ahead=10, normalize=TRUE, standardize=TRUE) {
  tvp.Phi = function (x, nstep = 10, ...) {
    nstep = abs(as.integer(nstep))
    K=nrow(x)
    p=floor(ncol(x)/K)
    A = array(0, c(K,K,nstep))
    for (i in 1:p){
      A[,,i]=x[,((i-1)*K+1):(i*K)]
    }
    
    Phi = array(0, dim = c(K, K, nstep + 1))
    Phi[, , 1] = diag(K)
    Phi[, , 2] = Phi[, , 1] %*% A[, , 1]
    if (nstep > 1) {
      for (i in 3:(nstep + 1)) {
        tmp1 = Phi[, , 1] %*% A[, , i - 1]
        tmp2 = matrix(0, nrow = K, ncol = K)
        idx = (i - 2):1
        for (j in 1:(i - 2)) {
          tmp2 = tmp2 + Phi[, , j + 1] %*% A[, , idx[j]]
        }
        Phi[, , i] = tmp1 + tmp2
      }
    }
    return(Phi)
  }
  A = tvp.Phi(Phi, n.ahead-1)
  gi = array(0, dim(A))
  sigmas = sqrt(diag(Sigma))
  for (j in 1:dim(A)[3]) {
    gi[,,j] = t(A[,,j] %*% Sigma %*% MASS::ginv(diag(sqrt(diag(Sigma)))))
  }
  if (standardize==TRUE){
    girf = array(NA, dim(gi))
    for (i in 1:dim(gi)[3]){
      girf[,,i] = gi[,,i] %*% MASS::ginv(diag(diag(gi[,,1])))
    }
    gi = girf
  }
  
  num = apply(gi^2, 1:2, sum)
  den = c(apply(num, 1, sum))
  fevd = t(num)/den
  nfevd = fevd
  if (normalize==TRUE) fevd = fevd/apply(fevd, 1, sum)
  return = list(GFEVD=fevd, GIRF=gi)
}
DCA = function(CV, digit=2){
  k = dim(CV)[1]
  CT = apply(CV,1:2,mean)*100 # spillover from others to one specific
  OWN = diag(diag(CT))
  TO = colSums(CT-OWN)
  FROM = rowSums(CT-OWN)
  NET = TO-FROM
  TCI = mean(TO)
  NPSO = CT-t(CT)
  NPDC = rowSums(NPSO>0)
  table = format(round(cbind(CT,FROM),digit),nsmall=digit)
  to = c(format(round(c(TO,sum(TO)),digit),nsmall=digit))
  net = c(format(round(c(NET, TCI),digit),nsmall=digit))
  npdc = c(format(round(NPDC,digit),nsmall=digit), "")
  inc = c(format(round(colSums(CT), digit),nsmall=digit), "TCI")
  TABLE = rbind(table,to,inc,net,npdc)
  colnames(TABLE) = c(rownames(CV),"FROM others")
  rownames(TABLE) = c(rownames(CV),"TO others","Inc. own","NET","NPDC")
  return = list(CT=CT,TCI=TCI,TCI_corrected=TCI*k/(k-1),TO=TO,FROM=FROM,NET=NET,NPSO=NPSO,NPDC=NPDC,TABLE=TABLE)
}
SummaryStatistics = function(data){
  data = as.matrix(data)
  moments = matrix(NA, ncol=ncol(data), nrow=14)
  colnames(moments)=colnames(data)
  rownames(moments)=c("Mean","Variance","Skewness","","Kurtosis","","JB","","ERS","","Q(20)","","Q2(20)","")
  for (i in 1:ncol(data)){
    moments[1,i] = mean(data[,i])
    moments[2,i] = var(data[,i])
    skew = moments::agostino.test(data[,i])
    moments[3,i] = skew$statistic[1]
    moments[4,i] = skew$p.value
    kurt = moments::anscombe.test(data[,i])
    moments[5,i] = kurt$statistic[1]-3
    moments[6,i] = kurt$p.value
    jb = moments::jarque.test(as.numeric(data[,i]))
    moments[7,i] = jb$statistic
    moments[8,i] = jb$p.value
    ers = urca::ur.ers(data[,i],type="DF-GLS",model="constant")
    moments[9,i] = ers@teststat
    moments[10,i]= ers@testreg$coefficients[1,4]
    bt = WeightedPortTest::Weighted.Box.test(data[,i], type="Ljung-Box", lag=20)
    moments[11,i] = bt$statistic
    moments[12,i] = bt$p.value
    bt2 = WeightedPortTest::Weighted.Box.test(data[,i], type="Ljung-Box", lag=20, sqrd.res=TRUE)
    moments[13,i] = bt2$statistic
    moments[14,i] = bt2$p.value
  }
  
  cc=c(4,6,8,10,12,14)
  moments = round(moments,3)
  MOMENTS = moments
  for (j in 1:k){
    for (i in 1:length(cc)){
      i = cc[i]
      if (moments[i,j]<=0.01) {
        MOMENTS[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"***",sep="")
        MOMENTS[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else if (moments[i,j]<=0.05) {
        MOMENTS[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"**",sep="")
        MOMENTS[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else if (moments[i,j]<=0.10) {
        MOMENTS[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"*",sep="")
        MOMENTS[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else {
        MOMENTS[(i-1),j] = format(round(moments[(i-1),j],3),nsmall=3)
        MOMENTS[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      }
    }
  }
  
  for (j in 1:k){
    i = 9
    if (moments[i,j]<=-2.57) {
      MOMENTS[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"***",sep="")
    } else if (moments[i,j]<=-1.96) {
      MOMENTS[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"**",sep="")
    } else if (moments[i,j]<=-1.62) {
      MOMENTS[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"*",sep="")
    } else {
      MOMENTS[i,j] = format(round(moments[i,j],3),nsmall=3)
    }
  }
  MOMENTS
}
net.contour = function(NET, selector, nlevels=20, threshold=0.01) {
  name = dimnames(NET)[[2]][selector]
  NET[which(NET<quantile(NET, threshold), arr.ind=TRUE)] = as.numeric(quantile(NET, threshold))
  NET[which(NET>quantile(NET, 1-threshold), arr.ind=TRUE)] = as.numeric(quantile(NET, 1-threshold))
  NET = NET[,selector,]
  lvls = seq(-ceiling(max(abs(NET))), ceiling(max(abs(NET))), length.out=nlevels)
  color.palette = function(n) {hcl.colors(n, "RdBu", rev=TRUE)}
  filled.contour(as.Date(rownames(NET)), as.numeric(colnames(NET)), NET,
                 col=color.palette(nlevels-1), levels=lvls, main=name, ylim=c(0,1))
}
