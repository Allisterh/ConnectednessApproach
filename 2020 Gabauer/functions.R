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
  inc = c(format(round(colSums(CT), digit),nsmall=digit), "TCI")
  npdc = c(format(round(NPDC,digit),nsmall=digit), TCI*k/(k-1))
  TABLE = rbind(table,to,inc,net,npdc)
  colnames(TABLE) = c(rownames(CV),"FROM others")
  rownames(TABLE) = c(rownames(CV),"TO others","Inc. own","NET","NPDC")
  PCI = matrix(NA, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      PCI[i,j] = 2*(CT[i,j]+CT[j,i])/(CT[i,i]+CT[i,j]+CT[j,i]+CT[j,j])
    }
  }
  return = list(CT=CT,TCI=TCI,TCI_corrected=TCI*k/(k-1),PCI=PCI,TO=TO,FROM=FROM,NET=NET,NPSO=NPSO,NPDC=NPDC,TABLE=TABLE)
}
VGFEVD = function(dcc.fit, h=100, standardize=FALSE) {
  NAMES = dcc.fit@model$modeldata$asset.names
  R = rcor(dcc.fit)
  R.bar = apply(R,1:2,mean)
  Q.bar = dcc.fit@mfit$Qbar
  H = rcov(dcc.fit)
  t = dim(H)[3]
  # alpha
  alpha = array(0,c(k,k,h))
  alpha[,,1] = diag(dcc.fit@mfit$matcoef[c(seq(3,(4*k),4)),1])
  # beta
  beta = diag(dcc.fit@mfit$matcoef[c(seq(4,(4*k),4)),1])
  # ALPHA
  ALPHA = dcc.fit@mfit$matcoef[(4*k+1),1]
  # BETA
  BETA = dcc.fit@mfit$matcoef[(4*k+2),1]

  H.hat = array(0,c(k,k,h+1))
  GVIRF = H.hat.shock = H.hat.no_shock = array(0,c(k,k,t,h+1))
  e = diag(k)
  for (i in 1:t) {
    H.hat[,,1] = H[,,i]
    Q.hat = H.hat
    Q.hat[,,1] = dcc.fit@mfit$Q[[i]]
    for (j in 1:h) {
      H.hat[,,j+1] = (alpha[,,j])%*%e^2 + beta%*%H.hat[,,j]
      D = diag(diag(H.hat[,,j+1])^0.5)
      u = D%*%e
      if (j==1) {
        Q.hat[,,2] = (1-ALPHA-BETA)*Q.bar + ALPHA*crossprod(u) + BETA*H.hat[,,1]
      } else {
        Q.hat[,,j+1] = (1-ALPHA-BETA)*Q.bar + (ALPHA+BETA)*Q.hat[,,j]
      }
      R.hat = diag(1/(diag(Q.hat[,,j+1])^0.5))%*%Q.hat[,,j+1]%*%(diag(1/diag(Q.hat[,,j+1])^0.5))
      H.hat[,,j+1] = D%*%R.hat%*%D
    }
    H.hat.shock[,,i,] = H.hat
  }
  if (standardize) {
    e = 0*diag(k)
    for (i in 1:t) {
      H.hat[,,1] = H[,,i]
      Q.hat = H.hat
      Q.hat[,,1] = dcc.fit@mfit$Q[[i]]
      for (j in 1:h) {
        H.hat[,,j+1] = beta%*%H.hat[,,j]
        D = diag(diag(H.hat[,,j+1])^0.5)
        if (j==1) {
          Q.hat[,,2] = (1-ALPHA-BETA)*Q.bar + BETA*H.hat[,,1]
        } else {
          Q.hat[,,j+1] = (1-ALPHA-BETA)*Q.bar+(ALPHA+BETA)*Q.hat[,,j]
        }
        R.hat = diag(1/(diag(Q.hat[,,j+1])^0.5))%*%Q.hat[,,j+1]%*%(diag(1/diag(Q.hat[,,j+1])^0.5))
        H.hat[,,j+1] = D%*%R.hat%*%D
      }
      H.hat.no_shock[,,i,] = H.hat
    }
    #for (i in 1:t) {
    #  GVIRF[,,i,] = H.hat.shock[,,i,] - H.hat.no_shock[,,i,]
    #}
  }# else {
  for (i in 1:t) {
    GVIRF[,,i,] = H.hat.shock[,,i,] - H.hat.no_shock[,,i,]
  }
  #GVIRF = H.hat.shock
  #}
  GVFEVD = array(NA, c(k,k,t))
  rownames(GVFEVD)=colnames(GVFEVD)=NAMES
  for (i in 1:t) {
    num = apply(GVIRF[,,i,]^2,1:2,sum)
    den = c(apply(num,1,sum))
    fevd = t(num)/den
    GVFEVD[,,i] = (fevd/apply(fevd, 1, sum))
  }
  return = list(GVFEVD=GVFEVD, GVIRF=GVIRF)
}
Moments = function(data){
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
    bt2 = WeightedPortTest::Weighted.Box.test(data[,i], type="Ljung-Box", lag=20, sqrd.res=T)
    moments[13,i] = bt2$statistic
    moments[14,i] = bt2$p.value
  }
  
  cc=c(4,6,8,10,12,14)
  moments = round(moments,3)
  moments1 = moments
  for (j in 1:k){
    for (i in 1:length(cc)){
      i = cc[i]
      if (moments[i,j]<=0.01) {
        moments1[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"***",sep="")
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else if (moments[i,j]<=0.05) {
        moments1[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"**",sep="")
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else if (moments[i,j]<=0.10) {
        moments1[(i-1),j] = paste(format(round(moments[(i-1),j],3),nsmall=3),"*",sep="")
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      } else {
        moments1[(i-1),j] = format(round(moments[(i-1),j],3),nsmall=3)
        moments1[i,j] = paste("(",format(round(moments[i,j],3),nsmall=3),")",sep="")
      }
    }
  }
  
  for (j in 1:k){
    i = 9
    if (moments[i,j]<=-2.57) {
      moments1[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"***",sep="")
    } else if (moments[i,j]<=-1.96) {
      moments1[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"**",sep="")
    } else if (moments[i,j]<=-1.62) {
      moments1[i,j] = paste(format(round(moments[i,j],3),nsmall=3),"*",sep="")
    } else {
      moments1[i,j] = format(round(moments[i,j],3),nsmall=3)
    }
  }
  moments1
}
