
library("MTS")
library("MASS")
SummaryStatistics = function(data){
  moments = matrix(NA, ncol=ncol(data), nrow=14)
  colnames(moments)=colnames(data)
  rownames(moments)=c("Mean","Variance","Skewness","","Kurtosis","","JB","","ERS","","Q2(20)","","ARCH(20)","")
  for (i in 1:ncol(data)){
    moments[1,i] = mean(data[,i])
    moments[2,i] = var(data[,i])
    skew = moments::agostino.test(data[,i])
    moments[3,i] = skew$statistic[1]
    moments[4,i] = skew$p.value
    kurt = moments::anscombe.test(data[,i])
    moments[5,i] = kurt$statistic[1]-3
    moments[6,i] = kurt$p.value
    jb = moments::jarque.test(data[,i])
    moments[7,i] = jb$statistic
    moments[8,i] = jb$p.value
    ers = urca::ur.ers(data[,i],type="DF-GLS",model="constant")
    moments[9,i] = ers@teststat
    moments[10,i]= ers@testreg$coefficients[1,4]
    bt = WeightedPortTest::Weighted.Box.test(data[,i], type="Ljung-Box", lag=20)
    bt2 = WeightedPortTest::Weighted.Box.test(data[,i], type="Ljung-Box", lag=20, sqrd.res=T)
    moments[11,i] = bt2$statistic
    moments[12,i] = bt2$p.value
    bt3 = WeightedPortTest::Weighted.LM.test(data[,i], h.t=c(rep(var(data[,i]),nrow(data))), type="partial", lag=20)
    moments[13,i] = bt3$statistic
    moments[14,i] = bt3$p.value
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
Table = function(p){
  para = as.matrix(p)
  if (sum(is.nan(p)==TRUE)==0) {
    k = ncol(para)
    kk = nrow(para)/2
    for (i in 1:kk){
      for (j in 1:k){
        if (para[i,j]==0){
          next
        } else if (abs((para[2*i,j]))<=0.01) {
          para[i,j] = paste(format(round(as.numeric(para[i,j]),3),nsmall=3),"***",sep="")
          para[(kk+i),j] = paste("(",format(round(as.numeric(para[(kk+i),j]),3),nsmall=3),")",sep="")
        } else if (abs((para[2*i,j]))<=0.05) {
          para[i,j] = paste(format(round(as.numeric(para[i,j]),3),nsmall=3),"**",sep="")
          para[(kk+i),j] = paste("(",format(round(as.numeric(para[(kk+i),j]),3),nsmall=3),")",sep="")
        } else if (abs((para[2*i,j]))<=0.10) {
          para[i,j] = paste(format(round(as.numeric(para[i,j]),3),nsmall=3),"*",sep="")
          para[(kk+i),j] = paste("(",format(round(as.numeric(para[(kk+i),j]),3),nsmall=3),")",sep="")
        } else if (abs((para[2*i,j]))>0.1) {
          para[i,j] = paste(format(round(as.numeric(para[i,j]),3),nsmall=3),sep="")
          para[(kk+i),j] = paste("(",format(round(as.numeric(para[(kk+i),j]),3),nsmall=3),")",sep="")
        }
      }
    }
  }
  para
}
SignBias_WARCH = function(ugarch.fit,lag=20) {
  sign.bias = signbias(ugarch.fit)[1,][1:2]
  warch = Weighted.LM.test(residuals(ugarch.fit), sigma(ugarch.fit)^2, lag=lag, type=c("correlation"), fitdf=2, weighted=TRUE)
  pval = c(sign.bias[2]$prob,warch$p.value)
  x = Table(matrix(as.matrix(sign.bias),ncol=1))
  y = Table(matrix(c(warch$statistic,warch$p.value),ncol=1))
  statistic = rbind(x,y)
  rownames(statistic) = c("Sign Bias","",paste0("WARCH(",lag,")"),"")
  return = list(statistic=statistic, pval=pval)
}
ValueAtRisk = function(ugarch.fit, ugarch.spec, prob=0.05, conf.level=0.90){
  setfixed(ugarch.spec) = as.list(coef(ugarch.fit))
  data = ugarch.fit@model$modeldata$data
  ugarch.filter = ugarchfilter(ugarch.spec, data, n.old=length(data))
  VaR = fitted(ugarch.filter) + sigma(ugarch.filter)*qdist(ugarch.fit@model$modeldesc$distribution, p=prob, mu=ugarch.fit@fit$matcoef[1,1], sigma=1,
                                         skew=ifelse(is.na(coef(ugarch.fit)["skew"]),0,coef(ugarch.fit)["skew"]), shape=coef(ugarch.fit)["shape"])
  var_test = VaRTest(prob, as.numeric(data), as.numeric(VaR), conf.level=conf.level)
  vardur_test = VaRDurTest(p, as.numeric(data), as.numeric(VaR),conf.level=conf.level)
  f = function(x) {
    qdist(ugarch.fit@model$modeldesc$distribution, p=x, mu=0,sigma=1,skew=ifelse(is.na(coef(ugarch.fit)["skew"]),0,coef(ugarch.fit)["skew"]), shape=coef(ugarch.fit)["shape"])
  }
  ES = fitted(ugarch.filter) + sigma(ugarch.filter)*integrate(f,0,prob)$value/prob
  ES = ESTest(prob, as.numeric(data), as.numeric(ES), VaR, boot=TRUE, n.boot=1000, conf.level=conf.level)
  decision = ifelse(ES$boot.p.value>0.10,"H0","H1")
  x = Table(c(var_test$uc.LRstat, var_test$uc.LRp))
  y = Table(c(vardur_test$rLL,vardur_test$LRp))
  z = Table(c(1, ES$boot.p.value))
  z[1,1] = decision
  statistic = rbind(x,y,z)
  pval = c(var_test$uc.LRp,round(ES$boot.p.value,3),vardur_test$LRp)
  rownames(statistic) = c("VaR","","CVaR","","VaR Dur.","")
  return = list(statistic=statistic,pval=pval)
}
InformationCriterion = function (ugarch.fit, ugarch.spec, prob=0.05, conf.level=0.90, lag=20) {
  qprob = qnorm(1-prob)
  loss = sum(abs(ugarch.fit@fit$robust.tval[-c(1:2)])<=qprob)
  if (is.na(ugarch.fit@fit$robust.matcoef[1,2])==FALSE) {
    if ("skew" %in% rownames(ugarch.fit@fit$robust.matcoef)) {
      upper = ugarch.fit@fit$robust.matcoef["skew",1] + qprob*ugarch.fit@fit$robust.matcoef["skew",2]
      lower = ugarch.fit@fit$robust.matcoef["skew",1] - qprob*ugarch.fit@fit$robust.matcoef["skew",2]
      if (upper>1 && lower<1) {
        loss = loss + 100
      }
    }
  }
  var = ValueAtRisk(ugarch.fit, ugarch.spec, prob, conf.level)$pval
  sbwarch = SignBias_WARCH(ugarch.fit, lag=lag)$pval
  
  t = length(ugarch.fit@fit$z)
  IC = -2*likelihood(ugarch.fit) + loss*log(t)
  IC = IC + sum(c(var,sbwarch)<0.10)*10^5
  IC = ifelse(is.na(IC), 10^8, IC)
  IC
}
BestGARCH = function(distr=c("norm","snorm","std","sstd","ged","sged"), models=c("sGARCH","eGARCH","gjrGARCH","iGARCH","TGARCH","AVGARCH","NGARCH","NAGARCH","APARCH","ALLGARCH"), data, ar=0, ma=0, prob=0.05, conf.level=0.90, lag=20){
  data = matrix(data, ncol=1)
  GARCH_IC = matrix(10^7, nrow=length(distr), ncol=length(models))
  colnames(GARCH_IC) = models
  rownames(GARCH_IC) = distr
  spec_list = list()
  for (i in 1:length(models)) {
    spec_list[[i]] = list()
  }
  for (j in 1:length(models)) {
    print(paste0("-",models[j]))
    for (i in 1:length(distr)) {
      print(paste0("--",distr[i]))
      if (models[j] %in% c("AVGARCH","TGARCH","APARCH","NAGARCH","NGARCH","ALLGARCH")) {
        #if (models[j]=="AVGARCH") {
        #  fixed.pars=list(gamma1=0,gamma2=0,delta=1)
        #} else if (models[j]=="TGARCH") {
        #  fixed.pars=list(delta=1)
        #} else {
        #  fixed.pars=list()
        #}
        ugarch.spec = ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                                 variance.model=list(model="fGARCH", submodel=models[j], garchOrder=c(1,1), variance.targeting=FALSE), 
                                 distribution.model=distr[i])
                                 #fixed.pars=fixed.pars)
      } else {
        ugarch.spec = ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                                 variance.model=list(model=models[j], garchOrder=c(1,1), variance.targeting=FALSE),
                                 distribution.model=distr[i])
      }
      ugarch.fit = ugarchfit(ugarch.spec, data, solver="hybrid", solver.list=list(outer.iter=10, inner.iter=1000))
      if (ugarch.fit@fit$convergence==0) {
        GARCH_IC[i,j] = InformationCriterion(ugarch.fit, ugarch.spec, prob=prob, conf.level=conf.level, lag=lag)
        spec_list[[i]][[j]] = ugarch.spec
      }
    }
  }
  return=list(GARCH_IC=GARCH_IC, spec_list=spec_list)
}
UninformativePrior = function(gamma, r, nlag){
  m = nlag*(r^2)
  A_prior = cbind(0*diag(r), matrix(0, ncol=(nlag-1)*r, nrow=r))
  aprior = c(A_prior)
  V_i = matrix(0, nrow=(m/r), ncol=r)
  for (i in 1:r){
    for (j in 1:(m/r)) {
      V_i[j,i] = gamma/(ceiling(j/r)^2)
    }
  }
  # Now V (MINNESOTA VARIANCE) is a diagonal matrix with diagonal elements the V_i'  
  V_i_T = t(V_i)
  Vprior = diag(c(V_i_T))
  diag(Vprior)
  return = list(aprior=aprior, Vprior=Vprior)
}
BayesPrior = function(Y, nlag){
  k = ncol(Y)
  vars = VAR(Y, p=nlag, include.mean=TRUE, output=FALSE)
  varcoef = vars$Phi
  SIGMA_OLS = vars$secoef
  Q_0 = vars$Sigma
  b_prior = 0*varcoef
  beta_0.var = diag(c(vars$secoef[-(k+1),]))^2
  return=list(aprior=b_prior,Vprior=beta_0.var,Q_0=Q_0)
}
TVPVAR = function(Y, l, nlag, prior, demean=TRUE){
  beta_0.mean = prior$aprior
  beta_0.var = prior$Vprior
  Q_0 = prior$Q_0
  if (is.null(Q_0)) {
    Q_0 = cov(Y)
  }
  
  create_RHS_NI = function(templag, r, nlag, t){
    K = nlag*(r^2)
    x_t = matrix(0, (t-nlag)*r, K)
    for (i in 1:(t-nlag)){
      ztemp=NULL
      for (j in 1:nlag){
        xtemp = templag[i,((j-1)*r+1):(j*r)]
        xtemp = t(kronecker(diag(r),xtemp))
        ztemp = cbind(ztemp, xtemp)
      }
      x_t[((i-1)*r+1):(i*r),] = ztemp
    }
    return=list(x_t=x_t, K=K)
  }
  Y = scale(Y,demean,FALSE)
  y_true = 0
  FPC = Y
  YX = cbind(Y,Y)
  nfac = 0
  p = n = ncol(Y)
  r = nfac + p
  m = nlag*(r^2)
  k = nlag*r
  t = nrow(FPC)
  q = n + p
  Q_0 = Q_0
  
  # Initialize matrices
  beta_0_prmean = beta_0.mean
  beta_0_prvar = beta_0.var
  
  beta_pred = matrix(0,m,t)
  beta_update = matrix(0,m,t)
  
  Rb_t = array(0,c(m,m,t))
  Sb_t = array(0,c(m,m,t))
  
  beta_t = array(0, c(k,k,t))
  Q_t = array(0, c(r,r,t))
  
  # Decay and forgetting factors
  l_2 = l[1]
  l_4 = l[2]
  
  # Define lags of the factors to be used in the state (VAR) equation         
  yy = FPC[(nlag+1):t,]      
  xx = embed(FPC,nlag+1)[,-c(1:ncol(FPC))]
  templag = embed(FPC,nlag+1)[,-c(1:ncol(FPC))]
  RHS1 = create_RHS_NI(templag,r,nlag,t);  
  Flagtemp = RHS1$x_t
  m = RHS1$K
  Flag = rbind(matrix(0, k,m), Flagtemp)
  
  ###-----| 1. KALMAN FILTER
  for (irep in 1:t){
    #-----| Update the state covariances
    # 1. Get the variance of the factor
    
    # Update Q[t]
    if (irep==1){
      Q_t[,,irep] = Q_0
    } else if (irep > 1) {
      if (irep <= (nlag+1)) { 
        Gf_t = 0.1*(t(matrix(FPC[irep,],nrow=1))%*%(FPC[irep,]))
      } else {
        Gf_t = t(yy[(irep-nlag),]-xx[(irep-nlag),]%*%t(B[1:r,1:k])) %*% (yy[(irep-nlag),]-xx[(irep-nlag),]%*%t(B[1:r,1:k]))
      }
      Q_t[,,irep] = l_2*Q_t[,,(irep-1)] + (1-l_2)*Gf_t[1:r,1:r]
    }
    # -for beta
    if (irep <= (nlag+1)) {
      beta_pred[,irep] = beta_0_prmean
      beta_update[,irep] = beta_pred[,irep]
      Rb_t[,,irep] = beta_0_prvar
    } else if (irep > (nlag+1)) {
      beta_pred[,irep] = beta_update[,(irep-1)]
      Rb_t[,,irep] = (1/l_4)*Sb_t[,,(irep-1)]
    }
    
    # -for beta
    if (irep >= (nlag+1)) {
      # 2/ Update VAR coefficients conditional on Principal Componets estimates
      Rx = Rb_t[,,irep]%*%t(Flag[((irep-1)*r+1):(irep*r),])
      KV_b = Q_t[,,irep] + Flag[((irep-1)*r+1):(irep*r),]%*%Rx
      KG = Rx%*%MASS::ginv(KV_b)
      beta_update[,irep] = matrix(beta_pred[,irep], ncol=1) + (KG%*%(t(matrix(FPC[irep,], nrow=1))-Flag[((irep-1)*r+1):(irep*r),]%*%matrix(beta_pred[,irep], ncol=1)) )
      Sb_t[,,irep] = Rb_t[,,irep] - KG%*%(Flag[((irep-1)*r+1):(irep*r),]%*%Rb_t[,,irep])
    }
    
    # Assign coefficients
    bb = matrix(beta_update[,irep], ncol=1)
    splace = 0
    biga = matrix(0, r,r*nlag)
    for (ii in 1:nlag) {                                          
      for (iii in 1:r) {           
        biga[iii,((ii-1)*r+1):(ii*r)] = t(bb[(splace+1):((splace+r)),1])
        splace = splace + r
      }
    }
    
    B = rbind(biga, cbind(diag(r*(nlag-1)), matrix(0, nrow=r*(nlag-1), ncol=r)))
    
    if ((max(abs(eigen(B)$values))<=1)||(irep==1)){
      beta_t[,,irep] = B
    } else {
      beta_t[,,irep] = beta_t[,,(irep-1)]
      beta_update[,irep] = 0.99*beta_update[,(irep-1)]
    }
  }
  
  return = list(beta_t=beta_t[1:ncol(Y),,], Q_t=Q_t)
}
GFEVD = function(Phi, Sigma, n.ahead=10,normalize=TRUE,standardize=TRUE) {
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
  A = tvp.Phi(Phi, (n.ahead-1))
  Sigma = Sigma
  gi = array(0, dim(A))
  sigmas = sqrt(diag(Sigma))
  for (j in 1:dim(A)[3]) {
    gi[,,j] = t(A[,,j]%*%Sigma%*%MASS::ginv(diag(sqrt(diag(Sigma)))))
  }
  if (standardize==TRUE){
    girf=array(NA, c(dim(gi)[1],dim(gi)[2], (dim(gi)[3])))
    for (i in 1:dim(gi)[3]){
      girf[,,i]=((gi[,,i])%*%MASS::ginv(diag(diag(gi[,,1]))))
    }
    gi=girf
  }
  
  num = apply(gi^2,1:2,sum)
  den = c(apply(num,1,sum))
  fevd = t(num)/den
  nfevd = fevd
  if (normalize==TRUE) {
    fevd=(fevd/apply(fevd, 1, sum))
  } else {
    fevd=(fevd)
  }
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
