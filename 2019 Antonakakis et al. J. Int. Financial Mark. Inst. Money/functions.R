
SummaryStatistics = function(data){
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
VAR = function (x, p=1, include.mean=T, fixed=NULL) {
  if (!is.matrix(x)) 
    x = as.matrix(x)
  Tn = dim(x)[1]
  k = dim(x)[2]
  if (p < 1) 
    p = 1
  idm = k * p
  ne = Tn - p
  ist = p + 1
  y = x[ist:Tn, ]
  if (include.mean) {
    idm = idm + 1
    xmtx = cbind(rep(1, ne), x[p:(Tn - 1), ])
  }
  else {
    xmtx = x[p:(Tn - 1), ]
  }
  if (p > 1) {
    for (i in 2:p) {
      xmtx = cbind(xmtx, x[(ist - i):(Tn - i), ])
    }
  }
  ndim = ncol(xmtx)
  if (length(fixed) == 0) {
    paridx = matrix(1, ndim, k)
  }
  else {
    paridx = fixed
  }
  res = NULL
  beta = matrix(0, ndim, k)
  sdbeta = matrix(0, ndim, k)
  npar = 0
  for (i in 1:k) {
    idx = c(1:ndim)[paridx[, i] == 1]
    resi = y[, i]
    if (length(idx) > 0) {
      xm = as.matrix(xmtx[, idx])
      npar = npar + dim(xm)[2]
      xpx = t(xm) %*% xm
      xpxinv = solve(xpx)
      xpy = t(xm) %*% as.matrix(y[, i], ne, 1)
      betai = xpxinv %*% xpy
      beta[idx, i] = betai
      resi = y[, i] - xm %*% betai
      nee = dim(xm)[2]
      sse = sum(resi * resi)/(Tn - p - nee)
      dd = diag(xpxinv)
      sdbeta[idx, i] = sqrt(dd * sse)
    }
    res = cbind(res, resi)
  }
  sse = t(res) %*% res/(Tn - p)
  aic = 0
  bic = 0
  hq = 0
  Phi = NULL
  Ph0 = NULL
  jst = 0
  if (include.mean) {
    Ph0 = beta[1, ]
    se = sdbeta[1, ]
    jst = 1
  }
  if (include.mean) {
    for (i in 1:k) {
      if (abs(Ph0[i]) > 1e-08) 
        npar = npar - 1
    }
  }
  for (i in 1:p) {
    phi = t(beta[(jst + 1):(jst + k), ])
    se = t(sdbeta[(jst + 1):(jst + k), ])
    jst = jst + k
    Phi = cbind(Phi, phi)
  }
  dd = det(sse)
  d1 = log(dd)
  aic = d1 + (2 * npar)/Tn
  bic = d1 + log(Tn) * npar/Tn
  hq = d1 + 2 * log(log(Tn)) * npar/Tn
  VAR <- list(data = x, cnst = include.mean, order = p, coef = beta, 
              aic = aic, bic = bic, hq = hq, residuals = res, secoef = sdbeta, 
              Sigma = sse, Phi = Phi, Ph0 = Ph0, fixed = fixed)
}
TVPVAR = function(Y, l, nlag, prior){
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
  Y = scale(Y,T,F)
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
  l_2 = l[2]
  l_4 = l[1]
  
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
BayesPrior = function(Y, nlag){
  k = ncol(Y)
  vars = MTS::VAR(Y, p=nlag, include.mean=TRUE, output=FALSE)
  varcoef = t(vars$Phi)
  SIGMA_OLS = vars$secoef
  Q_0 = vars$Sigma
  b_prior = varcoef
  beta_0.var = diag(c(vars$secoef[-(k+1),]))^2
  return=list(aprior=b_prior,Vprior=beta_0.var,Q_0=Q_0)
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
  warch = WeightedPortTest::Weighted.LM.test(residuals(ugarch.fit), sigma(ugarch.fit)^2, lag=lag, type=c("correlation"), fitdf=2, weighted=TRUE)
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
        ugarch.spec = ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                                 variance.model=list(model="fGARCH", submodel=models[j], garchOrder=c(1,1), variance.targeting=FALSE), 
                                 distribution.model=distr[i])
      } else {
        ugarch.spec = ugarchspec(mean.model=list(armaOrder=c(ar,ma)),
                                 variance.model=list(model=models[j], garchOrder=c(1,1), variance.targeting=FALSE),
                                 distribution.model=distr[i])
      }
      ugarch.fit = ugarchfit(ugarch.spec, data, solver="hybrid", solver.list=list(outer.iter=10, inner.iter=1000, eval.se=FALSE, tol=1e-12))
      if (ugarch.fit@fit$convergence==0) {
        GARCH_IC[i,j] = InformationCriterion(ugarch.fit, ugarch.spec, prob=prob, conf.level=conf.level, lag=lag)
        spec_list[[j]][[i]] = ugarch.spec
      }
    }
  }
  return=list(GARCH_IC=GARCH_IC, spec_list=spec_list)
}
values = function(x){
  x = as.matrix(x)
  res = matrix(NA, nrow=ncol(x), ncol=4)
  for (i in 1:ncol(x)){
    res[i,]=matrix(c(mean(x[,i]), sd(x[,i]), quantile(x[,i],0.05), quantile(x[,i],0.95)), nrow=1)
  }
  colnames(res)=c("Mean","Std.Dev.","5%","95%")
  res
}
HedgeRatio = function(data, H, method="cumsum") {
  k=ncol(data)
  NAMES = matrix(NA, ncol=k, nrow=k)
  for (i in 1:k){
    NAMES[i,]=paste0(colnames(data),"/",colnames(data)[i])
  }
  lowpre = NAMES[lower.tri(diag(k))]
  
  DM = array(NA, c(k, k, (nrow(data))))
  for (i in 1:k){
    for (j in 1:k){
      DM[i,j,] = H[i,j,] / H[j,j,]
    }
  }
  
  HR = array(1, c(k, 4, k))
  for (i in 1:k){
    for (j in 1:k){
      HR[i,,j] = values(DM[i,j,])
    }
  }
  colnames(HR) = c("Mean","Std.Dev.","5%","95%")
  dimnames(HR)[[1]] = colnames(data)
  
  for (i in 1:k){
    if (i==1){
      HRatio=HR[,,i]
    } else {
      HRatio=rbind(HRatio, HR[,,i])
    }
  }
  rownames(HRatio) = c(t(NAMES))
  
  t = nrow(data)
  pval = HE = matrix(NA,ncol=k,nrow=k)
  colnames(HE)=rownames(HE)=colnames(pval)=rownames(pval)=colnames(data)
  portfolio_return = cumulative_portfolio_return = cumulative_asset_return = array(NA,c(k,k,t))
  cumulative_ratio = array(NA,c(k,k))
  for (i in 1:k) {
    for (j in 1:k) {
      portfolio_return[i,j,] = data[,i] - DM[i,j,]*data[,j]
      HE[i,j] = 1 - var(portfolio_return[i,j,])/var(data[,i])
      pval[i,j] = var.test(x=portfolio_return[i,j,],y=data[,i],ratio=1)$p.value
      if (method=="cumsum") {
        cumulative_asset_return[i,j,] = cumsum(data[,i])
        cumulative_portfolio_return[i,j,] = cumsum(portfolio_return[i,j,])
      } else if (method=="cumprod") {
        cumulative_asset_return[i,j,] = cumprod(1+data[,i])-1
        cumulative_portfolio_return[i,j,] = cumprod(1+portfolio_return[i,j,])-1
      }
      cumulative_ratio[i,j] = cumulative_portfolio_return[i, j, dim(cumulative_portfolio_return)[3]] / cumulative_asset_return[i,j,dim(cumulative_portfolio_return)[3]]
    }
  }
  table=cbind(round(HRatio,2),c(t(HE)),c(t(pval)))
  table=table[-which(table[,1]==1),]
  colnames(table)=c("Mean","Std.Dev.","5%","95%","HE","p-value")
  
  return = list(hedge_ratio=DM, summary=format(round(table,2),nsmall=2), portfolio_return=portfolio_return, cumulative_portfolio_return=cumulative_portfolio_return, cumulative_asset_return=cumulative_asset_return, cumulative_ratio=cumulative_ratio)
}
PortfolioWeights = function(data, H, method="cumsum", type="long") {
  k=ncol(data)
  NAMES=matrix(NA, ncol=k, nrow=k)
  for (i in 1:k){
    NAMES[i,]=paste0(colnames(data),"/",colnames(data)[i])
  }
  lowpre=NAMES[lower.tri(diag(k))]
  col=sum(lower.tri(diag(k)))
  sel=which(lower.tri(diag(ncol(data)))==TRUE, arr.ind=T)
  
  summary = NULL
  portfolio_weights = array(1, c(k, k, (nrow(data))))
  for (i in 1:k){
    for (j in 1:k){
      if (i==j) {
        pwpre = c(10,10,10,10)
        summary = rbind(summary,pwpre)
      } else {
        RR = (H[j,j,]-H[i,j,])/(H[i,i,]-2*H[i,j,]+H[j,j,])
        if (type=="long") {
          RR = ifelse(RR>1,1,RR)
          RR = ifelse(RR<0,0,RR)
        }
        pwpre = values(RR)
        summary = rbind(summary,pwpre)
        colnames(summary)=c("Mean","Std.Dev.","5%","95%")
        portfolio_weights[i,j,] = RR
      }
    }
  }
  rownames(summary) = c(NAMES)
  
  t = nrow(data)
  pval = cumulative_ratio = HE = matrix(NA,ncol=k,nrow=k)
  rownames(cumulative_ratio)=colnames(cumulative_ratio)=rownames(HE)=colnames(HE)=colnames(data)
  portfolio_return = cumulative_portfolio_return = cumulative_asset_return = array(NA,c(k,k,t))
  for (i in 1:k) {
    for (j in 1:k) {
      portfolio_return[i,j,] = portfolio_weights[i,j,]*data[,i] + (1-portfolio_weights[i,j,])*data[,j]
      HE[i,j] = 1 - var(portfolio_return[i,j,])/var(data[,i])
      pval[i,j] = var.test(x=portfolio_return[i,j,],y=data[,i],ratio=1)$p.value
      if (method=="cumsum") {
        cumulative_asset_return[i,j,] = cumsum(data[,i])
        cumulative_portfolio_return[i,j,] = cumsum(portfolio_return[i,j,])
      } else if (method=="cumprod") {
        cumulative_asset_return[i,j,] = cumprod(1+data[,i])
        cumulative_portfolio_return[i,j,] = cumprod(1+portfolio_return[i,j,])
      }
      cumulative_ratio[i,j] = cumulative_portfolio_return[i,j,dim(cumulative_portfolio_return)[3]]/cumulative_asset_return[i,j,dim(cumulative_asset_return)[3]]
    }
  }
  table = round(cbind(summary,c(t(HE)),c(t(pval))),2)
  rownames(table)=c(NAMES)
  table = table[-which(table[,1]==10),]
  colnames(table)=c("Mean","Std.Dev.","5%","95%","HE","p-value")
  
  return = list(portfolio_weights=portfolio_weights, summary=format(round(table,2),nsmall=2), portfolio_return=portfolio_return, cumulative_portfolio_return=cumulative_portfolio_return, cumulative_asset_return=cumulative_asset_return, cumulative_ratio=cumulative_ratio)
}
