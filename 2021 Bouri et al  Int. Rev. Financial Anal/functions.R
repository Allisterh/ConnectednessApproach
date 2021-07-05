
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
