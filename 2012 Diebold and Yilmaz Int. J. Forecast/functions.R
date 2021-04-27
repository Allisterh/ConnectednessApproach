library("MASS")
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
