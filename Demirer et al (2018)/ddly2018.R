
### DEMIRER, M., DIEBOLD, FX., LIU, L., & YILMAZ, K. (2018)
### Estimating Global Bank Network Connectedness
### Journal Of Applied Econometrics
### replicated by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)

library("glmnet")
LASSO_RIDGE_VAR = function(y,p, alpha=0) {
  y = as.matrix(y)
  p = p
  res1 = coef1 = NULL
  k = ncol(y)
  for (i in 1:k) {
    yx = embed(y,p+1)
    x = yx[,-c(1:k)]
    z = y[-c(1:p),i]
    mod = cv.glmnet(x, z, alpha = alpha, type.measure="mae")
    pred = predict(mod, s = mod$lambda.min, newx=yx[,-c(1:k)])
    coef = predict(mod, type = "coefficients", s = mod$lambda.min)[-1]
    res = z - pred
    coef1 = rbind(coef1,coef) 
    res1 = cbind(res1,res)
  }  
  Q = (t(res1)%*%res1)/nrow(res1)
  results = list(B=coef1,Q=Q)
}
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
tvp.gfevd = function(model, Sigma, n.ahead=10,normalize=TRUE,standardize=TRUE) {
  A = tvp.Phi(model, (n.ahead-1))
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
  return = list(fevd=fevd, girf=gi, nfevd=nfevd)
}
DCA = function(CV){
  k = dim(CV)[1]
  SOFM = apply(CV,1:2,mean)*100 # spillover from others to one specific
  VSI = round(mean(100-diag(SOFM)),2)
  TO = colSums(SOFM-diag(diag(SOFM)))
  FROM = rowSums(SOFM-diag(diag(SOFM)))
  NET = TO-FROM
  NPSO = SOFM-t(SOFM)
  INC = rowSums(NPSO>0)
  ALL = rbind(format(round(cbind(SOFM,FROM),1),nsmall=1),c(format(round(TO,1),nsmall=1),format(round(sum(colSums(SOFM-diag(diag(SOFM)))),1),nsmall=1)),c(format(round(NET,1),nsmall=1),"TCI"),format(round(c(INC,VSI),1),nsmall=1))
  colnames(ALL) = c(rownames(CV),"FROM")
  rownames(ALL) = c(rownames(CV),"Contribution TO others","NET directional connectedness","NPDC transmitter")
  return = list(CT=SOFM,TCI=VSI,TO=TO,FROM=FROM,NET=NET,NPSO=NPSO,NPDC=INC,ALL=ALL)
}

path = file.path(file.choose()) # select dy2012.csv
DATA = read.csv(path)
DATE = as.Date(as.character(DATA[,1]))
Y = DATA[,-1]
k = ncol(Y)

#alpha = 0 # RIDGE VAR
alpha = 1 # LASSO VAR

### STATIC CONNECTEDNESS APPROACH
nlag = 4 # VAR(4)
nfore = 10 # 10-step ahead forecast
lasso_ridge_full = LASSO_RIDGE_VAR(Y, p=nlag, alpha=alpha)
CV_full = tvp.gfevd(lasso_ridge_full$B, lasso_ridge_full$Q, n.ahead=nfore)$fevd
rownames(CV_full)=colnames(CV_full)=colnames(Y)
print(DCA(CV_full))

### DYNAMIC CONNECTEDNESS APPROACH
t = nrow(Y)
space = 200 + nlag # 200 days rolling window estimation
CV = array(NA, c(k, k, (t-space)))
colnames(CV) = rownames(CV) = colnames(Y)
for (i in 1:dim(CV)[3]) {
  lasso_ridge_var = LASSO_RIDGE_VAR(Y[i:(space+i-1),], p=nlag, alpha=alpha)
  CV[,,i] = tvp.gfevd(lasso_ridge_var$B, lasso_ridge_var$Q, n.ahead=nfore)$fevd
  if (i%%100==0) {print(i)}
}

to = matrix(NA, ncol=k, nrow=t)
from = matrix(NA, ncol=k, nrow=t)
net = matrix(NA, ncol=k, nrow=t)
ct = npso = array(NA, c(k, k, t))
total = matrix(NA, ncol=1, nrow=t)
colnames(npso)=rownames(npso)=colnames(ct)=rownames(ct)=colnames(Y)
for (i in 1:t){
  CV = tvp.gfevd(B_t[,,i], Q_t[,,i], n.ahead=nfore)$fevd
  colnames(CV)=rownames(CV)=colnames(Y)
  vd = DCA(CV)
  ct[,,i] = vd$CT
  to[i,] = vd$TO/k
  from[i,] = vd$FROM/k
  net[i,] = vd$NET/k
  npso[,,i] = vd$NPSO/k
  ct[,,i]-t(ct[,,i])
  total[i,] = vd$TCI
}

nps = array(NA,c(t,k/2*(k-1)))
colnames(nps) = 1:ncol(nps)
jk = 1
for (i in 1:k) {
  for (j in 1:k) {
    if (j<=i) {
      next
    } else {
      nps[,jk] = npso[i,j,]
      colnames(nps)[jk] = paste0(colnames(Y)[i],"-",colnames(Y)[j])
      jk = jk + 1
    }
  }
}

### DYNAMIC TOTAL CONNECTEDNESS
date = DATE[-c(1:space)]
par(mfrow = c(1,1), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1,1) + .05, mgp = c(0, 0.1, 0))
plot(date,total, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(floor(min(total)),ceiling(max(total))),yaxs="i",xlab="",tck=0.01)
grid(NA,NULL,lty=1)
polygon(c(date,rev(date)),c(c(rep(0,nrow(total))),rev(total)),col="grey20", border="grey20")
box()

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
par(mfrow = c(k/2,2), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:k){
  plot(date,to[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"TO all others"),ylim=c(floor(min(to)),ceiling(max(to))),tck=0.01,yaxs="i")
  grid(NA,NULL,lty=1)
  polygon(c(date,rev(date)),c(c(rep(0,nrow(to))),rev(to[,i])),col="grey20", border="grey20")
  box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfrow = c(ceiling(k/2),2), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
for (i in 1:k){
  plot(date,from[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"FROM all others"),ylim=c(floor(min(from)),ceiling(max(from))),tck=0.01,yaxs="i")
  grid(NA,NULL,lty=1)
  polygon(c(date,rev(date)),c(c(rep(0,nrow(from))),rev(from[,i])),col="grey20", border="grey20")
  box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
par(mfrow = c(ceiling(k/2),2), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1,1) + .05, mgp = c(0, 0.1, 0))
for (i in 1:k){
  plot(date,net[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste("NET",colnames(Y)[i]),ylim=c(floor(min(net)),ceiling(max(net))),tck=0.01,yaxs="i")
  grid(NA,NULL,lty=1)
  polygon(c(date,rev(date)),c(c(rep(0,nrow(net))),rev(net[,i])),col="grey20", border="grey20")
  box()
}

### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
par(mfrow = c(ceiling(ncol(nps)/2),2), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1,1) + .05, mgp = c(0, 0.1, 0))
for (i in 1:ncol(nps)) {
  plot(date,nps[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=colnames(nps)[i],tck=0.02,yaxs="i",ylim=c(floor(min(nps)),ceiling(max(nps))))
  grid(NA,NULL,lty=1)
  polygon(c(date,rev(date)),c(c(rep(0,nrow(nps))),rev(nps[,i])),col="grey20", border="grey20")
  box()
}

### END
