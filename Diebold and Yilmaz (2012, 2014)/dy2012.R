
### DIEBOLD, FX., AND YILMAZ, K. (2012)
### BETTER TO GIVE THAN TO RECEIVE: PREDICTIVE DIRECTIONAL MEASUREMENT OF VOLATILITY SPILLOVERS
### International Journal of Forecasting
### replicated by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)

gfevd = function(model, n.ahead=10, normalize=TRUE, standardize=TRUE) {
   if (class(model) != "varest") {
      return("The model class is not varest!")
   }
   A = Phi(model, (n.ahead-1))
   epsilon = residuals(model)
   Sigma = t(epsilon)%*%epsilon / (model$obs)
   gi = array(0, dim(A))
   sigmas = sqrt(diag(Sigma))
   for (j in 1:dim(A)[3]) {
      gi[,,j] = t(A[,,j]%*%Sigma%*%solve(diag(sqrt(diag(Sigma)))))
   }
   
   if (standardize==TRUE){
      girf=array(NA, c(dim(gi)[1],dim(gi)[2], (dim(gi)[3])))
      for (i in 1:dim(gi)[3]){
         girf[,,i]=((gi[,,i])%*%solve(diag(diag(gi[,,1]))))
      }
      gi = girf
   }
   
   num = apply(gi^2,1:2,sum)
   den = c(apply(num,1,sum))
   fevd = t(num)/den
   if (normalize==TRUE) {
      fevd=(fevd/apply(fevd, 1, sum))
   } else {
      fevd=(fevd)
   }
   return = list(fevd=fevd, girf=gi)
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

### STATIC CONNECTEDNESS APPROACH
library("vars")
nlag = 4 # VAR(4)
nfore = 10 # 10-step ahead forecast
var_full = VAR(Y, p=nlag, type="const")
CV_full = gfevd(var_full, n.ahead=nfore)$fevd
rownames(CV_full)=colnames(CV_full)=colnames(Y)
print(DCA(CV_full))

### DYNAMIC CONNECTEDNESS APPROACH
t = nrow(Y)
space = 200 + nlag # 200 days rolling window estimation
CV = array(NA, c(k, k, (t-space)))
colnames(CV) = rownames(CV) = colnames(Y)
for (i in 1:dim(CV)[3]) {
  var = VAR(Y[i:(space+i-1),], p=nlag, type="const")
  CV[,,i] = gfevd(var, n.ahead=nfore)$fevd
  if (i%%500==0) {print(i)}
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
