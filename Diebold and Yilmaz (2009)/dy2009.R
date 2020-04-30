
### DIEBOLD, FX., AND YILMAZ, K. (2009)
### MEASURING FINANCIAL ASSET RETURN AND VOLATILITY SPILLOVERS, WITH APPLICATION TO GLOBAL EQUITY MARKETS
### Economic Journal
### by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)
library("vars")

ofevd = function(ofevd, n.ahead=10) {
  k = length(ofevd)
  fevd = matrix(NA,k,k)
  for (i in 1:k) {
    fevd[i,] = ofevd[[i]][n.ahead,]
  }
  result = list(fevd=fevd)
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
rownames(Y) = DATE
k = ncol(Y)

### STATIC CONNECTEDNESS APPROACH
nlag = 4 # VAR(4)
nfore = 10 # 10-step ahead forecast
var_full = VAR(Y, p=nlag, type="const")
CV_full = ofevd(fevd(var_full, n.ahead=nfore), n.ahead=nfore)$fevd
rownames(CV_full)=colnames(CV_full)=colnames(Y)
print(DCA(CV_full)$ALL)

### DYNAMIC CONNECTEDNESS APPROACH
t = nrow(Y)
space = 200 + nlag # 200 days rolling window estimation
CV = array(NA, c(k, k, (t-space)))
colnames(CV) = rownames(CV) = colnames(Y)
for (i in 1:(t-space)) {
  var = VAR(Y[i:(space+i-1),], p=nlag, type="const")
  CV[,,i] = ofevd(fevd(var, n.ahead=nfore), n.ahead=nfore)$fevd
  if (i%%500==0) {print(i)}
}

to = matrix(NA, ncol=k, nrow=(t-space))
from = matrix(NA, ncol=k, nrow=(t-space))
net = matrix(NA, ncol=k, nrow=(t-space))
ct = npso = array(NA, c(k, k, (t-space)))
total = matrix(NA, ncol=1, nrow=(t-space))
colnames(npso)=rownames(npso)=colnames(ct)=rownames(ct)=colnames(Y)
for (i in 1:(t-space)){
  vd = DCA(CV[,,i])
  ct[,,i] = vd$CT
  to[i,] = vd$TO/k
  from[i,] = vd$FROM/k
  net[i,] = vd$NET/k
  npso[,,i] = vd$NPSO/k
  total[i,] = vd$TCI
}

nps = array(NA,c((t-space),k/2*(k-1)))
colnames(nps) = 1:ncol(nps)
jk = 1
for (i in 1:k) {
  for (j in 1:k) {
    if (j<=i) {
      next
    } else {
      nps[,jk] = npso[j,i,]
      colnames(nps)[jk] = paste0(colnames(Y)[i],"-",colnames(Y)[j])
      jk = jk + 1
    }
  }
}

### DYNAMIC TOTAL CONNECTEDNESS
date = DATE[-c(1:space)]
par(mfrow = c(1,1), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1,1) + .05, mgp = c(0, 0.1, 0))
plot(date, total, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(floor(min(total)),ceiling(max(total))),yaxs="i",xlab="",tck=0.01)
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

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
print(DCA(ct/100)$ALL)

### END
