
### ANTONAKAKIS, N., CHATZIANTONIOU, I., AND GABAUER, D. (2020)
### REFINED MEASURES OF DYNAMIC CONNECTEDNESS BASED ON TIME-VARYING PARAMETERS VECTOR AUTREGRESSIONS
### Journal of Risk and Financial Management
### by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)
library("openxlsx")
library("parallel")
options(mc.cores=detectCores())
source("functions.R")

DATA = read.xlsx('./dy2012.xlsx')
date = as.Date(as.character(DATA[,1]))
Y = DATA[,-1]
k = ncol(Y)

### TVP-VAR
nlag = 4 # VAR(4)
nfore = 10 # 10-step ahead forecast
prior = UninformativePrior(0.1, k, nlag)
tvpvar = TVPVAR(Y, l=c(0.99, 0.99), nlag, prior)
B_t = tvpvar$beta_t
Q_t = tvpvar$Q_t

### DYNAMIC CONNECTEDNESS APPROACH
t = nrow(Y)
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
  to[i,] = vd$TO
  from[i,] = vd$FROM
  net[i,] = vd$NET
  npso[,,i] = vd$NPSO
  total[i,] = vd$TCI
  if (i%%100==0) print(paste0(round(100*i/t0,2),"%"))
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
      colnames(nps)[jk] = paste0(colnames(Y)[j],"-",colnames(Y)[i])
      jk = jk + 1
    }
  }
}

### DYNAMIC TOTAL CONNECTEDNESS
par(mfrow = c(1,1), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1,1) + .05, mgp = c(0, 0.1, 0))
plot(date,total, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(floor(min(total)),ceiling(max(total))),yaxs="i",xlab="",tck=0.01)
grid(NA,NULL,lty=1)
polygon(c(date,rev(date)),c(c(rep(0,nrow(total))),rev(total)),col="grey20", border="grey20")
box()

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
par(mfrow = c(ceiling(k/2),2), oma = c(0,1,0,0) + 0.02, mar = c(1,1,1,1) + .02, mgp = c(0, 0.1, 0))
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
