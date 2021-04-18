
### ANDO, T., GREENWOOD-NIMMO, M., AND SHIN, Y. (2020)
### QUANTILE CONNECTEDNESS: MODELLING TAIL BEHAVIOUR IN THE TOPOLOGY OF FINANCIAL NETWORKS
### Management Science
### replicated by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)
library("openxlsx")
library("parallel")
options(mc.cores=detectCores())
source("functions.R")

DATA = read.xlsx('./dy2012.xlsx')
DATE = as.Date(as.character(DATA[,1]))
Y = DATA[,-1]
k = ncol(Y)
NAMES = colnames(Y)

### STATIC CONNECTEDNESS APPROACH
tau = 0.5
nlag = 4 # VAR(4)
nfore = 10 # 10-step ahead forecast
qvar_full = QVAR(Y, p=nlag, tau=tau)
CV_full = GFEVD(qvar_full$B, qvar_full$Q, n.ahead=nfore)$GFEVD
rownames(CV_full)=colnames(CV_full)=NAMES
print(DCA(CV_full))

### DYNAMIC CONNECTEDNESS APPROACH
t = nrow(Y)
space = 200 + nlag # 200 days rolling window estimation
t0 = t-space

total = matrix(NA, ncol=1, nrow=t0)
gfevd = ct = npso = array(NA, c(k, k, t0))
net = from = to = matrix(NA, ncol=k, nrow=t0)
colnames(gfevd)=rownames(gfevd)=colnames(ct)=rownames(ct)=colnames(Y)
for (i in 1:(t-space)){
  qvar = QVAR(Y[i:(space+i-1),], p=nlag, tau=tau)
  gfevd[,,i] = GFEVD(qvar$B, qvar$Q, n.ahead=nfore)$GFEVD
  vd = DCA(gfevd[,,i])
  ct[,,i] = vd$CT
  to[i,] = vd$TO
  from[i,] = vd$FROM
  net[i,] = vd$NET
  npso[,,i] = vd$NPSO
  total[i,] = vd$TCI
  if (i%%100==0) print(paste0(round(100*i/t0,2),"%"))
}

### DYNAMIC TOTAL CONNECTEDNESS
date = DATE[-c(1:space)]
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date, total, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(floor(min(total)),ceiling(max(total))),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(total)),col="grey20", border="grey20")
box()

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
split = 2
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,to[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"TO all others"),ylim=c(floor(min(to)),ceiling(max(to))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(to[,i])),col="grey20", border="grey20")
  box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,from[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"FROM all others"),ylim=c(floor(min(from)),ceiling(max(from))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(from[,i])),col="grey20", border="grey20")
  box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,net[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste("NET",colnames(Y)[i]),ylim=c(floor(min(net)),ceiling(max(net))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(net[,i])),col="grey20", border="grey20")
  box()
}

### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
kk = k*(k-1)/2
lk = ceiling(sqrt(2))
par(mfcol=c(ceiling(kk/lk),lk), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k) {
  for (j in 1:k) {
    if (i<j) {
      plot(date,npso[j,i,], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso)),ceiling(max(npso))))
      grid(NA,NULL)
      polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(npso[j,i,])),col="grey20", border="grey20")
      box()
    }
  }
}

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
print(DCA(gfevd)$TABLE)

### END
