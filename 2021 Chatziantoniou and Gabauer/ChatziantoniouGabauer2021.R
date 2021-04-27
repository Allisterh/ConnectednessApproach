
### CHATZIANTONIOU, I., AND GABAUER, D. (2019)
### EMU-RISK SYNCHRONISATION AND FINANCIAL FRAGILITY THROUGH THE PRISM OF DYNAMIC CONNECTEDNESS
### by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)
library("openxlsx")
library("parallel")
options(mc.cores=detectCores())
source("functions.R")

DATA = read.xlsx('./dy2012.xlsx')
date = as.Date(as.character(DATA[,1]))
Y = DATA[,-1]
k = ncol(Y)
NAMES = colnames(Y)

### TVP-VAR
nlag = 4 # VAR(4)
nfore = 10 # 10-step ahead forecast
t = nrow(Y)

prior = MinnesotaPrior(0.1, k, nlag)
tvpvar = TVPVAR(Y, l=c(0.99, 0.99), nlag, prior)
B_t = tvpvar$beta_t
Q_t = tvpvar$Q_t

### DYNAMIC CONNECTEDNESS APPROACH
total = total_corr = matrix(NA, ncol=1, nrow=t)
gfevd = ct = npso = pci = array(NA, c(k, k, t))
net = to = from = matrix(NA, ncol=k, nrow=t)
colnames(gfevd)=rownames(gfevd)=colnames(ct)=rownames(ct)=colnames(Y)
for (i in 1:t){
  gfevd[,,i] = GFEVD(B_t[,,i], Q_t[,,i], n.ahead=nfore)$GFEVD
  CV = gfevd[,,i]
  dca = DCA(gfevd[,,i])
  ct[,,i] = dca$CT
  pci[,,i] = dca$PCI
  to[i,] = dca$TO
  from[i,] = dca$FROM
  net[i,] = dca$NET
  npso[,,i] = dca$NPSO
  pci[,,i] = dca$PCI
  total[i,] = dca$TCI
  total_corr[i,] = dca$TCI_corrected
  if (i%%100==0) print(paste0(round(100*i/t0,2),"%"))
}

### DYNAMIC TOTAL CONNECTEDNESS
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date,total_corr, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(floor(min(c(total,total_corr))),ceiling(max(c(total,total_corr)))),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(0,t)),rev(total_corr)),col="grey20", border="grey20")
lines(date,total,col="red")
box()
legend("topleft", c("TCI corrected","TCI"), fill=c("grey20", "red"))

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
split = 2
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,to[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"TO all others"),ylim=c(floor(min(to)),ceiling(max(to))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(to[,i])),col="grey20", border="grey20")
  box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,from[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"FROM all others"),ylim=c(floor(min(from)),ceiling(max(from))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(from[,i])),col="grey20", border="grey20")
  box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,net[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste("NET",colnames(Y)[i]),ylim=c(floor(min(net)),ceiling(max(net))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(net[,i])),col="grey20", border="grey20")
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
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(npso[j,i,])),col="grey20", border="grey20")
      box()
    }
  }
}

### PAIRWISE CONNECTEDNESS INDEX
par(mfcol=c(ceiling(kk/lk),lk), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k) {
  for (j in 1:k) {
    if (i<j) {
      plot(date,pci[j,i,], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(0,1))
      grid(NA,NULL)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(pci[j,i,])),col="grey20", border="grey20")
      box()
    }
  }
}

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
print(DCA(gfevd)$TABLE)

### END
