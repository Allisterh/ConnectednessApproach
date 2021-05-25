
library("openxlsx")
library("parallel")
library("RColorBrewer")
options(mc.cores=detectCores())
source("functions.R")

colors=c("black","grey30","grey20","grey80","springgreen4","springgreen2","steelblue4","steelblue1","maroon4","maroon1","orangered4","orangered","tan1","yellow3")
palette(colors)

DATA = read.xlsx('./ChatziantoniouGabauer2021.xlsx', detectDates=TRUE)
DATE = as.Date(DATA[,1],"%d.%m.%y")
RAW = na.omit(DATA[,-1])
NAMES = colnames(RAW)
k = ncol(RAW)

par(mfrow = c(1,1), oma = c(0.5,0.5,0,0) + 0.1, mar = c(0.5,1,1,0) + .5, mgp = c(3, 0.5, 0))
plot(DATE,RAW[,1], col=1,lwd=1, type="l", main="",xaxs="i", las=1,yaxs="i",ylim=c(min(RAW),max(RAW)))
grid()
for (i in 1:k){
  lines(DATE,RAW[,i],col=i)
}
legend("topleft", NAMES, fill=1:k)

Y = na.omit(diff(as.matrix(RAW)))
date = DATE[-1]

par(mfrow = c(4,3), oma = c(0.5,0.5,0,0) + 0.1, mar = c(0.5,1,1,0) + .5, mgp = c(3, 0.5, 0))
for (i in 1:k){
  plot(date,Y[,i], col="grey20",lwd=0.01, type="l", main=NAMES[i], ylim=c(-1.1,1.1),xaxs="i", las=1,yaxs="i")
  grid(NA,NULL,col="grey20")
  lines(date,Y[,i],col="steelblue4",lwd=0.5)
  box()
}

### TVP-VAR Estimation
nlag = 1
nfore = 20
prior = BayesPrior(Y[1:200,], nlag)
tvpvar = TVPVAR(Y, l=c(0.99, 0.99), nlag, prior)
B_t = tvpvar$beta_t
Q_t = tvpvar$Q_t

### DYNAMIC CONNECTEDNESS APPROACH
t = nrow(Y)
total = total_corr = matrix(NA, ncol=1, nrow=t)
gfevd = ct = npso = pci = influence = array(NA, c(k, k, t))
net = to = from = matrix(NA, ncol=k, nrow=t)
colnames(gfevd)=rownames(gfevd)=colnames(influence)=rownames(influence)=colnames(pci)=rownames(pci)=NAMES
for (i in 1:t){
  gfevd[,,i] = GFEVD(B_t[,,i], Q_t[,,i], n.ahead=nfore)$GFEVD
  dca = DCA(gfevd[,,i])
  ct[,,i] = dca$CT
  pci[,,i] = dca$PCI
  to[i,] = dca$TO
  from[i,] = dca$FROM
  net[i,] = dca$NET
  npso[,,i] = dca$NPSO
  influence[,,i] = dca$INFLUENCE
  pci[,,i] = dca$PCI
  total[i,] = dca$TCI
  total_corr[i,] = dca$TCI_corrected
  if (i%%100==0) print(paste0(round(100*i/t,2),"%"))
}

### DYNAMIC TOTAL CONNECTEDNESS
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date,total_corr, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(floor(min(c(total,total_corr))),ceiling(max(c(total,total_corr)))),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(0,t)),rev(total_corr)),col="grey20", border="grey20")
lines(date,total,col="red")
box()
legend("topright", c("TCI corrected","TCI"), fill=c("grey20", "red"))

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
split = 2
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,to[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(NAMES[i],"TO all others"),ylim=c(0,ceiling(max(to))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(to[,i])),col="grey20", border="grey20")
  box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,from[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(NAMES[i],"FROM all others"),ylim=c(0,100),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(from[,i])),col="grey20", border="grey20")
  box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,net[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste("NET",NAMES[i]),ylim=c(floor(min(net)),ceiling(max(net))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(net[,i])),col="grey20", border="grey20")
  box()
}

### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
i = 5
par(mfcol=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 1:k) {
  if (i!=j) {
    plot(date,npso[j,i,], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso)),ceiling(max(npso))))
    grid(NA,NULL)
    polygon(c(date,rev(date)),c(c(rep(0,t)),rev(npso[j,i,])),col="grey20", border="grey20")
    box()
  }
}

### PAIRWISE CONNECTEDNESS INDEX
par(mfrow=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 1:k) {
  if (i!=j) {
    plot(date,pci[j,i,], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(0,1))
    grid(NA,NULL)
    polygon(c(date,rev(date)),c(c(rep(0,t)),rev(pci[j,i,])),col="grey20", border="grey20")
    box()
  }
}

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
print(DCA(gfevd)$TABLE)


### RANKING OF PCI AND ABSOLUTE PAIRWISE INFLUENCE INDEX
l = 0
kk = k*(k-1)/2
influence_upper = pci_upper = array(NA, c(t,kk))
name_matrix = matrix(NA,k,k)
for (i in 1:k) {
  for (j in 1:k) {
    if (i>j) {
      l = l+1
      pci_upper[,l] = pci[j,i,]
      influence_upper[,l] = influence[j,i,]
      name_matrix[j,i] = paste0(NAMES[j],"-",NAMES[i])
    }
  }
}
lownames = name_matrix[upper.tri(name_matrix)]
colnames(influence_upper)=colnames(pci_upper)=lownames

rep = 1000
influence_mean = pci_mean = matrix(NA,nrow=rep,ncol=kk)
colnames(influence_mean)=colnames(pci_mean)=lownames
for (i in 1:rep) {
  pci_mean[i,] = apply(pci_upper[sample(1:t,replace=TRUE),],2,mean)
  influence_mean[i,] = apply(influence_upper[sample(1:t,replace=TRUE),],2,mean)
}

pci_ranking = names(sort(apply(pci_upper,2,mean),decreasing=T))
print(cbind(1:kk,pci_ranking))
par(mfrow = c(1,2), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1,1) + .05, mgp = c(0, 0.1, 0))
plot(sort(apply(pci_mean,2,mean), TRUE),xaxs="i",las=1,col="steelblue4",yaxs="i",xlab="",ylab="",tck=0.01,type="s",ylim=c(0,1))
grid()
lines(sort(apply(pci_mean,2,quantile,0.025), TRUE),type="s",pch="-",col="steelblue4",lty=3)
lines(sort(apply(pci_mean,2,quantile,0.975), TRUE),type="s",pch="-",col="steelblue4",lty=3)
box()
legend("topright", c("PCI"), fill=c("steelblue4"))

influence_ranking = names(sort(apply(influence_upper,2,mean)))
print(cbind(1:kk,influence_ranking))
plot(sort(apply(influence_mean,2,mean)),xaxs="i",las=1,col="steelblue4",yaxs="i",xlab="",ylab="",tck=0.01,type="s",ylim=c(0,0.5))
grid()
lines(sort(apply(influence_mean,2,quantile,0.025)),type="s",pch="-",col="steelblue4",lty=3)
lines(sort(apply(influence_mean,2,quantile,0.975)),type="s",pch="-",col="steelblue4",lty=3)
box()
legend("topright", c("Influence"), fill=c("steelblue4"))

# END
