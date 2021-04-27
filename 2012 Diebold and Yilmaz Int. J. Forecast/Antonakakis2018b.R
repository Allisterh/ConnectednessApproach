
### Dynamic Connectedness of UK Regional Property Returns. 
### ANTONAKAKIS, N., CHATZIANTONIOU, I., FLOROS, C., & GABAUER, D. (2018)
### URBAN STUIES

library("igraph")
library("openxlsx")
library("animation")
library("RColorBrewer")
library("parallel")
options(mc.cores=detectCores())
source("functions.R")

colors = c("black","tan4","springgreen4","springgreen2","steelblue4","steelblue1","maroon4",
           "maroon1","orangered4","orangered","tan1","yellow3","firebrick")
palette(colors)

DATA = read.xlsx("./Antonakakis2018b.xlsx", detectDates=TRUE)
RAW = DATA[,-c(1,ncol(DATA))]
DATE = as.Date(DATA[,1])
NAMES = colnames(RAW)
k = ncol(RAW)

# Figure 1: UK regional property price indices
par(mfcol=c(1,1), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(DATE, RAW[,1], type="l", xlab="", ylab="", main="", xaxs="i", yaxs="i", las=1, ylim=c(0,max(RAW)), tck=-0.01)
grid(NA, NULL)
for (i in 1:k) {
  lines(DATE, DATA[,i], col=i)
}
box()

delta = 4
Y = matrix(NA, nrow=(nrow(RAW)-delta), ncol=k)
for (i in 1:k) {
  Y[,i] = na.omit(diff(log(RAW[,i]), delta))
}
colnames(Y) = NAMES
date = DATE[-c(1:delta)]
t = nrow(Y)

# Figure 2: UK regional property returns
split = 3
par(mfrow=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 1:k) {
  plot(date, Y[,j], type="l", xlab="", ylab="", main=NAMES[j], xaxs="i", las=1, col="steelblue4", ylim=c(min(Y),max(Y)), tck=-0.05)
  grid(NA, NULL)
  lines(date, Y[,j], col="steelblue4")
  abline(h=0, lty=3)
  box()
}

#### ROLLING-WINDOW VAR ----
nlag = 2
nfore = 4
space = 60 + nlag

t0 = t-space
total = matrix(NA, ncol=1, nrow=t0)
gfevd = ct = npso = array(NA, c(k, k, t0))
net = from = to = matrix(NA, ncol=k, nrow=t0)
colnames(gfevd)=rownames(gfevd)=colnames(ct)=rownames(ct)=NAMES
for (i in 1:t0) {
  var = VAR(Y[i:(space+i-1),], p=nlag)
  gfevd[,,i] = GFEVD(var$Phi, var$Sigma, n.ahead=nfore)$GFEVD
  dca = DCA(gfevd[,,i])
  ct[,,i] = dca$CT
  to[i,] = dca$TO
  from[i,] = dca$FROM
  net[i,] = dca$NET
  npso[,,i] = dca$NPSO
  total[i,] = dca$TCI
  if (i%%100==0) print(paste0(round(100*i/t0,2),"%"))
}

### DYNAMIC TOTAL CONNECTEDNESS
date = DATE[-c(1:(space+delta))]
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date, total, type="h",col="steelblue4", las=1, main="",ylab="",ylim=c(70, 90),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)
lines(date,total,type="h",col="steelblue4", lwd=2)
box()

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
split = 4
par(mfrow=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date, to[,i], type="h", xlab="",ylab="",xaxs="i",col="steelblue4", las=1, main=paste(colnames(Y)[i],"TO all others"),ylim=c(0,ceiling(max(to))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  lines(date,to[,i],type="h",col="steelblue4", lwd=2)
  box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfrow=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,from[,i], xlab="",ylab="",type="h",xaxs="i",col="steelblue4", las=1, main=paste(colnames(Y)[i],"FROM all others"),ylim=c(0,ceiling(max(from))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  lines(date,from[,i],type="h",col="steelblue4")
  box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
par(mfrow=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,net[,i], xlab="",ylab="",type="h",xaxs="i",col="steelblue4", las=1, main=paste("NET",colnames(Y)[i]),ylim=c(floor(min(net)),ceiling(max(net))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  abline(h=0, lty=2)
  lines(date,net[,i],type="h",col="steelblue4")
  box()
}

### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
i = 9
print(NAMES[i])
par(mfcol=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 1:k) {
  if (i!=j) {
    plot(date,npso[j,i,], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso)),ceiling(max(npso))))
    grid(NA,NULL)
    polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(npso[j,i,])),col="steelblue4", border="steelblue4")
    box()
  }
}

### STATIC CONNECTEDNESS TABLE
var_sample = VAR(Y, p=nlag)
gfevd_sample = GFEVD(var_sample$Phi, var_sample$Sigma, n.ahead=nfore)$GFEVD
colnames(gfevd_sample)=rownames(gfevd_sample)=NAMES
print(DCA(gfevd_sample)$TABLE)

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
print(DCA(gfevd)$TABLE)

### END
