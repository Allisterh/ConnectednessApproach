
### GABAUER, D. (2020)
### Crude Oil Futures Contracts and Commodity Markets: New Evidence from a TVP-VAR Extended Joint Connectedness Approach
### JOURNAL OF FORECASTING
### replicated by David Gabauer
library("rmgarch")
source("functions.R")

DATA = read.xlsx('./gabauer2020.xlsx', detectDates=TRUE)
DATE = as.Date(as.character(DATA[,1]))
RAW = DATA[,-1]
NAMES = colnames(RAW)
t = nrow(RAW)
k = ncol(RAW)

split = 2
par(mfcol=c(1,1), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(DATE,scale(RAW[,1],TRUE,TRUE),type="l",xaxs="i",las=1,xlab="",ylab="",main="",tck=-0.02,col="grey20",ylim=c(min(scale(RAW)),max(scale(RAW))))
for (i in 1:k) {
  lines(DATE,scale(RAW[,i],TRUE,TRUE), col=i)
  grid(NA,NULL)
}
box()
legend("topleft", NAMES, fill=1:k)

Y = na.omit(100*diff(log(as.matrix(Y))))
date = DATE[-1]
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k) {
  plot(date,Y[,i],type="l",xaxs="i",las=1,xlab="",ylab="",main=NAMES[i],tck=-0.02,col="grey20",ylim=c(min(Y),max(Y)))
  grid(NA,NULL)
  lines(date,Y[,i], col="steelblue4")
  box()
}

summary_statistics = Moments(Y)
print(summary_statistics)
t = nrow(Y)

### Estimate DCC-GARCH(1,1) ----
garch11.spec = ugarchspec(mean.model=list(armaOrder=c(0,0)), 
                          variance.model=list(garchOrder=c(1,1), model="sGARCH"), 
                          distribution.model="norm")
dcc.garch11.spec = dccspec(uspec=multispec(replicate(k, garch11.spec)),
                           dccOrder=c(1,1), distribution="mvnorm")
dcc.fit = dccfit(dcc.garch11.spec, data=Y)
dcc.test = DCCtest(Y, garchOrder=c(1,1), solver="solnp", n.lags=5)

### GENERALISED VOLATILITY FORECAST ERROR VARIANCE DECOMPOSITION
nfore = 100
dcc_gfevd = VGFEVD(dcc.fit, nfore)
virf = dcc_gfevd$GVIRF

split = 2
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1)+0.5, mgp=c(0.5,0.5,0))
for (i in 1:k) {
  for (j in 1:k) {
    if (j==i) {
      plot(virf[i,j,1,-1],type="l",xaxs="i",las=1,xlab="",ylab="",col="grey20",tck=-0.02,main=paste(NAMES[i]),yaxs="i",ylim=c(floor(min(virf[i,j,1,])),ceiling(max(virf[i,j,1,]))))
      for (ij in 1:nfore) {
        lines(virf[i,j,ij,-1],col="grey20")
      }
      abline(h=0)
    }
  }
}

kk = k*(k-1)/2
split = 2
par(mfcol=c(ceiling(kk/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1)+0.5, mgp=c(0.5,0.5,0))
for (i in 1:k) {
  for (j in 1:k) {
    if (j>i) { 
      plot(virf[i,j,1,-1],type="l",xaxs="i",las=1,xlab="",ylab="",col="grey20",tck=-0.02,
           main=paste(NAMES[i],"-",NAMES[j]),yaxs="i",
           ylim=c(min(virf[i,j,1,-1]),max(virf[i,j,1,-1])))
      for (ij in 1:nfore) {
        lines(virf[i,j,ij,-1],col=paste0("grey",round((dim(virf)[4]-ij)/dim(virf)[4]*100)))
      }
      abline(h=0)
      box()
    }
  }
}

### CONNECTEDNESS APPROACH
net = to = from = matrix(NA,nrow=t,ncol=k)
total = matrix(NA,ncol=1,nrow=t)
ct = npso = array(NA,c(k,k,t))
for (i in 1:t) {
  dca = DCA(dcc_gfevd$GVFEVD[,,i])
  ct[,,i] = dca$CT
  total[i,] = dca$TCI_corrected
  to[i,] = dca$TO
  from[i,] = dca$FROM
  net[i,] = dca$NET
  npso[,,i] = dca$NPSO
}

### DYNAMIC TOTAL CONNECTEDNESS
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date, total, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(floor(min(total)),ceiling(max(total))),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(0,t)),rev(total)),col="grey20", border="grey20")
box()

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
split = 2
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,to[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"TO all others"),ylim=c(0,ceiling(max(to))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(to[,i])),col="grey20", border="grey20")
  box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,from[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"FROM all others"),ylim=c(0,ceiling(max(from))),tck=-0.02,yaxs="i")
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

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
print(DCA(dcc_gfevd$GVFEVD)$TABLE)

### END
