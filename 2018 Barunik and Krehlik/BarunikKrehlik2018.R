
### BARUNIK AND KREHLIK (2018)
### Measuring the Frequency Dynamics of Financial Connectedness and Systemic Risk 
### Journal of Financial Econometrics
### replicated by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details). Code derived from original frequencyConnectedness package

library("MTS")
library("frequencyConnectedness")
library("parallel")
library("openxlsx")
options(mc.cores=detectCores())
source("functions.R")

DATA = read.xlsx('./dy2012.xlsx')
DATE = as.Date(as.character(DATA[,1]))
Y = DATA[,-1]
k = ncol(Y)
NAMES = colnames(Y)

### STATIC CONNECTEDNESS APPROACH
nlag = 4 # VAR(4)
nfore = 100 # 100-step ahead forecast
partition = c(pi+0.00001, pi/5, pi/30, 0)

# Show that both functions are identical: 
var_full = vars::VAR(Y, p=nlag, type="const")
spilloverBK12(var_full, n.ahead=nfore, no.corr=F, partition)

# Advantages: Play with coefficients, use TVP-VAR, LASSO, Ridge, Elastic Net output, etc.
mts = MTSVAR(Y, p=nlag)
MTSspilloverBK12(mts$Phi, mts$Sigma, n.ahead=nfore, no.corr=F, partition)

### DYNAMIC CONNECTEDNESS APPROACH
t = nrow(Y)
space = 200+nlag
t0 = t-space
frequencies = length(partition)-1
net = from = to = array(NA, c(t0,k,frequencies))
npso = array(NA, c(t0, k*(k-1)/2, frequencies))
total = matrix(NA, nrow=t0, ncol=frequencies)
for (i in 1:t0) {
  var = MTSVAR(Y[i:(space+i-1),], p=nlag)
  DCA = MTSspilloverBK12(var$Phi, var$Sigma, n.ahead=nfore, no.corr=F, partition)
  TOTAL = overall(DCA, within=FALSE)
  TO = to(DCA, within=FALSE)
  FROM = from(DCA, within=FALSE)
  NET = net(DCA, within=FALSE)
  NPSO = pairwise(DCA, within=FALSE)
  for (j in 1:frequencies) {
    total[i,j] = TOTAL[[j]]
    to[i,,j] = TO[[j]]
    from[i,,j] = FROM[[j]]
    net[i,,j] = NET[[j]]
    npso[i,,j] = NPSO[[j]]
  }
  if (i%%100==0) print(paste0(round(100*i/t0,2),"%"))
}
date = DATE[-c(1:(space))]

sp = spilloverRollingBK12(Y, n.ahead=nfore, no.corr=F, func_est="VAR", params_est=list(p=nlag, type="const"), window=space, partition)
par(mfcol=c(frequencies,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plotOverall(sp)

### DYNAMIC TOTAL CONNECTEDNESS
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date, apply(total,1,sum),type="l",las=1,xaxs="i",ylim=c(0,max(apply(total,1,sum))),yaxs="i",tck=-0.02,main="",xlab="",ylab="")
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(apply(total,1,sum))),col="grey20", border="grey20")
for (i in 1:frequencies) {
  lines(date, total[,i],col=i+1)
}
box()

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
split = 2
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 1:k) {
  plot(date, apply(to[,j,],1,mean),type="l",las=1,xaxs="i",ylim=c(0,max(to)),yaxs="i",tck=-0.02,main=paste(NAMES[j],"TO Others"),xlab="",ylab="")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(apply(to[,j,],1,mean))),col="grey20", border="grey20")
  for (i in 1:frequencies) {
    lines(date, to[,j,i],col=i+1)
  }
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 1:k) {
  plot(date,apply(from[,j,],1,mean),type="h",las=1,xaxs="i",ylim=c(0,max(from)),yaxs="i",tck=-0.02,main=paste(NAMES[j],"FROM Others"),xlab="",ylab="")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(apply(to[,j,],1,mean))),col="grey20", border="grey20")
  for (i in 1:frequencies) {
    lines(date,from[,j,i],col=i+1)
  }
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 1:k) {
  plot(date,apply(net[,j,],1,mean),type="h",las=1,xaxs="i",ylim=c(min(net),max(net)),yaxs="i",tck=-0.02,main=paste("NET",NAMES[j]),xlab="",ylab="")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(apply(net[,j,],1,mean))),col="grey20", border="grey20")
  for (i in 1:frequencies) {
    lines(date,net[,j,i],col=i+1)
  }
}

### END
