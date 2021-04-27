
### GAMBA-SANTAMARIA, S., GOMEZ-GONZALEZ, J.E., HURTADO-GUARIN, J.L. AND MELO-VELANDIA, L.F. (2017)
### STOCK MARKET VOLATILITY SPILLOVERS: EVIDENCE FOR LATIN AMERICA
### Finance Research Letters
### by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)
library("MTS")
library("rmgarch")
library("openxlsx")
library("parallel")
options(mc.cores=detectCores())
source("functions.R")

DATA = read.xlsx('./dy2012.xlsx')
DATE = as.Date(as.character(DATA[,1]))
Y = DATA[,-1]
k = ncol(Y)
NAMES = colnames(Y)

### DYNAMIC CONNECTEDNESS APPROACH
nlag = 4   # VAR(4)
nfore = 10 # 10-step ahead forecast
var = VAR(Y, p=nlag)
B = var$Phi
t = nrow(var$residuals)
ugarch = ugarchspec(variance.model=list(garchOrder=c(1, 1), model="sGARCH"),
                    mean.model=list(armaOrder=c(0, 0)))
mgarch = multispec( replicate(k, ugarch) )
dccgarch_spec = cgarchspec(uspec = mgarch, dccOrder=c(1,1), asymmetric = FALSE,
                         distribution.model = list(copula = "mvt", method = "Kendall", time.varying = TRUE, transformation = "parametric"))
dcc_fit = cgarchfit(dccgarch_spec, data=var$residuals, solver="solnp", fit.control=list(eval.se = TRUE) )
Q_t = rcov(dcc_fit)

### DYNAMIC CONNECTEDNESS APPROACH
total = matrix(NA, ncol=1, nrow=t)
gfevd = ct = npso = array(NA, c(k, k, t))
net = from = to = matrix(NA, ncol=k, nrow=t)
colnames(gfevd)=rownames(gfevd)=colnames(ct)=rownames(ct)=NAMES
for (i in 1:t){
  gfevd[,,i] = GFEVD(B, Q_t[,,i], n.ahead=nfore)$GFEVD
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
date = DATE[-c(1:nlag)]
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date,total, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(floor(min(total)),ceiling(max(total))),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(0,t)),rev(total)),col="grey20", border="grey20")
box()

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

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
print(DCA(gfevd)$TABLE)

### END
