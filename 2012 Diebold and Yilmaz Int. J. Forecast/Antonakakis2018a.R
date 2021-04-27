
library("rmgarch")
library("moments")
library("openxlsx")
library("RColorBrewer")
source("functions.R")
colors=c("black","tan4","springgreen4","springgreen2","steelblue4","steelblue1","maroon4",
         "maroon1","orangered4","orangered","tan1","yellow3","firebrick")
palette(colors)

RAW = read.xlsx("./Antonakakis2018a.xlsx", detectDates=TRUE)
DATE = as.Date(RAW[,1])
DATA = RAW[,-1]
k = ncol(DATA)
NAMES = colnames(DATA)
head(DATA)

par(mfrow=c(1,1), oma=c(1,1,0,0) + 0.05, mar=c(1,1,1,1) + 0.3, mgp=c(0, 0.5, 0))
plot(DATE, scale(DATA[,1], TRUE, TRUE),type="l",ylim=c(-3, 4), las=1, xaxs="i", xlab="", ylab="", main="")
grid(NA, NULL)
for (i in k:1) {
  lines(DATE, scale(DATA[,i], TRUE, TRUE), col=i, lty=1)
}
box()
legend("topleft", legend=NAMES, fill=colors, ncol=1, cex=0.5, bty="n")

RETURN = DATA[-1,]
for (i in 1:k) {
  RETURN[,i] = 100*abs(diff(log(DATA[,i])))
}

# Figure 1: Oil and stock price volatility (absolute returns
date = DATE[-1]
par(mfrow=c(ceiling(k/3),3), oma=c(1,1,0,0) + 0.05, mar=c(1,1,1,1) + 0.3, mgp=c(0, 0.5, 0))
for (j in 1:k){
  plot(date, RETURN[,j], type="l", col="steelblue4", main=NAMES[j], xlab="", ylab="", las=1, xaxs="i",yaxs="i",tck=.01)
  grid(NA, NULL)
  lines(date, RETURN[,j], col="steelblue4")
  box()
}

#### DYNAMIC CONDITIONAL CORRELATION ----
u_spec = ugarchspec(variance.model=list(garchOrder=c(1,1), model="sGARCH"), 
                    mean.model = list(armaOrder=c(0, 0)),
                    distribution.model="norm")
m_spec = multispec(c( replicate(k, u_spec)))
dcc_spec = dccspec(uspec=m_spec, dccOrder=c(1, 1), model="DCC", distribution="mvt")
dcc_fit = dccfit(dcc_spec, data=RETURN, solver="solnp", fit.control=list(eval.se=TRUE))
H = rcov(dcc_fit)
R = rcor(dcc_fit)

# TESTING DYNAMIC CONDITIONAL CORRELATION
dcc_test = DCCtest(RETURN, garchOrder = c(1,1), solver = "solnp", n.lags=1); dcc_test

# VISUALIZE CONDITIONAL CORRELARIONS
split = 4
par(mfcol=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 2:k) {
  plot(date, R[1,j,], type="l", col="steelblue4", tck=-0.02,ylab="",xlab="",
       main=paste0("WTI", "/", NAMES[j]), las=1, ylim=c(floor(min(R[1,,])),max(R[1,-1,])), xaxs="i")
  grid(NA, NULL)
  lines(date, R[1,j,], col="steelblue4")
  box()
}

#### ROLLING-WINDOW VAR ----
nlag = 1
nfore = 30
space = 500 + nlag
t0 = nrow(RETURN) - space
gfevd = npso = array(NA, c(k, k, t0))
rownames(gfevd) = colnames(gfevd) = NAMES
total = matrix(NA, ncol=1, nrow=t0)
net = to = from = matrix(NA, ncol=k, nrow=t0)
for (i in 1:t0) {
  var = VAR(RETURN[i:(space+i-1),], p=nlag)
  gfevd[,,i] = GFEVD(var$Phi, var$Sigma, n.ahead=nfore)$GFEVD
  dca = DCA(gfevd[,,i])
  total[i,] = dca$TCI
  to[i,] = dca$TO
  from[i,] = dca$FROM
  net[i,] = dca$NET
  npso[,,i] = dca$NPSO
  if (i%%100==0) print(paste0(round(100*i/t0,2),"%"))
}

### DYNAMIC TOTAL CONNECTEDNESS
date = DATE[-c(1:(space+1))]
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date, total, type="l",xaxs="i",col="steelblue4", las=1, main="",ylab="",ylim=c(floor(min(total)),ceiling(max(total))),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(total)),col="steelblue4", border="steelblue4")
box()

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
split = 2
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,to[,i], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste(NAMES[i],"TO all others"),ylim=c(floor(min(to)),ceiling(max(to))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(to[,i])),col="steelblue4", border="steelblue4")
  box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,from[,i], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste(NAMES[i],"FROM all others"),ylim=c(floor(min(from)),ceiling(max(from))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(from[,i])),col="steelblue4", border="steelblue4")
  box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,net[,i], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=NAMES[i],ylim=c(floor(min(net)),ceiling(max(net))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(net[,i])),col="steelblue4", border="steelblue4")
  box()
}

### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
i = 1
par(mfcol=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 1:k) {
  if (i<j) {
    plot(date,npso[j,i,], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso)),ceiling(max(npso))))
    grid(NA,NULL)
    polygon(c(date,rev(date)),c(c(rep(0,t0)),rev(npso[j,i,])),col="steelblue4", border="steelblue4")
    abline(h=0,lty=2)
    box()
  }
}

### STATIC CONNECTEDNESS TABLE
var_sample = VAR(RETURN, p=nlag)
gfevd_sample = GFEVD(var_sample$Phi, var_sample$Sigma, n.ahead=nfore)$GFEVD
colnames(gfevd_sample)=rownames(gfevd_sample)=NAMES
print(DCA(gfevd_sample)$TABLE)

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
print(DCA(gfevd)$TABLE)

### END
