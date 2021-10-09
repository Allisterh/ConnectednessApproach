
### The Impact of Euro Through Time: Exchange Rate Dynamics Under Different Regimes
### ANTONAKAKIS, N., CHATZIANZTONIOU, I., & GABAUER, D. (2020)
### INTERNATIONAL JOURNAL OF FINANCE AND ECONOMICS
### replicated by David Gabauer

library("abind")
library("rmgarch")
library("openxlsx")
library("parallel")
library("RColorBrewer")
library("WeightedPortTest")
library("psych")

source("functions.R")

options(warn=-1)
options("mc.cores"=detectCores())
colors = c("black","tan4","springgreen4","springgreen2","steelblue4","steelblue1","maroon4","maroon1","orangered4","orangered","tan1","yellow3")
palette(colors)

DATA = read.xlsx("./Antonakakis2020.xlsx",detectDates=TRUE)
DATE = as.Date(DATA[,1],"%d-%m-%Y")
RAW = DATA[,-1]
k = ncol(RAW)

### BASE ALL ON OTHER CURRENCY
colnames(RAW)[1] = c("EUR(DM)")
NAMES=colnames(RAW)

ERMI = which(DATE=="1979-03-13")
ERMIadj = which(DATE=="1993-08-02")
ERMII = which(DATE=="1998-12-31")
GFC = which(DATE=="2008-09-10")
t = nrow(RAW)
IND = c(0,ERMI,ERMIadj,ERMII,GFC,t)
date = DATE

### PRINCIPAL COMPONENT ANALYSIS
fa = principal(RAW, nfactor=1, method="tenBerge")
Y = cbind(RAW,fa$scores)
NAMES = c(NAMES,"ME")
colnames(Y) = NAMES
k = ncol(Y)

split = 2
par(mfcol=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:(k-1)){
   plot(date, Y[,k], type="l", main=NAMES[i], xaxs="i", las=1, col="orange",ylim=c(-5,5),xlab="",ylab="")
   grid(NA,NULL)
   lines(date, Y[,i], col="steelblue4")
   abline(h=0)
}

### SPLIT DATASET INTO DIFFERENT PERIODS
data_list = date_list = list()
period_names = c("PRE ERM1","ERM1","PRE ERM2","PRE CRISIS","CRISIS")
for (l in 1:length(period_names)) {
  data_list[[l]] = Y[(IND[l]+1):IND[l+1],]
  date_list[[l]] = date[(IND[l]+1):IND[l+1]]
}
names(data_list) = period_names
periods = length(data_list)

digit = 2
SummaryStatistics_list = list()
for (i in 1:periods) {
   SummaryStatistics_list[[i]] = rbind(SummaryStatistics(data_list[[i]]), format(round(cor(data_list[[i]]),digit),nsmall=digit))
}


##### DCC-COPULA ANALYSIS
lag = 20
prob = 0.05
conf.level = 0.90
R = H = NULL
copula_fit = GARCH_spec = evaluation_list = list()
for (i in 1:periods) {
   print(period_names[i])
   spec = c()
   evaluation_matrix = matrix(NA, ncol=k, nrow=10)
   colnames(evaluation_matrix) = NAMES
   rownames(evaluation_matrix) = c("VaR","","CVaR","","VaR Dur.","","Sign Bias","",paste0("WARCH(",lag,")"),"")
   data = data_list[[i]]
   for (j in 1:k) {
      print(NAMES[j])
      bestgarch = BestGARCH(data=data[,j],prob=prob,conf.level=conf.level,lag=lag,
                            distr=c("norm","snorm","std","sstd","ged","sged"), 
                            models=c("sGARCH","iGARCH","eGARCH","gjrGARCH","AVGARCH","TGARCH"))
      GARCH_selection = which(bestgarch$GARCH_IC==min(bestgarch$GARCH_IC),arr.ind=TRUE)
      print(bestgarch$GARCH_IC)
      ugarch.spec = bestgarch[[2]][[GARCH_selection[2]]][[GARCH_selection[1]]]
      ugarch.fit = ugarchfit(ugarch.spec,data=Y[,j], solver="hybrid")
      evaluation_matrix[,j] = rbind(ValueAtRisk(ugarch.fit,ugarch.spec,prob=prob, conf.level=conf.level)$statistic,
                                    SignBias_WARCH(ugarch.fit,lag=lag)$statistic)
      spec = c(spec,ugarch.spec)
   }
   print(evaluation_matrix)
   evaluation_list[[i]] = evaluation_matrix
   GARCH_spec[[i]] = spec
   
   mgarch.spec = cgarchspec(uspec=multispec(spec), dccOrder=c(1,1), asymmetric=FALSE,  
                            distribution.model=list(copula="mvt", method="Kendall", time.varying=TRUE, transformation="parametric"))
   copula_fit[[i]] = cgarchfit(mgarch.spec, data=data, solver=c("hybrid","solnp"), fit.control=list(eval.se=TRUE))
   H = abind(H,rcov(copula_fit[[i]]))
   R = abind(R,rcor(copula_fit[[i]]))
}
names(evaluation_list) = period_names

### DYNAMIC CONDITIONAL CORRELATIONS
t = length(date)
kk = k*(k-1)/2
split = 3
par(mfcol=c(ceiling(kk/split),split), oma=c(0.5,2,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k) {
   for (j in 1:k) {
      if (i<j) {
         plot(date,R[j,i,], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(min(R),1))
         grid(NA,NULL)
         lines(date,R[j,i,],col="steelblue4")
         abline(v=date[IND],lty=2)
         abline(h=0,lty=2)
         box()
      }
   }
}


Y = matrix(NA,ncol=k, nrow=dim(H)[3])
for (i in 1:k) {
   Y[,i] = log(H[i,i,])
}
colnames(Y) = NAMES


### TVP-VAR Estimation
nlag  = 1
nfore = 5
prior = MinnesotaPrior(0.1, k, nlag)
tvp_var = TVPVAR(Y, l=c(0.99,0.99), nlag, prior)
B_t = tvp_var$beta_t
Q_t = tvp_var$Q_t

### DYNAMIC CONNECTEDNESS APPROACH
t = nrow(Y)
total = matrix(NA, ncol=1, nrow=t)
gfevd = ct = npso = array(NA, c(k, k, t))
net = to = from = matrix(NA, ncol=k, nrow=t)
colnames(gfevd)=rownames(gfevd)=NAMES
for (i in 1:t){
   gfevd[,,i] = GFEVD(B_t[,,i], Q_t[,,i], n.ahead=nfore)$GFEVD
   CV = gfevd[,,i]
   dca = DCA(gfevd[,,i])
   ct[,,i] = dca$CT
   to[i,] = dca$TO
   from[i,] = dca$FROM
   net[i,] = dca$NET
   npso[,,i] = dca$NPSO
   total[i,] = dca$TCI
   if (i%%100==0) print(paste0(round(100*i/t,2),"%"))
}

### DYNAMIC TOTAL CONNECTEDNESS
par(mfcol=c(1,1), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date,total, type="l",xaxs="i",col="steelblue4", las=1, main="",ylab="",ylim=c(0,100),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(0,t)),rev(total)),col="steelblue4", border="steelblue4")
abline(v=date[IND],lty=2)
box()

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
split = 1
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
   plot(date,to[,i], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste(NAMES[i],"TO all others"),ylim=c(0,ceiling(max(to))),tck=-0.02,yaxs="i")
   grid(NA,NULL)
   polygon(c(date,rev(date)),c(c(rep(0,t)),rev(to[,i])),col="steelblue4", border="steelblue4")
   abline(v=date[IND],lty=2)
   box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
   plot(date,from[,i], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste(NAMES[i],"FROM all others"),ylim=c(0,100),tck=-0.02,yaxs="i")
   grid(NA,NULL)
   polygon(c(date,rev(date)),c(c(rep(0,t)),rev(from[,i])),col="steelblue4", border="steelblue4")
   abline(v=date[IND],lty=2)
   box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
split=1
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
   plot(date,net[,i], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste("NET",NAMES[i]),ylim=c(floor(min(net)),ceiling(max(net))),tck=-0.02,yaxs="i")
   grid(NA,NULL)
   polygon(c(date,rev(date)),c(c(rep(0,t)),rev(net[,i])),col="steelblue4", border="steelblue4")
   abline(v=date[IND],lty=2)
   box()
}

### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
i = 7
print(NAMES[i])
par(mfcol=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 1:k) {
   if (i!=j) {
      plot(date,npso[j,i,], xlab="",ylab="",type="l",xaxs="i",col="steelblue4", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso)),ceiling(max(npso))))
      grid(NA,NULL)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(npso[j,i,])),col="steelblue4", border="steelblue4")
      abline(v=date[IND],lty=2)
      box()
   }
}

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
CT = list()
NET = OWN = matrix(NA, nrow=periods, ncol=k)
colnames(NET) = colnames(OWN) = NAMES
rownames(NET) = rownames(OWN) = period_names
for (i in 1:periods) {
  table = DCA(gfevd[,,(IND[i]+1):IND[i+1]])
  CT[[i]] = table$TABLE
  NET[i,] = table$NET
  OWN[i,] = diag(table$CT)
}
print(NET)

# END
