
### RETURN CONNECTEDNESS ACROSS ASSET CLASSES AROUND THE COVID-19 OUTBREAK
### BOURI, E., CEPNI, O., GABAUER, D., & GUPTA, R.
### INTERNATIONAL REVIEW OF FINANCIAL ANALYIS
### replicated by David Gabauer
### Dataset ends with 2020-04-03!

library("plyr")
library("openxlsx")
library("stargazer")
library("RColorBrewer")
colors=c("black","tan4","springgreen4","springgreen2","grey40","steelblue1","maroon4",
         "maroon1","orangered4","orangered","tan1","yellow3","firebrick")
palette(colors)
options(warn=-1)
options("mc.cores"=6)
source("functions.R")

#### DATA PREPARATION ----
RAW = read.xlsx("./data.xlsx", detectDates=TRUE)
RAW = na.omit(RAW)
DATE = as.Date(RAW[,1])

DATA = RAW[,-1]
k = ncol(DATA)
NAMES = colnames(DATA)

# FIGURE 1: Raw Series
par(mfcol=c(k, 1), oma = c(1,1,0,0) + 0.1, mar = c(1,0.5,0.5,0) + 1, mgp=c(0, .65, 0))
for (i in 1:k) {
   plot(DATE, DATA[,i], type="l", xaxs="i", las=1, xlab="", ylab="", main=NAMES[i], col="steelblue4")
   grid(NA,NULL)
   lines(DATE, DATA[,i], col="steelblue4")
   box()
}

#### DATA TRANSFORMATION ----
Y = DATA[-1,]
date = DATE[-1]
for (i in 1:k) {
   Y[,i] = 100*diff(log(DATA[,i]))
}

#FIGURE 2: Daily Percentage Changes
par(mfcol = c(k, 1), oma = c(1,1,0,0) + 0.1, mar = c(1,0.5,0.5,0) + 1, mgp=c(0, .65, 0))
for (i in 1:k) {
   plot(date,Y[,i],type="l",las=1,xlab="",ylab="",main=NAMES[i],ylim=c(-max(abs(Y[,i])),max(abs(Y[,i]))),col="steelblue4")
   grid(NA,NULL)
   lines(date,Y[,i],col="steelblue4")
   box()
}
# TABLE 1: Summary Statistics
print(SummaryStatistics(Y))
print(cor(Y))

### DYNAMIC CONNECTEDNESS APPROACH
nlag = 1
nfore = 20
t = nrow(Y)

# ROLLING WINDOW VAR
space = 200 + nlag
date_rw = date[-c(1:space)]
NET = TO = FROM =matrix(NA, ncol=k, nrow=(t-space))
CV = NPSO = array(NA, c(k, k, (t-space)))
TOTAL = matrix(NA, ncol=1, nrow=(t-space))
colnames(CV) = rownames(CV) = NAMES
for (i in 1:(t-space)) {
   var = VAR(Y[i:(space+i-1),], p=nlag)
   CV[,,i] = GFEVD(var$Phi, var$Sigma, n.ahead=nfore)$GFEVD
   vd = DCA(CV[,,i])
   TO[i,] = vd$TO
   FROM[i,] = vd$FROM
   NET[i,] = vd$NET
   NPSO[,,i] = vd$NPSO
   TOTAL[i,] = vd$TCI
}

# TVP-VAR
prior = BayesPrior(Y, nlag)
tvpvar = TVPVAR(Y, l=c(0.99, 0.99), nlag, prior)
B_t = tvpvar$beta_t
Q_t = tvpvar$Q_t

total = matrix(NA, ncol=1, nrow=t)
gfevd = npso = array(NA, c(k, k, t))
net = to = from = matrix(NA, ncol=k, nrow=t)
colnames(gfevd)=rownames(gfevd)=NAMES
for (i in 1:t) {
   gfevd[,,i] = GFEVD(B_t[,,i], Q_t[,,i], n.ahead=nfore)$GFEVD
   dca = DCA(gfevd[,,i])
   to[i,] = dca$TO
   from[i,] = dca$FROM
   net[i,] = dca$NET
   npso[,,i] = dca$NPSO
   total[i,] = dca$TCI
}


# FIGURE 3: Dynamic Total Connectedness
t = length(date)
par(mfcol = c(1, 1), oma = c(1,1,0,0) + 0.1, mar = c(1,0.5,0.5,0) + 1, mgp=c(0, .65, 0))
plot(date, total, type="l",col="grey20", las=1, main="",ylab="",yaxs="i",ylim=c(0,50),xlab="")
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(-1,t)),rev(total)),col="grey20", border="grey20")
lines(date_rw, TOTAL, col="red")
box()

# FIGURE 4: Net Total Directional Connectedness
par(mfcol = c(k, 1), oma = c(1,1,0,0) + 0.1, mar = c(1,0.5,0.5,0) + 1, mgp=c(0, .65, 0))
for (i in 1:k){
   plot(date,net[,i], xlab="",ylab="",type="l",col="grey20", las=1, main=paste("NET",NAMES[i]),ylim=c(-40,40),yaxs="i")
   grid(NA,NULL)
   abline(h=0, lty=3)
   polygon(c(date,rev(date)),c(rep(0,t),rev(net[,i])),col="grey20", border="grey20")
   lines(date_rw, NET[,i], col="red")
   box()
}

# FIGURE 5: Net Pairwise Directional Connectedness
kk = k*(k-1)/2
par(mfcol=c(ceiling(kk/2),2), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k) {
   for (j in 1:k) {
      if (i<j) {
         plot(date, npso[j,i,], xlab="", ylab="", type="l", col="grey20", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=0.02,yaxs="i",ylim=c(-15,5),cex.axis=0.75)
         grid(NA, NULL, lty=3)
         polygon(c(date, rev(date)), c(rep(0,t), rev(npso[j,i,])), col="grey20", border="grey20")
         lines(date_rw, NPSO[j,i,], col="red")
         box()
      }
   }
}

#### AVERAGE CONNECTEDNESS TABLE ----
ind = which(date=="2020-01-13")
print(DCA(gfevd[,,-c(ind:dim(gfevd)[3])])$TABLE)
print(DCA(gfevd[,,c(ind:dim(gfevd)[3])])$TABLE)

# END
