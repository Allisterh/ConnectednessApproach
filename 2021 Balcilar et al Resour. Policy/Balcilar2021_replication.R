
### CRUDE OIL FUTURES CONTRACTS AND COMMODITY MARKETS: 
### NEW EVIDENCE FROM A TVP-VAR EXTENDED JOINT CONNECTEDNESS APPROACH
### BALCILAR, M., GABAUER, D., & ZAGHUM, U. (2021)
### RESOURCES POLICY
### replicated by David Gabauer

### LOAD LIBRARIES ----
library("openxlsx")
library("stargazer")

### LOAD CODE ----
source("functions.R")

### DATA DESCRIPTION ----
RAW = read.xlsx("./data.xlsx", detectDates=TRUE, 1)
RAW = na.omit(RAW)
print(str(RAW))

DATA = RAW[,-1]
k = ncol(DATA)
print(paste("Using", k, "series, namely:"))
NAMES = colnames(DATA)
print(NAMES)
DATE = as.Date(RAW[,1], "%Y-%m-%d")
print(paste("From", DATE[1], "to", DATE[length(DATE)]))

### Calculate first log differences ----
date = DATE[-1]
Y = matrix(NA, ncol=k, nrow=nrow(DATA)-1)
for (i in 1:k) {
  Y[,i] = 100*(na.omit(diff(log(DATA[,i]))))
}

split = 2
par(mfrow=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k) {
  plot(date,Y[,i], type="l", las=1, xlab="", ylab="", main=NAMES[i], col="steelblue4", 
       xaxs="i", ylim=c(-12, 12),xaxt="n",cex.axis=0.75, tck=-0.025)
  grid(NA,NULL)
  lines(date,Y[,i],col="steelblue4")
  axis.Date(side=1, date, at=seq(date[1], tail(date, 1), by="years"), format="%Y", las=2, tck=-0.01, cex.axis=0.75)
  box()
}

print("Summary Statistics")
colnames(Y) = NAMES
summary_statistics = Moments(Y)
print(summary_statistics)

### TVP-VAR ----
p = 1  # lag length
H = 20 # forecast horizon
prior = UninformativePrior(0.1, k, p)
tvpvar = TVPVAR(Y, l=c(0.99, 0.99), p, prior)
B_t = tvpvar$beta_t
Q_t = tvpvar$Q_t


### DY2012 CONNECTEDNESS APPROACH ----
dy12 = DY12(B_t, Q_t, H, NAMES)
CT_dy12 = dy12$CT
TOTAL_dy12 = dy12$TOTAL
NET_dy12 = dy12$NET
NPSO_dy12 = dy12$NPSO
print("DY12; Averaged connectedness table")
print(dy12$TABLE)

### LW2020 JOINT CONNECTEDNESS APPROACH ----
lw20 = LW20(B_t, Q_t, H, NAMES)
CT_lw20 = lw20$CT
TOTAL_lw20 = lw20$TOTAL
NET_lw20 = lw20$NET
NPSO_lw20 = lw20$NPSO
print("DY12; Averaged connectedness table")
print(lw20$TABLE)

### CONNECTEDNESS PLOTS ----
# DYNAMIC TOTAL CONNECTEDNESS
t = length(TOTAL_lw20)
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date, TOTAL_lw20, type="l",xaxs="i",col="grey20", las=1, main="Dynamic Total Connectedness",ylab="",ylim=c(20,90),yaxs="i",xlab="",tck=-0.01,xaxt="n", cex.axis=0.75)
grid(NA,NULL,lty=3)
polygon(c(date, rev(date)), c(rep(0, t), rev(TOTAL_lw20)), col="grey20", border="grey20")
lines(date, TOTAL_dy12, col="red")
axis.Date(side=1, date, at=seq(date[1], tail(date, 1), by="years"), format="%Y", las=2, tck=-0.01, cex.axis=0.75)
box()

# NET TOTAL DIRECTIONAL CONNECTEDNESS
par(mfrow=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date, NET_lw20[,i], xlab="", ylab="", type="l", xaxs="i", col="grey20", las=1, main=paste("NET",NAMES[i]),ylim=c(-100,100),tck=0.01,yaxs="i",xaxt="n",cex.axis=0.75)
  grid(NA, NULL, lty=3)
  polygon(c(date, rev(date)), c(rep(0, t), rev(NET_lw20[,i])), col="grey20", border="grey20")
  lines(date, NET_dy12[,i], col="red")
  abline(h=0, lty=3)
  axis.Date(side=1, date, at=seq(date[1], tail(date, 1), by="years"), format="%Y", las=2, tck=-0.01, cex.axis=0.75)
  box()
}

# NET PAIRWISE DIRECTIONAL CONNECTEDNESS
j = 1
par(mfrow=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 2:k) {
  plot(date, NPSO_lw20[j,i,], xlab="", ylab="", type="l", xaxs="i", col="grey20", las=1, main=paste0(NAMES[j],"-",NAMES[i]),tck=0.02,yaxs="i",ylim=c(-10,10),xaxt="n",cex.axis=0.75)
  grid(NA, NULL, lty=3)
  polygon(c(date, rev(date)), c(rep(0,t), rev(NPSO_lw20[j,i,])), col="grey20", border="grey20")
  lines(date, NPSO_dy12[j,i,], col="red")
  axis.Date(side=1, date, at=seq(date[1], tail(date, 1), by="years"), format="%Y", las=2, tck=-0.01, cex.axis=0.75)
  box()
}

### END
