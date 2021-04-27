
### CRUDE OIL FUTURES CONTRACTS AND COMMODITY MARKETS: 
### NEW EVIDENCE FROM A TVP-VAR EXTENDED JOINT CONNECTEDNESS APPROACH
### BALCILAR, M., GABAUER, D., & ZAGHUM, U. (2021)
### RESOURCES POLICY
### replicated by David Gabauer
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
p = 1  # lag length
H = 20 # forecast horizon
prior = MinnesotaPrior(0.1, k, p)
tvpvar = TVPVAR(Y, l=c(0.99, 0.99), p, prior)
B_t = tvpvar$beta_t
Q_t = tvpvar$Q_t

dy12 = DY12(B_t, Q_t, H, NAMES)
ct_dy12 = dy12$CT
total_dy12 = dy12$TOTAL
to_dy12 = dy12$TO
from_dy12 = dy12$FROM
net_dy12 = dy12$NET
npso_dy12 = dy12$NPSO

### LW2020 JOINT CONNECTEDNESS APPROACH ----
lw20 = LW20(B_t, Q_t, H, NAMES)
ct_lw20 = lw20$CT
total_lw20 = lw20$TOTAL
to_lw20 = lw20$TO
from_lw20 = lw20$FROM
net_lw20 = lw20$NET
npso_lw20 = lw20$NPSO

### DYNAMIC TOTAL CONNECTEDNESS
date = DATE
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date, total_lw20, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(floor(min(c(total_dy12,total_lw20))),ceiling(max(c(total_dy12,total_lw20)))),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(0,t)),rev(total_lw20)),col="grey20", border="grey20")
lines(date, total_dy12, col="red")
box()

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
split = 2
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,to_lw20[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"TO all others"),ylim=c(0,ceiling(max(to_dy12))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(to_lw20[,i])),col="grey20", border="grey20")
  lines(date, to_dy12[,i], col="red")
  box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,from_lw20[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"FROM all others"),ylim=c(0,ceiling(max(from_dy12))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(from_lw20[,i])),col="grey20", border="grey20")
  lines(date, from_dy12[,i], col="red")
  box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,net_lw20[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste("NET",colnames(Y)[i]),ylim=c(floor(min(cbind(net_dy12,net_lw20))),ceiling(max(cbind(net_dy12,net_lw20)))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(net_lw20[,i])),col="grey20", border="grey20")
  lines(date, net_dy12[,i], col="red")
  box()
}

### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
kk = k*(k-1)/2
lk = ceiling(sqrt(2))
par(mfcol=c(ceiling(kk/lk),lk), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k) {
  for (j in 1:k) {
    if (i<j) {
      plot(date,npso_lw20[j,i,], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso_lw20)),ceiling(max(npso_lw20))))
      grid(NA,NULL)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(npso_lw20[j,i,])),col="grey20", border="grey20")
      lines(date,npso_dy12[j,i,],col="red")
      box()
    }
  }
}

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
print(dy12$TABLE)
print(lw20$TABLE)

### END
