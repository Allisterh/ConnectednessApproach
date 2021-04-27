
### CHATZIANTONIOU, I., GABAUER, D., & STENFORS, A. (2021)
### Interest Rate Swaps and the Transmission Mechanism of Monetary Policy:
### A Quantile Connectedness Approach
### ECONOMICS LETTERS
### by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)

library("openxlsx")
library("parallel")
options(warn=-1)
options(mc.cores=detectCores())
source("functions.R")

# LOAD DATA
RAW = read.xlsx("./Chatziantoniou2021.xlsx", detectDates=TRUE)
DATE = RAW[,1]
RAW = RAW[,-1]
NAMES = colnames(RAW) # variable names
k = ncol(RAW)         # number of variables

# FIGURE 1: 1-Year Interest Rate Swaps
par(mfcol=c(1,1), oma=c(1,1,0,0)+0.5, mar=c(1,2,1,1)+1.25, mgp=c(1.5,0.5, 0))
plot(DATE, RAW[,1], type="l", las=1, xaxs="i", tck=-0.02,ylim=c(min(RAW), max(RAW)),
     xlab="", ylab="1-year interest rate swaps (in %)", main="")
abline(h=0, lty=2)
grid(NA,NULL)
for (i in 1:k) {
  lines(DATE, RAW[,i],col=i,lty=1,lwd=1)
}
legend("topright", xpd=TRUE, NAMES, fill=c(1:k), bty='n', cex=1,ncol=1)


# CALCULATE FIRST DIFFERENCED 1Y IRS SERIES
Y = RAW[-1,]
colnames(Y)=NAMES
for (i in 1:k) {
  Y[,i] = diff(RAW[,i])
}

# FIGURE x: First differenced 1-Year Interest Rate Swaps
par(mfcol = c(ceiling(k/2),2), oma = c(.5,1,0,0) + 0.1, mar = c(.5,.5,0,0) + 1, mgp=c(.5, .5, 0))
for (i in 1:k) {
  plot(DATE[-1], Y[,i],type="l",las=1,xlab="",ylab="",main=NAMES[i],col="steelblue4", xaxs="i", ylim=c(-max(abs(Y)), max(abs(Y))))
  grid(NA,NULL)
  lines(DATE[-1], Y[,i],col="steelblue4")
}
# TABLE X: SUMMARY STATISTICS
print(SummaryStatistics(Y))


# QUANTILE VAR CONNECTEDNESS APPROACH
nlag = 1
nfore = 20
space = 200
date = DATE[-c(1:(space+nlag))]
t = length(date)

path_name = "./all_quantiles_IRS_1y.RData"
if (file.exists(path_name)) {
  load(path_name)
} else {
  # in the paper we have used n=100 - 100 quantile VARs are estimated - 
  # however this estimation takes a lot of time which is why we downsampled it to 10
  # if you want a single quantile VAR just exchange the seq() to the tau you want
  n = 10
  quantiles = seq(0.05, 0.95, 1/n) 
  print(quantiles)

  DCA_list = list()
  TCI = array(NA, c(t, length(quantiles)))
  NET = array(NA, c(t, k, length(quantiles)))
  dimnames(TCI) = list(date, quantiles)
  dimnames(NET) = list(as.character(date), NAMES, quantiles)
  for (j in 1:length(quantiles)) {
    gfevd = array(NA, c(k, k, t))
    total = matrix(NA, ncol=1, nrow=t)
    net = matrix(NA, ncol=k, nrow=t)
    colnames(gfevd) = rownames(gfevd) = NAMES
    for (i in 1:t) {
      qvar_est = QVAR(Y[i:(space+i-1),],p=nlag, tau=quantiles[j])
      gfevd[,,i] = GFEVD(qvar_est$B, qvar_est$Q, n.ahead=nfore)$GFEVD
      dca = DCA(gfevd[,,i]) 
      # dca stores many connectedness measures however for the replication only $NET and $TCI are needed
      net[i,] = dca$NET
      total[i,] = dca$TCI
    }
    TCI[,j] = total
    NET[,,j] = net
    DCA_list = c(DCA_list, list(DCA(gfevd)$TABLE))
    print(quantiles[j])
  }
  names(DCA_list) = quantiles
  save.image(path_name)
}

# FIGURE 2: 1-Year Dynamic Total Connectedness
dev.off()
filled.contour(date, quantiles, TCI, xlab="", ylab="", ylim=c(0,1))

# FIGURE 3-6: 1-Year Net Total Directional Connectedness (.)
for (i in 1:k) {
  net.contour(NET=NET, selector=i, nlevels=20)
}

### END
