
### On The Transmission Mechanism Of Country-Specific And International Economic Uncertainty Spillovers:
### Evidence From A TVP-VAR Connectedness Decomposition Approach
### GABAUER, D., and GUPTA, R. (2019)
### ECONOMICS LETTERS
### replicated by David Gabauer

library("openxlsx")
library("parallel")
options(mc.cores=detectCores())
source("functions.R")

DATA = read.xlsx("./dy2012.xlsx", detectDates=TRUE)
DATA = na.omit(DATA)
DATE = as.Date(DATA[,1])
Y = DATA[,-1]
NAMES = colnames(Y)
k = ncol(Y)
t = nrow(Y)

par(mfcol=c(1,1), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(DATE, Y[,1], type="l", ylim=c(min(Y), max(Y)), las=1, xlab="", ylab="", xaxs="i", yaxs="i")
grid(NA, NULL)
abline(h=0, lty=3)
for (i in 1:k) {
  lines(DATE, Y[,i], col=i)
}
box()
legend("topright", legend=NAMES, fill=1:k, ncol=1, cex=.75, bty="n")

#### TVP-VAR Estimation
nlag  = 1
nfore = 10
space = 200
prior = BayesPrior(Y[1:space,], nlag)
tvp_var = TVPVAR(Y, l=c(0.99,0.99), nlag, prior)
B_t = tvp_var$beta_t
Q_t = tvp_var$Q_t

t = nrow(Y)
total = matrix(NA,ncol=1,nrow=t)
net = to = from = matrix(NA, ncol=k, nrow=t)
npso = ct = gfevd = array(NA, c(k, k, t))
colnames(gfevd) = rownames(gfevd) = NAMES
for (i in 1:t) {
  gfevd[,,i] = GFEVD(B_t[,,i], Q_t[,,i], n.ahead=nfore, standardize=TRUE,normalize=TRUE)$GFEVD
  dca = DCA(gfevd[,,i])
  ct[,,i] = dca$CT
  total[i,] = dca$TCI
  to[i,] = dca$TO
  from[i,] = dca$FROM
  net[i,] = dca$NET
  npso[,,i] = dca$NPSO
}

# DECOMPOSITION BY GROUPS
ct_wo = ct
GROUPS = list(c(1,2), c(3,4))
for (i in 1:length(GROUPS)) {
  print(paste("Group",i,NAMES[GROUPS[i][[1]]]))
}
names(GROUPS) = c("SP500&R_10Y", "DJUBSCOM&USDX")
m = length(GROUPS)
npso_group = array(0, c(m, m, t))
dimnames(npso_group)[[1]] = dimnames(npso_group)[[2]] = names(GROUPS)
for (i in 1:m) {
  for (j in 1:m) {
    if (i!=j) {
      group_1 = GROUPS[i][[1]]
      group_2 = GROUPS[j][[1]]
      npso_group[i,j,] = apply(ct[group_1,group_2,],3,mean)
      npso_group[j,i,] = apply(ct[group_2,group_1,],3,mean)
      ct_wo[group_1,group_2,] = 0
      ct_wo[group_2,group_1,] = 0
    }
  }
}

npso_wo = array(NA, c(k, k, t))
to_wo = from_wo = net_wo = matrix(NA, ncol=k, nrow=t)
total_wo = matrix(NA, ncol=1, nrow=t)
for (i in 1:t) {
  to_wo[i,] = colSums(ct_wo[,,i] - diag(diag(ct_wo[,,i])))
  from_wo[i,] = rowSums(ct_wo[,,i] - diag(diag(ct_wo[,,i])))
  net_wo[i,] = to_wo[i,] - from_wo[i,]
  npso_wo[,,i] = ct_wo[,,i] - t(ct_wo[,,i])
  total_wo[i,] = mean(to_wo[i,])
}

total_group = list()
for (i in 1:m) {
  group = GROUPS[i][[1]]
  total_group[[i]] = rowSums(to_wo[,group])/k
}

# Figure 3: Dynamic total connectedness
date = DATE
par(mfcol=c(1,1), oma=c(1,1,0,0) + 0.05, mar=c(1,1,1,1) + .3, mgp=c(0, 0.5, 0))
plot(date, total, type="l",xaxs="i",ylim=c(0,max(total)), las=1, xlab="", ylab="", main="")
grid(NA, NULL)
lines(date, total_wo, col=m+2)
for (i in 1:length(GROUPS)) {
  lines(date, total_group[[i]], col=i+1)
}

# Figure 4: Net total directional connectedness
par(mfcol=c(ceiling(k/2),2), oma=c(1,1,0,0) + 0.05, mar=c(1,1,1,1) + .3, mgp=c(0, 0.5, 0))
for (i in 1:k) {
  plot(date, net[,i], type="l", xlab="", ylab="", las=1, main=paste("NET", NAMES[i]), ylim=c(min(net),max(net)), xaxs="i", yaxs="i")
  grid(NA, NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(net[,i])),col="grey20", border="grey20")
  lines(date, net_wo[,i], col=2)
  abline(h=0, lty=3)
  box()
}

# Figure 5: Total directional connectedness FROM others
par(mfcol=c(ceiling(k/2),2), oma=c(1,1,0,0) + 0.05, mar=c(1,1,1,1) + .3, mgp=c(0, 0.5, 0))
for (i in 1:k) {
  plot(date, from[,i], type="l", xlab="", ylab="", las=1, main=paste("FROM", NAMES[i]), ylim=c(0,max(from)), xaxs="i", yaxs="i")
  grid(NA, NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(from[,i])),col="grey20", border="grey20")
  lines(date, from_wo[,i], col=2)
  box()
}

# Figure 6: Total directional connectedness TO others
par(mfcol=c(ceiling(k/2),2), oma=c(1,1,0,0) + 0.05, mar=c(1,1,1,1) + .3, mgp=c(0, 0.5, 0))
for (i in 1:k) {
  plot(date, to[,i], type="l", xlab="", ylab="", las=1, main=paste("TO", NAMES[i]), ylim=c(0,max(to)), xaxs="i", yaxs="i")
  grid(NA, NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(to[,i])),col="grey20", border="grey20")
  lines(date, to_wo[,i], col=2)
  box()
}

# Figure 7 & 8: Internal & external net pairwise total directional connectedness
kk = k*(k-1)/2
split = 2
par(mfcol=c(ceiling(kk/split),split), oma=c(1,1,0,0) + 0.05, mar=c(1,1,1,1) + .3, mgp=c(0, 0.5, 0))
for (i in 1:k) {
  for (j in 1:k) {
    if (j<i) {
      plot(date, npso[i,j,], type="l", xlab="", ylab="", las=1, main=paste("NET", NAMES[j], '-', NAMES[i]), ylim=c(min(npso),max(npso)), xaxs="i", yaxs="i")
      grid(NA, NULL)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(npso[i,j,])),col="grey20", border="grey20")
      lines(date, npso_wo[i,j,], col=2)
      box()
    }
  }
}

# Figure 9: Market net total directional connectedness
mm = m*(m-1)/2
par(mfcol=c(mm,1), oma=c(1,1,0,0) + 0.05, mar=c(1,1,1,1) + .3, mgp=c(0, 0.5, 0))
for (i in 1:m) {
  for (j in 1:m) {
    if (i<j) {
      plot(date, npso_group[j,i,]-npso_group[i,j,], type="l", xlab="", ylab="", las=1, main=paste(names(GROUPS)[i], '-', names(GROUPS)[j]), ylim=c(-max(npso_group),max(npso_group)), xaxs="i", yaxs="i")
      grid(NA, NULL)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(npso_group[j,i,]-npso_group[i,j,])),col="grey20", border="grey20")
      abline(h=0, lty=3)
      lines(date, npso_group[i,j,], col=i+1)
      lines(date, npso_group[j,i,], col=i+2)
      box()
    }
  }
}

# Table 1 Connectedness table
dec_to = c(round(apply(to_wo, 2, mean), 1),"")
dec_net = c(round(apply(net_wo, 2, mean), 1), "")
dec_from = c(round(apply(from_wo, 2, mean), 1), rep("", 4), "TCI_I", round(mean(total_wo), 1))
CT = cbind(rbind(DCA(gfevd)$TABLE,dec_to,dec_net), dec_from)
print(CT)

### END
