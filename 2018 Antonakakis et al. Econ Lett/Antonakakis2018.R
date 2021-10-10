
library("xlsx")
library("waveslim")
source("functions.R")

RAW = read.xlsx("./data.xlsx", sheetIndex=1)
DATE = as.Date(RAW[,1])
DATA = RAW[,-1]
NAMES = colnames(DATA)
k = ncol(DATA)

Y = NULL
date = DATE[-1]
for (i in 1:k) {
  Y = cbind(Y, diff(DATA[,i]))
}
colnames(Y) = NAMES

waves = 8
Y_list = list()
for (i in 1:k) {
  fit = modwt(Y[,i], wf="la8", n.levels=waves)
  Y_list[[i]] = matrix(unlist(fit), ncol=waves+1)
}

final = list()
for (i in 1:waves) {
  data = NULL
  for (j in 1:k) {
    data = cbind(data, Y_list[[j]][,i])
  }
  colnames(data) = NAMES
  final[[i]] = data
}

split = 2
par(mfcol=c(ceiling(waves/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,2), mgp=c(0.5,0.5,0))
for (i in 1:length(final)) {
  plot(date, final[[i]][,1], type="l", col=j, las=1, xlab="", ylab="", main=paste("Level", i), xaxs="i", ylim=c(-max(abs(final[[i]])), max(abs(final[[i]]))))
  grid(NA, NULL, lty=3)
  for (j in 1:ncol(final[[i]])) {
    lines(date, final[[i]][,j], col=j)
  }
  box()
}

### TVP-VAR based Connectedness Approach
nlag  = 1
nfore = 10
space = 250
CT = TO = FROM = NET = NPSO = TOTAL = list()
for (j in 1:waves) {
  print(j)
  y = final[[j]]
  t = nrow(y)

  vars = MTS::VAR(y[1:space,], p=nlag, output=FALSE)
  prior = list()
  prior$aprior = vars$Phi
  prior$Vprior = 0.01*diag(k^2)
  tvpvar = TVPVAR(y, c(0.99,0.99), nlag, prior)#, Q)
  B_t = tvpvar$beta_t
  Q_t = tvpvar$Q_t
  date = DATE[-c(1:nlag)]

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
  CT[[j]] = DCA(gfevd)$TABLE
  TOTAL[[j]] = total
  TO[[j]] = to
  FROM[[j]] = from
  NET[[j]] = net
  NPSO[[j]] = npso
}

total_summary = NULL
to_summary = from_summary = net_summary = matrix(NA, nrow=waves, ncol=k)
colnames(to_summary) = colnames(from_summary) = colnames(net_summary) = NAMES
rownames(to_summary) = rownames(from_summary) = rownames(net_summary) = c(1:waves)
for (i in 1:waves) {
  total_summary[i] = mean(TOTAL[[i]])
  to_summary[i,] = apply(TO[[i]], 2, mean)
  from_summary[i,] = apply(FROM[[i]], 2, mean)
  net_summary[i,] = apply(NET[[i]], 2, mean)
}
names(total_summary) = rownames(to_summary)
print(to_summary)
print(from_summary)
print(net_summary)
print(total_summary)

### FIGURE 1: DYNAMIC TOTAL CONNECTEDNESS
par(mfrow = c(ceiling(waves/split),split), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1.2,1) + .05, mgp = c(0, 0.2, 0))
for (j in 1:waves) {
  plot(date, TOTAL[[j]], type="l",xaxs="i",col="grey20", las=1, main=paste("Components", j), xlab="", ylab="",ylim=c(0,80),yaxs="i",tck=-0.02)
  grid(NA, NULL, lty=3)
  polygon(c(date, rev(date)), c(c(rep(0, t)),rev(TOTAL[[j]])),col="grey20", border="grey20")
  box()
}

### AVERAGED DYNAMIC CONNECTEDNESS TABLES
wave = 1
print(CT[[wave]])

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
par(mfcol=c(5,3), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k) {
  plot(date,TO[[wave]][,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(NAMES[i],"TO all others"),ylim=c(0,quantile(TO[[j]],0.99)),tck=-0.02,yaxs="i")
  grid(NA,NULL,lty=3)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(TO[[wave]][,i])),col="grey20", border="grey20")
  box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
for (i in 1:k) {
  plot(date,FROM[[wave]][,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(NAMES[i],"FROM all others"),ylim=c(0,quantile(FROM[[j]],0.99)),tck=-0.02,yaxs="i")
  grid(NA,NULL,lty=3)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(FROM[[wave]][,i])),col="grey20", border="grey20")
  box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
for (i in 1:k) {
  plot(date,NET[[wave]][,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste("NET",NAMES[i]),tck=-0.02,yaxs="i", ylim=c(-quantile(abs(NET[[j]]),0.99), quantile(abs(NET[[j]]),0.99)))
  grid(NA,NULL,lty=3)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(NET[[wave]][,i])),col="grey20", border="grey20")
  box()
}

### FIGURE 3:
### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
kk = k*(k-1)/2
par(mfcol = c(ceiling(kk/split),split), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1.2,1) + .05, mgp = c(0, 0.2, 0))
for (i in 1:k) {
  for (l in 1:k) {
    if (i<l) {
      plot(date,NPSO[[wave]][l,i,], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste0(NAMES[i],"-",NAMES[l]),tck=-0.02,yaxs="i",ylim=c(-quantile(abs(NPSO[[j]]), 0.999), quantile(abs(NPSO[[j]]), 0.999)))
      grid(NA,NULL)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(NPSO[[wave]][l,i,])),col="grey20", border="grey20")
      box()
    }
  }
}

### END
