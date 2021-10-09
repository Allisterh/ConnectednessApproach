
### ANTONAKAKIS, N., CHATZIANTONIOU, I., AND GABAUER, D. (2019)
### CRYPTOCURRENCY MARKET CONTAGION: MARKET UNCERTAINTY, MARKET COMPLEXITY, AND DYNAMIC PORTFOLIOS (I)
### Journal of International Financial Markets, Institutions and Money
### by David Gabauer (https://sites.google.com/view/davidgabauer/contact-details)

# Information: As my laptop had a defect in 2019 I have lost the original code of this paper.
# Hence, this is an approximation of the original paper. Instead of the TVP-FAVAR method of Koop and Korobilis (2014)
# I am using the standard Principal Component Analysis and plug the factor scores into the TVP-VAR model of Koop and Korobilis (2014).
# Furthermore, I tried to downloaded all cryptocurrencies from CoinMarketCap, however, I could only retrieve a subsample of the original dataset.
# Thus, instead of 45 series that have been used for the PCA I am currently using only 9. Results, however, are qualitatively similar.
# Nevertheless, I wanted to make this code public as I received a lot of requests concerning the TVP-FAVAR model.
# Even though this is not the same it leads to similar results and I sincerely hope that it will be useful for future research.
# Kind regards, David

library("psych")
library("openxlsx")
library("parallel")
options(mc.cores=detectCores())
source("functions.R")

RAW = read.xlsx('./data.xlsx', detectDates=TRUE)
DATE = as.Date(RAW[,1])
DATA = RAW[,-1]
NAMES = colnames(DATA)
k = ncol(DATA)

date = DATE[-1]
Z = DATA[-1,]
for (i in 1:k) {
  x = embed(DATA[,i],2)
  Z[,i] = 100*(x[,1]-x[,2])/x[,1]
}
colnames(Z)=NAMES
y = Z[,1:9]
X = Z[,-c(1:9)]
PC = principal(X, nfactors=1)$scores

Y = data.frame(y, PC)
k = ncol(Y)
t = nrow(Y)
NAMES = colnames(Y)
NAMES[length(NAMES)] = "PC"

### TVP-VAR
nlag = 1
nfore = 20
prior = BayesPrior(Y[1:200,], nlag)
tvpvar = TVPVAR(Y, l=c(0.99, 0.99), nlag, prior)
B_t = tvpvar$beta_t
Q_t = tvpvar$Q_t

### DYNAMIC CONNECTEDNESS APPROACH
gfevd = npso = array(NA, c(k, k, t))
net = from = to = matrix(NA, ncol=k, nrow=t)
total = matrix(NA, ncol=1, nrow=t)
colnames(gfevd)=rownames(gfevd)=NAMES
for (i in 1:t) {
  gfevd[,,i] = GFEVD(B_t[,,i], Q_t[,,i], n.ahead=nfore)$GFEVD
  dca = DCA(gfevd[,,i])
  to[i,] = dca$TO
  from[i,] = dca$FROM
  net[i,] = dca$NET
  npso[,,i] = dca$NPSO
  total[i,] = dca$TCI
  if (i%%100==0) print(paste0(round(100*i/t,2),"%"))
}

### DYNAMIC TOTAL CONNECTEDNESS
par(mfcol = c(1,1), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1.2,1) + .05, mgp = c(0, 0.2, 0))
plot(date,total, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(0,100),yaxs="i",xlab="",tck=-0.01)
grid(NA,NULL,lty=3)
polygon(c(date,rev(date)),c(c(rep(0,nrow(total))),rev(total)),col="grey20", border="grey20")
box()

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
par(mfcol = c(ceiling(k/2),2), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1.2,1) + .05, mgp = c(0, 0.2, 0))
for (i in 1:k) {
  plot(date,to[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"TO all others"),ylim=c(floor(min(to)),ceiling(max(to))),tck=-0.02,yaxs="i")
  grid(NA,NULL,lty=3)
  polygon(c(date,rev(date)),c(c(rep(0,nrow(to))),rev(to[,i])),col="grey20", border="grey20")
  box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfcol = c(ceiling(k/2),2), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1.2,1) + .05, mgp = c(0, 0.2, 0))
for (i in 1:k) {
  plot(date,from[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"FROM all others"),ylim=c(floor(min(from)),ceiling(max(from))),tck=-0.02,yaxs="i")
  grid(NA,NULL,lty=3)
  polygon(c(date,rev(date)),c(c(rep(0,nrow(from))),rev(from[,i])),col="grey20", border="grey20")
  box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
par(mfcol = c(ceiling(k/2),2), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1.2,1) + .05, mgp = c(0, 0.2, 0))
for (i in 1:k) {
  plot(date,net[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste("NET",colnames(Y)[i]),ylim=c(floor(min(net)),ceiling(max(net))),tck=-0.02,yaxs="i")
  grid(NA,NULL,lty=3)
  polygon(c(date,rev(date)),c(c(rep(0,nrow(net))),rev(net[,i])),col="grey20", border="grey20")
  box()
}

### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
kk = k*(k-1)/2
par(mfcol=c(5,3), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:2) {
  for (j in 1:(k-1)) {
    if (i<j) {
      plot(date,npso[j,i,], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso)),ceiling(max(npso))))
      grid(NA,NULL)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(npso[j,i,])),col="grey20", border="grey20")
      box()
    }
  }
}

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
ind = which(date=="2017-07-20")
print(DCA(gfevd[,,c(1:ind)])$TABLE)
print(DCA(gfevd[,,-c(1:ind)])$TABLE)

### END OF CONNECTEDNESS ----

##### GARCH MODELLING
### DCC-GARCH t-COPULA ANALYSIS
k = ncol(y)
lag = 20
prob = 0.10
conf.level = 0.90

spec = c()
evaluation_matrix = matrix(NA, ncol=k, nrow=10)
colnames(evaluation_matrix) = NAMES
rownames(evaluation_matrix) = c("VaR","","CVaR","","VaR Dur.","","Sign Bias","",paste0("WARCH(",lag,")"),"")
for (j in 1:k) {
  print(NAMES[j])
  bestgarch = BestGARCH(data=y[,j],prob=prob,conf.level=conf.level,lag=lag,
                        distr=c("norm","snorm","std","sstd","ged","sged"), 
                        models=c("sGARCH","iGARCH","eGARCH","gjrGARCH", "AVGARCH", "TGARCH"))
  GARCH_selection = which(bestgarch$GARCH_IC==min(bestgarch$GARCH_IC),arr.ind=TRUE)
  print(bestgarch$GARCH_IC)
  ugarch.spec = bestgarch[[2]][[GARCH_selection[2]]][[GARCH_selection[1]]]
  ugarch.fit = ugarchfit(ugarch.spec,data=y[,j], solver="hybrid")
  evaluation_matrix[,j] = rbind(ValueAtRisk(ugarch.fit,ugarch.spec,prob=prob, conf.level=conf.level)$statistic,
                                SignBias_WARCH(ugarch.fit,lag=lag)$statistic)
  spec = c(spec, ugarch.spec)
}
print(evaluation_matrix)
mgarch.spec = cgarchspec(uspec=multispec(spec), dccOrder=c(1,1), asymmetric=FALSE,
                         distribution.model=list(copula="mvt", method="Kendall", time.varying=TRUE, transformation="parametric"))
copula_fit = cgarchfit(mgarch.spec, data=y, solver=c("hybrid","solnp"), fit.control=list(eval.se=TRUE))
H = rcov(copula_fit)
R = rcor(copula_fit)

### FIGURE4: DYNAMIC CONDITIONAL CORRELATIONS
par(mfcol = c(1,1), oma = c(.5,.5,0,0) + 0.5, mar = c(.5,.5,0.5,0) + 1, mgp=c(.5, .5, 0))
plot(date,R[1,2,], xlab="",ylab="",type="l",xaxs="i", las=1, main="",tck=-0.02,yaxs="i",ylim=c(min(R),1))
grid(NA,NULL)
for (i in 1:2) {
  for (j in 1:k) {
    if (i<j) {
      lines(date,R[j,i,], col="grey20")
    }
  }
}
lines(date, R[1,4,], col="steelblue4", lwd=2)
abline(h=0, lty=3)
abline(v=date[ind], col="steelblue4")
box()

#### HEDGING RATIOS AND PORTFOLIO WEIGHTS ----
PW = PortfolioWeights(y, H)$portfolio_weights
HR = HedgeRatio(y, H)$hedge_ratio

### FIGURE 5: DYNAMIC PORTFOLIO WEIGHTS - corrected
par(mfcol = c(5,3), oma = c(.5,.5,0,0) + 0.5, mar = c(.5,.5,0.5,0) + 1, mgp=c(.5, .5, 0))
for (i in 1:2) {
  for (j in 1:k) {
    if (i>=j) {
      next
    } else {
      plot(date, PW[i,j,], type="l", col="grey20", las=1, ylim=c(0,1), yaxs="i",
           main=paste0(NAMES[i],"/",NAMES[j]), xlab="", ylab="", xaxs="i")
      grid(NA, NULL)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(PW[i,j,])),col="grey20", border="grey20")
      abline(h=0)
      abline(v=date[ind], col="steelblue4")
      box()
    }
  }
}   

### FIGURE 6: DYNAMIC HEDGE RATIOS
par(mfcol = c(5,3), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1.2,1) + .05, mgp = c(0, 0.2, 0))
for (i in 1:2) {
  for (j in 1:k) {
    if (i>=j) {
      next
    } else {
      plot(date, HR[i,j,],type="l",xaxs="i",col="grey20",las=1,ylim=c(-1,2),
           main=paste0(NAMES[i],"/",NAMES[j]), xlab="",ylab="", tck=-0.02)
      grid(NA,NULL)
      polygon(c(date,rev(date)),c(c(rep(0,t)),rev(HR[i,j,])),col="grey20", border="grey20")
      abline(h=0)
      abline(v=date[ind], col="steelblue4")
      box()
    }
  }
}

# TABLE 3: DYNAMIC PORTFOLIO WEIGHTS AND HEDGING EFFECTIVENESS
print(PortfolioWeights(y[c(1:ind),], H[,,c(1:ind)])$summary)
print(PortfolioWeights(y[-c(1:ind),], H[,,-c(1:ind)])$summary)

# TABLE 4: DYNAMIC HEDGE RATIOS AND HEDGING EFFECTIVENESS - i and j needs to be changed
print(HedgeRatio(y[c(1:ind),], H[,,c(1:ind)])$summary)
print(HedgeRatio(y[-c(1:ind),], H[,,-c(1:ind)])$summary)

### END OF GARCH MODELLING ----
