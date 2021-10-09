
library("HH")
library("abind")
library("psych")
library("rmgarch")
library("openxlsx")
library("parallel")
library("stargazer")
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
  Z[,i] = (x[,1]-x[,2])/x[,1]
}
colnames(Z)=NAMES
Y = Z[,1:9]
k = ncol(Y)
NAMES = colnames(Y)

##### GARCH MODELLING
### DCC-GARCH t-COPULA ANALYSIS
lag = 20
prob = 0.10
conf.level = 0.95

spec = c()
copula_fit = GARCH_spec = evaluation_list = list()
evaluation_matrix = matrix(NA, ncol=k, nrow=10)
colnames(evaluation_matrix) = NAMES
rownames(evaluation_matrix) = c("VaR","","CVaR","","VaR Dur.","","Sign Bias","",paste0("WARCH(",lag,")"),"")
for (j in 1:k) {
  print(NAMES[j])
  #bestgarch = BestGARCH(data=Y[,j],prob=prob,conf.level=conf.level,lag=lag,
  #                      distr=c("norm","snorm","std","sstd","ged","sged"), 
  #                      models=c("sGARCH","iGARCH","eGARCH","gjrGARCH","AVGARCH","TGARCH"))
  bestgarch = BestGARCH(data=Y[,j],prob=prob,conf.level=conf.level,lag=lag,
                        distr=c("std"), 
                        models=c("gjrGARCH"))
  GARCH_selection = which(bestgarch$GARCH_IC==min(bestgarch$GARCH_IC),arr.ind=TRUE)
  print(paste(colnames(bestgarch$GARCH_IC)[GARCH_selection[2]], rownames(bestgarch$GARCH_IC)[GARCH_selection[1]]))
  print(bestgarch$GARCH_IC)
  ugarch.spec = bestgarch[[2]][[GARCH_selection[1]]][[GARCH_selection[2]]]
  ugarch.fit = ugarchfit(ugarch.spec,data=Y[,j])
  evaluation_matrix[,j] = rbind(ValueAtRisk(ugarch.fit,ugarch.spec,prob=prob, conf.level=conf.level)$statistic,
                                SignBias_WARCH(ugarch.fit,lag=lag)$statistic)
  spec = c(spec, ugarch.spec)
}
print(evaluation_matrix)
mgarch.spec = cgarchspec(uspec=multispec(spec), dccOrder=c(1,1), asymmetric=FALSE,  
                         distribution.model=list(copula="mvt", method="Kendall", time.varying=TRUE, transformation="parametric"))
copula_fit = cgarchfit(mgarch.spec, data=Y, solver=c("hybrid","solnp"), fit.control=list(eval.se=TRUE))
H = rcov(copula_fit)
R = rcor(copula_fit)

### FIGURE4: DYNAMIC CONDITIONAL CORRELATIONS
ind = which(date=="2017-07-20")
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
PW = PortfolioWeights(Y, H)$portfolio_weights
HR = HedgeRatio(Y, H)$hedge_ratio

### FIGURE 5: DYNAMIC PORTFOLIO WEIGHTS - corrected
par(mfcol = c(5,3), oma = c(.5,.5,0,0) + 0.5, mar = c(.5,.5,0.5,0) + 1, mgp=c(.5, .5, 0))
for (i in 1:k) {
  for (j in 1:k) {
    if (i>=j) {
      next
    } else {
      plot(date, PW[i,j,], type="h", col="grey20", las=1, ylim=c(0,1), yaxs="i",
           main=paste0(NAMES[i],"/",NAMES[j]), xlab="", ylab="", xaxs="i")
      grid(NA, NULL)
      #polygon(c(date,rev(date)),c(c(rep(0,t)),rev(PW[i,j,])),col="grey20", border="grey20")
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
print(PortfolioWeights(Y[c(1:ind),], H[,,c(1:ind)])$summary)
print(PortfolioWeights(Y[-c(1:ind),], H[,,-c(1:ind)])$summary)

# TABLE 4: DYNAMIC HEDGE RATIOS AND HEDGING EFFECTIVENESS
print(HedgeRatio(Y[c(1:ind),], H[,,c(1:ind)])$summary)
print(HedgeRatio(Y[-c(1:ind),], H[,,-c(1:ind)])$summary)

### END ----
