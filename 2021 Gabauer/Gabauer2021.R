
library("rmgarch")
library("openxlsx")
library("parallel")
library("RColorBrewer")
options(mc.cores=detectCores())
source("functions.R")

DATA = read.xlsx("./Gabauer2021.xlsx", detectDate=TRUE)
DATA = na.omit(DATA)
DATE = as.Date(DATA[,1],"%Y-%m-%d")
RAW = DATA[,-1]
NAMES = colnames(RAW)
k = ncol(RAW)

Y = na.omit(100*diff(log(as.matrix(RAW))))
date = DATE[-1]

start_date = which(date=="1979-03-13")
end_date = which(date=="1999-01-04")
ERM1 = rep(0,length(date))
ERM1[1:start_date] = 1
ERM1[end_date:length(ERM1)] = 1

par(mfrow = c(ceiling(k/2),2), oma = c(0.5,0.5,0,0) + 0.1, mar = c(0.5,1,1,0) + .5, mgp = c(3, 0.5, 0))
for (i in 1:k){
  plot(date,Y[,i], col="grey20",lwd=0.01, type="l", main=NAMES[i], ylim=c(-4.2,4.2),xaxs="i", las=1,yaxs="i")
  polygon(c(date,rev(date)),c(ERM1*c(rep(-100,length(date))),rev(ERM1*c(rep(100,length(date))))),col="grey85", border="grey85")
  grid(NA,NULL,col="grey20")
  lines(date,Y[,i],col="steelblue4",lwd=0.5)
  box()
}

MODEL = c("RETURN","VOLATILITY")
for (model in MODEL) {
  print(model)
  if (model=="RETURN" || model=="VOLATILITY") {
    if(model=="VOLATILITY") {
      Y = abs(Y)
    }
    
    print(paste0("-- TVP-VAR"))
    nlag = 1
    nfore = 10
    prior = BayesPrior(Y[1:200,], nlag)
    tvpvar = TVPVAR(Y, l=c(0.99, 0.99), nlag, prior)
    B_t = tvpvar$beta_t
    Q_t = tvpvar$Q_t

    print(paste0("-- GFEVD"))
    t = nrow(Y)
    k = ncol(Y)
    gfevd = array(NA, c(k, k, t))
    colnames(gfevd)=rownames(gfevd)=NAMES
    for (ik in 1:t){
      gfevd[,,ik] = GFEVD(B_t[,,ik], Q_t[,,ik], n.ahead=nfore, standardize=TRUE,normalize=TRUE)$GFEVD
    }
    gfevd = gfevd[,,start_date:end_date]

    print(paste0("-- CONNECTEDNESS APPROACH"))
    t = dim(gfevd)[3]
    total_corr = total = matrix(NA,ncol=1,nrow=t)
    influence = pci = npso = array(NA,c(k,k,t))
    net = to = from = matrix(NA, ncol=k, nrow=t)
    for (ik in 1:t){
      dca = DCA(gfevd[,,ik])
      total[ik,] = dca$TCI
      total_corr[ik,] = dca$TCI_corrected
      to[ik,] = dca$TO
      from[ik,] = dca$FROM
      net[ik,] = dca$NET
      npso[,,ik] = dca$NPSO
      influence[,,ik] = dca$INFLUENCE
      pci[,,ik] = dca$PCI
    }

    jk = 1
    kk = k/2*(k-1)
    influence_upper = pci_upper = array(NA,c(t,kk))
    colnames(pci_upper) = colnames(influence_upper) = 1:kk
    for (i in 1:k) {
      for (j in 1:k) {
        if (i>j) {
          pci_upper[,jk] = pci[j,i,]
          influence_upper[,jk] = influence[j,i,]
          colnames(pci_upper)[jk] = paste0(NAMES[j],"-",NAMES[i])
          colnames(influence_upper)[jk] = paste0(NAMES[j],"-",NAMES[i])
          jk = jk + 1
        }
      }
    }
  }
  rep = 1000
  pci_mean = influence_mean = matrix(NA,ncol=kk,nrow=rep)
  colnames(pci_mean) = colnames(influence_mean) = names(pci_upper)
  for (i in 1:rep) {
    pci_mean[i,] = apply(pci_upper[sample(1:t,replace=TRUE),],2,mean)
    influence_mean[i,] = apply(influence_upper[sample(1:t,replace=TRUE),],2,mean)
  }
  
  if (model=="RETURN") {
    results_return = list(net=net,to=to,from=from,total=total,total_corr=total_corr,npso=npso,
                          gfevd=gfevd,pci=pci,pci_mean=pci_mean,pci_upper=pci_upper,
                          influence_upper=influence_upper,influence_mean=influence_mean)
  } else if (model=="VOLATILITY") {
    results_volatility = list(net=net,to=to,from=from,total=total,total_corr=total_corr,npso=npso,
                              gfevd=gfevd,pci=pci,pci_mean=pci_mean,pci_upper=pci_upper,
                              influence_upper=influence_upper,influence_mean=influence_mean)
  }
}


### DYNAMIC TOTAL CONNECTEDNESS
date = DATE[start_date:end_date]
par(mfcol=c(1,1), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
plot(date,results_return$total_corr, type="l",xaxs="i",col="grey20", las=1, main="",ylab="",ylim=c(60,100),yaxs="i",xlab="",tck=-0.02)
grid(NA,NULL)
polygon(c(date,rev(date)),c(c(rep(0,t)),rev(results_return$total_corr)),col="grey20", border="grey20")
lines(date,results_volatility$total_corr,col="red")
box()
legend("bottomleft", c("TCI corrected","TCI"), fill=c("grey20", "red"))

### TOTAL DIRECTIONAL CONNECTEDNESS TO OTHERS
split = 2
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,results_return$to[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"TO all others"),ylim=c(floor(min(to)),ceiling(max(to))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(results_return$to[,i])),col="grey20", border="grey20")
  lines(date,results_volatility$to[,i],col="red")
  box()
}

### TOTAL DIRECTIONAL CONNECTEDNESS FROM OTHERS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,1.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,results_return$from[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste(colnames(Y)[i],"FROM all others"),ylim=c(floor(min(from)),ceiling(max(from))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(results_return$from[,i])),col="grey20", border="grey20")
  lines(date,results_volatility$from[,i],col="red")
  box()
}

### NET TOTAL DIRECTIONAL CONNECTEDNESS
par(mfcol=c(ceiling(k/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (i in 1:k){
  plot(date,results_return$net[,i], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste("NET",colnames(Y)[i]),ylim=c(floor(min(net)),ceiling(max(net))),tck=-0.02,yaxs="i")
  grid(NA,NULL)
  polygon(c(date,rev(date)),c(c(rep(0,t)),rev(results_return$net[,i])),col="grey20", border="grey20")
  lines(date,results_volatility$net[,i],col="red")
  box()
}

### NET PAIRWISE DIRECTIONAL CONNECTEDNESS
i = 4
split = 3
print(NAMES[i])
par(mfcol=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 1:k) {
  if (i!=j) {
    plot(date,results_return$npso[j,i,], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(floor(min(npso)),ceiling(max(npso))))
    grid(NA,NULL)
    polygon(c(date,rev(date)),c(c(rep(0,t)),rev(results_return$npso[j,i,])),col="grey20", border="grey20")
    lines(date,results_volatility$npso[j,i,],col="red")
    box()
  }
}

### PAIRWISE CONNECTEDNESS INDEX
par(mfcol=c(ceiling((k-1)/split),split), oma=c(0.5,0.5,0,0), mar=c(1.5,1,1.5,1), mgp=c(0.5,0.5,0))
for (j in 1:k) {
  if (i!=j) {
    plot(date,results_return$pci[j,i,], xlab="",ylab="",type="l",xaxs="i",col="grey20", las=1, main=paste0(NAMES[i],"-",NAMES[j]),tck=-0.02,yaxs="i",ylim=c(0,1))
    grid(NA,NULL)
    polygon(c(date,rev(date)),c(c(rep(0,t)),rev(results_return$pci[j,i,])),col="grey20", border="grey20")
    lines(date,results_volatility$pci[j,i,],col="gold4")
    box()
  }
}

### AVERAGE DYNAMIC CONNECTEDNESS TABLE
print(DCA(results_return$gfevd)$TABLE)
print(DCA(results_volatility$gfevd)$TABLE)

### RANKING OF PCI AND ABSOLUTE PAIRWISE INFLUENCE INDEX
par(mfrow = c(1,2), oma = c(0,1,0,0) + 0.05, mar = c(1,1,1,1) + .05, mgp = c(0, 0.1, 0))
plot(sort(apply(results_return$pci_mean,2,mean), TRUE),xaxs="i",las=1,col="steelblue4",yaxs="i",xlab="",ylab="",tck=0.01,type="s",ylim=c(0,1))
grid(NA,NULL)
lines(sort(apply(results_return$pci_mean,2,quantile,0.05), TRUE),type="s",pch="-",col="steelblue4",lty=3)
lines(sort(apply(results_return$pci_mean,2,quantile,0.95), TRUE),type="s",pch="-",col="steelblue4",lty=3)
lines(sort(apply(results_volatility$pci_mean,2,mean), TRUE),type="s",pch="-",col="gold3",lty=3)
lines(sort(apply(results_volatility$pci_mean,2,quantile,0.05), TRUE),type="s",pch="-",col="gold3",lty=3)
lines(sort(apply(results_volatility$pci_mean,2,quantile,0.95), TRUE),type="s",pch="-",col="gold3",lty=3)
abline(v=16,lty=3)
box()
legend("topright", c("RETURN","VOLATILITY"), fill=c("steelblue4","gold3","red"))

plot(sort(apply(results_return$influence_mean,2,mean)),xaxs="i",las=1,col="steelblue4",yaxs="i",xlab="",ylab="",tck=0.01,type="s",ylim=c(0,0.5))
grid(NA,NULL)
lines(sort(apply(results_return$influence_mean,2,quantile,0.05)),type="s",pch="-",col="steelblue4",lty=3)
lines(sort(apply(results_return$influence_mean,2,quantile,0.95)),type="s",pch="-",col="steelblue4",lty=3)
lines(sort(apply(results_volatility$influence_mean,2,mean)),type="s",pch="-",col="gold3",lty=3)
lines(sort(apply(results_volatility$influence_mean,2,quantile,0.05)),type="s",pch="-",col="gold3",lty=3)
lines(sort(apply(results_volatility$influence_mean,2,quantile,0.95)),type="s",pch="-",col="gold3",lty=3)
abline(v=16,lty=3)
box()
legend("topright", c("RETURN","VOLATILITY"), fill=c("steelblue4","gold3"))

# RANKING
kk = k*(k-1)/2
return_pci_ranking = names(sort(apply(results_return$pci_upper,2,mean),decreasing=TRUE))
volatility_pci_ranking = names(sort(apply(results_volatility$pci_upper,2,mean),decreasing=TRUE))
print(return_pci_ranking[1:20])
print(volatility_pci_ranking[1:20])

# END
