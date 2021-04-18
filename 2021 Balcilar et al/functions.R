tvp.Phi = function (x, nstep = 10, ...) {
  nstep = abs(as.integer(nstep))
  K=nrow(x)
  p=floor(ncol(x)/K)
  A = array(0, c(K,K,nstep))
  for (i in 1:p){
    A[,,i]=x[,((i-1)*K+1):(i*K)]
  }
  
  Phi = array(0, dim = c(K, K, nstep + 1))
  Phi[, , 1] = diag(K)
  Phi[, , 2] = Phi[, , 1] %*% A[, , 1]
  if (nstep > 1) {
    for (i in 3:(nstep + 1)) {
      tmp1 = Phi[, , 1] %*% A[, , i - 1]
      tmp2 = matrix(0, nrow = K, ncol = K)
      idx = (i - 2):1
      for (j in 1:(i - 2)) {
        tmp2 = tmp2 + Phi[, , j + 1] %*% A[, , idx[j]]
      }
      Phi[, , i] = tmp1 + tmp2
    }
  }
  return(Phi)
}
GFEVD = function(Phi, Sigma, n.ahead=10,normalize=TRUE,standardize=TRUE) {
  A = tvp.Phi(Phi, (n.ahead-1))
  Sigma = Sigma
  gi = array(0, dim(A))
  sigmas = sqrt(diag(Sigma))
  for (j in 1:dim(A)[3]) {
    gi[,,j] = t(A[,,j]%*%Sigma%*%MASS::ginv(diag(sqrt(diag(Sigma)))))
  }
  if (standardize==TRUE){
    girf=array(NA, c(dim(gi)[1],dim(gi)[2], (dim(gi)[3])))
    for (i in 1:dim(gi)[3]){
      girf[,,i]=((gi[,,i])%*%MASS::ginv(diag(diag(gi[,,1]))))
    }
    gi=girf
  }
  
  num = apply(gi^2,1:2,sum)
  den = c(apply(num,1,sum))
  fevd = t(num)/den
  nfevd = fevd
  if (normalize==TRUE) {
    fevd=(fevd/apply(fevd, 1, sum))
  } else {
    fevd=(fevd)
  }
  return = list(GFEVD=fevd, GIRF=gi)
}
DCA = function(CV, digit=2){
  k = dim(CV)[1]
  CT = apply(CV,1:2,mean)*100 # spillover from others to one specific
  OWN = diag(diag(CT))
  TO = colSums(CT-OWN)
  FROM = rowSums(CT-OWN)
  NET = TO-FROM
  TCI = mean(TO)
  NPSO = CT-t(CT)
  NPDC = rowSums(NPSO>0)
  table = format(round(cbind(CT,FROM),digit),nsmall=digit)
  to = c(format(round(c(TO,sum(TO)),digit),nsmall=digit))
  net = c(format(round(c(NET, TCI),digit),nsmall=digit))
  npdc = c(format(round(NPDC,digit),nsmall=digit), "")
  inc = c(format(round(colSums(CT), digit),nsmall=digit), "TCI")
  TABLE = rbind(table,to,inc,net,npdc)
  colnames(TABLE) = c(rownames(CV),"FROM others")
  rownames(TABLE) = c(rownames(CV),"TO others","Inc. own","NET","NPDC")
  return = list(CT=CT,TCI=TCI,TCI_corrected=TCI*k/(k-1),TO=TO,FROM=FROM,NET=NET,NPSO=NPSO,NPDC=NPDC,TABLE=TABLE)
}
DY12 = function(B_t, Q_t, nfore, NAMES) {
  k = dim(Q_t)[1]
  t = dim(Q_t)[3]
  CT_dy12 = array(NA, c(k, k, t))
  NPDC_dy12 = NET_dy12 = FROM_dy12 = TO_dy12 = matrix(NA, ncol=k, nrow=t)
  NPSO_dy12 = array(NA, c(k, k, t))
  TOTAL_dy12 = matrix(NA,ncol=1, nrow=t)
  colnames(CT_dy12) = rownames(CT_dy12) = NAMES
  for (i in 1:t) {
    CT_dy12[,,i] = GFEVD(B_t[,,i], Q_t[,,i], n.ahead=nfore, standardize=TRUE,normalize=TRUE)$GFEVD
    ct = DCA(CT_dy12[,,i])
    NPSO_dy12[,,i] = ct$NPSO
    TOTAL_dy12[i,] = ct$TCI
    TO_dy12[i,] = ct$TO
    FROM_dy12[i,] = ct$FROM
    NET_dy12[i,] = ct$NET
    NPDC_dy12[i,] = ct$NPDC
  }
  return = list(TABLE=DCA(CT_dy12)$TABLE, CT=CT_dy12, TOTAL=TOTAL_dy12, TO=TO_dy12, FROM=FROM_dy12, NET=NET_dy12, NPSO=NPSO_dy12)
}
LW20 = function(B_t, Q_t, nfore, NAMES) {
  k = dim(Q_t)[1]
  t = dim(Q_t)[3]
  TO_lw20 = matrix(NA, ncol=k, nrow=t)
  FROM_lw20 = matrix(NA, ncol=k, nrow=t)
  NET_lw20 = matrix(NA, ncol=k, nrow=t)
  CT_lw20 = NPSO_lw20 = array(NA, c(k, k, t))
  TOTAL_lw20 = matrix(NA, ncol=1, nrow=t)
  colnames(NPSO_lw20) = rownames(NPSO_lw20) = colnames(CT_lw20) = rownames(CT_lw20) = NAMES
  colnames(TO_lw20) = colnames(FROM_lw20) = colnames(NET_lw20) = NAMES
  colnames(TOTAL_lw20) = "TCI"
  for (ij in 1:t) {
    Phi = B_t[,,ij]         # the VAR coefficient matrices 
    A = tvp.Phi(Phi, nfore) # the VMA coefficient matrices
    sigma_eps = Q_t[,,ij]   # the shock covariance matrix
    K = dim(sigma_eps)[1]   # set the number of variables in the VAR
    
    # calculate Xi (the forecast error covariance matrix)
    Xi_h=array(0,dim=c(k,k,nfore))	  
    for (h in 1:nfore){
      Xi_h[,,h]=A[,,h]%*%sigma_eps%*%t(A[,,h]) # calculated Xi at each h
    }
    Xi=rowSums(Xi_h, dims=2) # sum them along THIRD dimension to form Xi  (note: because this is a row sum, dims=2, actually sums along the third dimension)
    
    # identity matrix
    I_K=diag(1,nrow=k,ncol=k)
    
    # calculate the gFEVD
    A_times_sigma_squar_h=array(0,dim=c(k,k,nfore))	
    for (h in 1:nfore){
      A_times_sigma_squar_h[,,h]=(A[,,h]%*%sigma_eps)*(A[,,h]%*%sigma_eps) # calculated A times Sigma squared (element by element squared) at each h
    }
    A_times_sigma_squar_sum=rowSums(A_times_sigma_squar_h, dims=2) # sum them along THIRD dimension (note: because this is a row sum, dims=2, actually sums along the third dimension)
    
    gFEVD=array(0,dim=c(k,k))
    for (i in 1:k){
      for (j in 1:k){
        gFEVD[i,j]=(A_times_sigma_squar_sum[i,j])/(sigma_eps[j,j]*Xi[i,i])
      }
    }
    
    # calculate gSOT
    gFEVD_row_sum = rowSums(gFEVD) # calculate gFEVD row sums
    gSOT = array(0,dim=c(k,k))
    for (i in 1:k){
      for (j in 1:k){
        gSOT[i,j] = (100*gFEVD[i,j])/gFEVD_row_sum[i] # calculate gSOT
      }
    }
    
    # calculate generalized total spillover from all others to variable i (S_gen)
    S_gen = array(0,dim=c(k))
    for (i in 1:k){
      S_gen[i] = sum(gSOT[i,])-gSOT[i,i] # generalized total spillover from others
    }
    
    # calculate the generalized spillover index (gSOI)
    gSOI = mean(S_gen)
    
    # Calculate the elimination matrices.
    # Usually denoted as a KxK-1 matrix M_i, here it is an array where M[,,1]=M_1, and in general M[,,i]=M_i.
    M = array(0,dim=c(k,k-1,k))
    for (i in 1:k){
      M[,,i] = I_K[,-i] # calculate the elimination matrices
    }
    
    # Calculate the joint total spillovers from all others to variable i (S_jnt)
    # calculate the numerator of S_jnt
    S_jnt_numerator_h = array(0,dim=c(k,nfore))
    for (i in 1:k){
      for (h in 1:nfore){
        S_jnt_numerator_h[i,h] = I_K[i,]%*%A[,,h]%*%sigma_eps%*%M[,,i]%*%(solve(t(M[,,i])%*%sigma_eps%*%M[,,i]))%*%t(M[,,i])%*%sigma_eps%*%t(A[,,h])%*%I_K[,i]     #calculate the numerator of S_jnt at each h
      }
    }
    
    S_jnt_numerator = array(0,dim=c(k))  
    for (i in 1:k){
      S_jnt_numerator[i] = sum(S_jnt_numerator_h[i,]) # calculate the numerator of j_jnt  (sum over h)
    }
    
    S_jnt=array(0,dim=c(k))
    for (i in 1:k){
      S_jnt[i]=(100*S_jnt_numerator[i])/Xi[i,i]
    }
    
    # calculate the joint spillover index (jSOI)
    jSOI = mean(S_jnt)
    gSOT_diag = gSOT
    diag(gSOT_diag) = 0
    lambda = S_jnt / apply(gSOT_diag, 1, sum)
    jSOT = gSOT
    for (i in 1:k) {
      jSOT[i,] = gSOT[i,]*lambda[i]
    }
    
    jSOT_diag = jSOT
    diag(jSOT_diag) = 0
    from_jnt = rowSums(jSOT_diag)
    to_jnt = colSums(jSOT_diag)
    jSOI = mean(to_jnt)
    diag(jSOT_diag) = 100-from_jnt
    
    CT_lw20[,,ij] = jSOT_diag
    TO_lw20[ij,] = to_jnt
    FROM_lw20[ij,] = from_jnt
    NET_lw20[ij,] = to_jnt-from_jnt
    NPSO_lw20[,,ij] = jSOT-t(jSOT)
    TOTAL_lw20[ij,] = jSOI
  }
  return = list(TABLE=DCA(CT_lw20/100)$TABLE, CT=CT_lw20, TOTAL=TOTAL_lw20, TO=TO_lw20, FROM=FROM_lw20, NET=NET_lw20, NPSO=NPSO_lw20)
}
UninformativePrior = function(gamma, r, nlag){
  m = nlag*(r^2)
  A_prior = cbind(0*diag(r), matrix(0, ncol=(nlag-1)*r, nrow=r))
  aprior = c(A_prior)
  V_i = matrix(0, nrow=(m/r), ncol=r)
  for (i in 1:r){
    for (j in 1:(m/r)) {
      V_i[j,i] = gamma/(ceiling(j/r)^2)
    }
  }
  # Now V (MINNESOTA VARIANCE) is a diagonal matrix with diagonal elements the V_i'  
  V_i_T = t(V_i)
  Vprior = diag(c(V_i_T))
  diag(Vprior)
  return = list(aprior=aprior, Vprior=Vprior)
}
TVPVAR = function(Y, l, nlag, prior){
  beta_0.mean = prior$aprior
  beta_0.var = prior$Vprior
  Q_0 = prior$Q_0
  if (is.null(Q_0)) {
    Q_0 = cov(Y)
  }
  
  create_RHS_NI = function(templag, r, nlag, t){
    K = nlag*(r^2)
    x_t = matrix(0, (t-nlag)*r, K)
    for (i in 1:(t-nlag)){
      ztemp=NULL
      for (j in 1:nlag){
        xtemp = templag[i,((j-1)*r+1):(j*r)]
        xtemp = t(kronecker(diag(r),xtemp))
        ztemp = cbind(ztemp, xtemp)
      }
      x_t[((i-1)*r+1):(i*r),] = ztemp
    }
    return=list(x_t=x_t, K=K)
  }
  Y = scale(Y,T,F)
  y_true = 0
  FPC = Y
  YX = cbind(Y,Y)
  nfac = 0
  p = n = ncol(Y)
  r = nfac + p
  m = nlag*(r^2)
  k = nlag*r
  t = nrow(FPC)
  q = n + p
  Q_0 = Q_0
  
  # Initialize matrices
  beta_0_prmean = beta_0.mean
  beta_0_prvar = beta_0.var
  
  beta_pred = matrix(0,m,t)
  beta_update = matrix(0,m,t)
  
  Rb_t = array(0,c(m,m,t))
  Sb_t = array(0,c(m,m,t))
  
  beta_t = array(0, c(k,k,t))
  Q_t = array(0, c(r,r,t))
  
  # Decay and forgetting factors
  l_2 = l[1]
  l_4 = l[2]
  
  # Define lags of the factors to be used in the state (VAR) equation         
  yy = FPC[(nlag+1):t,]      
  xx = embed(FPC,nlag+1)[,-c(1:ncol(FPC))]
  templag = embed(FPC,nlag+1)[,-c(1:ncol(FPC))]
  RHS1 = create_RHS_NI(templag,r,nlag,t);  
  Flagtemp = RHS1$x_t
  m = RHS1$K
  Flag = rbind(matrix(0, k,m), Flagtemp)
  
  ###-----| 1. KALMAN FILTER
  for (irep in 1:t){
    #-----| Update the state covariances
    # 1. Get the variance of the factor
    
    # Update Q[t]
    if (irep==1){
      Q_t[,,irep] = Q_0
    } else if (irep > 1) {
      if (irep <= (nlag+1)) { 
        Gf_t = 0.1*(t(matrix(FPC[irep,],nrow=1))%*%(FPC[irep,]))
      } else {
        Gf_t = t(yy[(irep-nlag),]-xx[(irep-nlag),]%*%t(B[1:r,1:k])) %*% (yy[(irep-nlag),]-xx[(irep-nlag),]%*%t(B[1:r,1:k]))
      }
      Q_t[,,irep] = l_2*Q_t[,,(irep-1)] + (1-l_2)*Gf_t[1:r,1:r]
    }
    # -for beta
    if (irep <= (nlag+1)) {
      beta_pred[,irep] = beta_0_prmean
      beta_update[,irep] = beta_pred[,irep]
      Rb_t[,,irep] = beta_0_prvar
    } else if (irep > (nlag+1)) {
      beta_pred[,irep] = beta_update[,(irep-1)]
      Rb_t[,,irep] = (1/l_4)*Sb_t[,,(irep-1)]
    }
    
    # -for beta
    if (irep >= (nlag+1)) {
      # 2/ Update VAR coefficients conditional on Principal Componets estimates
      Rx = Rb_t[,,irep]%*%t(Flag[((irep-1)*r+1):(irep*r),])
      KV_b = Q_t[,,irep] + Flag[((irep-1)*r+1):(irep*r),]%*%Rx
      KG = Rx%*%MASS::ginv(KV_b)
      beta_update[,irep] = matrix(beta_pred[,irep], ncol=1) + (KG%*%(t(matrix(FPC[irep,], nrow=1))-Flag[((irep-1)*r+1):(irep*r),]%*%matrix(beta_pred[,irep], ncol=1)) )
      Sb_t[,,irep] = Rb_t[,,irep] - KG%*%(Flag[((irep-1)*r+1):(irep*r),]%*%Rb_t[,,irep])
    }
    
    # Assign coefficients
    bb = matrix(beta_update[,irep], ncol=1)
    splace = 0
    biga = matrix(0, r,r*nlag)
    for (ii in 1:nlag) {                                          
      for (iii in 1:r) {           
        biga[iii,((ii-1)*r+1):(ii*r)] = t(bb[(splace+1):((splace+r)),1])
        splace = splace + r
      }
    }
    
    B = rbind(biga, cbind(diag(r*(nlag-1)), matrix(0, nrow=r*(nlag-1), ncol=r)))
    
    if ((max(abs(eigen(B)$values))<=1)||(irep==1)){
      beta_t[,,irep] = B
    } else {
      beta_t[,,irep] = beta_t[,,(irep-1)]
      beta_update[,irep] = 0.99*beta_update[,(irep-1)]
    }
  }
  
  return = list(beta_t=beta_t[1:ncol(Y),,], Q_t=Q_t)
}
