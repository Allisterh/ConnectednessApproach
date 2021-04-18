MTSVAR = function (x, p = 1, include.mean = T, fixed = NULL) {
  if (!is.matrix(x)) 
    x = as.matrix(x)
  Tn = dim(x)[1]
  k = dim(x)[2]
  if (p < 1) 
    p = 1
  idm = k * p
  ne = Tn - p
  ist = p + 1
  y = x[ist:Tn, ]
  if (include.mean) {
    idm = idm + 1
    xmtx = cbind(rep(1, ne), x[p:(Tn - 1), ])
  }
  else {
    xmtx = x[p:(Tn - 1), ]
  }
  if (p > 1) {
    for (i in 2:p) {
      xmtx = cbind(xmtx, x[(ist - i):(Tn - i), ])
    }
  }
  ndim = ncol(xmtx)
  if (length(fixed) == 0) {
    paridx = matrix(1, ndim, k)
  }
  else {
    paridx = fixed
  }
  res = NULL
  beta = matrix(0, ndim, k)
  sdbeta = matrix(0, ndim, k)
  npar = 0
  for (i in 1:k) {
    idx = c(1:ndim)[paridx[, i] == 1]
    resi = y[, i]
    if (length(idx) > 0) {
      xm = as.matrix(xmtx[, idx])
      npar = npar + dim(xm)[2]
      xpx = t(xm) %*% xm
      xpxinv = solve(xpx)
      xpy = t(xm) %*% as.matrix(y[, i], ne, 1)
      betai = xpxinv %*% xpy
      beta[idx, i] = betai
      resi = y[, i] - xm %*% betai
      nee = dim(xm)[2]
      sse = sum(resi * resi)/(Tn - p - nee)
      dd = diag(xpxinv)
      sdbeta[idx, i] = sqrt(dd * sse)
    }
    res = cbind(res, resi)
  }
  sse = t(res) %*% res/(Tn - p)
  aic = 0
  bic = 0
  hq = 0
  Phi = NULL
  Ph0 = NULL
  jst = 0
  if (include.mean) {
    Ph0 = beta[1, ]
    se = sdbeta[1, ]
    jst = 1
  }
  if (include.mean) {
    for (i in 1:k) {
      if (abs(Ph0[i]) > 1e-08) 
        npar = npar - 1
    }
  }
  for (i in 1:p) {
    phi = t(beta[(jst + 1):(jst + k), ])
    se = t(sdbeta[(jst + 1):(jst + k), ])
    jst = jst + k
    Phi = cbind(Phi, phi)
  }
  dd = det(sse)
  d1 = log(dd)
  aic = d1 + (2 * npar)/Tn
  bic = d1 + log(Tn) * npar/Tn
  hq = d1 + 2 * log(log(Tn)) * npar/Tn
  VAR <- list(data = x, cnst = include.mean, order = p, coef = beta, 
              aic = aic, bic = bic, hq = hq, residuals = res, secoef = sdbeta, 
              Sigma = sse, Phi = Phi, Ph0 = Ph0)
}
MTSirf = function (Phi, Sig, lag = 10, orth = FALSE) {
  if (!is.matrix(Phi)) 
    Phi = as.matrix(Phi)
  if (!is.matrix(Sig)) 
    Sig = as.matrix(Sig)
  k = nrow(Phi)
  m = ncol(Phi)
  p = floor(m/k)
  Si = diag(rep(1, k))
  wk = c(Si)
  awk = c(wk)
  acuwk = c(awk)
  if (p < 1) 
    p = 1
  if (lag < 1) 
    lag = 1
  for (i in 1:lag) {
    if (i <= p) {
      idx = (i - 1) * k
      tmp = Phi[, (idx + 1):(idx + k)]
    }
    else {
      tmp = matrix(0, k, k)
    }
    jj = i - 1
    jp = min(jj, p)
    if (jp > 0) {
      for (j in 1:jp) {
        jdx = (j - 1) * k
        idx = (i - j) * k
        w1 = Phi[, (jdx + 1):(jdx + k)]
        w2 = Si[, (idx + 1):(idx + k)]
        tmp = tmp + w1 %*% w2
      }
    }
    Si = cbind(Si, tmp)
    wk = cbind(wk, c(tmp))
    awk = awk + c(tmp)
    acuwk = cbind(acuwk, awk)
  }
  orSi = NULL
  wk1 = NULL
  awk1 = NULL
  acuwk1 = NULL
  if (orth) {
    m1 = chol(Sig)
    P = t(m1)
    wk1 = cbind(wk1, c(P))
    awk1 = wk1
    acuwk1 = wk1
    orSi = cbind(orSi, P)
    for (i in 1:lag) {
      idx = i * k
      w1 = Si[, (idx + 1):(idx + k)]
      w2 = w1 %*% P
      orSi = cbind(orSi, w2)
      wk1 = cbind(wk1, c(w2))
      awk1 = awk1 + c(w2)
      acuwk1 = cbind(acuwk1, awk1)
    }
  }
  tdx = c(1:(lag + 1)) - 1
  par(mfcol = c(k, k), mai = c(0.3, 0.3, 0.3, 0.3))
  if (orth) {
    gmax = max(wk1)
    gmin = min(wk1)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
    gmax = max(acuwk1)
    gmin = min(acuwk1)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
  } else {
    gmax = max(wk)
    gmin = min(wk)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
    gmax = max(acuwk)
    gmin = min(acuwk)
    cx = (gmax - gmin)/10
    gmax = gmax + cx
    gmin = gmin - cx
  }
  PHI1 = t(Si)
  k = ncol(PHI1)
  phi1 = list()
  for (i in 1:k) {
    ind = seq(i,nrow(PHI1),k)
    phi1[[i]] = PHI1[ind,]
  }
  
  return(phi1)
}
MTSfftGenFEVD = function(Phi,Sigma, n.ahead = 10, no.corr = F, range) {
  # Warn if the n.ahead is too low.
  if (n.ahead < 100) {
    warning("The frequency decomposition works with unconditional IRF. You have opted for 
            IRF with horizon lower than 100 periods. This might cause trouble, some frequencies
            might not be estimable depending on the bounds settings.")
  }
  # Get the unorthogonalized impulse responses (essentially Wold decomposition
  # coefficients thats why the name Phi.)
  
  Phi1 = MTSirf(Phi,Sig=Sigma,lag=n.ahead,orth=FALSE)
  # Get the Fourier transform of the impulse responses
  fftir1 <- lapply(Phi1, function(i) apply(i, 2, fft))
  # Transform them into shape we work with
  fftir1 <- lapply(1:(n.ahead+1), function(j) sapply(fftir1, function(i) i[j,]))
  if (no.corr) {
    Sigma <- diag(diag(Sigma))
  }
  denom1 <- diag(
    Re(
      Reduce('+', lapply(fftir1, function(i) 
        i %*% Sigma %*% t( Conj(i) ) / (n.ahead + 1)
      )[range]
      )
    )
  )
  # Compute the enumerator of the equation
  enum1 <- lapply(fftir1, function(i) 
    ( abs( i %*% Sigma ) )^2 / (n.ahead+1)
  )
  # Compute the fevd table be dividing the individual elements
  tab1 <- lapply(enum1, function(i) 
    sapply(1:nrow(i), function(j) 
      i[j, ] / ( denom1[j] * diag(Sigma) ) 
    ) 
  )
  # Compute the totals over the range for standardization
  tot <- apply(Reduce('+', tab1[range]), 2, sum)	
  # Standardize so that it sums up to one row-wise
  tab1 <- lapply(tab1, function(i) t(i)/tot)
  return(tab1)
}
MTSspilloverFft = function(func, Phi,Sigma, n.ahead, partition, no.corr = F) {
  f <- get( func )
  new_p <- getPartition( partition, n.ahead )
  range <- sort( unique( do.call(c, new_p) ) )
  decomp <- f(Phi,Sigma, n.ahead, no.corr = no.corr, range = range)
  for (i in 1:length(decomp)) {
    rownames(decomp[[i]]) <- colnames(decomp[[i]]) <- 1:ncol(Sigma)
  }
  tables <- lapply(new_p, function(j) Reduce('+', decomp[j]))
  
  return(structure(list(tables = tables, bounds = partition, date = NULL), class = "spillover_table"))
}
MTSspilloverBK12 = function(Phi,Sigma, n.ahead = 100, no.corr, partition) {
  return(MTSspilloverFft("MTSfftGenFEVD", Phi,Sigma, n.ahead, partition, no.corr = no.corr))
}
