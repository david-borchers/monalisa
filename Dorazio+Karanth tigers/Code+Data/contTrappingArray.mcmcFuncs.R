library(raster)
library(fields)
library(mvtnorm)


### Compute Markov chain used to estimate summaries of the posterior distribution

fit.contTrappingArray = function(niter=100, niterInterval=10, sgrid, indForNA, traploc, y, X, W, Z, T.day, T.nite,  adaptiveMHoutput=TRUE) {

    ## ... define utility function for computing index of non-missing values in layer of raster object indexed by indForNA
    ComputeIndexForNonmissingValues = function(rasterObj, indForNA) {
        rasterObj.valuesVector = NA
        if (nlayers(rasterObj)==1)  rasterObj.valuesVector = values(rasterObj)
        if (nlayers(rasterObj)>1) rasterObj.valuesVector = values(rasterObj)[,indForNA]
        !is.na(rasterObj.valuesVector)
    }

    
    ## ... compute variables used in model
    n = nrow(y)
    K = ncol(y)
    sgrid.loc = xyFromCell(sgrid, 1:ncell(sgrid))

    ind = ComputeIndexForNonmissingValues(sgrid, indForNA)
    sgrid.loc = sgrid.loc[ind, ]
    sgrid.ncell = nrow(sgrid.loc)
    sgrid.cellsize = prod(res(sgrid))
    sgrid.area = sgrid.ncell * sgrid.cellsize
    distmat = rdist(sgrid.loc, traploc)   # matrix of distances between trap locations and sgrid locations


    

## ... define utility functions
    
xyLimitsOfRasterCells = function(r, indForNA) {
    ## returns coordinates of corners of all cells of raster r whose value is not NA
    ind = ComputeIndexForNonmissingValues(r, indForNA)
    loc = xyFromCell(r, 1:ncell(r))
    loc = loc[ind, ]
    xyres = res(r)
    xyLimits = cbind(xmin=loc[,1]-xyres[1]/2, xmax=loc[,1]+xyres[1]/2, ymin=loc[,2]-xyres[2]/2, ymax=loc[,2]+xyres[2]/2)
    xyLimits
}

    
xyInRaster = function(loc, xyLimits) {
    ## returns TRUE if loc corresponds to any location included within the set of xyLimits
    inRaster = loc[1]>=xyLimits[,'xmin'] & loc[1]<=xyLimits[,'xmax'] & loc[2]>=xyLimits[,'ymin'] & loc[2]<=xyLimits[,'ymax']
    any(inRaster)
}

    
    
## ... define functions used in Metropolis-Hastings sampling

Lambda = function(beta, X, sgrid.cellsize) {
  lambda = exp(X %*% beta)
  sgrid.cellsize*sum(lambda)
}


piZeroLambda = function(alpha, sigma, xi, beta, X, sgrid.cellsize, W, T.nite, T.day, distmat) {  
  lambda = exp(X %*% beta)
  sgrid.ncell = dim(distmat)[1]
  K = dim(distmat)[2]
  logPsi = W %*% alpha
  logPsiMat = matrix(rep(logPsi, each=sgrid.ncell), ncol=K)
  logPhiMat = logPsiMat - 0.5*(distmat/sigma)^2
  Tmat = matrix(rep(T.nite + T.day*exp(xi), each=sgrid.ncell), ncol=K)
  PhiMat = exp(logPhiMat) * Tmat
  retVal = sgrid.cellsize * sum(lambda*exp(apply(-PhiMat,1,sum)))
  retVal
}



logDensity.beta = function(beta, alpha, sigma, xi, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, X.obs, beta.mu, beta.sigma) {
  logDensity.prior = sum(dnorm(beta, mean=beta.mu, sd=beta.sigma, log=TRUE))
  logDensity = -Lambda(beta, X, sgrid.cellsize) + n0*log(piZeroLambda(alpha, sigma, xi, beta, X, sgrid.cellsize, W, T.nite, T.day, distmat)) + sum(X.obs%*%beta)
  logDensity + logDensity.prior
}


neglogDensity.beta = function(beta, alpha, sigma, xi, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, X.obs, beta.mu, beta.sigma) {
   (-1)*logDensity.beta(beta, alpha, sigma, xi, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, X.obs, beta.mu, beta.sigma)
}

negGrad.beta = function(beta, alpha, sigma, xi, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, X.obs, beta.mu, beta.sigma) {
  lambda = exp(X %*% beta)
  sgrid.ncell = dim(distmat)[1]
  K = dim(distmat)[2]
  logPsi = W %*% alpha
  logPsiMat = matrix(rep(logPsi, each=sgrid.ncell), ncol=K)
  logPhiMat = logPsiMat - 0.5*(distmat/sigma)^2
  Tmat = matrix(rep(T.nite + T.day*exp(xi), each=sgrid.ncell), ncol=K)
  PhiMat = exp(logPhiMat) * Tmat
  piZeroLambda = sgrid.cellsize * sum(lambda*exp(apply(-PhiMat,1,sum)))
  cterm = exp(apply(-PhiMat,1,sum))
  temp = as.vector(lambda) * (cterm*n0/piZeroLambda - 1)
  grad = -(beta-beta.mu)/(beta.sigma^2) + apply(X.obs,2,sum) + sgrid.cellsize * apply(X*temp,2,sum)
  (-1) * grad
}



logDensity.alphaANDlogsigma = function(param, beta, y, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, distmat.obs, Z, alpha.mu, alpha.sigma, sigma.nu, sigma.scale, xi.mu, xi.sigma) {
    sigma = exp(param[1])
    xi = param[2]
    alpha = param[-(1:2)]
    logDensity.alphaprior = sum(dnorm(alpha, mean=alpha.mu, sd=alpha.sigma, log=TRUE))
    logDensity.xiprior = dnorm(xi, mean=xi.mu, sd=xi.sigma, log=TRUE)
    logDensity.sigmaprior = -0.5 * (sigma.nu+1) * log(1 + (1/sigma.nu)*(sigma/sigma.scale)^2)
    logDensity.prior = logDensity.alphaprior + logDensity.sigmaprior + logDensity.xiprior + log(sigma)
    logPsi = W %*% alpha
    n = dim(y)[1]
    K = dim(distmat)[2]
    logPsiMat = matrix(rep(logPsi, each=n), ncol=K)
    logPhiMat = logPsiMat - 0.5*(distmat.obs/sigma)^2
    Tmat = matrix(rep(T.nite + T.day*exp(xi), each=n), ncol=K)
    PhiMat = exp(logPhiMat) * Tmat
    logDensity = n0*log(piZeroLambda(alpha, sigma, xi, beta, X, sgrid.cellsize, W, T.nite, T.day, distmat)) - sum(PhiMat)
    for (i in 1:n) {
        for (k in 1:K) {
            if (y[i,k]>0) {
                logDensity = logDensity + y[i,k]*logPhiMat[i,k] + sum(xi*Z[i,k, 1:y[i,k]])
            }
        }
    }
    logDensity + logDensity.prior
}

neglogDensity.alphaANDlogsigma = function(param, beta, y, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, distmat.obs, Z, alpha.mu, alpha.sigma, sigma.nu, sigma.scale, xi.mu, xi.sigma) {
   (-1)*logDensity.alphaANDlogsigma(param, beta, y, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, distmat.obs, Z, alpha.mu, alpha.sigma, sigma.nu, sigma.scale, xi.mu, xi.sigma)
}



## this conditional is used for predicting locations of n0 individuals not observed at any of K traps
logDensity.spred = function(s, Xloc, alpha, sigma, xi, beta, X.pred, W, T.nite, T.day) {
    logLambda = sum(X.pred*beta)
    smat = matrix(s, nrow=1)
    distVec = as.vector(rdist(Xloc, smat))
    logOfPhi = W %*% alpha - 0.5*(distVec/sigma)^2
    T = T.nite + T.day*exp(xi)
    Phi = T*exp(logOfPhi)
    retVal = logLambda - sum(Phi)
    retVal
}




logDensity.s = function(s, yvec, Xloc, alpha, sigma, xi, beta, X.obs, W, T.nite, T.day, Zmat) {
    logLambda = sum(X.obs*beta)
    smat = matrix(s, nrow=1)
    distVec = as.vector(rdist(Xloc, smat))
    logOfPhi = W %*% alpha - 0.5*(distVec/sigma)^2
    T = T.nite + T.day*exp(xi)
    Phi = T*exp(logOfPhi)
    retVal = logLambda - sum(Phi)
    for (k in 1:length(yvec)) {
        if (yvec[k]>0) {
            retVal = retVal + yvec[k]*logOfPhi[k] + sum(xi*Zmat[k,1:yvec[k]])
        }
    }       
    retVal
}

neglogDensity.s = function(s, yvec, Xloc, alpha, sigma, xi, beta, X.obs, W, T.nite, T.day, Zmat) {
    (-1)*logDensity.s(s, yvec, Xloc, alpha, sigma, xi, beta, X.obs, W, T.nite, T.day, Zmat)
}

negGrad.s = function(s, yvec, Xloc, alpha, sigma, xi, beta, X.obs, W, T.nite, T.day, Zmat) {
    ## this code assumes that X.obs does not include s as a regressor
    smat = matrix(s, nrow=1)
    distVec = as.vector(rdist(Xloc, smat))
    logOfPhi = W %*% alpha - 0.5*(distVec/sigma)^2
    T = T.nite + T.day*exp(xi)
    cterm = yvec - T*exp(logOfPhi)
    g1 =  sum(cterm * (smat[,1] - Xloc[,1]) )
    g2 =  sum(cterm * (smat[,2] - Xloc[,2]) )
    g = matrix(c(g1, g2), ncol=1) / (sigma*sigma)
    g
}

negHess.s = function(s, yvec, Xloc, alpha, sigma, xi, beta, X.obs, W, T.nite, T.day, Zmat) {
    ## this code assumes that X.obs does not include s as a regressor
    smat = matrix(s, nrow=1)
    distVec = as.vector(rdist(Xloc, smat))
    logOfPhi = W %*% alpha - 0.5*(distVec/sigma)^2
    T = T.nite + T.day*exp(xi)
    cterm = yvec - T*exp(logOfPhi)

    h11 =  sum( (yvec-cterm) * (smat[,1] - Xloc[,1]) * (smat[,1] - Xloc[,1]) / (sigma*sigma) + cterm )
    h22 =  sum( (yvec-cterm) * (smat[,2] - Xloc[,2]) * (smat[,2] - Xloc[,2]) / (sigma*sigma) + cterm )
    h12 =  sum( (yvec-cterm) * (smat[,1] - Xloc[,1]) * (smat[,2] - Xloc[,2]) / (sigma*sigma)         )
    h21 = h12
    H = matrix(c(h11, h12, h12, h22), ncol=2) / (sigma*sigma)
    H
}


##  Begin MCMC
    
## .... assign prior hyperparameter values
beta.mu = rep(0, ncol(X))
beta.sigma = rep(10, ncol(X))

alpha.mu = rep(0, ncol(W))
alpha.sigma = rep(10, ncol(W))

sigma.nu = 2
sigma.scale = 10

xi.mu = 0
xi.sigma = 10


## .... initialize Gibbs sampler for model that retains activity centers of observed individuals

## .... initialize Markov chain
n0 = 1.0*n

beta = rep(0, ncol(X))
beta0 = log((n+n0) / sgrid.area)
beta[1] = beta0

alpha = rep(0, ncol(W))
alpha0 = log( 10 * mean(apply(y,2,mean)/(T.nite+T.day)) )
alpha[1] = alpha0
xi = 0

distmat.traps = rdist(traploc, traploc)
distvec.traps = distmat.traps[col(distmat.traps) < row(distmat.traps)]
sigma = quantile(distvec.traps, probs=0.25)

relCapFreq = y/apply(y,1,sum)
loc = relCapFreq %*% traploc  # activity centers of observed individuals
X.obs = matrix(nrow=n, ncol=dim(X)[2]) # values of X at activity centers
distmat.obs = rdist(loc, traploc)
for (i in 1:n) {
    ## .... find location on sgrid that is closest to loc[i,] and assign its value of X
    distVec = as.vector(rdist(sgrid.loc, matrix(loc[i,],nrow=1)))
    indVec = sort.int(distVec, method='quick', index.return=TRUE)$ix
    X.obs[i,] = X[indVec[1], ]
}


## .... assign values for adaptive, random-walk Metropolis sampling
batch = 0
batchsize = 50
mindeltaLogSigmaForRW = 0.01
delta = mindeltaLogSigmaForRW


sigmaOfProposal.s = rep(0.1*sigma, n)
naccept.s = rep(0, n)
nattempt.s = rep(0, n)
sgrid.xylimits = xyLimitsOfRasterCells(sgrid, indForNA)


sigmaOfProposal.beta = 0.1
naccept.beta = 0
nattempt.beta = 0

sigmaOfProposal.alpha = 0.1
naccept.alpha = 0
nattempt.alpha = 0

## .... compute Gibbs draws

continueGibbs = TRUE
draw = 0
ndraws = niter

CR = '\n'
cat('Begin Gibbs sampling:', CR, CR)


alpha.names = paste('alpha', 1:dim(W)[2], sep='')
beta.names = paste('beta', 1:dim(X)[2], sep='')

cat(c('n0','sigma', 'xi', alpha.names, beta.names), sep=',', file='mc.csv')
cat(CR, file='mc.csv', append=TRUE)

locx.names = paste('s.x', 1:n, sep='')
cat(locx.names, sep=',', file='mcForSlocx.csv')
cat(CR, file='mcForSlocx.csv', append=TRUE)
locx.names = paste('s.y', 1:n, sep='')
cat(locx.names, sep=',', file='mcForSlocy.csv')
cat(CR, file='mcForSlocy.csv', append=TRUE)


start.time = Sys.time()

while(continueGibbs) {
  
  draw = draw + 1
  drawinterval = niterInterval
  if (draw == round(draw/drawinterval)*drawinterval)  {
      end.time = Sys.time()
      elapsed.time = difftime(end.time, start.time, units='mins')
      cat('..... drawing sample #', draw, ' after ', elapsed.time, ' minutes', CR)
  }

  ## update the increment/decrement for adaptive Metropolis-Hastings samplers
  if (floor(draw/batchsize)*batchsize == draw) {
    batch = batch + 1
    if (1/sqrt(batch) < mindeltaLogSigmaForRW)  delta = 1/sqrt(batch)
  }
  
  
  ##  draw n0
  n0 = rpois(1, piZeroLambda(alpha, sigma, xi, beta, X, sgrid.cellsize, W, T.nite, T.day, distmat))
  

  ##  draw loc (= activity center of each observed individual)
  for (i in 1:n) {
      ## ... find mode and hessian of unnormalized conditional density function
      fit = optim(par=loc[i,], fn=neglogDensity.s, gr=negGrad.s, method='BFGS', hessian=FALSE, yvec=y[i,], Xloc=traploc, alpha=alpha, sigma=sigma, xi=xi, beta=beta, X.obs=X.obs[i,], W=W, T.nite=T.nite, T.day=T.day, Zmat=Z[i,,])

      fit.gradient = negGrad.s(fit$par, yvec=y[i,], Xloc=traploc, alpha=alpha, sigma=sigma, xi=xi, beta=beta, X.obs=X.obs[i,], W=W, T.nite=T.nite, T.day=T.day, Zmat=Z[i,,])
      fit.hessian = negHess.s(fit$par, yvec=y[i,], Xloc=traploc, alpha=alpha, sigma=sigma, xi=xi, beta=beta, X.obs=X.obs[i,], W=W, T.nite=T.nite, T.day=T.day, Zmat=Z[i,,])
      
      if (fit$conv==0 & sum(abs(fit.gradient))<0.0001*abs(fit$value) & all(diag(fit.hessian)>0) )  {  # test validity of MH proposal

          ## ... draw candidate for Metropolis-Hastings sampler
          meanOfProposal = fit$par
          ev = eigen(fit.hessian)$values
          if (any(ev<0) | ev[2]/ev[1]<0.0001) {
              fit.hessian[1,2] = 0
              fit.hessian[2,1] = 0
          }
          SigmaOfProposal = chol2inv(chol(fit.hessian))
          loc.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))
          ## ... find X value for candidate
          distVec = as.vector(rdist(sgrid.loc, matrix(loc.cand,nrow=1)))
          indVec = sort.int(distVec, method='quick', index.return=TRUE)$ix
          X.obs.cand = X[indVec[1], ]

          ## ... only use candidate if it lies within the spatial domain
          if (xyInRaster(loc.cand, sgrid.xylimits)) {
              logL.loc.cand = logDensity.s(loc.cand, yvec=y[i,], Xloc=traploc, alpha=alpha, sigma=sigma, xi=xi, beta=beta, X.obs=X.obs.cand, W=W, T.nite=T.nite, T.day=T.day, Zmat=Z[i,,])
              logL.loc = logDensity.s(loc[i,], yvec=y[i,], Xloc=traploc, alpha=alpha, sigma=sigma, xi=xi, beta=beta, X.obs=X.obs[i,], W=W, T.nite=T.nite, T.day=T.day, Zmat=Z[i,,])
              logQ.loc.cand = dmvnorm(loc.cand, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
              logQ.loc = dmvnorm(loc[i,], mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
              logR = logL.loc.cand - logL.loc - logQ.loc.cand + logQ.loc
              if (logR >= 0 | runif(1,0,1) <= exp(logR)) {
                  loc[i,] = loc.cand
                  X.obs[i,] = X.obs.cand
              }
              else {
                  ## ... draw second candidate for delayed random-walk Metropolis sampler
                  SigmaOfProposal2 = diag(rep(sigmaOfProposal.s[i]^2, 2))
                  meanOfProposal2 = loc[i,]
                  loc.cand2 = as.vector(rmvnorm(1, mean=meanOfProposal2, sigma=SigmaOfProposal2))
                  ## ... find X value for candidate
                  distVec = as.vector(rdist(sgrid.loc, matrix(loc.cand2,nrow=1)))
                  indVec = sort.int(distVec, method='quick', index.return=TRUE)$ix
                  X.obs.cand2 = X[indVec[1], ]

                  ## ... only use second candidate if it lies within the spatial domain
                  if (xyInRaster(loc.cand2, sgrid.xylimits)) {
                      logL.loc.cand2 = logDensity.s(loc.cand2, yvec=y[i,], Xloc=traploc, alpha=alpha, sigma=sigma, xi=xi, beta=beta, X.obs=X.obs.cand2, W=W, T.nite=T.nite, T.day=T.day, Zmat=Z[i,,])
                      logQ.loc.cand2 = dmvnorm(loc.cand2, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
                      logR2 = logL.loc.cand - logL.loc.cand2 - logQ.loc.cand + logQ.loc.cand2
                      if (logR2 < 0) {
                          logAlpha2 = logL.loc.cand2 - logL.loc - log(1-exp(logR)) + log(1-exp(logR2))
                          if (logAlpha2 >=0 | runif(1,0,1) <= exp(logAlpha2)) {
                              loc[i,] = loc.cand2
                              X.obs[i,] = X.obs.cand2
                              naccept.s[i] = naccept.s[i] + 1
                          }
                      }
                  }
                  nattempt.s[i] = nattempt.s[i] + 1
              }
          }
         
      } else {
          ## ... draw candidate for random-walk Metropolis sampler
          SigmaOfProposal = diag(rep(sigmaOfProposal.s[i]^2, 2))
          meanOfProposal = loc[i,]
          loc.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))
          ## ... find X value for candidate
          distVec = as.vector(rdist(sgrid.loc, matrix(loc.cand,nrow=1)))
          indVec = sort.int(distVec, method='quick', index.return=TRUE)$ix
          X.obs.cand = X[indVec[1], ]

          ## ... only use candidate if it lies within the spatial domain
          if (xyInRaster(loc.cand, sgrid.xylimits)) {
              logL.loc.cand = logDensity.s(loc.cand, yvec=y[i,], Xloc=traploc, alpha=alpha, sigma=sigma, xi=xi, beta=beta, X.obs=X.obs.cand, W=W, T.nite=T.nite, T.day=T.day, Zmat=Z[i,,])
              logL.loc = logDensity.s(loc[i,], yvec=y[i,], Xloc=traploc, alpha=alpha, sigma=sigma, xi=xi, beta=beta, X.obs=X.obs[i,], W=W, T.nite=T.nite, T.day=T.day, Zmat=Z[i,,])
              logR = logL.loc.cand - logL.loc
              if (logR >= 0 | runif(1,0,1) <= exp(logR)) {
                  loc[i,] = loc.cand
                  X.obs[i,] = X.obs.cand
                  naccept.s[i] = naccept.s[i] + 1
              }
          }
          nattempt.s[i] = nattempt.s[i] + 1
      }
      
      ## ... update proposal variance for random-walk Metropolis sampler
      if (floor(draw/batchsize)*batchsize == draw) {
          if (nattempt.s[i]>0 & adaptiveMHoutput) cat("i = ", i, " s(i) acceptRate = ", naccept.s[i]/nattempt.s[i], " in ", nattempt.s[i], " attempts", CR)
          SigmaDiff = ifelse(naccept.s[i] > 0.234*nattempt.s[i], exp(2*delta), exp(-2*delta))
          sigmaOfProposal.s[i] = sigmaOfProposal.s[i] * SigmaDiff
          nattempt.s[i] = 0  # reset counter for next batch
          naccept.s[i] = 0   # reset counter for next batch
      }
  }
  distmat.obs = rdist(loc, traploc)



  ## draw beta
  ## ... find mode and hessian of unnormalized conditional density function
  fit = optim(par=beta, fn=neglogDensity.beta, gr=negGrad.beta, method='BFGS', hessian=TRUE, alpha=alpha, sigma=sigma, xi=xi, X=X, sgrid.cellsize=sgrid.cellsize, W=W, T.nite=T.nite, T.day=T.day, distmat=distmat, n0=n0, X.obs=X.obs, beta.mu=beta.mu, beta.sigma=beta.sigma)

  fit.gradient = negGrad.beta(fit$par, alpha, sigma, xi, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, X.obs, beta.mu, beta.sigma)
  fit.hessian = fit$hessian

  ## ... draw candidate using multivariate normal distribution as proposal  
  if (fit$conv==0 & sum(abs(fit.gradient))<0.0001*abs(fit$value) & all(diag(fit.hessian)>0) )  {  # test validity of proposal
      meanOfProposal = fit$par
      ev = eigen(fit.hessian)$values
      if (any(ev<0) | (length(ev)>1 & ev[2]/ev[1]<0.0001)) { fit.hessian = diag(diag(fit.hessian)) }
      SigmaOfProposal = chol2inv(chol(fit.hessian))
      beta.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))

      logL.beta.cand = logDensity.beta(beta.cand, alpha, sigma, xi, X, sgrid.cellsize, W, T.nite=T.nite, T.day=T.day, distmat, n0, X.obs, beta.mu, beta.sigma)
      logL.beta = logDensity.beta(beta, alpha, sigma, xi, X, sgrid.cellsize, W, T.nite=T.nite, T.day=T.day, distmat, n0, X.obs, beta.mu, beta.sigma)
      logQ.beta.cand = dmvnorm(beta.cand, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
      logQ.beta =  dmvnorm(beta, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
      logR = logL.beta.cand - logL.beta - logQ.beta.cand + logQ.beta
      if (logR >= 0 | runif(1,0,1) <= exp(logR)) {
          beta = beta.cand
      }
      else {
          ## ... draw second candidate for delayed random-walk Metropolis sampler
          SigmaOfProposal2 = diag(rep(sigmaOfProposal.beta^2, length(beta)), nrow=length(beta), ncol=length(beta))
          meanOfProposal2 = beta
          beta.cand2 = as.vector(rmvnorm(1, mean=meanOfProposal2, sigma=SigmaOfProposal2))
          logL.beta.cand2 = logDensity.beta(beta.cand2, alpha, sigma, xi, X, sgrid.cellsize, W, T.nite=T.nite, T.day=T.day, distmat, n0, X.obs, beta.mu, beta.sigma)
          logQ.beta.cand2 = dmvnorm(beta.cand2, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
          logR2 = logL.beta.cand - logL.beta.cand2 - logQ.beta.cand + logQ.beta.cand2
          if (logR2 < 0) {
              logAlpha2 = logL.beta.cand2 - logL.beta - log(1-exp(logR)) + log(1-exp(logR2))
              if (logAlpha2 >=0 | runif(1,0,1) <= exp(logAlpha2)) {
                  beta = beta.cand2
                  naccept.beta = naccept.beta + 1
              }
          }
          nattempt.beta = nattempt.beta + 1
      }
  }
  else {
      ## ... draw candidate for random-walk Metropolis sampler
      SigmaOfProposal = diag(rep(sigmaOfProposal.beta^2, length(beta)), nrow=length(beta), ncol=length(beta))
      meanOfProposal = beta
      beta.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))
      logL.beta.cand = logDensity.beta(beta.cand, alpha, sigma, xi, X, sgrid.cellsize, W, T.nite=T.nite, T.day=T.day, distmat, n0, X.obs, beta.mu, beta.sigma)
      logL.beta = logDensity.beta(beta, alpha, sigma, xi, X, sgrid.cellsize, W, T.nite=T.nite, T.day=T.day, distmat, n0, X.obs, beta.mu, beta.sigma)
      logR = logL.beta.cand - logL.beta
      if (logR >= 0 | runif(1,0,1) <= exp(logR)) {
          beta = beta.cand
          naccept.beta = naccept.beta + 1
      }
      nattempt.beta = nattempt.beta + 1
  }
  ## ... update proposal variance for random-walk Metropolis sampler
  if (floor(draw/batchsize)*batchsize == draw) {
      if (nattempt.beta>0 & adaptiveMHoutput) cat(" beta acceptRate = ", naccept.beta/nattempt.beta, " in ", nattempt.beta, " attempts", CR)
      SigmaDiff = ifelse(naccept.beta > 0.234*nattempt.beta, exp(2*delta), exp(-2*delta))
      sigmaOfProposal.beta = sigmaOfProposal.beta * SigmaDiff
      nattempt.beta = 0  # reset counter for next batch
      naccept.beta = 0   # reset counter for next batch
  }




   
  ## draw alpha and sigma and xi
  ## ... find mode and hessian of unnormalized conditional density function
  param = c(log(sigma), xi, alpha)
  fit = optim(par=param, fn=neglogDensity.alphaANDlogsigma, method='BFGS', hessian=TRUE, beta=beta, y=y, X=X, sgrid.cellsize=sgrid.cellsize, W=W, T.nite=T.nite, T.day=T.day, distmat=distmat, n0=n0, distmat.obs=distmat.obs, Z=Z, alpha.mu=alpha.mu, alpha.sigma=alpha.sigma, sigma.nu=sigma.nu, sigma.scale=sigma.scale, xi.mu=xi.mu, xi.sigma=xi.sigma)
  fit.hessian = fit$hessian

  ## ... draw candidate using  normal distribution as proposal  
  if (fit$conv==0 & all(diag(fit.hessian)>0) )  {  # test validity of proposal
      meanOfProposal = fit$par
      ev = eigen(fit.hessian)$values
      if (any(ev<0) | ev[2]/ev[1]<0.0001) { fit.hessian = diag(diag(fit.hessian)) }
      SigmaOfProposal = chol2inv(chol(fit.hessian))
      param.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))

      logL.param.cand = logDensity.alphaANDlogsigma(param.cand, beta, y, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, distmat.obs, Z, alpha.mu, alpha.sigma, sigma.nu, sigma.scale, xi.mu, xi.sigma)
      logL.param = logDensity.alphaANDlogsigma(param, beta, y, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, distmat.obs, Z, alpha.mu, alpha.sigma, sigma.nu, sigma.scale, xi.mu, xi.sigma)
      logQ.param.cand = dmvnorm(param.cand, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
      logQ.param =  dmvnorm(param, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
      logR = logL.param.cand - logL.param - logQ.param.cand + logQ.param
      if (logR >= 0 | runif(1,0,1) <= exp(logR)) {
          sigma = exp(param.cand[1])
          xi = param.cand[2]
          alpha = param.cand[-(1:2)]
      }
      else {
          ## ... draw second candidate for delayed random-walk Metropolis sampler
          SigmaOfProposal2 = diag(rep(sigmaOfProposal.alpha^2, length(param)))
          meanOfProposal2 = param
          param.cand2 = as.vector(rmvnorm(1, mean=meanOfProposal2, sigma=SigmaOfProposal2))
          logL.param.cand2 = logDensity.alphaANDlogsigma(param.cand2, beta, y, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, distmat.obs, Z, alpha.mu, alpha.sigma, sigma.nu, sigma.scale, xi.mu, xi.sigma)
          logQ.param.cand2 = dmvnorm(param.cand2, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
          logR2 = logL.param.cand - logL.param.cand2 - logQ.param.cand + logQ.param.cand2
          if (logR2 < 0) {
              logAlpha2 = logL.param.cand2 - logL.param - log(1-exp(logR)) + log(1-exp(logR2))
              if (logAlpha2 >=0 | runif(1,0,1) <= exp(logAlpha2)) {
                  sigma = exp(param.cand[1])
                  xi = param.cand[2]
                  alpha = param.cand[-(1:2)]
                  naccept.alpha = naccept.alpha + 1
              }
          }
          nattempt.alpha = nattempt.alpha + 1
      }
   }
  else {
      ## ... draw candidate for random-walk Metropolis sampler
      meanOfProposal = param
      SigmaOfProposal = diag(rep(sigmaOfProposal.alpha^2, length(param)))
      param.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))
      logL.param.cand = logDensity.alphaANDlogsigma(param.cand, beta, y, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, distmat.obs, Z, alpha.mu, alpha.sigma, sigma.nu, sigma.scale, xi.mu, xi.sigma)
      logL.param = logDensity.alphaANDlogsigma(param, beta, y, X, sgrid.cellsize, W, T.nite, T.day, distmat, n0, distmat.obs, Z, alpha.mu, alpha.sigma, sigma.nu, sigma.scale, xi.mu, xi.sigma)
      logR = logL.param.cand - logL.param
      if (logR >= 0 | runif(1,0,1) <= exp(logR)) {
          sigma = exp(param.cand[1])
          xi = param.cand[2]
          alpha = param.cand[-(1:2)]
          naccept.alpha = naccept.alpha + 1
      }
      nattempt.alpha = nattempt.alpha + 1
  }
  ## ... update proposal variance for random-walk Metropolis sampler
  if (floor(draw/batchsize)*batchsize == draw) {
      if (nattempt.alpha>0 & adaptiveMHoutput) cat(" alpha acceptRate = ", naccept.alpha/nattempt.alpha, " in ", nattempt.alpha, " attempts", CR)
      SigmaDiff = ifelse(naccept.alpha > 0.234*nattempt.alpha, exp(2*delta), exp(-2*delta))
      sigmaOfProposal.alpha = sigmaOfProposal.alpha * SigmaDiff
      nattempt.alpha = 0  # reset counter for next batch
      naccept.alpha = 0   # reset counter for next batch
  }
  
  

  
  ## output results of Gibbs iteration to files
  cat(c(n0, sigma, xi, alpha, beta), sep=',', file='mc.csv', append=TRUE)
  cat(CR, file='mc.csv', append=TRUE)
  cat(loc[,1], sep=',', file='mcForSlocx.csv', append=TRUE)
  cat(CR, file='mcForSlocx.csv', append=TRUE)
  cat(loc[,2], sep=',', file='mcForSlocy.csv', append=TRUE)
  cat(CR, file='mcForSlocy.csv', append=TRUE)

  if (draw == ndraws) {
    cat('Completed ', ndraws, ' draws of MCMC algorithm', CR)
    numOfDraws = as.integer(readline(prompt='Enter additional number of MCMC draws -> '))

    if (numOfDraws == 0) {
      continueGibbs = FALSE
    }
    else {
      ndraws = ndraws + numOfDraws
    }
  }

}  # end of while loop

invisible()
}  # end of fitting model

