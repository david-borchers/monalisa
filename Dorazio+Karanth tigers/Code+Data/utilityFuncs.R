# Define utility functions


logit = function(xx) { log(xx) - log(1-xx) }

expit = function(xx) { ifelse(xx<0, exp(xx) /(1 + exp(xx)),  1/(1 + exp(-xx))) }

lowerQuantile = function(z) {  quantile(z, probs=.025)  }

upperQuantile = function(z) {  quantile(z, probs=.975)  }



# compute Monte Carlo standard error of a functional (FUNC) of a Markov chain (x) using the Subsampling Bootstrap Method (SBM)
mcse = function(x, FUNC,  batchsize=floor(sqrt(length(x)))) {
  n = length(x)
  nb = n - batchsize + 1
  bval = rep(NA,nb)
  for (j in 1:nb) {
    ind = j:(batchsize+j-1)
    bval[j] = FUNC(x[ind])
  }
  ##  var.bval = var(bval)*(nb-1) * n * batchsize / ( (n-batchsize) * nb )    #  OLBM method
  var.bval = var(bval)*(nb-1) * batchsize / nb
  
  list(se=sqrt(var.bval / n), batchsize=batchsize)
}


# compute Monte Carlo standard error of Pr(X=xValue) from a Markov chain (x) of a discrete random variable using the Subsampling Bootstrap Method (SBM)
mcse.ProbDiscreteRV = function(x, xValue,  batchsize=floor(sqrt(length(x)))) {
 
  n = length(x)
  nb = n - batchsize + 1
  bval = rep(NA,nb)
  for (j in 1:nb) {
    ind = j:(batchsize+j-1)
    bval[j] = mean(x[ind]==xValue)
  }
  ##  var.bval = var(bval)*(nb-1) * n * batchsize / ( (n-batchsize) * nb )    #  OLBM method
  var.bval = var(bval)*(nb-1) * batchsize / nb
  
  list(se=sqrt(var.bval / n), batchsize=batchsize)
}



# compute Monte Carlo estimators of mean, variance, and standard deviation (and MCSE of each) for n independent draws of x, where distribution of x is not necessarily known
MonteCarloEsts = function(x) {
    n = length(x)
    xsq = x*x
    xmat = cbind(x, xsq)

    mu = apply(xmat, 2, mean)
    Sigma = var(xmat)

    xmean.est = mu[1]
    xmean.mcse = sqrt(Sigma[1,1]/n)

    xvar.est = mu[2] - mu[1]^2
    xvar.mcse =  sqrt( (Sigma[2,2] + 4*mu[1]*(mu[1]*Sigma[1,1] - Sigma[1,2]))/n )

    xsd.est = sqrt(xvar.est)
    xsd.mcse = sqrt( (0.25*Sigma[2,2] + mu[1]*mu[1]*Sigma[1,1] - mu[1]*Sigma[1,2]) / (n * (mu[2] - mu[1]^2)) )

    retval = c(as.vector(xmean.est), as.vector(xmean.mcse), as.vector(xvar.est), as.vector(xvar.mcse),
        as.vector(xsd.est), as.vector(xsd.mcse))
    names(retval) = c('mean.est', 'mean.mcse', 'var.est', 'var.mcse', 'sd.est', 'sd.mcse')
    retval
}
