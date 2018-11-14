library(raster)
library(fields)
library(mvtnorm)

source('utilityFuncs.R')
source('contTrappingArray.simData.R')
source('contTrappingArray.mcmcFuncs.R')
source('contTrappingArrayOfIPP.mcmcFuncs.R')



### Write column headings of simulation output
CR = '\n'
names.Full = c('n0','sigma', 'xi', 'alpha1', 'alpha2', 'beta1', 'beta2')
cat(names.Full, sep=',', file='postMean.Full.csv')
cat(CR, file='postMean.Full.csv', append=TRUE)
cat(names.Full, sep=',', file='postLowerLimit.Full.csv')
cat(CR, file='postLowerLimit.Full.csv', append=TRUE)
cat(names.Full, sep=',', file='postUpperLimit.Full.csv')
cat(CR, file='postUpperLimit.Full.csv', append=TRUE)

names.Res = c('n0','sigma', 'alpha1', 'alpha2', 'beta1', 'beta2')
cat(names.Res, sep=',', file='postMean.Restricted.csv')
cat(CR, file='postMean.Restricted.csv', append=TRUE)
cat(names.Res, sep=',', file='postLowerLimit.Restricted.csv')
cat(CR, file='postLowerLimit.Restricted.csv', append=TRUE)
cat(names.Res, sep=',', file='postUpperLimit.Restricted.csv')
cat(CR, file='postUpperLimit.Restricted.csv', append=TRUE)

cat(c('N.ipp', 'n'), sep=',', file='sizesOfPopAndSample.csv')
cat(CR, file='sizesOfPopAndSample.csv', append=TRUE)


set.seed(11)
nsim = 9
for (sim in 1:nsim) {


###  Simulate camera-trap data
d = simData(beta=c(1.4, 0.8), alpha=c(-0.7, 1.0), sigma=0.4, xi=-1)

sgrid = d$sgrid
traploc = d$traploc
y = d$y
X = d$X
W = d$W
Z = d$Z
T.day = d$T.day
T.nite = d$T.nite
indForNA = 1  #  field in shape file used to determine cells that contain missing values

cat(c(d$N.ipp, nrow(y)), sep=',', file='sizesOfPopAndSample.csv', append=TRUE)
cat(CR, file='sizesOfPopAndSample.csv', append=TRUE)

    

### Fit model to data
fit.contTrappingArray(niter=1200, niterInterval=100, sgrid, indForNA, traploc, y, X, W, Z, T.day, T.nite, FALSE)


### Compute summaries of posterior

## ... read Markov chains from files
mc = as.matrix(read.csv('mc.csv'))

## ...estimate posterior means and quantiles
burnin = 200
out = -(1:burnin)
mc = mc[out, ]

postMean = colMeans(mc)
postLowerLimit = apply(mc, 2, quantile, probs=0.025)
postUpperLimit = apply(mc, 2, quantile, probs=0.975)


### Write estimates to output files
cat(postMean, sep=',', file='postMean.Full.csv', append=TRUE)
cat(CR, file='postMean.Full.csv', append=TRUE)
cat(postLowerLimit, sep=',', file='postLowerLimit.Full.csv', append=TRUE)
cat(CR, file='postLowerLimit.Full.csv', append=TRUE)
cat(postUpperLimit, sep=',', file='postUpperLimit.Full.csv', append=TRUE)
cat(CR, file='postUpperLimit.Full.csv', append=TRUE)





### Fit restricted model to data
fit.contTrappingArrayOfIPP(niter=1200, niterInterval=100, sgrid, indForNA, traploc, y, X, W, T.day+T.nite, FALSE)


### Compute summaries of posterior

## ... read Markov chains from files
mc = as.matrix(read.csv('mc.csv'))

## ...estimate posterior means and quantiles
burnin = 200
out = -(1:burnin)
mc = mc[out, ]

postMean = colMeans(mc)
postLowerLimit = apply(mc, 2, quantile, probs=0.025)
postUpperLimit = apply(mc, 2, quantile, probs=0.975)


### Write estimates to output files
cat(postMean, sep=',', file='postMean.Restricted.csv', append=TRUE)
cat(CR, file='postMean.Restricted.csv', append=TRUE)
cat(postLowerLimit, sep=',', file='postLowerLimit.Restricted.csv', append=TRUE)
cat(CR, file='postLowerLimit.Restricted.csv', append=TRUE)
cat(postUpperLimit, sep=',', file='postUpperLimit.Restricted.csv', append=TRUE)
cat(CR, file='postUpperLimit.Restricted.csv', append=TRUE)

}  # end of sim loop
