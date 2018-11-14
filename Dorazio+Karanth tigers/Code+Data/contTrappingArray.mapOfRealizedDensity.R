library(raster)
library(fields)


EstimateMapOfRealizedDensity = function(burnin, sgrid, indForNA, traploc, y, X, W, Z, T.day, T.nite) {

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

    
xyCountInRaster = function(loc, xyLimits) {
    ## returns a vector of counts of the locations in loc that occur within each value of xyLimits
    retVal = rep(0, nrow(xyLimits))
    for (i in 1:nrow(xyLimits)) {
        inRaster = loc[,1]>=xyLimits[i,'xmin'] & loc[,1]<=xyLimits[i,'xmax'] & loc[,2]>=xyLimits[i,'ymin'] & loc[,2]<=xyLimits[i,'ymax']
        retVal[i] = sum(inRaster)
    }
    retVal
}

    
    
## ... define functions used in Metropolis-Hastings sampling

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

    

                                        
### Compute posterior of map of realized abundance for visualizing map of realized density

## ... read Markov chains from files
out = -(1:burnin)
mc = as.matrix(read.csv('mc.csv'))
mc = mc[out, ]

## ... read Markov chains of locations from files
mc.sx = as.matrix(read.csv('mcForSlocx.csv'))
mc.sy = as.matrix(read.csv('mcForSlocy.csv'))
mc.sx = mc.sx[out, ]
mc.sy = mc.sy[out, ]


    
## ... form discrete approximation of spatial domain for computing maps of predictions
spred = sgrid
ind = ComputeIndexForNonmissingValues(spred, indForNA)
    
## ... for each draw of Markov chain, tabulate number of individual activity centers located within each raster cell for the n individuals detected at least once
spred.xylimits = xyLimitsOfRasterCells(spred, indForNA)
ndraws = nrow(mc.sx)
mc.n = matrix(nrow=ndraws, ncol=nrow(spred.xylimits))
for (draw in 1:ndraws) {
    loc = cbind(mc.sx[draw, ], mc.sy[draw, ])
    mc.n[draw, ] = xyCountInRaster(loc, spred.xylimits)
}

## ... now calculate posterior summary for map of realized density of n individuals
npost.mean = apply(mc.n, 2, mean)
temp = raster(spred)
values(temp)[ind] = npost.mean
names(temp) = 'DetectedIndivs'
spred = addLayer(spred, temp)



## ... for each draw of Markov chain, draw locations of n0 activity centers of individuals not detected by any of K traps
## (note:  these locations are drawn from a single conditional distribution with density function [s|y=0] )
spred.loc = xyFromCell(spred, 1:ncell(spred))
spred.loc = spred.loc[ind, ]
spred.ncell = sum(ind)
spred.cellsize = prod(res(spred))
spred.distmat = rdist(spred.loc, traploc)
## .... assign covariates at each location on discrete prediction grid
X.pred = X
alpha.names = paste('alpha', 1:ncol(W), sep='')
beta.names = paste('beta', 1:ncol(X), sep='')


ndraws = nrow(mc)
mc.n0 = matrix(nrow=ndraws, ncol=nrow(spred.xylimits))

for (draw in 1:ndraws) {
    alpha = mc[draw, alpha.names]
    beta = mc[draw, beta.names]
    sigma = mc[draw, 'sigma']
    xi = mc[draw, 'xi']
    piZeroLambdaValue = piZeroLambda(alpha, sigma, xi, beta, X.pred, spred.cellsize, W, T.nite, T.day, spred.distmat)
    probS = rep(NA, spred.ncell)
    for (i in 1:spred.ncell) {
        probS[i] = spred.cellsize * exp(logDensity.spred(spred.loc[i,], traploc, alpha, sigma, xi, beta, X.pred[i,], W, T.nite, T.day)) / piZeroLambdaValue
    }
    
    ## draw n0 locations from same conditional distribution of [s | y=0]
    ind.spred = sample(1:spred.ncell, size=mc[,'n0'], replace=TRUE, prob=probS)
    loc  = spred.loc[ind.spred, ]
    mc.n0[draw, ] = xyCountInRaster(loc, spred.xylimits)
}


## ... now calculate posterior summary for map of realized density of n0 individuals
n0post.mean = apply(mc.n0, 2, mean)
temp = raster(spred)
values(temp)[ind] = n0post.mean
names(temp) = 'UndetectedIndivs'
spred = addLayer(spred, temp)



## ... now calculate posterior summary for map of realized density of n+n0 individuals
Npost.mean = apply(mc.n+mc.n0, 2, mean)
temp = raster(spred)
values(temp)[ind] = Npost.mean
names(temp) = 'AllIndivs'
spred = addLayer(spred, temp)

Npost.var = apply((mc.n+mc.n0)^2, 2, mean) - Npost.mean^2
temp = raster(spred)
values(temp)[ind] = sqrt(Npost.var)
names(temp) = 'SDofAllIndivs'
spred = addLayer(spred, temp)
    
spred
    }
