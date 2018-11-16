library(raster)
library(fields)
library(mvtnorm)
library(RColorBrewer)
library(maptools)

source('utilityFuncs.R')
source('ReadTigerData.R')
source('contTrappingArray.mcmcFuncs.R')
source('contTrappingArray.mapOfRealizedDensity.R')



## ... define utility function for computing index of non-missing values in layer of raster object indexed by indForNA
ComputeIndexForNonmissingValues = function(rasterObj, indForNA) {
    rasterObj.valuesVector = NA
    if (nlayers(rasterObj)==1)  rasterObj.valuesVector = values(rasterObj)
    if (nlayers(rasterObj)>1) rasterObj.valuesVector = values(rasterObj)[,indForNA]
    !is.na(rasterObj.valuesVector)
}


###  Read in camera-trap data
resolution = 500  #  resolution (meters) used to create discrete grid of spatial domain
shapefile = 'NHstatespace-utm'
trapLocFile = 'trap_2015_coord.csv'
trapOpFile = 'trap_op_ind.csv'
sunriseANDsunsetFile = 'sunriseANDsunsetData.txt'
detectionFile = 'Capt_date_time_noHrep.csv'

d = ReadCameraTrapData(resolution, shapefile, trapLocFile, trapOpFile, sunriseANDsunsetFile, detectionFile)


#> names(d)
#[1] "sgrid"                     "traploc"                   "trapcov"                   "trapDaytimeEffort"        
#[5] "trapNitetimeEffort"        "captureFreqMatrix"         "captureTimeCovariateArray"
#
# $sgrid is "RasterStack" object with a bunch of "RasterLayer" objects (e.g. d$sgrid$TLS)
# $traploc is coordinates of 162 traps
# $trapcov is NULL
# $trapDaytimeEffort seems to be time on effort in day
# $trapNitetimeEffort seems to be time on effort at night
# $captureFreqMatrix is 86x162 matrix (presumably 86 individuals' counts on 162 traps)
# $captureTimeCovariateArray is 86x162x5 matrix ... looks like last dimension represents day=1/night=0, and 
#                            it is 5 long as max(d$captureFreqMatrix)==5.

plot(d$sgrid$TLS)
points(d$traploc,pch="+")









### Assign observed data to variables used in MCMC algorithm

## ... rescale units of spatial domain from meters to kilometers
distanceScaler = 1000
sgrid.extent = extent(as.matrix(extent(d$sgrid)) / distanceScaler)
sgrid = setExtent(d$sgrid, sgrid.extent)

## ... rescale units of trap locations from meters to kilometers
traploc = d$traploc / distanceScaler

y = d$captureFreqMatrix
Z = d$captureTimeCovariateArray
T.day = d$trapDaytimeEffort
T.nite = d$trapNitetimeEffort

parkIndicator = 'TLS'  # Tiger Land Status:  binary indicator equaling 1 if pixel is w/i park
ind = ComputeIndexForNonmissingValues(sgrid, parkIndicator)
sgrid.ncell = sum(ind)

X = matrix(1, nrow=sgrid.ncell)
W = matrix(1, nrow=nrow(traploc))



### Fit model to data
fit.contTrappingArray(niter=2500, niterInterval=50, sgrid, parkIndicator, traploc, y, X, W, Z, T.day, T.nite, FALSE)



### Compute summaries of posterior

## ... read Markov chains from files
mc = as.matrix(read.csv('mc.csv'))

## ...compute posterior means and quantiles
burnin = 500
out = -(1:burnin)
mc = mc[out, ]

psi.nite = exp(mc[,'alpha1'])
psi.day = exp(mc[,'alpha1'] + mc[,'xi'])
lambda0 = exp(mc[,'beta1'])

mc = cbind(mc, lambda0=lambda0, psiNite=psi.nite, psiDay=psi.day)

prob.quantiles = c(.50, .025, .975)  # for credible limits
post.stats = cbind(apply(mc,2,mean),  t(apply(mc, 2, quantile, probs=prob.quantiles)) )

prob.names = paste(as.character(100*prob.quantiles), '%', sep='')
post.names = dimnames(post.stats)[[1]]
dimnames(post.stats) = list(post.names, c('Mean', prob.names))

## .... calculate Monte Carlo standard errors
post.stats.MCSE = post.stats
for (i in 1:dim(mc)[2]) {
  postVec = as.vector(mc[,i])
  post.stats.MCSE[i,1] = mcse(postVec, mean)$se
  post.stats.MCSE[i,2] = mcse(postVec, median)$se
  post.stats.MCSE[i,3] = mcse(postVec, lowerQuantile)$se
  post.stats.MCSE[i,4] = mcse(postVec, upperQuantile)$se
}



### Print summaries of posterior
CR = '\n'
cat (CR, CR, 'Bayesian estimates of model parameters', CR, CR)
print(round(post.stats, 3))

cat (CR, CR, 'Monte Carlo SE of Bayesian estimates', CR, CR)
print(round(post.stats.MCSE,4))



### Compute posterior of map of mean realized abundance for visualizing map of mean realized density

spred = EstimateMapOfRealizedDensity(burnin, sgrid, parkIndicator, traploc, y, X, W, Z, T.day, T.nite)

## ... now compute map of individual activity centers within Nagarahole only
temp=raster(spred)
ind.nag = !(values(spred)[,parkIndicator]!=1 | is.na(values(spred)[,parkIndicator]))
values(temp)[ind.nag] = values(spred)[ind.nag, 'AllIndivs']
names(temp) = 'AllTigersInNagarahole'
spred = addLayer(spred, temp)

temp=raster(spred)
values(temp)[ind.nag] = values(spred)[ind.nag, 'SDofAllIndivs']
names(temp) = 'SDofAllTigersInNagarahole'
spred = addLayer(spred, temp)


pdf(file='MapOfEstimatedTigerAbund.pdf', onefile=FALSE)
par(mar=c(5,6,1,1))
##pal = brewer.pal(9, 'Blues')[-(1:3)]
pal = brewer.pal(9, 'YlOrRd')[-1]
plot(spred, 'AllTigersInNagarahole', col=pal, cex.axis=1.4, cex.lab=1.6, xlab='Easting (km)', ylab='Northing (km)', main='')
dev.off()

postscript(file='MapOfEstimatedTigerAbund.eps', onefile=FALSE)
par(mar=c(5,6,1,1))
##pal = brewer.pal(9, 'Blues')[-(1:3)]
pal = brewer.pal(9, 'YlOrRd')[-1]
plot(spred, 'AllTigersInNagarahole', col=pal, cex.axis=1.4, cex.lab=1.6, xlab='Easting (km)', ylab='Northing (km)', main='')
dev.off()

pdf(file='MapOfEstimatedSDofTigerAbund.pdf', onefile=FALSE)
par(mar=c(5,6,1,1))
##pal = brewer.pal(9, 'Blues')[-(1:3)]
pal = brewer.pal(9, 'YlOrRd')[-1]
plot(spred, 'SDofAllTigersInNagarahole', col=pal, cex.axis=1.4, cex.lab=1.6, xlab='Easting (km)', ylab='Northing (km)', main='')
dev.off()

## ... compute prediction of mean tiger abundance within Nagarahole only
spred.TigerAbundOfNag = sum(values(spred)[ind.nag, 'AllTigersInNagarahole'])
spred.areaOfNag = sum(ind.nag) * prod(res(spred))
cat('Total Abundance = ', spred.TigerAbundOfNag, '   Area = ', spred.areaOfNag, '\n')
