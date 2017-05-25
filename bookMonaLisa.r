load("/Users/dlb/Research/SCR book/bookV1/monalisamesh.RData")

plotcovariate(mlmesh,covariate="D",contour=FALSE,col=gray.colors(50))

# "habitat" covariate
cD = covariates(mlmesh)$D-mean(covariates(mlmesh)$D)
habitat = (log((cD+1)))
plot(cD,habitat)
covariates(mlmesh)$habitat = habitat
plotcovariate(mlmesh,covariate="habitat",contour=FALSE,col=gray.colors(50))

# Generate a population
Dcov = 0.9-covariates(mlmesh)$D # invert to make dark low density (works better with estimate plot)

cDcov = cos(4*pi*Dcov/max(Dcov))
plot(Dcov,cDcov)
hist(cDcov)
mlDbase = exp(max(cDcov)-cDcov)*20 - 20
#mlDbase = exp(max(Dcov)-Dcov)*20 - 20
scaleup = 25
shiftup =750
mlD = shiftup + mlDbase*scaleup
pop=sim.popn(D=mlD, core=mlmesh, model2D="IHP", seed=12345)
N=dim(pop)[1];N # check simulated population size
# plot mesh with individuals' locations and detectors overlaid
plot(mlmesh, border=1,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.05)


# make a grid of detectors
dets=make.grid(nx=10,ny=10,spacex=8,spacey=6,originxy=c(22,13),detector="count")
dets1=make.grid(nx=7,ny=7,spacex=5,spacey=4,originxy=c(55,33),detector="count")
#plot(mlmesh, border=1,dots=FALSE, col="white", meshcol="gray")
plot(dets1,add=TRUE)

# Make mask for this grid
mask1=make.mask(dets1,buffer=12,type="trapbuffer")
names(covariates(mask1))
mask1 = addCovariates(mask1,mlmesh,"habitat") # add habitat as covariate

## Generate capture histories
## first set number of occasions
lambda0=0.5;sigma=2.5
nt <- 1
capthist=sim.capthist(dets1,popn=pop, detectfn="HHN",detectpar=list(lambda0=lambda0,sigma=sigma), noccasions=nt, nsessions=1,seed=12345)
summary(capthist)
n=dim(capthist)[1];n
plot(mlmesh, border=5,dots=FALSE, col="white", meshcol="gray")
points(pop$x,pop$y,pch=19,cex=0.25)
plot(dets1,add=TRUE)
plot(capthist, border=sigma, tracks=TRUE, varycol=FALSE,gridlines=FALSE,rad=3,add=TRUE)

# Constant D estimate:
smfit0=secr.fit(capthist,mask=mask1)
smfit0.Dhat = predictDsurface(smfit0)
names(covariates(smfit0.Dhat))
plot(smfit0.Dhat,contour=FALSE,col=gray.colors(40))
plot(dets1,add=TRUE)
region.N(smfit0)

# Habitat D estimate:
smfit.hab=secr.fit(capthist,mask=mask1,model=D~habitat)
smfit.hab.Dhat = predictDsurface(smfit.hab)
plotcovariate(mask1,covariate="habitat",col=gray.colors(40),contour=FALSE)
plot(smfit.hab.Dhat,contour=FALSE,col=gray.colors(40))
plot(dets1,add=TRUE)
region.N(smfit.hab)


# estimate:
smfit1=secr.fit(capthist,model=list(D~s(x,y,k=10)),mask=mask1)
smfit1.Dhat = predictDsurface(smfit1)
names(covariates(smfit1.Dhat))
plotcovariate(smfit1.Dhat,covariate="D.0",contour=FALSE,col=gray.colors(40))
plot(dets1,add=TRUE)
region.N(smfit1)
