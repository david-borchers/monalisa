# Code to produce Figure 5

# Libraries we need
library(nimble)
library(coda)
library(nimbleSCR)
library(spatstat)
# Objects we need
load("../output/revision/mona-inputs.RData")
# Functions we need
source("Functions.R")

## ---------------------------------------------------------------------------------------
# Creating the data objects we need for all MCMC samples
## ---------------------------------------------------------------------------------------

## Function we will use
# The only argument we provide is 'j': this will be a value in {1, 2, 3} and represents the index of the objects we want to work with from the RData objects we have loaded in. If want objects generated using 18, 52 or 111 sampling occasions, j=1,2,3 respectively.
organise.data = function(j) {
  # Number of sampling occasions used for simulated data
  nocc  <- capthists_few_alloccs_3x3$noccasions[j]

  # Summing capture histories over all of the simulated sampling occasions
  all.dat <- capthists_few_alloccs_3x3$capthist[[j]]
  all.mat <- matrix(0, nrow=nrow(all.dat[,1,]), ncol=ncol(all.dat[,1,]))
  for (i in 1:nocc) {
    all.mat <- all.mat + all.dat[,i,]
  }

  # Trap locations
  trap.loc <- attributes(all.dat)$traps

  # xlim, ylim (for our map area)
  xlim <- c(0.5, 50.5)
  ylim <- c(0.5, 50.5)

  # Creating the data object
  data <- list(encounter.data = all.mat, trap.loc = trap.loc, xlim = xlim, ylim = ylim, n.occasions = nocc)
  data
}

# Data object for 18 sampling occasions
data.18occ <- organise.data(1)

# Data object for 52 sampling occasins
data.52occ <- organise.data(2)

# Data object for 111 sampling occasions
data.111occ <- organise.data(3)

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need for the plots in Row 1
## ---------------------------------------------------------------------------------------

##### Running the MCMC #####

## Running MCMC for simulated data from 18, 52 and 111 sampling occasions.
results.18occ <- run.MCMC(data=data.18occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
save(results.18occ, file="MCMC_Results/Figure5/HomPP_18occ.RData")

results.52occ <- run.MCMC(data=data.52occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
save(results.52occ, file="MCMC_Results/Figure5/HomPP_52occ.RData")

results.111occ <- run.MCMC(data=data.111occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
save(results.111occ, file="MCMC_Results/Figure5/HomPP_111occ.RData")

## Loading in the RData files -- uncomment this if don't want to run the MCMC
#load("MCMC_Results/Figure5/HomPP_18occ.RData")
#load("MCMC_Results/Figure5/HomPP_52occ.RData")
#load("MCMC_Results/Figure5/HomPP_111occ.RData")

## Note that the burn-in iterations for these samples are discarded automatically, so we don't need to worry about this.
## Checking trace plots to make sure everything looks okay
# 18 sampling occasions -- looks good
check.trace.plots(results.18occ)
# 52 sampling occasions -- looks good
check.trace.plots(results.52occ)
# 111 sampling occasions -- looks good
check.trace.plots(results.111occ)

##### Creating the objects we need for the RACD plots #####

## Row 1 consists of RACD maps. So, we will create vectors that contain the posterior mean of the number of activity centres in each pixel -- these are the density values for each pixel in RACD maps (created using MCMC).
racd.18occ <- no.movement.density.vector(results=results.18occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
racd.52occ <-  no.movement.density.vector(results=results.52occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
racd.111occ <- no.movement.density.vector(results=results.111occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need for the plots in Row 2
## ---------------------------------------------------------------------------------------

##### Covariate value for each pixel #####
# Subsetting the covariate values from the data we loaded in above
mona.densities <-  small_blurry_mona_df[,c("x", "y", "Dblur")]
# Re-ordering 'mona.densities', so order of pixels matches order of pixels in 'pixel.centres' object
split <-  split(mona.densities, mona.densities$y)
mona.densities <-  do.call("rbind", split)
rownames(mona.densities) = NULL
# Now, subsetting covariate vector only so is in corresponding order to centres in 'pixel.centres'
dblur <-  mona.densities[,"Dblur"]
# Logging the covariate
dblur <-  log(dblur)

##### 'pixel.info' object needed for MCMC #####
pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
pixel.info <- cbind(pixel.centres, dblur)

##### Running the MCMC #####

# 18 sampling occasions, saving the results (takes under 10 minutes to run)
inhom.results.18occ <- run.MCMC.inhom(data=data.18occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
save(inhom.results.18occ, file="MCMC_Results/Figure5/InhomPP_18occ.RData")

# 52 sampling occasions, saving the results (takes under 10 minutes to run)
inhom.results.52occ <- run.MCMC.inhom(data=data.52occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
save(inhom.results.18occ, file="MCMC_Results/Figure5/InhomPP_52occ.RData")

# 111 sampling occasions, saving the results (takes under 10 minutes to run)
inhom.results.52occ <- run.MCMC.inhom(data=data.111occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
save(inhom.results.18occ, file="MCMC_Results/Figure5/InhomPP_111occ.RData")

## Loading in the RData files -- uncomment this if don't want to run the MCMC
#load("MCMC_Results/Figure5/InhomPP_18occ.RData")
#load("MCMC_Results/Figure5/InhomPP_52occ.RData")
#load("MCMC_Results/Figure5/InhomPP_111occ.RData")

## Note that with these MCMC samples, the burn-in isn't discarded automatically. So, each object contains data from 101,000 MCMC iterations
## Checking trace plots to decide how many iterations to discard as burn-in
# 18 sampling occasions
system.time(check.trace.plots(inhom.results.18occ))
# 52 sampling occasions

# 111 sampling occasions

## ---------------------------------------------------------------------------------------
# Putting together Figure 5
## ---------------------------------------------------------------------------------------
