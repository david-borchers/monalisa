## Code to produce Figure 5
# Working directory should be the 'bayesian_code' folder

## In this figure, the first column uses simulated data w/ 7 sampling occasions. The second column uses simulated data w/ 25 sampling occasions, and the third column uses simulated data w/ 55 sampling occasions. 

## Libraries we need
library(nimble)
library(coda)
library(nimbleSCR)
library(spatstat)
library(ggplot2)
library(dplyr)
library(stringr)
library(purrr)
library(secr)
library(patchwork)
library(ggpubr)
## Objects we need
load("../output/revision/mona-inputs.RData")
load("../output/revision/mona-results.RData")
## Functions we need
source("Functions.R")

## ---------------------------------------------------------------------------------------
# Creating the data objects we need for all MCMC samples
## ---------------------------------------------------------------------------------------

## Function we will use
# The only argument we provide is 'j': this will be a value in {1, 2, 3} and represents the index of the objects we want to work with from the RData objects we have loaded in. If want objects generated using 7, 25 or 55 sampling occasions, j=1,2,3 respectively.
organise.data = function(j) {
  # Number of sampling occasions used for simulated data
  nocc  <- capthists_few_alloccs_7x7$noccasions[j]

  # Summing capture histories over all of the simulated sampling occasions
  all.dat <- capthists_few_alloccs_7x7$capthist[[j]]
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

# Data object for 7 sampling occasions
data.7occ <- organise.data(1)

# Data object for 25 sampling occasions
data.25occ <- organise.data(2)

# Data object for 55 sampling occasions
data.55occ <- organise.data(3)

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need for the plots in Row 1
## ---------------------------------------------------------------------------------------

## Row 1 consists of RACD maps. The SCR models that we fit to create these maps assume that the state process (the random process governing the distribution of the activity centres) is a homogeneous Poisson process.

##### Running the MCMC #####

## Uncomment the lines below if want to run the MCMC 
## Running MCMC for simulated data from 7, 25 and 55 sampling occasions.
#results.7occ <- run.MCMC(data=data.7occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.7occ, file="MCMC_Results/Figure5/HomPP_7occ.RData")

#results.25occ <- run.MCMC(data=data.25occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.25occ, file="MCMC_Results/Figure5/HomPP_25occ.RData")

#results.55occ <- run.MCMC(data=data.55occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.55occ, file="MCMC_Results/Figure5/HomPP_55occ.RData")

## Loading in the RData files containing the MCMC results 
load("MCMC_Results/Figure5/HomPP_7occ.RData")
load("MCMC_Results/Figure5/HomPP_25occ.RData")
load("MCMC_Results/Figure5/HomPP_55occ.RData")

## Note that the burn-in iterations for these samples are discarded automatically, so we don't need to worry about this.
## Checking trace plots to make sure everything looks okay
# 7 sampling occasions -- looks good
check.trace.plots(results.7occ)
# 25 sampling occasions
check.trace.plots(results.25occ)
# 55 sampling occasions
check.trace.plots(results.55occ)
# For 25 and 55 sampling occasions, we have very restricted mixing of N (N only ranges over a small range of values). As D is calculated directly from N, this means we also have restricted mixing of D. Note that for 25 and 55 sampling occasions, the lower limit of N is equal to the total number of observed animals, which seems sensible. 
# Changing the value of M or the prior for sigma doesn't change this. Lambda0 and sigma seem to be mixing well, which seems to further reinforce that changing the priors for lambda0 and sigma would not help.
# It seems likely that in these cases, having many traps (49 traps) and many sampling occasions (25, 55 sampling occasions) means that we have enough information that we become fairly sure of the true value of N in the region of interest, resulting in the limited mixing we see. This seems to be supported by the fact that if we run: 'table(results.55occ[,"N"]); table(results.25occ[,"N"])', we can see that with 55 sampling occasions, the mixing for N is considerably more limited than for 25 occasions (we are even more sure about the value of N, due to the increase in sampling occasions). So, we will continue.


##### Creating the objects we need for the RACD plots #####

## Row 1 consists of RACD maps. So, we will create vectors that contain the posterior mean of the number of activity centres in each pixel -- these are the density values for each pixel in RACD maps (created using MCMC).
racd.7occ <- no.movement.density.vector(results=results.7occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
racd.25occ <-  no.movement.density.vector(results=results.25occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
racd.55occ <- no.movement.density.vector(results=results.55occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))

## ---------------------------------------------------------------------------------------
# Creating the objects we specifically need for the plots in Row 2
## ---------------------------------------------------------------------------------------

## Row 2 consists of EACD maps. The SCR models that we fit to create these maps assume that the state process (the random process governing the distribution of the activity centres) is an inhomogeneous Poisson process.

##### Covariate value for each pixel #####

# Note that we use the same covariate as in 'bayesian_code/Figure4_code.R'
# Subsetting the covariate values from the data we loaded in above
mona.densities <-  small_blurry_mona_df[,c("x", "y", "Dblur")]
# Re-ordering 'mona.densities', so order of pixels matches order of pixels in 'pixel.centres' object
split <-  split(mona.densities, mona.densities$y)
mona.densities <-  do.call("rbind", split)
rownames(mona.densities) = NULL
# Now, subsetting covariate vector only so is in corresponding order to centres in 'pixel.centres'
dblur <-  mona.densities[,"Dblur"]
# Logging the covariate, so we have the values of log(Dblur) (this is the covariate we will use to fit our SCR models)
log.dblur <-  log(dblur)

##### 'pixel.info' object needed for MCMC #####

## Uncomment if want to run MCMC 
#pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
#pixel.info <- cbind(pixel.centres, log.dblur)

##### Running the MCMC #####

## Uncomment the lines below if want to run the MCMC. Note that all of the MCMC chains below take below 10 minutes to run. 
# 7 sampling occasions, saving the results 
#inhom.results.7occ <- run.MCMC.inhom(data=data.7occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.7occ, file="MCMC_Results/Figure5/InhomPP_7occ.RData")

# 25 sampling occasions, saving the results
#inhom.results.25occ <- run.MCMC.inhom(data=data.25occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.25occ, file="MCMC_Results/Figure5/InhomPP_25occ.RData")

# 55 sampling occasions, saving the results
#inhom.results.55occ <- run.MCMC.inhom(data=data.55occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.55occ, file="MCMC_Results/Figure5/InhomPP_55occ.RData")

## Loading in the RData files containing the MCMC results 
load("MCMC_Results/Figure5/InhomPP_7occ.RData")
load("MCMC_Results/Figure5/InhomPP_25occ.RData")
load("MCMC_Results/Figure5/InhomPP_55occ.RData")

## Note that with these MCMC samples, the burn-in isn't discarded automatically. So, each object contains data from 101,000 MCMC iterations
## Checking trace plots to decide how many iterations to discard as burn-in -- if we don't discard any iterations, we clearly see some burn-in on the trace plots. If we discard 1000 iterations as burn-in, the trace plots look good (see below)
# 18 sampling occasions
check.trace.plots(inhom.results.7occ[-c(1:1000),], inhom=T)
# 52 sampling occasions
check.trace.plots(inhom.results.25occ[-c(1:1000),], inhom=T)
# 111 sampling occasions
check.trace.plots(inhom.results.55occ[-c(1:1000),], inhom=T)
# Overall, things look good. We once again see restricted mixing in the trace plots for N and D as we did above. Once again, we believe this is due to an increase in information meaning that we have increased confidence in the value of N (and therefore D), so we will continue.

## So, discarding 1000 iterations as burn-in for our 3 MCMC samples
inhom.results.7occ <- inhom.results.7occ[-c(1:1000),]
inhom.results.25occ <- inhom.results.25occ[-c(1:1000),]
inhom.results.55occ <- inhom.results.55occ[-c(1:1000),]

##### Creating the objects we need for the EACD plots #####

## Creating vectors containing density values for each pixel when working with 7/25/55 sampling occasions
eacd.7occ <- eacd.density.vector(results=inhom.results.7occ, covariate=log.dblur, nPix=2500)
eacd.25occ <- eacd.density.vector(results=inhom.results.25occ, covariate=log.dblur, nPix=2500)
eacd.55occ <- eacd.density.vector(results=inhom.results.55occ, covariate=log.dblur, nPix=2500)

## ---------------------------------------------------------------------------------------
# Putting together Figure 5 (based on the code in 'code/revision/mona-plots.R')
## ---------------------------------------------------------------------------------------

## Creating an object (labelled 'predicted_densities_all') that summarises all of the information that we will use to create the plots included in Figure 5
pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50) # Pixel centres that we are working with
# Information we will use in Row 1
row1.1 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~1", 2500), occasions=rep(7, 2500), array_size=rep("7x7", 2500), value=racd.7occ) # Plot 1 of row 1
row1.2 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~1", 2500), occasions=rep(25, 2500), array_size=rep("7x7", 2500), value=racd.25occ) # Plot 2 of row 1
row1.3 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~1", 2500), occasions=rep(55, 2500), array_size=rep("7x7", 2500), value=racd.55occ) # Plot 3 of row 1
# Information we will use in row 2
row2.1 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~log(Dblur)", 2500), occasions=rep(7, 2500), array_size=rep("7x7", 2500), value=eacd.7occ) # Plot 1 of row 2
row2.2 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~log(Dblur)", 2500), occasions=rep(25, 2500), array_size=rep("7x7", 2500), value=eacd.25occ) # Plot 2 of row 2
row2.3 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~log(Dblur)", 2500), occasions=rep(55, 2500), array_size=rep("7x7", 2500), value=eacd.55occ) # Plot 3 of row 2
# Combining all of the information into one object 
predicted_densities_all <- rbind(row1.1, row1.2, row1.3, row2.1, row2.2, row2.3)

## Object that contains information on detectors 
detectors_df_all <- res_acd %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

nn <- 3
occ <-capthists_few_alloccs_7x7$noccasions
asz <- c("7x7")

chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_7x7$capthist, summary, terse = TRUE)))
chs <- chs %>% dplyr::filter(Occasions %in% occ)
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist()

predicted_densities_all$occasions2 <- factor(predicted_densities_all$occasions, 
                                             levels = occ,
                                             labels = capthist_labels)

detectors_df_all$occasions2 <- factor(detectors_df_all$occasions, 
                                      levels = occ,
                                      labels = capthist_labels)

predicted_densities_all$covtype2 <- factor(predicted_densities_all$covtype, 
                                           levels = unique(predicted_densities_all$covtype),
                                           labels = c("Realised AC", "Expected AC"))

detectors_df_all$covtype2 <- factor(detectors_df_all$covtype, 
                                    levels = unique(detectors_df_all$covtype),
                                    labels = c("Realised AC", "Expected AC"))

## Setting the max value of the colour scale (making it the same for all plots)
xx <- predicted_densities_all %>% filter(array_size == "7x7", occasions %in% capthists_few_alloccs_7x7$noccasions[1:nn])
maxval <- max(xx$value)

p2b <- predicted_densities_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Light grays"), limits = c(0,maxval)) +
  facet_grid(covtype2 ~ occasions2) +
  geom_point(data = detectors_df_all %>% filter(occasions %in% occ[1:nn], array_size %in% asz), inherit.aes = T,
             colour = "gray80", pch = 4, size = 2) +
  geom_point(data = simulated_points, inherit.aes = F, aes(x=x,y=y),
             colour = "darkorange", pch = 16, size = 1, alpha = 0.5) +
  theme_bw(base_size = 14) +
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.spacing=unit(-1, "lines"),
        strip.background = element_rect(fill=NA, colour = NA), 
        legend.position="none", legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2b

