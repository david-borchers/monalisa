## Code to produce Figure 4
# Working directory should be the 'bayesian_code' folder

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

## Uncomment the lines below if want to run the MCMC 
## Running MCMC for simulated data from 18, 52 and 111 sampling occasions.
#results.18occ <- run.MCMC(data=data.18occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.18occ, file="MCMC_Results/Figure4/HomPP_18occ.RData")

#results.52occ <- run.MCMC(data=data.52occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.52occ, file="MCMC_Results/Figure4/HomPP_52occ.RData")

#results.111occ <- run.MCMC(data=data.111occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
#save(results.111occ, file="MCMC_Results/Figure4/HomPP_111occ.RData")

## Loading in the RData files containing the MCMC results 
load("MCMC_Results/Figure4/HomPP_18occ.RData")
load("MCMC_Results/Figure4/HomPP_52occ.RData")
load("MCMC_Results/Figure4/HomPP_111occ.RData")

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
# Logging the covariate, so we have the values of log(Dblur) (this is the covariate we will use to fit our SCR models)
log.dblur <-  log(dblur)

##### 'pixel.info' object needed for MCMC #####

## Uncomment if want to run M]CMC 
#pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50)
#pixel.info <- cbind(pixel.centres, log.dblur)

##### Running the MCMC #####

## Uncomment the lines below if want to run the MCMC. Note that depending on the computer, all of the MCMC chains below can take below 5 or 10 minutes to run. 
# 18 sampling occasions, saving the results 
#inhom.results.18occ <- run.MCMC.inhom(data=data.18occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.18occ, file="MCMC_Results/Figure4/InhomPP_18occ.RData")

# 52 sampling occasions, saving the results
#inhom.results.52occ <- run.MCMC.inhom(data=data.52occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.52occ, file="MCMC_Results/Figure4/InhomPP_52occ.RData")

# 111 sampling occasions, saving the results
#inhom.results.111occ <- run.MCMC.inhom(data=data.111occ, pixel.info=pixel.info, M=300, inits.vec=c(10, 4, 0), n.iter=100000, n.burn=1000)
#save(inhom.results.111occ, file="MCMC_Results/Figure4/InhomPP_111occ.RData")

## Loading in the RData files containing the MCMC results 
load("MCMC_Results/Figure4/InhomPP_18occ.RData")
load("MCMC_Results/Figure4/InhomPP_52occ.RData")
load("MCMC_Results/Figure4/InhomPP_111occ.RData")

## Note that with these MCMC samples, the burn-in isn't discarded automatically. So, each object contains data from 101,000 MCMC iterations
## Checking trace plots to decide how many iterations to discard as burn-in -- if we don't discard any iterations, we clearly see some burn-in on the trace plots. If we discard 1000 iterations as burn-in, the trace plots look good (see below)
# 18 sampling occasions'
check.trace.plots(inhom.results.18occ[-c(1:1000),], inhom=T)
# 52 sampling occasions
check.trace.plots(inhom.results.52occ[-c(1:1000),], inhom=T)
# 111 sampling occasions
check.trace.plots(inhom.results.111occ[-c(1:1000),], inhom=T)

## So, discarding 1000 iterations as burn-in for our 3 MCMC samples
inhom.results.18occ <- inhom.results.18occ[-c(1:1000),]
inhom.results.52occ <- inhom.results.52occ[-c(1:1000),]
inhom.results.111occ <- inhom.results.111occ[-c(1:1000),]

##### Creating the objects we need for the EACD plots #####

## Creating vectors containing density values for each pixel when working with 18/52/111 sampling occasions
eacd.18occ <- eacd.density.vector(results=inhom.results.18occ, covariate=dblur, nPix=2500)
eacd.52occ <- eacd.density.vector(results=inhom.results.52occ, covariate=dblur, nPix=2500)
eacd.111occ <- eacd.density.vector(results=inhom.results.111occ, covariate=dblur, nPix=2500)

## ---------------------------------------------------------------------------------------
# Putting together Figure 4 (based on the code in 'code/revision/mona-plots.R')
## ---------------------------------------------------------------------------------------

## Creating an object (labelled 'predicted_densities_all') that summarises all of the information that we will use to create the plots included in Figure 4
pixel.centres <- centres(xlim=c(0.5,50.5), ylim=c(0.5,50.5), x.pixels=50, y.pixels=50) # Pixel centres that we are working with
# Information we will use in Row 1
row1.1 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~1", 2500), occasions=rep(18, 2500), array_size=rep("3x3", 2500), value=racd.18occ) # Plot 1 of row 1
row1.2 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~1", 2500), occasions=rep(52, 2500), array_size=rep("3x3", 2500), value=racd.52occ) # Plot 2 of row 1
row1.3 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~1", 2500), occasions=rep(111, 2500), array_size=rep("3x3", 2500), value=racd.111occ) # Plot 3 of row 1
# Information we will use in row 2
row2.1 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~log(Dblur)", 2500), occasions=rep(18, 2500), array_size=rep("3x3", 2500), value=eacd.18occ) # Plot 1 of row 2
row2.2 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~log(Dblur)", 2500), occasions=rep(52, 2500), array_size=rep("3x3", 2500), value=eacd.52occ) # Plot 2 of row 2
row2.3 <- data.frame(x=pixel.centres[,1], y=pixel.centres[,2], covtype=rep("D~log(Dblur)", 2500), occasions=rep(111, 2500), array_size=rep("3x3", 2500), value=eacd.111occ) # Plot 3 of row 2
# Combining all of the information into one object 
predicted_densities_all <- rbind(row1.1, row1.2, row1.3, row2.1, row2.2, row2.3)

## Object that contains information on detectors 
detectors_df_all <- res_acd %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()
head(detectors_df_all)

## Headings for the columns of Figure 4
nn <- 3
occ <- capthists_few_alloccs_3x3$noccasions
asz <- c("3x3")
chs <- data.frame(do.call(rbind, lapply(capthists_few_alloccs_3x3$capthist, summary, terse = TRUE)))
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = chs$Detections, .y = chs$Animals, .f = paster) %>% unlist() 

## Adding these headings to the 'predicted_densities_all' and 'detectors_df_all' objects
predicted_densities_all$occasions2 <- factor(predicted_densities_all$occasions, 
                                             levels = occ,
                                             labels = capthist_labels)
detectors_df_all$occasions2 <- factor(detectors_df_all$occasions, 
                                      levels = occ,
                                      labels = capthist_labels)

## Adding row headings to the 'predicted_densities_all' and 'detectors_df_all' objects
predicted_densities_all$covtype2 <- factor(predicted_densities_all$covtype, 
                                           levels = unique(predicted_densities_all$covtype),
                                           labels = c("Realised AC", "Expected AC"))
detectors_df_all$covtype2 <- factor(detectors_df_all$covtype, 
                                    levels = unique(detectors_df_all$covtype),
                                    labels = c("Realised AC", "Expected AC"))

## Setting the max value of the colour scale (making it the same for all plots)
xx <- predicted_densities_all %>% filter(array_size == "3x3", occasions %in% capthists_few_alloccs_3x3$noccasions[1:nn])
maxval <- max(xx$value)

## Creating Figure 4
p2a <- predicted_densities_all %>%
  filter(occasions %in% occ[1:nn], array_size %in% asz) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Light grays"), limits = c(0,maxval)) +
  facet_grid(covtype2 ~ occasions2) +
  geom_point(data = detectors_df_all %>% filter(occasions %in% occ[1:nn], array_size %in% asz), inherit.aes = T,
             colour = "gray80", pch = 4, size = 2) +
  geom_point(data = simulated_points, inherit.aes = F, aes(x=x,y=y),
             colour = "darkorange", pch = 16, size = 1, alpha = 0.5) +
  coord_equal() +
  theme_classic(base_size = 14) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        panel.spacing=unit(-1, "lines"),
        strip.background = element_rect(fill=NA, colour = NA), 
        legend.position="none", legend.key.width = unit(3, "cm"),
        legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

## Saving Figure 4
ggsave("mona_3x3.png", p2a, width=8, height=6, dpi=600)
