## Code to create Figure 9 of the paper

## ---------------------------------------------------------------------------------------

# Creating the necessary data objects -- need to run whether or not need to run MCMC

## ---------------------------------------------------------------------------------------

## First, sourcing in capture histories to use
load("../output/capthists.RData")

# In Figure 9, we have columns with: 31 inividuals, 40 individuals and 44 individuals. We are working with 3/10/20 occasions, one set of simulated data for each number of occasions.

## 3 sampling occasions (first column)
first.col <- capthists_realised_and_expected_acd_few$capthist[[1]]
# Summing capture histories over all of the 3 sampling occasions
encounterdat.3occ <- matrix(0, nrow=nrow(first.col[,1,]), ncol=ncol(first.col[,1,]))
for (i in 1:3) {
  encounterdat.3occ <- encounterdat.3occ + first.col[,i,]
}
# Trap locations
trap.loc <- attributes(first.col)$traps
# xlim, ylim (we know these)
xlim <- c(0.5, 50.5)
ylim <- c(0.5, 50.5)
# Creating the data object for Figure 9, 3 sampling occasions
data.3occ <- list(encounter.data = encounterdat.3occ, trap.loc = trap.loc, xlim = xlim, ylim = ylim, n.occasions = 3)
#sum(encounterdat.3occ) # 85 detections

## 10 sampling occasions (second column)
second.col <- capthists_realised_and_expected_acd_few$capthist[[101]]
# Summing the capture histories over all 10 sampling occasions
encounterdat.10occ <- matrix(0, nrow=nrow(second.col[,1,]), ncol=ncol(second.col[,1,]))
for (i in 1:10) {
  encounterdat.10occ <- encounterdat.10occ + second.col[,i,]
}
# Creating the data object (uses same trap locs, xlim, ylim as above)
data.10occ <- list(encounter.data = encounterdat.10occ, trap.loc = trap.loc, xlim = xlim, ylim = ylim, n.occasions = 10)
#sum(encounterdat.10occ) # 293 detections

## 20 sampling occasions (third column)
third.col <- capthists_realised_and_expected_acd_few$capthist[[201]]
# Summing the capture histories over all 20 sampling occasions
encounterdat.20occ <- matrix(0, nrow=nrow(third.col[,1,]), ncol=ncol(third.col[,1,]))
for (i in 1:20) {
  encounterdat.20occ <- encounterdat.20occ + third.col[,i,]
}
# Creating the data object
data.20occ <- list(encounter.data = encounterdat.20occ, trap.loc = trap.loc, xlim = xlim, ylim = ylim, n.occasions = 20)
#sum(encounterdat.20occ) # 536 detections

## ---------------------------------------------------------------------------------------

# Running the MCMC

## ---------------------------------------------------------------------------------------

# Libraries we need
library("nimble")

# Function we need
source("MCMC_Function.R")

## Running MCMC for simulated data from 3, 10 and 20 sampling occasions.
results.3occ <- run.MCMC(data=data.3occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
save(results.3occ, file="Figure 9/Fig9_MCMC_3occ.RData")

results.10occ <- run.MCMC(data=data.10occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
save(results.10occ, file="Figure 9/Fig9_MCMC_10occ.RData")

results.20occ <- run.MCMC(data=data.20occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
save(results.20occ, file="Figure 9/Fig9_MCMC_20occ.RData")

# Loading in the RData files -- uncomment this if don't want to run the MCMC
#load("Figure 9/Fig9_MCMC_3occ.RData")
#load("Figure 9/Fig9_MCMC_10occ.RData")
#load("Figure 9/Fig9_MCMC_20occ.RData")

# Burn-in
results.3occ <- results.3occ[-c(1:500),]
results.10occ <- results.10occ[-c(1:500),]
results.20occ <- results.20occ[-c(1:500),]

# Checking the trace plots
# 3 sampling occasions -- all looking v good
plot(results.3occ[,"lambda0"], type='l')
plot(results.3occ[,"sigma"], type='l')
plot(results.3occ[,"N"], type='l')
plot(results.3occ[,"D"], type='l')
# 10 sampling occasions -- v good
plot(results.10occ[,"lambda0"], type='l')
plot(results.10occ[,"sigma"], type='l')
plot(results.10occ[,"N"], type='l')
plot(results.10occ[,"D"], type='l')
# 20 sampling occaions -- looking v good
plot(results.20occ[,"lambda0"], type='l')
plot(results.20occ[,"sigma"], type='l')
plot(results.20occ[,"N"], type='l')
plot(results.20occ[,"D"], type='l')

## ---------------------------------------------------------------------------------------

# Creating the objects we need

## ---------------------------------------------------------------------------------------

## Row 1: RACD maps. Are creating vectors that contain the density values for each pixel. Need to run this section to create first row of plots:
source("DensityVectorFunction_RACDMaps.R")

density.3occ.all <- no.movement.density.vector(results=results.3occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
density.10occ.all <-  no.movement.density.vector(results=results.10occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
density.20occ.all <- no.movement.density.vector(results=results.20occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))

# Coordinates of pixel centres we are working with
source("RUDMaps_Functions.R")

pixel.centres <- centres(xrange=c(0.5,50.5), yrange=c(0.5,50.5), x.pixels=50, y.pixels=50)

## ---------------------------------------------------------------------------------------

# Running code to create maps

## ---------------------------------------------------------------------------------------

# Libraries needed
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(purrr)

# Loading RData object that we need
load("../output/mona_raw_outputs.RData")

# Process the outputs
detectors_df_all <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

# ---------------
## Creating objects that contain the density values for each map, along with other required information
# 3 sampling occasions
df.3occ.all <- data.frame(pixel.centres, density.3occ.all, as.factor(rep("D~1", 2500)), as.factor(rep("None", 2500)), as.factor(rep(3,2500)))
names(df.3occ.all) <- c("x", "y", "value", "covtype", "movetype", "occasions")
# 10 sampling occasions
df.10occ.all <- data.frame(pixel.centres, density.10occ.all, as.factor(rep("D~1", 2500)),  as.factor(rep("None", 2500)), as.factor(rep(10,2500)))
names(df.10occ.all) <- c("x", "y", "value", "covtype", "movetype", "occasions")
# 20 sampling occasions
df.20occ.all <- data.frame(pixel.centres, density.20occ.all, as.factor(rep("D~1", 2500)), as.factor(rep("None", 2500)), as.factor(rep(20,2500)))
names(df.20occ.all) <- c("x", "y", "value", "covtype", "movetype", "occasions")
## Combining these objects into one data frame
ac_densities <- rbind.data.frame(df.3occ.all, df.10occ.all, df.20occ.all)
# ---------------

# Detectors are the same for all plots so just extract unique combos of (x,y)
detectors <- detectors_df_all %>% group_by(x,y) %>% count()

# Column labels for plots
capthist_labels <-  paste(c(sum(encounterdat.3occ), sum(encounterdat.10occ), sum(encounterdat.20occ)), "detections\n", paste("(", c(nrow(encounterdat.3occ), nrow(encounterdat.10occ), nrow(encounterdat.20occ)), sep=""),  "individuals)")
capthist_labels <-  c(capthist_labels) # Adding this 'extra' level for Ian's code below to work

# Relabel factor levels for occasion variable
ac_densities$occasions <- factor(ac_densities$occasions,
                                            levels = c(3,10,20),
                                            labels = capthist_labels)

# Scale the plots to have min 0 and max 1
# (need to think about best way to scale things for visualisation)
#ac_densities_with_movement2 <- ac_densities_with_movement %>%
#  group_by(occasions, movetype) %>%
#  mutate(value = (value - min(value)) / (max(value) - min(value))) %>%
#  ungroup()

#pal <- wes_palette("Zissou1", 100, type = "continuous")

p2a <- ac_densities %>%
  filter(movetype == "None") %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = pal, limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3)) +
  facet_grid(movetype ~ occasions) +
  geom_point(data = detectors, aes(x,y),
             colour = "black", pch = 4, size = 2) +
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right", legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank(),
        strip.text.y=element_blank())

p2a

ggsave("Figure9.png", p2a, width=8, height=5, dpi = 600, bg="white")
