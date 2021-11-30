## Code to create Figure 9 of the paper

## ---------------------------------------------------------------------------------------

# Creating the necessary data objects

## ---------------------------------------------------------------------------------------

## First, sourcing in capture histories to use
load("../output/capthists.RData")

# In Figure 9, we have columns with: 31 inividuals, 40 individuals and 44 individuals. We are working with 3/10/20 occasions, one set of simulated data for each number of occasions.

capthists_realised_and_expected_acd_few$noccasions[c(1,101,201)]
capthists_realised_and_expected_acd_few$secr.fitformula[c(1,101,201)]
# So, we want to use the 1st list element for the 1st column, 101th list element for 2nd column and 201st list element for 3rd column

# Column 1, 3 sampling occasions
first.col = capthists_realised_and_expected_acd_few$capthist[[1]]
# Summing capture histories over all of the 3 sampling occasions
encounterdat.3occ = matrix(0, nrow=nrow(first.col[,1,]), ncol=ncol(first.col[,1,]))
for (i in 1:3) {
  encounterdat.3occ = encounterdat.3occ + first.col[,i,]
}
# Trap locations
trap.loc = attributes(first.col)$traps
# xlim, ylim (we know these)
xlim = c(0.5, 50.5)
ylim = c(0.5, 50.5)
# Creating the data object for Figure 9, 3 sampling occasions
data.3occ = list(encounter.data = encounterdat.3occ, trap.loc = trap.loc, xlim = xlim, ylim = ylim, n.occasions = 3)
sum(encounterdat.3occ) # 85 detections

## 10 sampling occasions (second column)
second.col = capthists_realised_and_expected_acd_few$capthist[[101]]
# Summing the capture histories over all 10 sampling occasions
encounterdat.10occ = matrix(0, nrow=nrow(second.col[,1,]), ncol=ncol(second.col[,1,]))
for (i in 1:10) {
  encounterdat.10occ = encounterdat.10occ + second.col[,i,]
}
# Creating the data object (uses same trap locs, xlim, ylim as above)
data.10occ = list(encounter.data = encounterdat.10occ, trap.loc = trap.loc, xlim = xlim, ylim = ylim, n.occasions = 10)
sum(encounterdat.10occ) # 293 detections

## 20 sampling occasions (third column)
third.col = capthists_realised_and_expected_acd_few$capthist[[201]]
# Summing the capture histories over all 20 sampling occasions
encounterdat.20occ = matrix(0, nrow=nrow(third.col[,1,]), ncol=ncol(third.col[,1,]))
for (i in 1:20) {
  encounterdat.20occ = encounterdat.20occ + third.col[,i,]
}
sum(encounterdat.20occ)
# Creating the data object
data.20occ = list(encounter.data = encounterdat.20occ, trap.loc = trap.loc, xlim = xlim, ylim = ylim, n.occasions = 20)
sum(encounterdat.20occ) # 536 detections

## NEED to confirm with Ian that use 'capthists.RData' here? First few lines of plot_realised_usage_few.R don't seem to use this RData file!

## ---------------------------------------------------------------------------------------

# Running the MCMC

## ---------------------------------------------------------------------------------------

# Libraries we need
library("nimble")

# Function we need
source("MCMC_Function.R")

# Running MCMC for simulated data from 3, 10 and 20 sampling occasions
results.3occ = run.MCMC(data=data.3occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
save(results.3occ, file="Fig9_MCMC_3occ.RData")

results.10occ = run.MCMC(data=data.10occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
save(results.10occ, file="Fig9_MCMC_10occ.RData")

results.20occ = run.MCMC(data=data.20occ, M=300, parameters=c("lambda0", "coeff", "sigma", "N", "D", "z", "s"), n.iter=10000, n.burn=1000, lambda0.start=1, log_coeff.start=-5)
save(results.20occ, file="Fig9_MCMC_20occ.RData")

# Uncomment if want to load in the RData files instead:
# load("Figure 9/Fig9_MCMC_3occ.RData")
# load("Figure 9/Fig9_MCMC_10occ.RData")
# load("Figure 9/Fig9_MCMC_20occ.RData")

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

## Row 1: RACD maps. Are creating vectors that contain the density values for each pixel

source("DensityVectorFunction_RACDMaps.R")

density.3occ.all = no.movement.density.vector(results=results.3occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
density.10occ.all = no.movement.density.vector(results=results.10occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
density.20occ.all = no.movement.density.vector(results=results.20occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))

## Row 2: RUD maps. Once again, creating vectors that contain the density values for each pixel (we use more steps/create more objects in the process to do so than above)

source("RUDMaps_Functions.R")

# Coordinates for all pixel centres
pixel.centres = centres(xrange=c(0.5,50.5), yrange=c(0.5,50.5), x.pixels=50, y.pixels=50)

# Creating activity centre matrices (for explanation, see 'RUDMaps_Functions.R')
activity.centres.3occ = activity.matrices(results=results.3occ, M=300)
activity.centres.10occ = activity.matrices(results=results.10occ, M=300)
activity.centres.20occ = activity.matrices(results=results.20occ, M=300)

# Matrix containing z-values for each MCMC iteration (for each number of samp occ)
z.values.3occ = extract.z.values(results.3occ)
z.values.10occ =  extract.z.values(results.10occ)
z.values.20occ = extract.z.values(results.20occ)

# Creating the density vectors
density.3occ.all.withmov = density.vector(results=results.3occ, activity.centres=activity.centres.3occ, pixel.centres=pixel.centres, z.values=z.values.3occ)
save(density.3occ.all.withmov, file="Fig9.3occ.RData")
density.10occ.all.withmov = density.vector(results=results.10occ, activity.centres=activity.centres.10occ, pixel.centres=pixel.centres, z.values=z.values.10occ)
save(density.10occ.all.withmov, file="Fig9.10occ.RData")
density.20occ.all.withmov = density.vector(results=results.20occ, activity.centres=activity.centres.20occ, pixel.centres=pixel.centres, z.values=z.values.20occ)
save(density.20occ.all.withmov, file="Fig9.20occ.RData")


## ---------------------------------------------------------------------------------------

# Running Ian's code

## ---------------------------------------------------------------------------------------

# Libraries needed
library(tidyverse)
library(viridis)
library(patchwork)
library(scales)

# Loading RData object that we need
load("../output/capthist_summaries_100sim.RData")

# process the outputs
detectors_df_all <- fig67_results_100sim %>% purrr::map_depth(2, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

## Creating objects that contain the density values for each map, along with other required information
# 3 sampling occasions
df.3occ.all <- data.frame(pixel.centres, density.3occ.all, as.factor(rep(3,2500)),
                         as.factor(rep("None", 2500)))
names(df.3occ.all) <- c("x", "y", "value", "occasions", "movetype")
# 10 sampling occasions
df.10occ.all <- data.frame(pixel.centres, density.10occ.all, as.factor(rep(10,2500)),
                          as.factor(rep("None", 2500)))
names(df.10occ.all) <- c("x", "y", "value", "occasions", "movetype")
# 20 sampling occasions
df.20occ.all <- data.frame(pixel.centres, density.20occ.all, as.factor(rep(20,2500)),
                          as.factor(rep("None", 2500)))
names(df.20occ.all) <- c("x", "y", "value", "occasions", "movetype")
## Combining these objects into one data frame
ac_densities_without_movement <- rbind.data.frame(df.3occ.all, df.10occ.all, df.20occ.all)

# 3 sampling occasions
df.3occ.all.withmov <- data.frame(pixel.centres, density.3occ.all.withmov, as.factor(rep(3,2500)),
                                 as.factor(rep("With movement", 2500)))
names(df.3occ.all.withmov) <- c("x", "y", "value", "occasions", "movetype")
# 10 sampling occasions
df.10occ.all.withmov <- data.frame(pixel.centres, density.10occ.all.withmov, as.factor(rep(10,2500)),
                                  as.factor(rep("With movement", 2500)))
names(df.10occ.all.withmov) <- c("x", "y", "value", "occasions", "movetype")
# 20 sampling occasions
df.20occ.all.withmov <- data.frame(pixel.centres, density.20occ.all.withmov, as.factor(rep(20,2500)),
                                  as.factor(rep("With movement", 2500)))
names(df.20occ.all.withmov) <- c("x", "y", "value", "occasions", "movetype")
## Combining these objects into one data frame
ac_densities_with_movement <- rbind.data.frame(df.3occ.all.withmov, df.10occ.all.withmov, df.20occ.all.withmov)

## Combining both data frames now
ac_densities_with_movement = rbind.data.frame(ac_densities_without_movement, ac_densities_with_movement)

# detectors are the same for all plots so just extract unique combos of (x,y)
detectors <- detectors_df_all %>% group_by(x,y) %>% count()

# Column labels for plots
capthist_labels = paste(c(sum(encounterdat.3occ), sum(encounterdat.10occ), sum(encounterdat.20occ)), "detections\n", paste("(", c(nrow(encounterdat.3occ), nrow(encounterdat.10occ), nrow(encounterdat.20occ)), sep=""),  "individuals)")

# relabel factor levels for occasion variable
ac_densities_with_movement$occasions <- factor(ac_densities_with_movement$occasions,
                                            levels = c(3,10,20),
                                            labels = capthist_labels)

# scale the plots to have min 0 and max 1
# (need to think about best way to scale things for visualisation)
ac_densities_with_movement2 <- ac_densities_with_movement %>%
  group_by(occasions, movetype) %>%
  mutate(value = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

#pal <- wes_palette("Zissou1", 100, type = "continuous")
p1 <- ac_densities_with_movement %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis") +
  facet_grid(movetype ~ occasions) +
  geom_point(data = detectors, aes(x,y),
             colour = "black", pch = 4, size = 2) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p1

p2a <- ac_densities_with_movement %>%
  filter(movetype == "None") %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = pal, limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3)) +
  facet_grid(movetype ~ occasions) +
  geom_point(data = detectors, aes(x,y),
             colour = "black", pch = 4, size = 2) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right", legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2b <- ac_densities_with_movement %>%
  filter(movetype == "With movement") %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = pal, limits = c(0,0.1), breaks = c(0,0.05,0.1)) +
  #scale_fill_viridis(direction = 1, option = "viridis", trans = "log10", limits = c(1,103), oob = squish) +
  scale_fill_viridis(direction = 1, option = "viridis", trans = "log10") +
  facet_grid(movetype ~ occasions) +
  geom_point(data = detectors, aes(x,y),
             colour = "black", pch = 4, size = 2) +
  theme(strip.text.x = element_blank(),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right", legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


p2 <- p2a / p2b

p2

ggsave("Figure9.png", p2, width=8, height=6, dpi = 600)
