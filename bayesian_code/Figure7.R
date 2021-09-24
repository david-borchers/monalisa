## Code to create Figure 7 of the paper

## ---------------------------------------------------------------------------------------

# Creating the necessary data objects

## ---------------------------------------------------------------------------------------

## First, sourcing in 'posthoc_extract_chs.R'
load("../output/mona_raw_outputs_100sim.RData")
source("../code/posthoc_extract_chs.R")
#load("mona_raw_outputs_100sim.RData") # for NeSI
#source("posthoc_extract_chs.R") # for NeSI

## So, it appears that:
# ch7b is the plot from traps placed at the top right
# ch7c is the plot from traps placed at the bottom left
# ch7e is the plot from traps placed at the top right
# ch7f is the plot from traps placed at the bottom left

## Are going to create 4 EACD maps, two sets.
# One set uses the 'Dgood' covariate, and uses the trap array at the top right and the bottom left.
# The other uses the 'Dblur' covariate, and also uses the trap array at the top right and the bottom left.

## So, we need to use MCMC to fit an inhomogeneous density model -- i.e. to generate MCMC samples for beta0 and beta1, which we will then use to estimate the inhomogeneous PP, and therefore to create the maps.

## We begin by creating the necessary data objects. Here, we are working with one set of simulated data only, which has 20 sampling occasions. So, for each map, we need to create data objects that contain an encounter data matrix found by summing over these 20 sampling occasions

## ch7b
# Encounter matrix
encounterdat.ch7b = matrix(0, nrow=nrow(ch7b[,1,]), ncol=ncol(ch7b[,1,]))
for (i in 1:20) {
  encounterdat.ch7b = encounterdat.ch7b + ch7b[,i,]
}
## Trap locations
ch7b.traploc = attributes(ch7b)$traps
# xlim, ylim
xlim = c(0.5, 50.5)
ylim = c(0.5, 50.5)
# Data object
data.ch7b = list(encounter.data = encounterdat.ch7b, trap.loc = ch7b.traploc, xlim = xlim, ylim = ylim, n.occasions = 20)

## ch7c
encounterdat.ch7c = matrix(0, nrow=nrow(ch7c[,1,]), ncol=ncol(ch7c[,1,]))
for (i in 1:20) {
  encounterdat.ch7c = encounterdat.ch7c + ch7c[,i,]
}
ch7c.traploc = attributes(ch7c)$traps
data.ch7c = list(encounter.data = encounterdat.ch7c, trap.loc = ch7c.traploc, xlim = xlim, ylim = ylim, n.occasions = 20)

## ch7e
encounterdat.ch7e = matrix(0, nrow=nrow(ch7e[,1,]), ncol=ncol(ch7e[,1,]))
for (i in 1:20) {
  encounterdat.ch7e = encounterdat.ch7e + ch7e[,i,]
}
ch7e.traploc = attributes(ch7e)$traps
data.ch7e = list(encounter.data = encounterdat.ch7e, trap.loc = ch7e.traploc, xlim = xlim, ylim = ylim, n.occasions = 20)

## ch7f
encounterdat.ch7f = matrix(0, nrow=nrow(ch7f[,1,]), ncol=ncol(ch7f[,1,]))
for (i in 1:20) {
  encounterdat.ch7f = encounterdat.ch7f + ch7f[,i,]
}
ch7f.traploc = attributes(ch7f)$traps
data.ch7f = list(encounter.data = encounterdat.ch7f, trap.loc = ch7f.traploc, xlim = xlim, ylim = ylim, n.occasions = 20)

## ---------------------------------------------------------------------------------------

# Running the MCMC

## ---------------------------------------------------------------------------------------

# Libraries we need
library(nimble)
library(coda)
library(nimbleSCR)

# Sourcing in the dpoisLocal_normal_2() function that we will need to run the MCMC
#source("dpoisLocal_normal.R")

# Sourcing in the function that we'll need to run the MCMC
source("MCMC_Function_Inhomogeneous.R")
library("spatstat")

## ------------------------------

# ch7b

## ------------------------------

# We use the R files in the 'Figure 7' folder to run the MCMC using NeSI, not the lines below!
#ch7b.sample = run.MCMC.inhom(data=data.ch7b, M=9000, mona.column="Dgood", n.iter=2000)
#save(ch7b.sample, file="ch7b.RData")
#rm(ch7b.sample)
#ch7c.sample = run.MCMC.inhom(data.ch7c, M=9000, mona.column="Dgood", n.iter=11000)
#save(ch7c.sample, file="ch7c.RData")
#rm(ch7c.sample)
#ch7e.sample = run.MCMC.inhom(data.ch7e, M=9000, mona.column="Dblur", n.iter=11000)
#save(ch7e.sample, file="ch7e.RData")
#rm(ch7e.sample)
#ch7f.sample = run.MCMC.inhom(data.ch7f, M=9000, mona.column="Dblur", n.iter=1000)
#save(ch7f.sample, file="ch7f.RData")
#rm(ch7f.sample)

## Loading in what we get from NeSI from the separate R files in the 'Figure 7' folder
load("Figure 7/ch7b.RData")
#ch7b.sample = sample
load("Figure 7/ch7c.RData")
#ch7c.sample = sample
load("Figure 7/ch7e.RData")
#ch7e.sample = sample
load("Figure 7/ch7f.RData")
#ch7f.sample = sample

# All have about 120,000 iterations. Can shave down later to be the same amount
dim(ch7b.sample)
dim(ch7c.sample)
dim(ch7e.sample)
dim(ch7f.sample)

# Removing 1000 iterations as burn-in
ch7b.sample = ch7b.sample[-c(1:1000),]
ch7c.sample = ch7c.sample[-c(1:1000),]
ch7e.sample = ch7e.sample[-c(1:1000),]
ch7f.sample = ch7f.sample[-c(1:1000),]

# Looking at trace plots
# ch7b
plot(ch7b.sample[,"lambda0"], type="l")
plot(ch7b.sample[,"sigma"], type="l")
plot(ch7b.sample[,"beta0"], type="l")
plot(ch7b.sample[,"beta1"], type="l")
plot(ch7b.sample[,"N"], type="l")
plot(ch7b.sample[,"D"], type="l")
# Looking very good!

# ch7c
plot(ch7c.sample[,"lambda0"], type="l")
plot(ch7c.sample[,"sigma"], type="l")
plot(ch7c.sample[,"beta0"], type="l")
plot(ch7c.sample[,"beta1"], type="l")
plot(ch7c.sample[,"N"], type="l")
plot(ch7c.sample[,"D"], type="l")
# Why no mixing of lambda0 and sigma? Starting value for sigma seems bad?

# ch7e
plot(ch7e.sample[,"lambda0"], type="l")
plot(ch7e.sample[,"sigma"], type="l")
plot(ch7e.sample[,"beta0"], type="l")
plot(ch7e.sample[,"beta1"], type="l")
plot(ch7e.sample[,"N"], type="l")
plot(ch7e.sample[,"D"], type="l")
# Why no mixing of lambda0 and sigma? Starting value for sigma seems bad?

# ch7f
plot(ch7f.sample[,"lambda0"], type="l")
plot(ch7f.sample[,"sigma"], type="l")
plot(ch7f.sample[,"beta0"], type="l")
plot(ch7f.sample[,"beta1"], type="l")
plot(ch7f.sample[,"N"], type="l")
plot(ch7f.sample[,"D"], type="l")
#  Why no mixing of lambda0 and sigma? Starting value for sigma seems bad?

## ---------------------------------------------------------------------------------------

# Necessary objects

## ---------------------------------------------------------------------------------------

## We want to create a density vector for each map, using the covariate values and posterior estimates for beta0 and beta1

## Covariate values
# Dgood
load("../output/mona_inputs.RData")
# Extracting x,y-coordinates and corresponding Dgood values (we use the x,y-coordinates to make sure the vector is ordered correctly, see below)
dgood.df = mona_df[,c("x", "y", "Dgood")]
# Removing duplicate rows
sequence = seq(1, 7500, by=3)
dgood.df = dgood.df[sequence,]
# Reordering the data frame in the same way as we do in the run.MCMC.inhom() function, so the Dgood values are in the desired order (so they correspond with the order of the pixel centres that we will generate using the centres() function)
split = split(dgood.df, dgood.df$y)
dgood.df = do.call("rbind", split)
# Now, extracting the 'Dgood' values
dgood = dgood.df$Dgood

# Dblur
# Doing the same thing as above, but to extract the Dblur values
dblur.df = mona_df[,c("x", "y", "Dblur")]
dblur.df = dblur.df[sequence,]
split = split(dblur.df, dblur.df$y)
dblur.df = do.call("rbind", split)
dblur = dblur.df$Dblur

## Posterior estimates for beta0 and beta1. We will use the posterior mean.
# ch7b
ch7b.beta0 = mean(ch7b.sample[,"beta0"])
ch7b.beta1 = mean(ch7b.sample[,"beta1"])
# Previously, the means were higher than they are -- 2.28 and 1.14 vs 1.38 and 0.20 now

# ch7c
ch7c.beta0 = mean(ch7c.sample[,"beta0"])
ch7c.beta1 = mean(ch7c.sample[,"beta1"])
# Kind of used to be higher before -- 1.25 and -0.017 vs 1.14 and 0.14 now

# ch7e
ch7e.beta0 = mean(ch7e.sample[,"beta0"])
ch7e.beta1 = mean(ch7e.sample[,"beta1"])
# Mean used to be higher -- 3.35 and 2.24 vs 1.50 and 0.32 now

# ch7f
ch7f.beta0 = mean(ch7f.sample[,"beta0"])
ch7f.beta1 = mean(ch7f.sample[,"beta1"])
# Mean used to be higher -- 3.22 and 1.97 vs 0.85 and -0.13 now!

## Checking if the other posterior means make sense!
# True N is 7451
mean(ch7b.sample[,"N"]) # Now is 7722, was 7364
mean(ch7c.sample[,"N"]) # Now is 6580, was 8992 (but with wrong covariate, I think)
mean(ch7e.sample[,"N"]) # Now is 7801, was 7552
mean(ch7f.sample[,"N"]) # Now is 6891, was 8431 (I think also with wrong covariate)
# Doesn't seem so bad?

# True sigma is 2
mean(ch7b.sample[,'sigma']) # Now is 1.99, was 2.00
mean(ch7c.sample[,'sigma']) # Now is 2.02, was 0.008 (wrong covariate)
mean(ch7e.sample[,'sigma']) # Now is 1.99, was 1.99
mean(ch7f.sample[,'sigma']) # Now is 2.02, was 2.02
# Looks pretty good!

# True lambda0 is 0.69*20 = 13.8
mean(ch7b.sample[,'lambda0']) # Now is 13.2, was 13.5
mean(ch7c.sample[,'lambda0']) # Now is 13.3, was 0.5 (wrong covariate)
mean(ch7e.sample[,'lambda0']) # Now is 13.2, was 13.5
mean(ch7f.sample[,'lambda0']) # Now is 13.3, was 13.4
# Also looking really good!

# Doesn't really seem like estimation isn't going wrong here, even compared to before? So v different means for beta0 and beta1 isn't too worrying? Is okay that plots are 'dim' then?

## Therefore, creating the density vector for each plot
# ch7b
ch7b.density = exp(ch7b.beta0 + (ch7b.beta1 * log(dgood)))

# ch7c
ch7c.density = exp(ch7c.beta0 + (ch7c.beta1 * log(dgood)))

# ch7e
ch7e.density = exp(ch7e.beta0 + (ch7e.beta1 * log(dblur)))

# ch7f
ch7f.density = exp(ch7f.beta0 + (ch7f.beta1 * log(dblur)))

## ---------------------------------------------------------------------------------------

# Using Ian's code to create plots

## ---------------------------------------------------------------------------------------

library(tidyverse)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(purrr)

mona_df <- mona_df %>% dplyr::select(-cc) %>% distinct()

# process the covariates - contains duplicated pixel centres, each row of duplicate gives Dgood/Dblur value
predicted_densities_covs <- mona_df %>% select(-D) %>%
  gather(grid, value, -x, -y) %>% arrange(x,y) %>%
  mutate(covtype = grid,
         array_origin = "none") %>%
  select(-grid) %>%
  filter(covtype %in% c("Dgood", "Dblur"))

# process the outputs
# Appears to contain density values for 3 different covariates, 2 different trap arrays and 1 or 20 sampling occasions -- so contains 12 sets of 2500 values = 30000 rows
predicted_densities_all <- fig5_results %>% purrr::map("predicted_densities") %>% map_df(bind_rows)
# Contains coordinates of detectors for same 12 sets as above
detectors_df_all <- fig5_results %>% purrr::map("detectors_df") %>% map_df(bind_rows)

# change covariate variable names to be consistent with mona_inputs (removing D~ from data frames)
predicted_densities_all <- predicted_densities_all %>% mutate(covtype = str_remove(covtype, "D~"))
detectors_df_all <- detectors_df_all %>% mutate(covtype = str_remove(covtype, "D~"))

# only want the results for 1 occasion (change if desired) -- i.e. only results for 20 sampling occasions, so now 15000 rows (dealing with 6 sets now)
predicted_densities_all <- predicted_densities_all %>% filter(occasions == 20)

# choose covariates we want - only want to work with Dgood and Dblur, so now 10000 rows (4 sets now)
predicted_densities_all <- predicted_densities_all %>% filter(covtype %in% c("Dgood", "Dblur"))
detectors_df_all <- detectors_df_all %>% filter(covtype %in% c("Dgood", "Dblur"))

# choose variables we need
predicted_densities_all <- predicted_densities_all %>%
  select(x, y, value = prob_ac, covtype, occasions, array_origin) %>%
  mutate(value = value / 10000)

predicted_densities_all %>% group_by(covtype, array_origin) %>% summarize(mv = mean(value))


## -----------------------------------------------------------------------------

# So now, we want to compare the density values. First, from the preceding line, we can see that the mean density values are:
# 3.62 for Dblur, 15_15 (ch7e)
# 3.30 for Dblur, 27_31 (ch7f)
# 3.04 for Dgood, 15_15 (ch7b)
# 2.98 for Dgood, 27_31 (ch7c)
# For us, the mean density values are:
mean(ch7e.density) # 3.013 vs 3.62
mean(ch7f.density) # 3.386 vs 3.30
mean(ch7b.density) # 2.943 vs 3.04
mean(ch7c.density) # 3.074 vs 2.98
# Mmm, okay?

# Now, creating scatterplots of my density values vs Ian's
# First, we want to extract the values from Ian's data frame and put them in the same order as ours
head(dgood.df)
head(predicted_densities_all)
# Ian's density values for ch7b (Dgood, 15_15), and just extracting pixel centres and density values
ian.ch7b = predicted_densities_all[1:2500,]
ian.ch7b = ian.ch7b[,c("x", "y", "value")]
# Ian's density values for ch7e (Dblur, 15_15)
ian.ch7e = predicted_densities_all[2501:5000,]
ian.ch7e = ian.ch7e[,c("x", "y", "value")]
# Ian's density values for ch7c (Dgood, 27_31)
ian.ch7c = predicted_densities_all[5001:7500,]
ian.ch7c = ian.ch7c[,c("x", "y", "value")]
# Ian's density values for ch7f (Dblur, 27_31)
ian.ch7f = predicted_densities_all[7501:10000,]
ian.ch7f = ian.ch7f[,c("x", "y", "value")]
# Reordering these in the same way we did the mona_df object when finding our density values
# ch7b
split.ch7b = split(ian.ch7b, ian.ch7b$y)
ian.ch7b = do.call("rbind", split.ch7b)
# ch7e
split.ch7e = split(ian.ch7e, ian.ch7e$y)
ian.ch7e = do.call("rbind", split.ch7e)
# ch7c
split.ch7c = split(ian.ch7c, ian.ch7c$y)
ian.ch7c = do.call("rbind", split.ch7c)
# ch7f
split.ch7f = split(ian.ch7f, ian.ch7f$y)
ian.ch7f = do.call("rbind", split.ch7f)
# Checking the pixel ordering matches what we used for our density values
head(dgood.df)
head(ian.ch7b)
head(ian.ch7e)
head(ian.ch7c)
head(ian.ch7f)
all.equal(dgood.df[,"x"], ian.ch7b[,"x"], ian.ch7e[,"x"], ian.ch7c[,"x"], ian.ch7f[,"x"])
all.equal(dgood.df[,"y"], ian.ch7b[,"y"], ian.ch7e[,"y"], ian.ch7c[,"y"], ian.ch7f[,"y"])
# Looks like all density vectors are in the same order!

# So now, adding our densities to these data frames
ian.ch7b[,'rishika'] = ch7b.density
ian.ch7e[,'rishika'] = ch7e.density
ian.ch7c[,'rishika'] = ch7c.density
ian.ch7f[,'rishika'] = ch7f.density

# And now creating a scatterplot of my estimates against Ian's, so that each point's x-value should be my density value and its y-value should be Ian's density value. Saving these in a pdf to send to Ben.
pdf("Density scatterplots.pdf", height=20, width=20)
# ch7b
plot(ian.ch7b[,'rishika'], ian.ch7b[,'value'], main="ch7b", xlab="Rishika", ylab="Ian")
# Adding a line through intercept with gradient 1
abline(a=0, b=1)
# ch7c
plot(ian.ch7c[,'rishika'], ian.ch7c[,'value'], main="ch7c", xlab="Rishika", ylab="Ian")
abline(a=0, b=1)
# ch7e
plot(ian.ch7e[,'rishika'], ian.ch7e[,'value'], main="ch7e", xlab="Rishika", ylab="Ian")
abline(a=0, b=1)
# ch7f
plot(ian.ch7f[,'rishika'], ian.ch7f[,'value'], main="ch7f", xlab="Rishika", ylab="Ian")
abline(a=0, b=1)
dev.off()

## And now, fitting a model using secr.fit() to our data objects, comparing the estimates to the posterior means we are using to estimate beta0 and beta1
library("secr")
# ch7b
ch7b.beta0
ch7b.beta1
# Posterior mean of beta0 is 2.28 and beta1 is 1.14
# Mask, I think
mlmesh = read.mask(data=mona_df)
# Using secr.fit
ch7b.fit = secr.fit(capthist=ch7b, model=D~log(Dgood), mask=mlmesh, detectfn="HHN")
ch7b.fit
# I think that beta1 is 1.14! :) However, I'm not too sure what beta0 is here

# ch7c
ch7c.beta0
ch7c.beta1
# Posterior mean of beta0 is 2.33 and beta1 is 1.14
# secr.fit
ch7c.fit = secr.fit(capthist=ch7c, model=D~log(Dgood), mask=mlmesh, detectfn="HHN")
ch7c.fit
# Looks like beta1 is 1.14 again!!! Unsure of beta0, unfortunately

# ch7e
ch7e.beta0
ch7e.beta1
# Posterior mean of beta0 is 3.33 and beta1 is 2.23
# secr.fit
ch7e.fit = secr.fit(capthist=ch7e, model=D~log(Dblur), mask=mlmesh, detectfn="HHN")
ch7e.fit
# BETA1 IS 2.22!!!!!!!

# ch7f
ch7f.beta0
ch7f.beta1
# Posterior mean of beta0 is 3.23 and beta1 is 1.97
ch7f.fit = secr.fit(capthist=ch7f, model=D~log(Dblur), mask=mlmesh, detectfn="HHN")
ch7f.fit
# BETA1 IS 1.98 I MIGHT CRY

## -----------------------------------------------------------------------------

## Seems like now is a good time to replace predicted_densities_all so that it now contains our density values!!
# Matrix of pixel centres
source("RUDMaps_Functions.R")
pixel.centres = centres(xrange=c(0.5,50.5), yrange=c(0.5,50.5), x.pixels=50, y.pixels=50)

# ch7b
# Here, covtype=Dgood, occasions=20, array_origin=15_15
ch7b.df = data.frame(pixel.centres, ch7b.density, rep("Dgood", 2500), rep(20, 2500), rep("15_15", 2500))
names(ch7b.df) = c("x", "y", "value", "covtype", "occasions", "array_origin")

# ch7c
# covtype=Dgood, occasions=20, array_origin=27_31
ch7c.df = data.frame(pixel.centres, ch7c.density, rep("Dgood", 2500), rep(20, 2500), rep("27_31", 2500))
names(ch7c.df) = c("x", "y", "value", "covtype", "occasions", "array_origin")

# ch7e
# covtype=Dblur, occasions=20, array_origin=15_15
ch7e.df = data.frame(pixel.centres, ch7e.density, rep("Dblur", 2500), rep(20, 2500), rep("15_15", 2500))
names(ch7e.df) = c("x", "y", "value", "covtype", "occasions", "array_origin")

# ch7f
# covtype=Dblur, occasions=20, array_origin=27_31
ch7f.df = data.frame(pixel.centres, ch7f.density, rep("Dblur", 2500), rep(20, 2500), rep("27_31", 2500))
names(ch7f.df) = c("x", "y", "value", "covtype", "occasions", "array_origin")

# And now, combining into one data frame
predicted_densities_all = rbind(ch7b.df, ch7c.df, ch7e.df, ch7f.df)

# scale the covariate plots to have the same mean as the density plots
predicted_densities_covs <- predicted_densities_covs %>%
  group_by(covtype, array_origin) %>%
  mutate(value = value * mean(predicted_densities_all$value) / mean(predicted_densities_covs$value)) %>%
  ungroup()

# # scale the plots to all sum to one
# predicted_densities_all <- predicted_densities_all %>%
#   group_by(covtype, array_origin) %>%
#   mutate(value = value / sum(value)) %>%
#   ungroup()

sigma <- detectors_df_all$sigma[1]
buffersigmas <- 4
orgn <- data.frame(x = c(15, 15, 27, 27), y = c(15, 31, 15, 31))
yy <- data.frame(array_origin = paste0(orgn$x,"_",orgn$y), xmin = orgn$x - buffersigmas*sigma, xmax = orgn$x + (4+buffersigmas)*sigma, ymin = orgn$y-buffersigmas*sigma, ymax = pmin(50, orgn$y+(6+buffersigmas)*sigma))
yy <- yy %>% filter(array_origin %in% c("15_15", "27_31"))
all_segs <- data.frame(array_origin = as.character(), xmin = as.numeric(), xmax = as.numeric(), ymin = as.numeric(), ymax = as.numeric())
segl = 2
for(i in 1:nrow(yy)){
  array_origin = yy$array_origin[i]
  xmin = yy$xmin[i]
  xmax = yy$xmax[i]
  ymin = yy$ymin[i]
  ymax = yy$ymax[i]
  segs <- data.frame(array_origin = array_origin,
                     xmin = c(xmin, xmax-segl, xmin, xmax-segl, xmin, xmin, xmax, xmax),
                     xmax = c(xmin + segl, xmax, xmin + segl, xmax, xmin, xmin, xmax, xmax),
                     ymin = c(ymin, ymin, ymax, ymax, ymin, ymax-segl, ymin, ymax-segl),
                     ymax = c(ymin, ymin, ymax, ymax, ymin + segl, ymax, ymin + segl, ymax))
  all_segs <- rbind(all_segs, segs)
}
common_area <- data.frame(xmin = max(orgn$x) - buffersigmas*sigma,
                          xmax = min(orgn$x) + (4+buffersigmas)*sigma,
                          ymin = max(orgn$y) - buffersigmas*sigma,
                          ymax = min(orgn$y) + (6+buffersigmas)*sigma)

brew.cols <- brewer.pal(6, "Accent")[-c(1,5)]
#state.cols <- c(mycols, "#B8DE29FF")

orgn <- c(15, 15)
orgn_str <- paste0(orgn[1],"_",orgn[2])

plot_mona <- function(orgn, densities = predicted_densities_all){
  orgn_str <- paste0(orgn[1],"_",orgn[2])
  p <- densities %>%
    filter(array_origin == orgn_str, occasions == 20) %>%
    filter(x >= orgn[1] - buffersigmas*sigma, x <= orgn[1] + (4+buffersigmas)*sigma, y >= orgn[2]-buffersigmas*sigma, y <= orgn[2]+(6+buffersigmas)*sigma) %>%
    ggplot(aes(x, y)) +
    geom_raster(aes(fill = value)) +
    geom_point(data = detectors_df_all %>% filter(array_origin == orgn_str, occasions == 20),
               inherit.aes = T, aes(colour = array_origin), pch = 4, size = 3) +
    geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 1, linetype = 2,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    geom_segment(data = all_segs %>% filter(array_origin == orgn_str), inherit.aes = F, fill = NA, size = 2,
                 aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin)) +
    scale_fill_viridis(direction = 1, limits = c(0,15)) +
    scale_colour_manual(name = "",
                        values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                        breaks=c("15_15", "15_31", "27_15", "27_31")) +
    scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title=element_blank(),legend.position="right", legend.key.height=unit(0.7,"cm"), legend.title=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  return(p)
}

p1good <- plot_mona(orgn = c(15,15), densities = predicted_densities_all %>% filter(covtype == "Dgood"))
p2good <- plot_mona(orgn = c(27,31), densities = predicted_densities_all %>% filter(covtype == "Dgood"))

p1blur <- plot_mona(orgn = c(15,15), densities = predicted_densities_all %>% filter(covtype == "Dblur"))
p2blur <- plot_mona(orgn = c(27,31), densities = predicted_densities_all %>% filter(covtype == "Dblur"))

trap_labels <- detectors_df_all %>% group_by(array_origin) %>%
  summarize(x = mean(x), y = mean(y), array_origin = first(array_origin)) %>%
  mutate(label = c("(c)", "(b)"))

pgoodcov <- predicted_densities_covs %>%
  filter(array_origin == "none", covtype == "Dgood") %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value * 1.5)) +
  # geom_point(data = detectors_df_all %>% filter(occasions == 20),
  #            inherit.aes = T, aes(colour = array_origin, shape = array_origin), size = 2) +
  geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 1, linetype = 2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_segment(data = all_segs, inherit.aes = F, fill = NA, size = 1,
               aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin)) +
  geom_text(data = trap_labels,
            inherit.aes = T, aes(colour = array_origin, label = label), size = 6) +
  scale_fill_viridis(direction = 1, limits = c(0,15)) +
  scale_colour_manual(name="",
                      values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                      breaks=c("15_15", "15_31", "27_15", "27_31")) +
  scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title=element_blank(),legend.position="right", legend.key.height=unit(0.7,"cm"), legend.title=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

pgoodcov

pblurcov <- predicted_densities_covs %>%
  filter(array_origin == "none", covtype == "Dblur") %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = 1.8 * value)) +
  # geom_point(data = detectors_df_all %>% filter(occasions == 20),
  #            inherit.aes = T, aes(colour = array_origin, shape = array_origin), size = 2) +
  geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 1, linetype = 2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_segment(data = all_segs, inherit.aes = F, fill = NA, size = 1,
               aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin)) +
  geom_text(data = trap_labels %>% mutate(label = c("(f)", "(e)")),
            inherit.aes = T, aes(colour = array_origin, label = label), size = 6) +
  scale_fill_viridis(direction = 1, limits = c(0,15)) +
  scale_colour_manual(name="",
                      values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                      breaks=c("15_15", "15_31", "27_15", "27_31")) +
  scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title=element_blank(),legend.position="right", legend.key.height=unit(0.7,"cm"), legend.title=element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


pp <- (pgoodcov | p2good | p1good) / (pblurcov | p2blur | p1blur) + plot_layout(nrow = 2) +
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(',
                  tag_sep = '', tag_suffix = ')') &
  theme(plot.tag = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        plot.tag.position = "bottom")

pp

ggsave("Figure7.png", pp, width=7.5, height=5.8, dpi = 600)

## Isn't bright enough?!
