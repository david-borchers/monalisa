## Code to put together uncertainty plots of Figures 7 and 9: each pixel will be coloured based on the posterior standard error of its estimated density

# Libraries we need
library(nimble)
library(coda)
library(nimbleSCR)

## ----------------------- Uncertainty plot for Figure 7 ----------------------------------

## ----------------------------------------------------------------------------------------

## Retrieving the necessary objects from the MCMC results

## ----------------------------------------------------------------------------------------

## First, retrieving the necessary MCMC objects:
load("Figure 7/ch7b.RData")
load("Figure 7/ch7c.RData")
load("Figure 7/ch7e.RData")
load("Figure 7/ch7f.RData")

# Removing burn-in -- ran 101,000 iterations, discaring 2500 iterations as burn-in. Can see Figure7.R for checks of the MCMC results, including trace plots
ch7b.sample <-  ch7b.sample[-c(1:2500),]
ch7c.sample <-  ch7c.sample[-c(1:2500),]
ch7e.sample <- ch7e.sample[-c(1:2500),]
ch7f.sample <- ch7f.sample[-c(1:2500),]

## ----------------------------------------------------------------------------------------

## Finding the standard error of posterior distribution of density in each pixel

## ----------------------------------------------------------------------------------------

load("../output/mona_inputs.RData")

# Subsetting the covariates that we will use (we use a covariate we label 'Dgood' for plots (b) and (c), and a covariate labelled 'Dblur' for plots (e) and (f))
# dgood
dgood.df <- mona_df[,c("x", "y", "Dgood_bigD")]
split <- split(dgood.df, dgood.df$y)
dgood.df <- do.call("rbind", split)
dgood <- dgood.df[,"Dgood_bigD"]
# dblur
dblur.df <- mona_df[,c("x", "y", "Dblur_bigD")]
split <- split(dblur.df, dblur.df$y)
dblur.df <- do.call("rbind", split)
dblur <- dblur.df[,"Dblur_bigD"]

# Finding the posterior distribution for the density of each pixel in plots (b), (c), (e) and (f), and saving the standard error of this distribution for each pixel:
# ch7b
ch7b.cv <- numeric()
for (i in 1:2500) {
  # Posterior distribution for the density of the ith pixel
  density.vec <- exp(ch7b.sample[,'beta0'] + ch7b.sample[,'beta1']*(log(dgood)[i]))
  # Saving density value for ith pixel. This should be the posterior sd for the density for the ith pixel
  ch7b.cv[i] <- sd(density.vec)/mean(density.vec)
}

# ch7c
ch7c.cv <- numeric()
for (i in 1:2500) {
  density.vec <- exp(ch7c.sample[,'beta0'] + ch7c.sample[,'beta1']*(log(dgood)[i]))
  ch7c.cv[i] <- sd(density.vec)/mean(density.vec)
}

# ch7e
ch7e.cv <- numeric()
for (i in 1:2500) {
  density.vec <- exp(ch7e.sample[,'beta0'] + ch7e.sample[,'beta1']*(log(dblur)[i]))
  ch7e.cv[i] <- sd(density.vec)/mean(density.vec)
}

# ch7f
ch7f.cv <- numeric()
for (i in 1:2500) {
  density.vec <- exp(ch7f.sample[,'beta0'] + ch7f.sample[,'beta1']*(log(dblur)[i]))
  ch7f.cv[i] <- sd(density.vec)/mean(density.vec)
}

## ----------------------------------------------------------------------------------------

## Creating the plot

## ----------------------------------------------------------------------------------------

# Libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(purrr)

# We need to load in these files
load("../output/mona_raw_outputs.RData")

# Process the covariates -- this object contains the covariate values used to create the first plots in both rows
predicted_densities_covs <- mona_df %>% dplyr::select(x,y,Dgood_bigD,Dblur_bigD) %>%
  pivot_longer(cols=c(Dgood_bigD,Dblur_bigD), names_to = "covtype") %>% arrange(x,y) %>%
  mutate(value = value/10000,
         array_origin = "none") %>%
  mutate(covtype = str_remove(covtype, "_bigD"))

# Process the outputs -- predicted_densities_all is what we want to replace with our MCMC results, as it contains density values for each pixel
predicted_densities_all <- res_expected_acd_many %>% purrr::map("predicted_densities") %>% map_df(bind_rows)
detectors_df_all <- res_expected_acd_many %>% purrr::map("detectors_df") %>% map_df(bind_rows)

# Change covariate variable names to be consistent with mona_inputs
predicted_densities_all <- predicted_densities_all %>% mutate(covtype = str_remove(covtype, "D~log\\(")) %>% mutate(covtype = str_remove(covtype, "_bigD\\)"))
detectors_df_all <- detectors_df_all %>% mutate(covtype = str_remove(covtype, "D~log\\(")) %>% mutate(covtype = str_remove(covtype, "_bigD\\)"))

# Choose covariates we want
predicted_densities_all <- predicted_densities_all %>% filter(covtype %in% c("Dgood", "Dblur"))
detectors_df_all <- detectors_df_all %>% filter(covtype %in% c("Dgood", "Dblur"))

# Choose variables we need
predicted_densities_all <- predicted_densities_all %>%
  select(x, y, value = prob_ac, covtype, occasions, array_origin) %>%
  mutate(value = value / 10000)

# ---------------
# Now, we replace predicted_densities_all so that it contains our density values
# Matrix of pixel centres
source("RUDMaps_Functions.R")
pixel.centres <- centres(xrange=c(0.5,50.5), yrange=c(0.5,50.5), x.pixels=50, y.pixels=50)

# ch7b
# Here, covtype=Dgood, occasions=20, array_origin=27_31
ch7b.df <- data.frame(pixel.centres, ch7b.cv, rep("Dgood", 2500), rep(20, 2500), rep("27_31", 2500))
names(ch7b.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
# Reordering to match pixel order
split <-  split(ch7b.df, ch7b.df$y)
ch7b.df <-  do.call("rbind", rev(split))

# ch7c
# covtype=Dgood, occasions=20, array_origin=15_15
ch7c.df <- data.frame(pixel.centres, ch7c.cv, rep("Dgood", 2500), rep(20, 2500), rep("15_15", 2500))
names(ch7c.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
split <-  split(ch7c.df, ch7c.df$y)
ch7c.df <-  do.call("rbind", rev(split))

# ch7e
# covtype=Dblur, occasions=20, array_origin=27_31
ch7e.df <- data.frame(pixel.centres, ch7e.cv, rep("Dblur", 2500), rep(20, 2500), rep("27_31", 2500))
names(ch7e.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
split <-  split(ch7e.df, ch7e.df$y)
ch7e.df <-  do.call("rbind", rev(split))

# ch7f
# covtype=Dblur, occasions=20, array_origin=15_15
ch7f.df <- data.frame(pixel.centres, ch7f.cv, rep("Dblur", 2500), rep(20, 2500), rep("15_15", 2500))
names(ch7f.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
split <-  split(ch7f.df, ch7f.df$y)
ch7f.df <-  do.call("rbind", rev(split))

# And now, combining into one data frame
predicted_densities_all <- rbind(ch7b.df, ch7e.df, ch7c.df, ch7f.df)
# ---------------

# Scale the covariate plots to have the same mean as the density plots
predicted_densities_covs <- predicted_densities_covs %>%
  group_by(covtype, array_origin) %>%
  mutate(value = value * mean(predicted_densities_all$value) / mean(predicted_densities_covs$value)) %>%
  ungroup()

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

fill_max <- 0.25 # Max standard deviation value is 0.541, so setting fill_max to 1 should make colours easy to see (e.g. if used 'fill_max=15' like for Figure 7, will just have solid colours for all plots, since we are working with values in our pixels that are much smaller than 15)

plot_mona <- function(orgn, densities = predicted_densities_all, legend=FALSE){
  orgn_str <- paste0(orgn[1],"_",orgn[2])
  p <- densities %>%
    filter(array_origin == orgn_str, occasions == 20) %>%
    filter(x >= orgn[1] - buffersigmas*sigma, x <= orgn[1] + (4+buffersigmas)*sigma, y >= orgn[2]-buffersigmas*sigma, y <= orgn[2]+(6+buffersigmas)*sigma) %>%
    ggplot(aes(x, y)) +
    geom_raster(aes(fill = value)) +
    geom_point(data = detectors_df_all %>% filter(array_origin == orgn_str, occasions == 20), inherit.aes = T, aes(colour = array_origin), pch = 4, size = 3, show.legend=F) +
    geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 1, linetype = 2,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    geom_segment(data = all_segs %>% filter(array_origin == orgn_str), inherit.aes = F, fill = NA, size = 2,
                 aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin), show.legend=F) +
    scale_fill_viridis(direction = 1, limits = c(0,fill_max)) +
    scale_colour_manual(name = "",
                        values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                        breaks=c("15_15", "15_31", "27_15", "27_31")) +
    scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title=element_blank(),
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())
  if (legend) {
    p  <- p + theme(legend.position="right", legend.key.height=unit(0.7,"cm"))
  } else {
        p  <- p + theme(legend.position="none")
    }
  return(p)
}

p1good <- plot_mona(orgn = c(15,15), densities = predicted_densities_all %>% filter(covtype == "Dgood"), legend=TRUE)
p2good <- plot_mona(orgn = c(27,31), densities = predicted_densities_all %>% filter(covtype == "Dgood"))

p1blur <- plot_mona(orgn = c(15,15), densities = predicted_densities_all %>% filter(covtype == "Dblur"), legend=TRUE)
p2blur <- plot_mona(orgn = c(27,31), densities = predicted_densities_all %>% filter(covtype == "Dblur"))

trap_labels <- detectors_df_all %>% group_by(array_origin) %>%
  summarize(x = mean(x), y = mean(y), array_origin = first(array_origin)) %>%
  mutate(label = c("(c)", "(b)"))

pgoodcov <- predicted_densities_covs %>%
  filter(array_origin == "none", covtype == "Dgood") %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  # geom_point(data = detectors_df_all %>% filter(occasions == 20),
  #            inherit.aes = T, aes(colour = array_origin, shape = array_origin), size = 2) +
  geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 1, linetype = 2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_segment(data = all_segs, inherit.aes = F, fill = NA, size = 1,
               aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin)) +
  geom_text(data = trap_labels,
            inherit.aes = T, aes(colour = array_origin, label = label), size = 6) +
  scale_fill_viridis(direction = 1, limits = c(0,fill_max)) +
  scale_colour_manual(name="",
                      values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                      breaks=c("15_15", "15_31", "27_15", "27_31")) +
  scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

pgoodcov

pblurcov <- predicted_densities_covs %>%
  filter(array_origin == "none", covtype == "Dblur") %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  # geom_point(data = detectors_df_all %>% filter(occasions == 20),
  #            inherit.aes = T, aes(colour = array_origin, shape = array_origin), size = 2) +
  geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 1, linetype = 2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_segment(data = all_segs, inherit.aes = F, fill = NA, size = 1,
               aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin)) +
  geom_text(data = trap_labels %>% mutate(label = c("(f)", "(e)")),
            inherit.aes = T, aes(colour = array_origin, label = label), size = 6) +
  scale_fill_viridis(direction = 1, limits = c(0,fill_max)) +
  scale_colour_manual(name="",
                      values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                      breaks=c("15_15", "15_31", "27_15", "27_31")) +
  scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title=element_blank(),legend.position="none",
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

ggsave("Figure7_UncertaintyPlot.png", pp, width=7.5, height=5.8, dpi = 600)

# Probably need to do something to fix the alignment of the (c) and (f) captions -- these captions come from 'patchwork' package, so either try to fiddle w/ code or use photo editor later as Ian did.

## ----------------------- Uncertainty plot for Figure 9 ----------------------------------

# Basically, what we are already doing with the RACD maps is finding the posterior mean of the number of activity centres in a given pixel.
# To create our uncertainty plots, we want to find the posterior standard deviation of the number of activity centres in each pixel, and divide this by the posterior mean of the number of activity centres in each pixel (therefore finding the CV for each pixel).

# Function to find the posterior standard deviation for the number of activity centres in each pixel, modelled after no.movement.density.vector(). Posterior sd's will be in corresponding order to posterior means found using no.movement.density.vector().
sd.RACD <- function(xlim, ylim, results, M) {

  ## Pixel centres we are working with
  xg <- seq(xlim[1], xlim[2], by=1)
  yg <- seq(ylim[1], ylim[2], by=1)

  ## Extracting z-values
  # Names of variables that have been monitored
  names <- names(results[1,])
  # Extracting "z" values from MCMC results
  Z <- results[,grep("z", names)]

  ## Extracting activity centres
  # Extracting "s" values from MCMC results (i.e. extracting sampled activity centres)
  S <- results[,grep("s[^i]", names)]
  # x-coordinates of all activity centres
  Sx <- S[,1:M]
  # y-coordinates of all activity centres
  Sy <- S[,-(1:M)]

  ## For each MCMC iteration, storing the number of animals alive and with their activity centres in each cell
  # Number of pixel centres
  npix <-  (length(xg) - 1) * (length(yg) - 1)
  Dn.vals <-  matrix(0, nrow=nrow(results), ncol=npix)
  for (i in 1:nrow(results)) {
    if ((i %% 100) == 0) print(i) # Track progress
    Sxout <-  Sx[i,][Z[i,] == 1]
    Sxout <-  cut(Sxout, breaks=xg, include.lowest=TRUE)
    Syout <-  Sy[i,][Z[i,] == 1]
    Syout <-  cut(Syout, breaks=yg, include.lowest=TRUE)
    Dn.vals[i,] <-  as.vector(table(Sxout, Syout))
  }

  ## Posterior standard deviation for number of activity centres in each pixel
  density.vector <- apply(Dn.vals, 2, sd)
}

# Loading in MCMC results
load("Figure 9/Fig9_MCMC_3occ.RData")
load("Figure 9/Fig9_MCMC_10occ.RData")
load("Figure 9/Fig9_MCMC_20occ.RData")

# Calculating CV for each pixel
source("DensityVectorFunction_RACDMaps.R")
# 3 sampling occasions
cv.3occ <-  sd.RACD(results=results.3occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))/no.movement.density.vector(results=results.3occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
# 10 sampling occasions
cv.10occ <-  sd.RACD(results=results.10occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))/no.movement.density.vector(results=results.10occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
# 20 sampling occasions
cv.20occ <-  sd.RACD(results=results.20occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))/no.movement.density.vector(results=results.20occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))

summary(cv.3occ); summary(cv.10occ); summary(cv.20occ)
# Why do we get NA's for 10 occ and 20 occ? Looking into this:
sd.10occ <-  sd.RACD(results=results.10occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
summary(sd.10occ) # No NA's here
mean.10occ <-  no.movement.density.vector(results=results.10occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
summary(mean.10occ) # No NA's here
0 / 0 # Maybe this is what's happening? Taking a closer look:
cbind(sd.10occ[which(sd.10occ==0)], mean.10occ[which(sd.10occ==0)]) # Yes, and there are 18 of these
summary(cv.10occ) # And we have 18 NA's here... so it looks like we have NA's when the posterior mean for a cell is 0. Seems okay? Checking that the same is occurring for 20 occ...
sd.20occ <-  sd.RACD(results=results.20occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
mean.20occ <-  no.movement.density.vector(results=results.20occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
cbind(sd.20occ[which(sd.20occ==0)], mean.20occ[which(sd.20occ==0)]) # 25 rows of matching 0's
summary(cv.20occ) # 25 NA's. So yes, looks like NA's are occurring due to posterior means of 0. I think is okay...

# Now, wanting to take a look at when we have CV values of 100 (which seems to also occur for 10 occ and 20 occ):
cbind(sd.10occ[which(cv.10occ==100)], mean.10occ[which(cv.10occ==100)]) # So these occur when we have a v v small posterior mean (0.0001) and a sd of 0.01. Seems not crazy?
cbind(sd.20occ[which(cv.20occ==100)], mean.20occ[which(cv.20occ==100)]) # Same thing here. Seems okay?

# The fact that we only get v low posterior means/posterior means of 0 when we increase sampling occasions seems to be what we would expect. Looking at the RACDs for 10 occ and 20 occ, the nature of these maps means that as we get more information the surface becomes more 'spiky' around the detectors. We expect to get lower posterior means in this area than if we had fewer sampling occ (e.g. 3 samp occ). And otherwise, we don't expect density in the areas away from the detector to change. Our results seem to make sense so far! This means that we see a greater range of posterior means as we move from 3, to 10, to 20 occ and we also see a greater range of CV values as we move from 3 to 10 to 20 samp occ. So the values we have *seem* to make sense/have a reasonable explanation.

# Maybe we can replace the NA's with 0's? Just to represent the fact that the posterior means are 0 when we see CV = NA?
cv.10occ[which(is.na(cv.10occ))] <-  0
summary(cv.10occ) # Cool, have replaced the NA's with 0's here
cv.20occ[which(is.na(cv.20occ))] <-  0
summary(cv.20occ)

## ----------------------------------------------------------------------------------------

## Summarising the simulated data (needed for column headings in plot)

## ----------------------------------------------------------------------------------------

## First, sourcing in capture histories to use
load("../output/capthists.RData")

# 3 sampling occasions (first column)
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

## 10 sampling occasions (second column)
second.col <- capthists_realised_and_expected_acd_few$capthist[[101]]
# Summing the capture histories over all 10 sampling occasions
encounterdat.10occ <- matrix(0, nrow=nrow(second.col[,1,]), ncol=ncol(second.col[,1,]))
for (i in 1:10) {
  encounterdat.10occ <- encounterdat.10occ + second.col[,i,]
}
# Creating the data object (uses same trap locs, xlim, ylim as above)
data.10occ <- list(encounter.data = encounterdat.10occ, trap.loc = trap.loc, xlim = xlim, ylim = ylim, n.occasions = 10)

## 20 sampling occasions (third column)
third.col <- capthists_realised_and_expected_acd_few$capthist[[201]]
# Summing the capture histories over all 20 sampling occasions
encounterdat.20occ <- matrix(0, nrow=nrow(third.col[,1,]), ncol=ncol(third.col[,1,]))
for (i in 1:20) {
  encounterdat.20occ <- encounterdat.20occ + third.col[,i,]
}
# Creating the data object
data.20occ <- list(encounter.data = encounterdat.20occ, trap.loc = trap.loc, xlim = xlim, ylim = ylim, n.occasions = 20)

## ----------------------------------------------------------------------------------------

## Creating the plot

## ----------------------------------------------------------------------------------------

# Libraries needed
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(purrr)

# Loading RData object that we need
load("../output/mona_raw_outputs.RData")

# Coordinates of pixel centres we are working with
source("RUDMaps_Functions.R")
pixel.centres <- centres(xrange=c(0.5,50.5), yrange=c(0.5,50.5), x.pixels=50, y.pixels=50)

# Process the outputs
detectors_df_all <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

# ---------------
## Creating objects that contain the density values for each map, along with other required information
# 3 sampling occasions
df.3occ.all <- data.frame(pixel.centres, cv.3occ, as.factor(rep("D~1", 2500)), as.factor(rep("None", 2500)), as.factor(rep(3,2500)))
names(df.3occ.all) <- c("x", "y", "value", "covtype", "movetype", "occasions")
# 10 sampling occasions
df.10occ.all <- data.frame(pixel.centres, cv.10occ, as.factor(rep("D~1", 2500)),  as.factor(rep("None", 2500)), as.factor(rep(10,2500)))
names(df.10occ.all) <- c("x", "y", "value", "covtype", "movetype", "occasions")
# 20 sampling occasions
df.20occ.all <- data.frame(pixel.centres, cv.20occ, as.factor(rep("D~1", 2500)), as.factor(rep("None", 2500)), as.factor(rep(20,2500)))
names(df.20occ.all) <- c("x", "y", "value", "covtype", "movetype", "occasions")
## Combining these objects into one data frame
ac_densities <- rbind.data.frame(df.3occ.all, df.10occ.all, df.20occ.all)
# ---------------

# Detectors are the same for all plots so just extract unique combos of (x,y)
detectors <- detectors_df_all %>% group_by(x,y) %>% count()

# Column labels for plots
capthist_labels <-  paste(c(sum(encounterdat.3occ), sum(encounterdat.10occ), sum(encounterdat.20occ)), "detections\n", paste("(", c(nrow(encounterdat.3occ), nrow(encounterdat.10occ), nrow(encounterdat.20occ)), sep=""),  "individuals)")

# Relabel factor levels for occasion variable
ac_densities$occasions <- factor(ac_densities$occasions,
                                            levels = c(3,10,20),
                                            labels = capthist_labels)

p2a <- ac_densities %>%
  filter(movetype == "None") %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = pal, limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,100), breaks = seq(0, 100, by=20)) +
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

ggsave("Figure9_UncertaintyPlot.png", p2a, width=8, height=2.5, dpi = 600, bg="white")
