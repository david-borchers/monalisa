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
  # Saving density value for ith pixel. This should be the posterior mean for the density for the ith pixel
  ch7b.cv[i] <- sd(density.vec)/mean(density.vec)
}

# ch7c
ch7c.cv <- numeric()
for (i in 1:2500) {
  # Posterior distribution for the density of the ith pixel
  density.vec <- exp(ch7c.sample[,'beta0'] + ch7c.sample[,'beta1']*(log(dgood)[i]))
  # Saving density value for ith pixel. This should be the posterior mean for the density for the ith pixel
  ch7c.cv[i] <- sd(density.vec)/mean(density.vec)
}

# ch7e
ch7e.cv <- numeric()
for (i in 1:2500) {
  # Posterior distribution for the density of the ith pixel
  density.vec <- exp(ch7e.sample[,'beta0'] + ch7e.sample[,'beta1']*(log(dblur)[i]))
  # Saving density value for ith pixel. This should be the posterior mean for the density for the ith pixel
  ch7e.cv[i] <- sd(density.vec)/mean(density.vec)
}

# ch7f
ch7f.cv <- numeric()
for (i in 1:2500) {
  # Posterior distribution for the density of the ith pixel
  density.vec <- exp(ch7f.sample[,'beta0'] + ch7f.sample[,'beta1']*(log(dblur)[i]))
  # Saving density value for ith pixel. This should be the posterior mean for the density for the ith pixel
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
split = split(ch7b.df, ch7b.df$y)
ch7b.df = do.call("rbind", rev(split))

# ch7c
# covtype=Dgood, occasions=20, array_origin=15_15
ch7c.df <- data.frame(pixel.centres, ch7c.cv, rep("Dgood", 2500), rep(20, 2500), rep("15_15", 2500))
names(ch7c.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
split = split(ch7c.df, ch7c.df$y)
ch7c.df = do.call("rbind", rev(split))

# ch7e
# covtype=Dblur, occasions=20, array_origin=27_31
ch7e.df <- data.frame(pixel.centres, ch7e.cv, rep("Dblur", 2500), rep(20, 2500), rep("27_31", 2500))
names(ch7e.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
split = split(ch7e.df, ch7e.df$y)
ch7e.df = do.call("rbind", rev(split))

# ch7f
# covtype=Dblur, occasions=20, array_origin=15_15
ch7f.df <- data.frame(pixel.centres, ch7f.cv, rep("Dblur", 2500), rep(20, 2500), rep("15_15", 2500))
names(ch7f.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
split = split(ch7f.df, ch7f.df$y)
ch7f.df = do.call("rbind", rev(split))

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

# Basically, what we are already doing with the RACD maps is finding the posterior mean of the number of activity centres in a given pixel
# So, what we want is to find the posterior standard deviation of the number of activity centres in a pixel, and divide it by the posterior mean to give the CV (coefficient of variation), which we will then plot in our uncertainty plots!

load("Figure 9/Fig9_MCMC_3occ.RData")
results = results.3occ
xlim = c(0.5, 50.5)
ylim = c(0.5, 50.5)
M = 300

# Points at which local density will be estimated
xg <- seq(xlim[1], xlim[2], by=1)
yg <- seq(ylim[1], ylim[2], by=1)

# Extracting z-values
# Names of variables that have been monitored
names <- names(results[1,])
# Extracting "z" values from MCMC results
Z <- results[,grep("z", names)]

# Extracting activity centres
# Extracting "s" values from MCMC results (i.e. extracting sampled activity centres)
S <- results[,grep("s[^i]", names)]
# x-coordinates of all activity centres
Sx <- S[,1:M]
# y-coordinates of all activity centres
Sy <- S[,-(1:M)]

# For each iteration, want number of animals alive and in each cell -- this is what we store in the 'Dn.vals' matrix
head(Sx)
Dn.vals = matrix(0, nrow=10000, ncol=2500)
for (i in 1:10000) {
  if ((i %% 100) == 0) print(i)
  Sxout = Sx[i,][Z[i,] == 1]
  Sxout = cut(Sxout, breaks=xg, include.lowest=TRUE)
  Syout = Sy[i,][Z[i,] == 1]
  Syout = cut(Syout, breaks=yg, include.lowest=TRUE)
  Dn.vals[i,] = as.vector(table(Sxout, Syout))
}

# Posterior mean of number of animals in each pixel
range(apply(Dn.vals, 2, mean))
posterior.means = apply(Dn.vals, 2, mean)
# WANT TO check that this is the same as the density values used in our RACD map -- if so, can confirm that with an RACD map, we are indeed finding the posterior mean, and we just need to find the posterior sd to find the CV values we want to use in our uncertainty plots

# Density vector we use in our RACD maps:
source("DensityVectorFunction_RACDMaps.R")
density.3occ.all <- no.movement.density.vector(results=results.3occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
range(density.3occ.all)

head(Dn.vals)
head(posterior.means); head(density.3occ.all)
all.equal(posterior.means, density.3occ.all) # Cool! :)
# SO what we are doing with an RACD map IS just finding a posterior mean for each cell!

## Later, can amend no.movement.density.vector() function so that uses a for loop similar to the above -- easier to amend later to work with posterior distn's than the current approach.
## Make sure that somewhere in code, e.g. maybe in no.movement.density.vector(), have some comment about what would do if pixel area was not 1!!! (e.g. would divide mean by pixel area, I think??)
