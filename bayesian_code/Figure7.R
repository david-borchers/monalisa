## Code to create Figure 7 of the paper

## The files 'ch7b.R', 'ch7c.R', 'ch7e.R' and 'ch7f.R' in the 'Figure 7' folder contain the code required to run the MCMC for Figure 7.

# Libraries we need
library(nimble)
library(coda)
library(nimbleSCR)

## ----------------------------------------------------------------------------------------

# Checking the MCMC objects

## ----------------------------------------------------------------------------------------

load("Figure 7/ch7b.RData")
load("Figure 7/ch7c.RData")
load("Figure 7/ch7e.RData")
load("Figure 7/ch7f.RData")

# Removing burn-in -- so ran 101,000 iterations, discaring 2500 iterations as burn-in

ch7b.sample <-  ch7b.sample[-c(1:2500),]
ch7c.sample <-  ch7c.sample[-c(1:2500),]
ch7e.sample <- ch7e.sample[-c(1:2500),]
ch7f.sample <- ch7f.sample[-c(1:2500),]


# Checking and saving trace plots
pdf("Figure 7/Figure7_TracePlots.pdf", height=20, width=20)
# ch7b
par(mfrow=c(3,2), mar=c(1.8, 1.8, 1.8, 1.8), oma=c(0, 0, 4, 0))
plot(ch7b.sample[,'lambda0'], type="l", main="lambda0", ylab="lambda0 value")
plot(ch7b.sample[,'sigma'], type="l", main="sigma", ylab="sigma value")
plot(ch7b.sample[,'beta0'], type="l", main="beta0", ylab="beta0 value")
plot(ch7b.sample[,'beta1'], type="l", main="beta1", ylab="beta1 value")
plot(ch7b.sample[,'N'], type="l", main="N", ylab="N value")
plot(ch7b.sample[,'D'], type="l", main="D", ylab="D value")
mtext("ch7b trace plots", outer=TRUE, cex=1.5)
# ch7c
plot(ch7c.sample[,'lambda0'], type="l", main="lambda0", ylab="lambda0 value")
plot(ch7c.sample[,'sigma'], type="l", main="sigma", ylab="sigma value")
plot(ch7c.sample[,'beta0'], type="l", main="beta0", ylab="beta0 value")
plot(ch7c.sample[,'beta1'], type="l", main="beta1", ylab="beta1 value")
plot(ch7c.sample[,'N'], type="l", main="N", ylab="N value")
plot(ch7c.sample[,'D'], type="l", main="D", ylab="D value")
mtext("ch7c trace plots", outer=TRUE, cex=1.5)
# ch7e
plot(ch7e.sample[,'lambda0'], type="l", main="lambda0", ylab="lambda0 value")
plot(ch7e.sample[,'sigma'], type="l", main="sigma", ylab="sigma value")
plot(ch7e.sample[,'beta0'], type="l", main="beta0", ylab="beta0 value")
plot(ch7e.sample[,'beta1'], type="l", main="beta1", ylab="beta1 value")
plot(ch7e.sample[,'N'], type="l", main="N", ylab="N value")
plot(ch7e.sample[,'D'], type="l", main="D", ylab="D value")
mtext("ch7e trace plots", outer=TRUE, cex=1.5)
# ch7f
plot(ch7f.sample[,'lambda0'], type="l", main="lambda0", ylab="lambda0 value")
plot(ch7f.sample[,'sigma'], type="l", main="sigma", ylab="sigma value")
plot(ch7f.sample[,'beta0'], type="l", main="beta0", ylab="beta0 value")
plot(ch7f.sample[,'beta1'], type="l", main="beta1", ylab="beta1 value")
plot(ch7f.sample[,'N'], type="l", main="N", ylab="N value")
plot(ch7f.sample[,'D'], type="l", main="D", ylab="D value")
mtext("ch7f trace plots", outer=TRUE, cex=1.5)
dev.off()
# ch7e trace plots for beta0 and beta1 look like they are hitting some sort of limit, but we haven't placed a limit on beta0 and beta1?


# Checking posterior means for beta0 and beta1
# ch7b
ch7b.beta0 <- mean(ch7b.sample[,"beta0"])
ch7b.beta1 <- mean(ch7b.sample[,"beta1"])
# ch7c
ch7c.beta0 <- mean(ch7c.sample[,"beta0"])
ch7c.beta1 <- mean(ch7c.sample[,"beta1"])
# ch7e
ch7e.beta0 <- mean(ch7e.sample[,"beta0"])
ch7e.beta1 <- mean(ch7e.sample[,"beta1"])
# ch7f
ch7f.beta0 <- mean(ch7f.sample[,"beta0"])
ch7f.beta1 <- mean(ch7f.sample[,"beta1"])



# Checking posterior means for N. True N is 7451 - posterior means seem fine, ch7f being about 1000 above likely due to 'bad' covariate
mean(ch7b.sample[,"N"])
mean(ch7c.sample[,"N"])
mean(ch7e.sample[,"N"])
mean(ch7f.sample[,"N"])

# Checking posterior means for sigma. True sigma is 2 - looks pretty good
mean(ch7b.sample[,'sigma'])
mean(ch7c.sample[,'sigma'])
mean(ch7e.sample[,'sigma'])
mean(ch7f.sample[,'sigma'])

# Checking posterior means for lambda0. True lambda0 is 0.69*20 = 13.8 - looks pretty good
mean(ch7b.sample[,'lambda0'])
mean(ch7c.sample[,'lambda0'])
mean(ch7e.sample[,'lambda0'])
mean(ch7f.sample[,'lambda0'])



# Comparing MCMC results to secr.fit -- overall, looking really good! So looks like the MCMC ran 'correctly'
library("secr")
# Mask
mlmesh <- read.mask(data=mona_df)
# Loading in capture histories
load("../output/capthists.RData")
ch7b.ch7e.capthist <- capthists_expected_acd_many$capthist[1][[1]]
ch7c.ch7f.capthist <- capthists_expected_acd_many$capthist[2][[1]]

# ch7b - looking good!
ch7b.fit <- secr.fit(capthist=ch7b.ch7e.capthist, model=D~log(Dgood_bigD), mask=mlmesh, detectfn="HHN") # logged covariate
## Looking at estimates of beta0 -- looking pretty similar
ch7b.beta0 # beta0 is -10.6741 from MCMC
log(exp(-1.4489734)/10000) # beta0 is -10.65931 from secr.fit
## Comparing to the MCMC credible intervals. They look pretty good:
quantile(ch7b.sample[,'beta1'], probs=c(0.025, 0.5, 0.975))
# secr: (1.002, 1.270) and now is (1.007, 1.272)
quantile(ch7b.sample[,'lambda0']/20, probs=c(0.025, 0.5, 0.975))
# secr: is (0.653, 0.697), and now is (0.653, 0.697)!!!
quantile(ch7b.sample[,'sigma'], probs=c(0.025, 0.5, 0.975))
# secr: is (1.967, 2.024) and now is (1.967, 2.024)!!!

# ch7c - looking good!
ch7c.fit <- secr.fit(capthist=ch7c.ch7f.capthist, model=D~log(Dgood_bigD), mask=mlmesh, detectfn="HHN") # logged covariate
## Looking at estimates of beta0 -- looking pretty similar
ch7c.beta0 # -10.6641 from MCMC
log(exp(-1.4462406)/10000) # -10.65658 from secr.fit
## Comparing to the MCMC credible intervals. They look pretty good:
quantile(ch7c.sample[,'beta1'], probs=c(0.025, 0.5, 0.975))
# secr: (1.023, 1.257) and now is (1.026, 1.261)
quantile(ch7c.sample[,'lambda0']/20, probs=c(0.025, 0.5, 0.975))
# secr: is (0.652, 0.696), and now is (0.652, 0.696)
quantile(ch7c.sample[,'sigma'], probs=c(0.025, 0.5, 0.975))
# secr: is (1.992, 2.051) and now is (1.992, 2.052)

# ch7e - looking good! :)
ch7e.fit <- secr.fit(capthist=ch7b.ch7e.capthist, model=D~log(Dblur_bigD), mask=mlmesh, detectfn="HHN") # logged covariate
## Looking at estimates of beta0 -- looking pretty similar
ch7e.beta0 # -21.8867 from MCMC
log(exp(-12.8007829)/10000) # -22.01112 from secr.fit
## Comparing to the MCMC credible intervals. They look pretty good:
quantile(ch7e.sample[,'beta1'], probs=c(0.025, 0.5, 0.975))
# secr: (1.928, 2.518) and now is (1.942, 2.464)
quantile(ch7e.sample[,'lambda0']/20, probs=c(0.025, 0.5, 0.975))
# secr: is (0.651, 0.695), and now is (0.652, 0.695)
quantile(ch7e.sample[,'sigma'], probs=c(0.025, 0.5, 0.975))
# secr: is (1.962, 2.019) and now is (1.961, 2.019)

# ch7e - looking good!!!
ch7f.fit <- secr.fit(capthist=ch7c.ch7f.capthist, model=D~log(Dblur_bigD), mask=mlmesh, detectfn="HHN") # logged covariate
## Looking at estimates of beta0 -- looking pretty similar
ch7f.beta0 # -19.25819 from MCMC
log(exp(-10.1519810)/10000) # -19.362 from secr.fit
## Comparing to the MCMC credible intervals. They look pretty good:
quantile(ch7f.sample[,'beta1'], probs=c(0.025, 0.5, 0.975))
# secr: (1.762, 2.203) and now is (1.758, 2.185)
quantile(ch7f.sample[,'lambda0']/20, probs=c(0.025, 0.5, 0.975))
# secr: is (0.650, 0.694), and now is (0.650, 0.694)
quantile(ch7f.sample[,'sigma'], probs=c(0.025, 0.5, 0.975))
# secr: is (1.989, 2.048) and now is (1.990, 2.049)


## ---------------------------------------------------------------------------------------

# Necessary objects

## ---------------------------------------------------------------------------------------

# Creating a density vector (vector containing density values for each pixel) for each figure

load("../output/mona_inputs.RData")

# Subsetting the dgood and dblur values we will use
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

# ch7b
ch7b.density <- numeric()
for (i in 1:2500) {
  # Posterior distribution for the density of the ith pixel
  density.vec <- exp(ch7b.sample[,'beta0'] + ch7b.sample[,'beta1']*(log(dgood)[i]))
  # Saving density value for ith pixel. This should be the posterior mean for the density for the ith pixel
  ch7b.density[i] <- mean(density.vec)
}

# ch7c
ch7c.density <- numeric()
for (i in 1:2500) {
  # Posterior distribution for the density of the ith pixel
  density.vec <- exp(ch7c.sample[,'beta0'] + ch7c.sample[,'beta1']*(log(dgood)[i]))
  # Saving density value for ith pixel. This should be the posterior mean for the density for the ith pixel
  ch7c.density[i] <- mean(density.vec)
}

# ch7e
ch7e.density <- numeric()
for (i in 1:2500) {
  # Posterior distribution for the density of the ith pixel
  density.vec <- exp(ch7e.sample[,'beta0'] + ch7e.sample[,'beta1']*(log(dblur)[i]))
  # Saving density value for ith pixel. This should be the posterior mean for the density for the ith pixel
  ch7e.density[i] <- mean(density.vec)
}

# ch7f
ch7f.density <- numeric()
for (i in 1:2500) {
  # Posterior distribution for the density of the ith pixel
  density.vec <- exp(ch7f.sample[,'beta0'] + ch7f.sample[,'beta1']*(log(dblur)[i]))
  # Saving density value for ith pixel. This should be the posterior mean for the density for the ith pixel
  ch7f.density[i] <- mean(density.vec)
}

## ---------------------------------------------------------------------------------------

# Creating the plots

## ---------------------------------------------------------------------------------------

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
load("../output/mona_inputs.RData")

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

# Re-saving frequentist density values for a quick comparison later:
prev_predicted_densities_all <- predicted_densities_all
prev_predicted_densities_all %>% group_by(covtype, array_origin) %>% summarize(mv = mean(value))
# Frequentist mean density values seem to match nicely with MCMC mean density values for ch7b, ch7c, ch7e, ch7f. Can frequentist means above to following means from MCMC:
c(mean(ch7f.density), mean(ch7e.density), mean(ch7c.density), mean(ch7b.density))

# ---------------
# Now, we replace predicted_densities_all so that it contains our density values
# Matrix of pixel centres
source("RUDMaps_Functions.R")
pixel.centres <- centres(xrange=c(0.5,50.5), yrange=c(0.5,50.5), x.pixels=50, y.pixels=50)

# ch7b
# Here, covtype=Dgood, occasions=20, array_origin=15_15
ch7b.df <- data.frame(pixel.centres, ch7b.density, rep("Dgood", 2500), rep(20, 2500), rep("15_15", 2500))
names(ch7b.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
# Reordering to match pixel order
split = split(ch7b.df, ch7b.df$y)
ch7b.df = do.call("rbind", rev(split))

# ch7c
# covtype=Dgood, occasions=20, array_origin=27_31
ch7c.df <- data.frame(pixel.centres, ch7c.density, rep("Dgood", 2500), rep(20, 2500), rep("27_31", 2500))
names(ch7c.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
split = split(ch7c.df, ch7c.df$y)
ch7c.df = do.call("rbind", rev(split))

# ch7e
# covtype=Dblur, occasions=20, array_origin=15_15
ch7e.df <- data.frame(pixel.centres, ch7e.density, rep("Dblur", 2500), rep(20, 2500), rep("15_15", 2500))
names(ch7e.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
split = split(ch7e.df, ch7e.df$y)
ch7e.df = do.call("rbind", rev(split))

# ch7f
# covtype=Dblur, occasions=20, array_origin=27_31
ch7f.df <- data.frame(pixel.centres, ch7f.density, rep("Dblur", 2500), rep(20, 2500), rep("27_31", 2500))
names(ch7f.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
split = split(ch7f.df, ch7f.df$y)
ch7f.df = do.call("rbind", rev(split))

# And now, combining into one data frame
predicted_densities_all <- rbind(ch7b.df, ch7e.df, ch7c.df, ch7f.df)
# ---------------

## Quick check: can plot frequentist vs Bayesian densities
#plot(predicted_densities_all$value, prev_predicted_densities_all$value)
#abline(0, 1, col="goldenrod", lwd=2)
## Looks pretty good

# Scale the covariate plots to have the same mean as the density plots
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

fill_max <- 15 # 9.1 is max in exp data, or 15 to make fill colour scale same as realised_ac_many plot

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
    scale_fill_viridis(direction = 1, limits = c(0,fill_max)) +
    scale_colour_manual(name = "",
                        values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                        breaks=c("15_15", "15_31", "27_15", "27_31")) +
    scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title=element_blank(),legend.position="none",
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

ggsave("Figure7.png", pp, width=7.5, height=5.8, dpi = 600)
