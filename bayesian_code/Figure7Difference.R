## Code to generate a plot showing the different density values across pixels for both versions of Figure 7.
# Will use the same MCMC objects as in Figure7.R, so will assume that all necessary MCMC checks have already been conducted.

# ----------------

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(purrr)

# Matrix of density values used to make frequentist version of Figure 7
load("../output/mona_raw_outputs.RData")
load("../output/mona_inputs.RData")

# process the outputs -- I think predicted_densities_all is what we want to replace, seems to contain density values for each pixel
predicted_densities_all <- res_expected_acd_many %>% purrr::map("predicted_densities") %>% map_df(bind_rows)

# change covariate variable names to be consistent with mona_inputs
predicted_densities_all <- predicted_densities_all %>% mutate(covtype = str_remove(covtype, "D~log\\(")) %>% mutate(covtype = str_remove(covtype, "_bigD\\)"))

# choose covariates we want
predicted_densities_all <- predicted_densities_all %>% filter(covtype %in% c("Dgood", "Dblur"))

# choose variables we need
predicted_densities_all <- predicted_densities_all %>%
  select(x, y, value = prob_ac, covtype, occasions, array_origin) %>%
  mutate(value = value / 10000)

## For plotting later:
# process the covariates
predicted_densities_covs <- mona_df %>% dplyr::select(x,y,Dgood_bigD,Dblur_bigD) %>%
  pivot_longer(cols=c(Dgood_bigD,Dblur_bigD), names_to = "covtype") %>% arrange(x,y) %>%
  mutate(value = value/10000,
         array_origin = "none") %>%
  mutate(covtype = str_remove(covtype, "_bigD"))
# scale the covariate plots to have the same mean as the density plots
predicted_densities_covs <- predicted_densities_covs %>%
  group_by(covtype, array_origin) %>%
  mutate(value = value * mean(predicted_densities_all$value) / mean(predicted_densities_covs$value)) %>%
  ungroup()

# Re-saving the matrix with a different name:
freq_predicted_densities_all = predicted_densities_all


# ----------------

# Loading in Figure 7 MCMC results
load("Figure 7/ch7b.RData")
load("Figure 7/ch7c.RData")
load("Figure 7/ch7e.RData")
load("Figure 7/ch7f.RData")
# Removing burn-in -- ran 101,000 iterations, discarding 2500 iterations as burn-in
ch7b.sample <-  ch7b.sample[-c(1:2500),]
ch7c.sample <-  ch7c.sample[-c(1:2500),]
ch7e.sample <- ch7e.sample[-c(1:2500),]
ch7f.sample <- ch7f.sample[-c(1:2500),]

# ----------------

# Creating density vectors for each Bayesian plot

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

# ----------------

# Matrix of density values used to make Bayesian version of Figure 7:
# Matrix of pixel centres
source("RUDMaps_Functions.R")
pixel.centres <- centres(xrange=c(0.5,50.5), yrange=c(0.5,50.5), x.pixels=50, y.pixels=50)

# ch7b
# Here, covtype=Dgood, occasions=20, array_origin=15_15
ch7b.df <- data.frame(pixel.centres, ch7b.density, rep("Dgood", 2500), rep(20, 2500), rep("15_15", 2500))
names(ch7b.df) <- c("x", "y", "value", "covtype", "occasions", "array_origin")
# Reordering to match Ian's pixel order
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

# ----------------

# Vector of percentage differences for densities across each plot.

# Checking that both matrices have everything in the same order
# x-values
all.equal(predicted_densities_all[,1], freq_predicted_densities_all[,1])
# y-values
all.equal(predicted_densities_all[,2], freq_predicted_densities_all[,2])
# Covariate type
all.equal(as.character(predicted_densities_all[,4]), freq_predicted_densities_all[,4])
# Occasions
all.equal(predicted_densities_all[,5], freq_predicted_densities_all[,5])
# Array origin
all.equal(as.character(predicted_densities_all[,6]), freq_predicted_densities_all[,6])
# Looking good, so we continue:


# Percentage differences, as a percentage of the densities in the frequentist plot
perc.diff = (predicted_densities_all$value - freq_predicted_densities_all$value)/(freq_predicted_densities_all$value) * 100


# Trying out a plot using Ian's code:
# Replacing the 'value' column in predicted_densities_all with these percentage differences
predicted_densities_all$value = perc.diff

# processing output for detectors
detectors_df_all <- res_expected_acd_many %>% purrr::map("detectors_df") %>% map_df(bind_rows)
# change covariate variable names to be consistent with mona_inputs
detectors_df_all <- detectors_df_all %>% mutate(covtype = str_remove(covtype, "D~log\\(")) %>% mutate(covtype = str_remove(covtype, "_bigD\\)"))
# choose covariates we want
detectors_df_all <- detectors_df_all %>% filter(covtype %in% c("Dgood", "Dblur"))

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
    scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red", space="Lab", limits=c(-100,100)) +
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

# Both of these look good:
p1good <- plot_mona(orgn = c(15,15), densities = predicted_densities_all %>% filter(covtype == "Dgood"))
# Looking at data used here
predicted_densities_all %>% filter(covtype=="Dgood") %>% filter(array_origin=="15_15")
p2good <- plot_mona(orgn = c(27,31), densities = predicted_densities_all %>% filter(covtype == "Dgood"))
predicted_densities_all %>% filter(covtype=="Dgood") %>% filter(array_origin=="27_31")

p1blur <- plot_mona(orgn = c(15,15), densities = predicted_densities_all %>% filter(covtype == "Dblur"))
predicted_densities_all %>% filter(covtype=="Dblur") %>% filter(array_origin=="15_15")
summary((predicted_densities_all %>% filter(covtype=="Dblur") %>% filter(array_origin=="15_15"))$value)
p2blur <- plot_mona(orgn = c(27,31), densities = predicted_densities_all %>% filter(covtype == "Dblur"))
predicted_densities_all %>% filter(covtype=="Dblur") %>% filter(array_origin=="27_31")
summary((predicted_densities_all %>% filter(covtype=="Dblur") %>% filter(array_origin=="27_31"))$value)

trap_labels <- detectors_df_all %>% group_by(array_origin) %>%
  summarize(x = mean(x), y = mean(y), array_origin = first(array_origin)) %>%
  mutate(label = c("(c)", "(b)"))

summary((predicted_densities_covs %>% filter(array_origin == "none", covtype == "Dgood"))$value)

fill_max <- 15

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

ggsave("Figure7Difference.png", pp, width=7.5, height=5.8, dpi=600)

# So now, need to:
# * Confirm colour palette
# * Add legend to interpret percentage differences
# * Also add similar legend to Figure 7, just to interpret densities
# (Maybe ask Ian about the last two)
