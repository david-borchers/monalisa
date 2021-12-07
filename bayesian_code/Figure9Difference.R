## Plot to visualise differences in density for both versions of Figure 9

# Libraries needed
library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(purrr)
library(pals)

# ----------------

# Creating matrix of densities with frequentist method

source("../code/add_movement_to_acs.R")
load("../output/mona_raw_outputs.RData")
load("../output/capthist_summaries_100sim.RData")

# process the outputs -- note that some of the below code is also needed for plotting later!
# average over 100 sim
predicted_densities_all <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "predicted_densities") %>% map_df(bind_rows)
predicted_densities_all$sim_id <- rep(rep(1:100, each = 2500),9)
# REMOVE COMMENT IF YOU WANT just one sim (first one)
predicted_densities_all <- predicted_densities_all %>% filter(sim_id == 1)
detectors_df_all <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()
estimated_sigma <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "estimated_sigma") %>% map_df(bind_rows)
estimated_sigma <- estimated_sigma %>% mutate(sim_id = rep(1:100, 9))
# REMOVE COMMENT IF YOU WANT just one sim (first one)
estimated_sigma <- estimated_sigma %>% filter(sim_id == 1)

# only want the "no covariate" runs
predicted_densities_all <- predicted_densities_all %>% filter(covtype == "D~1") %>% droplevels()
detectors_df_all <- detectors_df_all %>% filter(covtype == "D~1") %>% droplevels()
estimated_sigma <- estimated_sigma %>% filter(covtype == "D~1") %>% droplevels()

# add "movetype" variable to say whether plot is for no move, move, or move (captures only)
predicted_densities_all$movetype <- "None"
detectors_df_all$movetype <- "None"
estimated_sigma$movetype <- "None"

# set up empty dataframe to store results
ac_densities_with_movement <- data.frame(x = as.integer(),
                                         y = as.integer(),
                                         value = as.numeric(),
                                         covtype = as.character(),
                                         movetype = as.character(),
                                         occasions = as.integer())

occ_list <- unique(predicted_densities_all$occasions)

for(i in occ_list){

  # extract the AC density to work with in this iteration
  this_ac_density <- predicted_densities_all %>% filter(occasions == i)

  # extract the appropriate sigma
  this_sigma <- as.numeric(estimated_sigma[estimated_sigma$occasions == i &
                                             estimated_sigma$movetype == "None",
                                           "sigma"])

  # for each cell in this_ac_density surface, redistribute the probability mass in this
  # cell across neighbouring cells according to the movement kernel
  densities_per_cell <- map2(this_ac_density$x, this_ac_density$y,
                             .f = add_movement_to_acs,
                             ac_densities = this_ac_density, sigma = this_sigma, named_density = "expnumber_ac")

  # combine list elements
  this_ac_density <- this_ac_density %>% mutate(value = densities_per_cell %>% purrr::reduce(`+`))

  # put back into main data frame storing results
  ac_densities_with_movement <- rbind.data.frame(ac_densities_with_movement, this_ac_density)

}

ac_densities_with_movement$movetype <- "With movement"

## now do exactly the same thing for captured animals

# process the outputs for all animals

predicted_densities_capt <- predicted_densities_all
detectors_df_capt <- detectors_df_all
estimated_sigma_capt <- estimated_sigma

# add "movetype" variable to say whether plot is for no move, move, or move (captures only)
predicted_densities_capt$movetype <- "With movement (seen only)"
detectors_df_capt$movetype <- "With movement (seen only)"
estimated_sigma_capt$movetype <- "With movement (seen only)"

# set up empty dataframe to store results
ac_densities_with_movement_capt <- data.frame(x = as.integer(),
                                         y = as.integer(),
                                         value = as.numeric(),
                                         covtype = as.character(),
                                         movetype = as.character(),
                                         occasions = as.integer())

occ_list <- unique(predicted_densities_capt$occasions)

for(i in occ_list){

  # extract the AC density to work with in this iteration
  this_ac_density <- predicted_densities_capt %>% filter(occasions == i)

  # extract the appropriate sigma
  this_sigma <- as.numeric(estimated_sigma_capt[estimated_sigma_capt$occasions == i,
                                           "sigma"])

  # for each cell in this_ac_density surface, redistribute the probability mass in this
  # cell across neighbouring cells according to the movement kernel
  densities_per_cell <- map2(this_ac_density$x, this_ac_density$y,
                             .f = add_movement_to_acs,
                             ac_densities = this_ac_density, sigma = this_sigma, named_density = "expnumber_ac_seenonly")

  # combine list elements
  this_ac_density <- this_ac_density %>% mutate(value = densities_per_cell %>% purrr::reduce(`+`))

  # put back into main data frame storing results
  ac_densities_with_movement_capt <- rbind.data.frame(ac_densities_with_movement_capt, this_ac_density)

}

ac_densities_with_movement_capt$movetype <- "With movement (seen only)"

# combine no move, move, and move (capt only) results
predicted_densities_all <- predicted_densities_all %>% mutate(value = expnumber_ac)
ac_densities_with_movement <- rbind.data.frame(predicted_densities_all,
                                               ac_densities_with_movement)

# Renaming final data frame
freq_ac_densities_with_movement = ac_densities_with_movement

# ----------------

# Using same MCMC objects as in Figure9.R. Creating matrix of densities for MCMC version of Figure 9

# Loading in MCMC results
load("Figure 9/Fig9_MCMC_3occ.RData")
load("Figure 9/Fig9_MCMC_10occ.RData")
load("Figure 9/Fig9_MCMC_20occ.RData")
# Burn-in
results.3occ = results.3occ[-c(1:500),]
results.10occ = results.10occ[-c(1:500),]
results.20occ = results.20occ[-c(1:500),]

# Density vectors for row 1
source("DensityVectorFunction_RACDMaps.R")
density.3occ.all <- no.movement.density.vector(results=results.3occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
density.10occ.all <-  no.movement.density.vector(results=results.10occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))
density.20occ.all <- no.movement.density.vector(results=results.20occ, M=300, xlim=c(0.5, 50.5), ylim=c(0.5, 50.5))

# Density vectors for row 2
load("Figure 9/Fig9_3occ.RData")
load("Figure 9/Fig9_10occ.RData")
load("Figure 9/Fig9_20occ.RData")

source("RUDMaps_Functions.R")

# Coordinates for all pixel centres
pixel.centres <- centres(xrange=c(0.5,50.5), yrange=c(0.5,50.5), x.pixels=50, y.pixels=50)

## Creating the matrix
# 3 sampling occasions
df.3occ.all <- data.frame(pixel.centres, density.3occ.all, as.factor(rep("D~1", 2500)), as.factor(rep("None", 2500)), as.factor(rep(3,2500)))
names(df.3occ.all) <- c("x", "y", "value", "covtype", "movetype", "occasions")
# Reordering to match Ian's pixel order
split = split(df.3occ.all, df.3occ.all$y)
df.3occ.all = do.call("rbind", rev(split))
# 10 sampling occasions
df.10occ.all <- data.frame(pixel.centres, density.10occ.all, as.factor(rep("D~1", 2500)),  as.factor(rep("None", 2500)), as.factor(rep(10,2500)))
names(df.10occ.all) <- c("x", "y", "value", "covtype", "movetype", "occasions")
split = split(df.10occ.all, df.10occ.all$y)
df.10occ.all = do.call("rbind", rev(split))
# 20 sampling occasions
df.20occ.all <- data.frame(pixel.centres, density.20occ.all, as.factor(rep("D~1", 2500)), as.factor(rep("None", 2500)), as.factor(rep(20,2500)))
names(df.20occ.all) <- c("x", "y", "value", "covtype", "movetype", "occasions")
split = split(df.20occ.all, df.20occ.all$y)
df.20occ.all = do.call("rbind", rev(split))
# Combining these objects into one data frame
ac_densities_without_movement <- rbind.data.frame(df.3occ.all, df.10occ.all, df.20occ.all)

# 3 sampling occasions
df.3occ.all.withmov <- data.frame(pixel.centres, density.3occ.all.withmov, as.factor(rep("D~1", 2500)), as.factor(rep("With movement", 2500)), as.factor(rep(3,2500)))
names(df.3occ.all.withmov) <- c("x", "y", "value", "covtype", "movetype", "occasions")
split = split(df.3occ.all.withmov, df.3occ.all.withmov$y)
df.3occ.all.withmov = do.call("rbind", rev(split))
# 10 sampling occasions
df.10occ.all.withmov <- data.frame(pixel.centres, density.10occ.all.withmov, as.factor(rep("D~1", 2500)), as.factor(rep("With movement", 2500)), as.factor(rep(10,2500)))
names(df.10occ.all.withmov) <- c("x", "y", "value", "covtype", "movetype", "occasions")
split = split(df.10occ.all.withmov, df.10occ.all.withmov$y)
df.10occ.all.withmov = do.call("rbind", rev(split))
# 20 sampling occasions
df.20occ.all.withmov <- data.frame(pixel.centres, density.20occ.all.withmov, as.factor(rep("D~1", 2500)), as.factor(rep("With movement", 2500)), as.factor(rep(20,2500)))
names(df.20occ.all.withmov) <- c("x", "y", "value", "covtype", "movetype", "occasions")
split = split(df.20occ.all.withmov, df.20occ.all.withmov$y)
df.20occ.all.withmov = do.call("rbind", rev(split))
# Combining these objects into one data frame
ac_densities_with_movement <- rbind.data.frame(df.3occ.all.withmov, df.10occ.all.withmov, df.20occ.all.withmov)

## Combining both data frames now
ac_densities_with_movement  <- rbind.data.frame(ac_densities_without_movement, ac_densities_with_movement)

# ----------------

# Before we continue, checking that both matrices contain everything in the same order:
all.equal(ac_densities_with_movement$x, freq_ac_densities_with_movement$x)
all.equal(ac_densities_with_movement$y, freq_ac_densities_with_movement$y)
all.equal(as.character(ac_densities_with_movement$covtype), freq_ac_densities_with_movement$covtype)
all.equal(as.character(ac_densities_with_movement$movetype), freq_ac_densities_with_movement$movetype)
all.equal(ac_densities_with_movement$occasions, as.factor(freq_ac_densities_with_movement$occasions))
# All looking good! Continuing:

# So, we want to replace the 'value' column in ac_densities_with_movement with the percentage difference between both value columns, as a percentage of the frequentist densities
perc.diff = (ac_densities_with_movement$value - freq_ac_densities_with_movement$value)/(freq_ac_densities_with_movement$value) * 100
# Replacing the values
ac_densities_with_movement$value = perc.diff

# ----------------

# Continuing with the plot:

# detectors are the same for all plots so just extract unique combos of (x,y)
detectors <- detectors_df_all %>% group_by(x,y) %>% count()

# capture histories from get_capthist_summaries.r
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = ch_fig6_firstsim$n_detections, .y = ch_fig6_firstsim$n_animals, .f = paster) %>% unlist()

# relabel factor levels for occasion variable
ac_densities_with_movement$occasions <- factor(ac_densities_with_movement$occasions,
                                            levels = c(1,3,10,20),
                                            labels = capthist_labels)

# # scale the plots to have min 0 and max 1
# # (need to think about best way to scale things for visualisation)
# ac_densities_with_movement2 <- ac_densities_with_movement %>%
#   group_by(occasions, movetype) %>%
#   mutate(value = (value - min(value)) / (max(value) - min(value))) %>%
#   ungroup()

#pal <- wes_palette("Zissou1", 100, type = "continuous")
p1 <- ac_densities_with_movement %>%
  filter(occasions != capthist_labels[1]) %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red", space="Lab", limits=c(-100, 100)) +
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
  filter(occasions != capthist_labels[1], movetype == "None") %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = pal, limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3)) +
  #scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red", space="Lab", limits=c(-100, 100)) +
  #scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3)) +
  scale_fill_gradientn(name = "% Change", colours = coolwarm(22), limits = c(-100,100),
                       breaks = c(-100,-50,0,50,100), labels = c("-100", "-50", "0", "50", "100")) +
  facet_grid(movetype ~ occasions) +
  geom_point(data = detectors, aes(x,y),
             colour = "black", pch = 4, size = 2) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right", legend.key.height = unit(0.7,"cm"), #legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2b <- ac_densities_with_movement %>%
  filter(occasions != capthist_labels[1], movetype == "With movement") %>%
  ggplot(aes(x, y)) +
  geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = pal, limits = c(0,0.1), breaks = c(0,0.05,0.1)) +
  #scale_fill_viridis(direction = 1, option = "viridis", trans = "log10", limits = c(1,103), oob = squish) +
  #scale_fill_viridis(direction = 1, option = "viridis", trans = "log10") +
  #scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red", space="Lab", limits=c(-100, 100)) +
  scale_fill_gradientn(name = "% Change", colours = coolwarm(22), limits = c(-100,100),
                       breaks = c(-100,-50,0,50,100), labels = c("-100", "-50", "0", "50", "100")) +
  facet_grid(movetype ~ occasions) +
  geom_point(data = detectors, aes(x,y),
             colour = "black", pch = 4, size = 2) +
  theme(strip.text.x = element_blank(),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right", legend.key.height = unit(0.7,"cm"), #legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


#p2 <- p2a / p2b / 2c
p2 <- p2a / p2b

p2

ggsave("Figure9Difference.png", p2, width=8, height=6, dpi=600)

# NOTE that only in 0.0128 of the pixels (1.3%) do we see a percentage difference greater than 100%. No pixels have a difference lower than -100%. Therefore, for interpretability, we have coloured these 1.3% of pixels as if they have a percentage difference of 100% (otherwise, we will see a large decrease in the interpretability of the plots).
