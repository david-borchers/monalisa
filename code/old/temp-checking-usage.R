# master file

library(dplyr)
library(stringr)
library(purrr)
library(secr)
library(doParallel)
library(ggplot2)

source("code/predicted_densities_for_D0.R")
source("code/run_secr.R")
source("code/add_movement_to_acs.R")

load("output/mona_inputs.RData")
load("output/capthists.RData")

#mona_df <- mona_df %>% select(-cc) %>% distinct()
mlmesh <- read.mask(data = mona_df)

# I load the activity centers generated for the paper, comment out as desired
load("output/simulated_densities.Rdata")

# example of single run (not used in paper) 
x <- run_secr(simulated_points = simulated_points_few,
              secr.fitformula = "D~1", 
              dx = 8, dy = 8, nx = 3, ny = 3, 
              xorig = 10, yorig = 26, 
              sigma = 4, lambda0 = 0.69, noccasions = c(20), 
              capthist = capthists_realised_and_expected_acd_few$capthist[[201]])

i = 601
x2 <- run_secr(simulated_points = simulated_points_few,
              secr.fitformula = parlist3$secr.fitformula[i], 
              dx = parlist3$dx[i], dy = parlist3$dy[i], 
              nx = parlist3$nx[i], ny = parlist3$ny[i], 
              xorig = parlist3$xorig[i], yorig = parlist3$yorig[i], 
              sigma = parlist3$sigma[i], lambda0 = parlist3$lambda0[i], 
              noccasions = parlist3$noccasions[i], 
              capthist = parlist3$capthist[[i]])

# fitted secr model and correct realised AC density from fx.total
m2 <- x$cfit
alld <- fx.total(m2)
alld_df <- data.frame(x = alld$x, y = alld$y, 
                      d_seen = covariates(alld)[,"D.fx"] * attr(mlmesh, "area"), 
                      d_unseen = covariates(alld)[,"D.nc"] * attr(mlmesh, "area"), 
                      d_all = covariates(alld)[,"D.sum"] * attr(mlmesh, "area"))

# realised AC density from my results object
x_eac <- x$predicted_densities$expnumber_ac

# check differences
max(abs(alld_df$d_all-x_eac)) # realised AC density
max(abs(alld_df$d_seen-x$predicted_densities$expnumber_ac_seenonly)) # realised AC density for detected animals
## THERE IS A DIFFERENCE, DRIVEN BY NON-DET ANIMALS

# realised AC density from function used to generate these (run inside main run_secr function)
xd <- predicted_densities_for_D0(m2, traps(m2$capthist), mlmesh, captures.only = FALSE)
max(abs(alld_df$d_all-xd$expnumber_ac))
max(abs(alld_df$d_seen-xd$expnumber_ac_seenonly))
## HERE NO DIFFERENCE

# run_secr must be doing something weird

# run secr line by line
max(abs(alld_df$d_all-predsD0$expnumber_ac))
max(abs(alld_df$d_seen-predsD0$expnumber_ac_seenonly))
expnumber_ac





p1 <- ggplot(alld_df, aes(x = x, y = y)) + 
  geom_raster(aes(fill = d_all)) +
  scale_fill_viridis(direction = 1, option = "viridis") + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="right",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

n <-  dim(m2$capthist)[1]
prob_activity_center_per_seen_indiv <- fxi.secr(m2, i=1:n)
# add up the densities over captured individuals to get expected activity centers
# per cell for captured individuals, divide by n to get density
expnumber_activity_center_all_seen_indiv <- purrr::reduce(prob_activity_center_per_seen_indiv, `+`) 
prob_activity_center_all_seen_indiv <- expnumber_activity_center_all_seen_indiv / n ## this isn't right but isn't used in anything important

p_nd.ac <- 1 - pdot(mlmesh, traps(m2$capthist), detectfn="HHN", detectpar=detectpar(m2), noccasions = 1)
p_ac <- 1 / nrow(mlmesh)
p_nd.nac <- sum(p_nd.ac * 1/(nrow(mlmesh) - 1)) - ((1/(nrow(mlmesh)-1)) * p_nd.ac)
p_nd <- p_nd.ac * p_ac + p_nd.nac * (1 - p_ac)

prob_activity_center_all_unseen_indiv <- p_nd.ac * p_ac / p_nd
estimated_N <- round(region.N(m2)["R.N","estimate"])
unseen_n <- estimated_N - n
exp_ac_unseen <- prob_activity_center_all_unseen_indiv * unseen_n

max(abs(alld_df$d_unseen-exp_ac_unseen))
max(abs(alld_df$d_seen-expnumber_activity_center_all_seen_indiv))

expnumber_activity_center <- n * prob_activity_center_all_seen_indiv +
  unseen_n * prob_activity_center_all_unseen_indiv
max(abs(alld_df$d_all-expnumber_activity_center))


# extract the AC density to work with in this iteration
this_ac_density <- alld_df

# extract the appropriate sigma
this_sigma <- 4

# for each cell in this_ac_density surface, redistribute the probability mass in this
# cell across neighbouring cells according to the movement kernel
densities_per_cell <- map2(this_ac_density$x, this_ac_density$y, 
                           .f = add_movement_to_acs, 
                           ac_densities = this_ac_density, sigma = this_sigma, named_density = "d_all")

# combine list elements
this_ac_density <- this_ac_density %>% mutate(value = densities_per_cell %>% purrr::reduce(`+`))

# put back into main data frame storing results
ac_densities_with_movement <- rbind.data.frame(ac_densities_with_movement, this_ac_density)

p2 <- ggplot(this_ac_density, aes(x = x, y = y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis") + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="right",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

library(patchwork)
p1 + p2



prob_activity_center_all_seen_indiv <- expnumber_activity_center_all_seen_indiv / n
#### MY WAY

# average over 100 sim
predicted_densities_all <- x$predicted_densities
detectors_df_all <- x$detectors_df
estimated_sigma <- x$estimated_sigma

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

# combine no move, move
predicted_densities_all <- predicted_densities_all %>% mutate(value = expnumber_ac)
ac_densities_with_movement <- rbind.data.frame(predicted_densities_all, 
                                               ac_densities_with_movement)


# detectors are the same for all plots so just extract unique combos of (x,y)
detectors <- detectors_df_all %>% group_by(x,y) %>% count()

p1 <- ac_densities_with_movement %>% 
  filter(occasions != capthist_labels[1]) %>%
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
  filter(occasions != capthist_labels[1], movetype == "None") %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = pal, limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3)) + 
  scale_fill_viridis(direction = 1, option = "viridis") + 
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
  filter(occasions != capthist_labels[1], movetype == "With movement") %>%
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

#p2 <- p2a / p2b / 2c
p2 <- p2a / p2b

p2


### BR WAY

xt <- fxi.secr(m0)
length(xt[[1]])

xt2 <- 1-pdot(m0)
