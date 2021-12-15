# adds the movement component to an existing map of the activity center densities

library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(scales)
library(purrr)

source("code/add_movement_to_acs.R")

load("output/mona_raw_outputs.RData")
load("output/capthist_summaries_100sim.RData")

# process the outputs
# average over 100 sim
predicted_densities_all <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "predicted_densities") %>% map_df(bind_rows)
predicted_densities_all$sim_id <- rep(rep(1:100, each = 2500),9)
# REMOVE COMMENT IF YOU WANT just one sim (first one)
predicted_densities_all <- predicted_densities_all %>% filter(sim_id == 1)
# # average over sims
# predicted_densities_all <- predicted_densities_all %>% 
#   group_by(x, y, covtype, occasions, array_size, array_spacing, array_origin, lambda0, sigma, n_pts) %>%
#   summarize(prob_ac = mean(prob_ac),
#             expnumber_ac = mean(expnumber_ac),
#             prob_ac_seenonly = mean(prob_ac_seenonly),
#             expnumber_ac_seenonly = mean(expnumber_ac_seenonly)) %>% ungroup()

# predicted_densities_all <- fig6_results %>% 
#   purrr::map("predicted_densities") %>% map_df(bind_rows) %>%
#   select(x, y, value, covtype, occasions)
detectors_df_all <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()
estimated_sigma <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "estimated_sigma") %>% map_df(bind_rows)
estimated_sigma <- estimated_sigma %>% mutate(sim_id = rep(1:100, 9))
# REMOVE COMMENT IF YOU WANT just one sim (first one)
estimated_sigma <- estimated_sigma %>% filter(sim_id == 1)
# # average over sims
# estimated_sigma <- estimated_sigma %>% 
#   group_by(covtype, occasions, array_size, array_spacing, array_origin, lambda0, sigma, n_pts) %>%
#   summarize(est_sigma = mean(est_sigma)) %>% ungroup()

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
                                               ac_densities_with_movement, 
                                               ac_densities_with_movement_capt)


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
  scale_fill_viridis(direction = 1, option = "viridis") + 
  facet_grid(movetype ~ occasions) +
  geom_point(data = detectors, aes(x,y), 
             colour = "black", pch = 4, size = 2) + 
  coord_equal() +
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
  coord_equal() +
  theme(strip.text.x = element_blank(), 
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(), 
        legend.position="right", legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2c <- ac_densities_with_movement %>% 
  filter(occasions != capthist_labels[1], movetype == "With movement (seen only)") %>%
 # mutate(value = pmax(0.01,value)) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  #scale_fill_gradientn(colours = pal, limits = c(0,0.1), breaks = c(0,0.05,0.1)) + 
  #scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,0.1), breaks = c(0,0.05,0.1)) + 
  #scale_fill_viridis(direction = 1, option = "viridis", trans = "log10", limits = c(1,103), oob = squish) + 
  scale_fill_viridis(direction = 1, option = "viridis", trans = "log10", limits = c(0.001,103), oob = squish,
                     breaks = c(0.01,0.1,1,10,100), labels = c("0.01", "0.1", "1", "10", "100")) + 
  #scale_fill_viridis(direction = 1, option = "viridis", trans = "log10") + 
  facet_grid(movetype ~ occasions) +
  geom_point(data = detectors, aes(x,y), 
             colour = "black", pch = 4, size = 2) + 
  coord_equal() +
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

#ggsave("paper/mona_with_movement_unstd.png", p1, width=8, height=6, dpi = 600)
#ggsave("paper/mona_with_movement.png", p2, width=8, height=6, dpi = 600)
ggsave("paper/mona_with_movement.png", p2, width=8, height=5, dpi = 600)

