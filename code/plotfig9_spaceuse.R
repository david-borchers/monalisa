# adds the movement component to an existing map of the activity center densities

library(tidyverse)
library(viridis)

source("code/add_movement_to_acs.R")

load("output/mona_raw_outputs.RData")
load("output/capthist_summaries.RData")

# process the outputs for all animals
predicted_densities_all <- fig6_results %>% 
  purrr::map("predicted_densities") %>% map_df(bind_rows) %>%
  select(x, y, value, covtype, occasions)
detectors_df_all <- fig6_results %>% purrr::map("detectors_df") %>% map_df(bind_rows)
estimated_sigma <- fig6_results %>% purrr::map("estimated_sigma") %>% map_df(bind_rows)

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
                             ac_densities = this_ac_density, sigma = this_sigma)

  # combine list elements
  this_ac_density <- this_ac_density %>% mutate(value = densities_per_cell %>% purrr::reduce(`+`))
  
  # put back into main data frame storing results
  ac_densities_with_movement <- rbind.data.frame(ac_densities_with_movement, this_ac_density)
  
}

ac_densities_with_movement$movetype <- "With movement"

## now do exactly the same thing for captured animals

# process the outputs for all animals
predicted_densities_capt <- fig7_results %>% 
  purrr::map("predicted_densities") %>% map_df(bind_rows) %>%
  select(x, y, value, covtype, occasions)
detectors_df_capt <- fig7_results %>% purrr::map("detectors_df") %>% map_df(bind_rows)
estimated_sigma_capt <- fig7_results %>% purrr::map("estimated_sigma") %>% map_df(bind_rows)

# only want the "no covariate" runs
predicted_densities_capt <- predicted_densities_capt %>% filter(covtype == "D~1") %>% droplevels()
detectors_df_capt <- detectors_df_capt %>% filter(covtype == "D~1") %>% droplevels()
estimated_sigma_capt <- estimated_sigma_capt %>% filter(covtype == "D~1") %>% droplevels()

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
                             ac_densities = this_ac_density, sigma = this_sigma)
  
  # combine list elements
  this_ac_density <- this_ac_density %>% mutate(value = densities_per_cell %>% purrr::reduce(`+`))
  
  # put back into main data frame storing results
  ac_densities_with_movement_capt <- rbind.data.frame(ac_densities_with_movement_capt, this_ac_density)
  
}

ac_densities_with_movement_capt$movetype <- "With movement (seen only)"

# combine no move, move, and move (capt only) results
ac_densities_with_movement <- rbind.data.frame(predicted_densities_all, 
                                               ac_densities_with_movement, 
                                               ac_densities_with_movement_capt)

# scale the plots to have min 0 and max 1 
# (need to think about best way to scale things for visualisation)
ac_densities_with_movement <- ac_densities_with_movement %>% 
  group_by(occasions, movetype) %>%
  mutate(value = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

# detectors are the same for all plots so just extract unique combos of (x,y)
detectors <- detectors_df_all %>% group_by(x,y) %>% count()

# capture histories from get_capthist_summaries.r
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = ch_fig8$n_detections, .y = ch_fig8$n_animals, .f = paster) %>% unlist() 

# relabel factor levels for occasion variable
ac_densities_with_movement$occasions <- factor(ac_densities_with_movement$occasions, 
                                            levels = c(1,3,10,20),
                                            labels = capthist_labels)


p1 <- ac_densities_with_movement %>% 
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis() + 
  facet_grid(movetype ~ occasions) +
  geom_point(data = detectors, aes(x,y), 
             colour = "red", pch = 4, alpha = 0.2, size = 1) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("paper/mona_with_movement.png", p1, width=8, height=6, dpi = 600)

