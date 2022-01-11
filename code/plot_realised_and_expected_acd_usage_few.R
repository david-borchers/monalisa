## plotting figure 6: Intensity surfaces .

library(dplyr)
library(ggplot2)
library(stringr)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(purrr)

source("code/add_movement_to_acs.R")

load("output/mona_raw_outputs.RData")
load("output/capthist_summaries_100sim.RData")

# controls if results should be averaged over simulations or just shown for first sim (see last few lines)
average_over_sims <- TRUE

# process the outputs
predicted_densities_all <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "predicted_densities") %>% map_df(bind_rows)
predicted_densities_all$sim_id <- rep(rep(1:100, each = 2500),9)
if(average_over_sims){
  # average over 100 sim
  predicted_densities_all <- predicted_densities_all %>% 
    group_by(x, y, covtype, occasions, array_size, array_spacing, array_origin, lambda0, sigma, n_pts) %>%
    summarize(prob_ac = mean(prob_ac),
              expnumber_ac = mean(expnumber_ac)) %>% ungroup()
  } else {
  # just one sim (first one)
  predicted_densities_all <- predicted_densities_all %>% filter(sim_id == 1)
}

# predicted densities are calculated in different ways for uniform and non-uniform cases. The latter uses an 
# secr function, the former 'by hand'. 
# For D~1 models, expected numbers of ac's are in the expnumber_ac column, as you'd expect
# For D~covariate models, expected numbers of ac's must be calculated as predict density (in the prob_ac column, 
# confusingly) / 10000 (the area of a grid cell in ha)
predicted_densities_all$value <- predicted_densities_all$expnumber_ac * str_detect(predicted_densities_all$covtype, "~1") +
  predicted_densities_all$prob_ac/ 10000 * !str_detect(predicted_densities_all$covtype, "~1") 

# check that its worked
predicted_densities_all %>% group_by(occasions, covtype) %>% summarize(mp = mean(prob_ac),
                                                                       me = mean(expnumber_ac),
                                                                       mv = mean(value))

detectors_df_all <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

# change covariate variable names to be consistent with mona_inputs
predicted_densities_all <- predicted_densities_all %>% mutate(covtype = str_remove(covtype, "D~")) %>% mutate(covtype = str_remove(covtype, "log\\(")) %>% mutate(covtype = str_remove(covtype, "_smallD\\)"))
detectors_df_all <- detectors_df_all %>% mutate(covtype = str_remove(covtype, "D~")) %>% mutate(covtype = str_remove(covtype, "log\\(")) %>% mutate(covtype = str_remove(covtype, "_smallD\\)"))

# recode and reorder factor levels
predicted_densities_all$covtype <- factor(predicted_densities_all$covtype, 
                                          levels = c("1", "Dgood", "Dblur"),
                                          labels = c("Realised AC", "Expected AC (S)", "Expected AC (M)"))

detectors_df_all$covtype <- factor(detectors_df_all$covtype, 
                                   levels = c("1", "Dgood", "Dblur"),
                                   labels = c("Realised AC", "Expected AC (S)", "Expected AC (M)"))

# capture histories from get_capthist_summaries.r
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
if(average_over_sims){
  # average over 100 sim
  capthist_labels <- map2(.x = ch_fig6_sim100$n_detections, .y = ch_fig6_sim100$n_animals, .f = paster) %>% unlist() 
} else {
  # just one sim (first one)
  capthist_labels <- map2(.x = ch_fig6_firstsim$n_detections, .y = ch_fig6_firstsim$n_animals, .f = paster) %>% unlist() 
}
#capthist_labels <- capthist_labels[-1] # didn't use 1 occasion in the end

predicted_densities_all$occasions <- factor(predicted_densities_all$occasions, 
                                            levels = c(1,3,10,20),
                                            labels = capthist_labels)

detectors_df_all$occasions <- factor(detectors_df_all$occasions, 
                                     levels = c(1,3,10,20),
                                     labels = capthist_labels)

predicted_densities_all %>% group_by(occasions, covtype) %>% summarize(minv = min(value), maxv = max(value))

p2a <- predicted_densities_all %>%
  filter(occasions != capthist_labels[1], covtype == "Realised AC") %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3)) +
  facet_grid(covtype ~ occasions) +
  geom_point(data = detectors_df_all %>% filter(occasions != capthist_labels[1], covtype == "None"), inherit.aes = T, 
             colour = "black", pch = 4, size = 2) + 
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right", legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2b <- predicted_densities_all %>%
  filter(occasions != capthist_labels[1], covtype == "Expected AC (S)") %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,0.12), breaks = c(0,0.04,0.08,0.12)) +
  facet_grid(covtype ~ occasions) +
  geom_point(data = detectors_df_all %>% filter(occasions != capthist_labels[1], covtype == "Expected AC (S)"), inherit.aes = T, 
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

p2c <- predicted_densities_all %>%
  filter(occasions != capthist_labels[1], covtype == "Expected AC (M)") %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,0.12), breaks = c(0,0.04,0.08,0.12)) +
  facet_grid(covtype ~ occasions) +
  geom_point(data = detectors_df_all %>% filter(occasions != capthist_labels[1], covtype == "Expected AC (M)"), inherit.aes = T, 
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

p2 <- p2a / p2b / p2c

p2

# adding realised usage

detectors_df_all <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()
estimated_sigma <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "estimated_sigma") %>% map_df(bind_rows)
estimated_sigma <- estimated_sigma %>% mutate(sim_id = rep(1:100, 9))
if(average_over_sims){
  # average over 100 sim
  estimated_sigma <- estimated_sigma %>%
    group_by(covtype, occasions, array_size, array_spacing, array_origin, lambda0, sigma, n_pts) %>%
    summarize(est_sigma = mean(est_sigma)) %>% ungroup()
} else {
  # just one sim (first one)
  estimated_sigma <- estimated_sigma %>% filter(sim_id == 1)
}

# only want the "no covariate" runs
predicted_densities_all <- predicted_densities_all %>% filter(covtype == "Realised AC") %>% droplevels()
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
occ_list_sigma <- unique(estimated_sigma$occasions)

for(i in 1:length(occ_list)){
  
  # extract the AC density to work with in this iteration
  this_ac_density <- predicted_densities_all %>% filter(occasions == occ_list[i])
  
  # extract the appropriate sigma
  this_sigma <- as.numeric(estimated_sigma[estimated_sigma$occasions == occ_list_sigma[i] & 
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

ac_densities_with_movement$movetype <- "Realised usage"

# detectors are the same for all plots so just extract unique combos of (x,y)
detectors <- detectors_df_all %>% group_by(x,y) %>% count()

p2d <- ac_densities_with_movement %>% 
  filter(occasions != capthist_labels[1]) %>%
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


p2d

p2 <- p2a / p2b / p2c / p2d

p2

ggsave("paper/mona_peaky_avgd.png", p2, width=7, height=8, dpi = 600)

# # here only realised AC and realised usage, used with average_over_sims = FALSE for comparison with Bayesian
# p2 <- p2a / p2d
# p2
# ggsave("paper/appendix/figure/fig9-1-mle-onesim.png", p2, width=9, height=6, dpi = 600)

