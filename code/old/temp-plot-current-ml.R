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
load("output/simulated_densities.Rdata")

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


detectors_df_all <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

p2a <- predicted_densities_all %>%
  filter(occasions == 20, covtype == "D~1") %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3)) +
  facet_grid(covtype ~ occasions) +
  geom_point(data = detectors_df_all %>% filter(occasions == 20, covtype == "D~1"), inherit.aes = T, 
             colour = "white", pch = 4, size = 2) + 
  geom_point(data = simulated_points_few, inherit.aes = F, aes(x=x,y=y),
              colour = "red", pch = 1, size = 1.5) +
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right", legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2a
