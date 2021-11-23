## plotting figure 6: Intensity surfaces .

library(dplyr)
library(ggplot2)
library(stringr)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(purrr)

load("output/mona_raw_outputs.RData")
load("output/capthist_summaries_100sim.RData")

# process the outputs
# average over 100 sim
predicted_densities_all <- res_realised_and_expected_acd_few %>% purrr::map_depth(1, "predicted_densities") %>% map_df(bind_rows)
# REMOVE COMMENT IF YOU WANT just one sim (first one)
#predicted_densities_all <- predicted_densities_all %>% filter(sim_id == 1)
# average over sims
predicted_densities_all <- predicted_densities_all %>% 
  group_by(x, y, covtype, occasions, array_size, array_spacing, array_origin, lambda0, sigma, n_pts) %>%
  summarize(prob_ac = mean(prob_ac),
            expnumber_ac = mean(expnumber_ac)) %>% ungroup()
# predicted densities are calculated in different ways for uniform and non-uniform cases. The latter uses an 
# secr function, the former 'by hand'. Things are in a bit of a mess: I output prob_ac and expnumber_ac, and plot
# expected numbers of ac's in the density plots for the paper.
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
                                          labels = c("None", "Strong", "Moderate"))

detectors_df_all$covtype <- factor(detectors_df_all$covtype, 
                                   levels = c("1", "Dgood", "Dblur"),
                                   labels = c("None", "Strong", "Moderate"))

# capture histories from get_capthist_summaries.r
paster <- function(nd,na){
  paste0(nd," detections\n(",na, " individuals)")
}
capthist_labels <- map2(.x = ch_fig6_sim100$n_detections, .y = ch_fig6_sim100$n_animals, .f = paster) %>% unlist() 
#capthist_labels <- capthist_labels[-1] # didn't use 1 occasion in the end

predicted_densities_all$occasions <- factor(predicted_densities_all$occasions, 
                                            levels = c(1,3,10,20),
                                            labels = capthist_labels)

detectors_df_all$occasions <- factor(detectors_df_all$occasions, 
                                     levels = c(1,3,10,20),
                                     labels = capthist_labels)

predicted_densities_all %>% group_by(occasions, covtype) %>% summarize(minv = min(value), maxv = max(value))

p1 <- predicted_densities_all %>%
  filter(occasions != capthist_labels[1]) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis") +
  facet_grid(covtype ~ occasions) +
  geom_point(data = detectors_df_all %>% filter(occasions != capthist_labels[1]), inherit.aes = T, 
             colour = "red", pch = 4, alpha = 0.5, size = 1) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p1

p2a <- predicted_densities_all %>%
  filter(occasions != capthist_labels[1], covtype == "None") %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,0.3), breaks = c(0,0.1,0.2,0.3)) +
  facet_grid(covtype ~ occasions) +
  geom_point(data = detectors_df_all %>% filter(occasions != capthist_labels[1], covtype == "None"), inherit.aes = T, 
             colour = "black", pch = 4, size = 2) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right", legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2b <- predicted_densities_all %>%
  filter(occasions != capthist_labels[1], covtype == "Strong") %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,0.12), breaks = c(0,0.04,0.08,0.12)) +
  facet_grid(covtype ~ occasions) +
  geom_point(data = detectors_df_all %>% filter(occasions != capthist_labels[1], covtype == "Strong"), inherit.aes = T, 
             colour = "black", pch = 4, size = 2) + 
  theme(strip.text.x = element_blank(), 
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="right", legend.key.height = unit(0.7,"cm"), legend.title = element_blank(),
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2c <- predicted_densities_all %>%
  filter(occasions != capthist_labels[1], covtype == "Moderate") %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(direction = 1, option = "viridis", limits = c(0,0.12), breaks = c(0,0.04,0.08,0.12)) +
  facet_grid(covtype ~ occasions) +
  geom_point(data = detectors_df_all %>% filter(occasions != capthist_labels[1], covtype == "Moderate"), inherit.aes = T, 
             colour = "black", pch = 4, size = 2) + 
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

#ggsave("paper/mona_peaky_unstd.png", p1, width=8, height=6, dpi = 600)
#ggsave("paper/mona_peaky.png", p2, width=8, height=6, dpi = 600)
ggsave("paper/mona_peaky_avgd.png", p2, width=7, height=7, dpi = 600)

