## plotting figure 6: Intensity surfaces .

library(tidyverse)
library(viridis)

load("output/mona_raw_outputs.RData")

# process the outputs
predicted_densities_all <- fig6_results %>% purrr::map("predicted_densities") %>% map_df(bind_rows)
detectors_df_all <- fig6_results %>% purrr::map("detectors_df") %>% map_df(bind_rows)

# change covariate variable names
predicted_densities_all <- predicted_densities_all %>% mutate(covtype = str_remove(covtype, "D~"))
detectors_df_all <- detectors_df_all %>% mutate(covtype = str_remove(covtype, "D~"))

# scale the plots to have min 0 and max 1 
# (need to think about best way to scale things for visualisation)
predicted_densities_all <- predicted_densities_all %>% 
  select(x, y, value, covtype, occasions) %>%
  group_by(covtype, occasions) %>% 
  mutate(value = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

# recode and reorder factor levels
predicted_densities_all$covtype <- factor(predicted_densities_all$covtype, 
                                          levels = c("1", "Dgood", "Dblur"),
                                          labels = c("None", "Strong", "Moderate"))

detectors_df_all$covtype <- factor(detectors_df_all$covtype, 
                                   levels = c("1", "Dgood", "Dblur"),
                                   labels = c("None", "Strong", "Moderate"))

predicted_densities_all$occasions <- factor(predicted_densities_all$occasions, 
                                          levels = c(1,3,10,20),
                                          labels = c("1 occasion", "3 occasions", "10 occasions", "20 occasions"))

detectors_df_all$occasions <- factor(detectors_df_all$occasions, 
                                     levels = c(1,3,10,20),
                                     labels = c("1 occasion", "3 occasions", "10 occasions", "20 occasions"))

p1 <- predicted_densities_all %>%
  ggplot(aes(x, y)) + 
  geom_raster(data = predicted_densities_all, aes(fill = value)) +
  scale_fill_viridis() +
  facet_grid(covtype ~ occasions) +
  geom_point(data = detectors_df_all, inherit.aes = T, 
             colour = "red", pch = 4, alpha = 0.5, size = 1) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("mona_results/mona_peaky.png", p1, width=8, height=6, dpi = 600)
