## plotting figure 5: Intensity surfaces estimated using a model with density a 
# function on one of four simulated spatially-varying covariates.

library(tidyverse)
library(viridis)

load("output/mona_raw_outputs.RData")
load("output/mona_inputs.RData")

# process the covariates
predicted_densities_covs <- mona_df %>% select(-D) %>% 
  gather(grid, value, -x, -y) %>% arrange(x,y) %>% 
  mutate(covtype = grid,
         array_origin = "none") %>%
  select(-grid) %>% 
  filter(covtype %in% c("Dgood", "Dblur", "Dldv"))

# scale the plots to all sum to one
predicted_densities_covs <- predicted_densities_covs %>% 
  group_by(covtype, array_origin) %>% 
  mutate(value = value / sum(value)) %>%
  ungroup()

# process the outputs
predicted_densities_all <- fig5_results %>% purrr::map("predicted_densities") %>% map_df(bind_rows)
detectors_df_all <- fig5_results %>% purrr::map("detectors_df") %>% map_df(bind_rows)

# change covariate variable names to be consistent with mona_inputs
predicted_densities_all <- predicted_densities_all %>% mutate(covtype = str_remove(covtype, "D~"))
detectors_df_all <- detectors_df_all %>% mutate(covtype = str_remove(covtype, "D~"))

# only want the results for 1 occasion (change if desired)
predicted_densities_all <- predicted_densities_all %>% filter(occasions == 1)

# scale the plots to all sum to one
predicted_densities_all <- predicted_densities_all %>% 
  select(x, y, value, covtype, array_origin) %>%
  group_by(covtype, array_origin) %>% 
  mutate(value = value / sum(value)) %>%
  ungroup()

# combine covariates and results
predicted_densities_all <- rbind(predicted_densities_all, predicted_densities_covs)


# recode and reorder factor levels
predicted_densities_all$covtype <- factor(predicted_densities_all$covtype, 
                                          levels = c("Dgood", "Dblur", "Dldv"),
                                          labels = c("Strong", "Moderate", "Weak"))

detectors_df_all$covtype <- factor(detectors_df_all$covtype, 
                                   levels = c("Dgood", "Dblur", "Dldv"),
                                   labels = c("Strong", "Moderate", "Weak"))

predicted_densities_all$array_origin <- factor(predicted_densities_all$array_origin,
                                         levels = c("none", "27_31", "15_15"),
                                         labels = c("Covariate surface", "Estimated density (Array #1) ", "Estimated density (Array #2)"))

detectors_df_all$array_origin <- factor(detectors_df_all$array_origin,
                                  levels = c("none", "27_31", "15_15"),
                                  labels = c("Covariate surface", "Estimated density (Array #1) ", "Estimated density (Array #2)"))

p1 <- predicted_densities_all %>% 
  ggplot(aes(x, y)) + 
  geom_raster(data = predicted_densities_all, aes(fill = value)) +
  scale_fill_viridis() +
  facet_grid(array_origin ~ covtype) +
  geom_point(data = detectors_df_all, inherit.aes = T, 
             colour = "red", pch = 4, alpha = 0.5, size = 1) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("paper/mona_covariates.png", p1, width=6, height=6, dpi = 600)
