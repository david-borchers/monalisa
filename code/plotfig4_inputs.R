## plotting part of tidy_mona.R 

library(tidyverse)
library(viridis)

resized_x <- 50

load("output/mona_inputs.RData")
load("output/simulated_densities_for_paper.Rdata")

actual_densities <- data.frame(x = mona_df$x, y = mona_df$y, value = mona_df$D)

predicted_densities_all <- original_densities %>% mutate(grid = "original")
predicted_densities_all <- rbind(predicted_densities_all,
                                 actual_densities %>% mutate(grid = "lowres"),
                                 simulated_densities_df %>% mutate(grid = "simulated"),
                                 simulated_densities_small_df %>% mutate(grid = "simulated small")) 
predicted_densities_all <- predicted_densities_all %>% group_by(grid) %>%
  mutate(value = (value / sum(value)) * (n() / resized_x^2)) %>% ungroup()

detectors_df_all <- data.frame(x = as.integer(), y = as.integer(), grid = as.character())

rm(original_densities, actual_densities, simulated_densities_df, simulated_densities_small_df)

minvalue <- predicted_densities_all %>% filter(grid == "simulated") %>% select(value) %>% min() %>% as.numeric()
maxvalue <- predicted_densities_all %>% filter(grid == "simulated") %>% select(value) %>% max() %>% as.numeric()

predicted_densities_all$grid <- factor(predicted_densities_all$grid, 
                                           levels = c("original", "lowres", "simulated", "simulated small"),
                                           labels = c("Original", "Low Res", "Simulated", "Simulated Small"))

detectors_df_all$grid <- factor(detectors_df_all$grid, 
                                levels = c("original", "lowres", "simulated", "simulated small"),
                                labels = c("Original", "Low Res", "Simulated", "Simulated Small"))

grids_for_rest <- c("Low Res", "Simulated", "Simulated Small")

predicted_densities_orig <- predicted_densities_all %>% filter(grid == "Original") %>% droplevels()
predicted_densities_rest <- predicted_densities_all %>% filter(grid %in% grids_for_rest) %>% droplevels()

# p1 <- predicted_densities_all %>% 
#   ggplot(aes(x, y)) + 
#   geom_raster(data = predicted_densities_orig, aes(fill = value)) +
#   geom_raster(data = predicted_densities_rest, aes(fill = value)) +
#   scale_fill_viridis(limits = c(minvalue, maxvalue)) +
#   facet_wrap(~ grid, nrow = 1) +
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank(),legend.position="none",
#         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),plot.background=element_blank())

p1 <- predicted_densities_rest %>% 
  mutate(grid = recode_factor(grid, "Low Res" = "True Density",
                              "Simulated" = "Realisation 1",
                              "Simulated Small" = "Realisation 2")) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(limits = c(minvalue, maxvalue)) +
  facet_wrap(~ grid, nrow = 1) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("paper/mona_inputdata.png", p1, width=7, height=3, dpi = 300)
