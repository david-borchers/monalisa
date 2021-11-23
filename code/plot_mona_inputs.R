## plotting part of tidy_mona.R 

library(dplyr)
library(ggplot2)
library(viridis)

resized_x <- 50

load("output/mona_inputs.RData")
load("output/simulated_densities.Rdata")

actual_densities <- data.frame(x = mona_df$x, y = mona_df$y, value = mona_df$D_unstd)

true_density <- mona_df %>% dplyr::select(x,y,value = D_bigD) %>% mutate(value = value/10000) %>% mutate(grid = "True Density")
realisation_1 <- simulated_densities_df %>% dplyr::select(x,y,value = dvalue) %>% mutate(value = value * 7451) %>% mutate(grid = "Realisation 1")
maxD <- max(true_density$value, realisation_1$value)
realisation_2 <- simulated_densities_small_df %>% dplyr::select(x,y,value = value) %>% mutate(value = ifelse(value==0,0,maxD)) %>% mutate(grid = "Realisation 2")

df <- rbind(true_density, realisation_1, realisation_2)

p1 <- df %>% 
  mutate(grid = factor(grid, levels = c("True Density", "Realisation 1", "Realisation 2"))) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(limits = c(0, maxD)) +
  facet_wrap(~ grid, nrow = 1) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p1

ggsave("paper/mona_inputdata.png", p1, width=7, height=3, dpi = 300)
