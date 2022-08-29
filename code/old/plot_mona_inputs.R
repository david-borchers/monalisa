## plotting part of tidy_mona.R 

library(dplyr)
library(ggplot2)
library(viridis)
library(colorspace)
library(patchwork)

resized_x <- 50

load("output/mona_inputs.RData")
load("output/simulated_densities.Rdata")

actual_densities <- data.frame(x = mona_df$x, y = mona_df$y, value = mona_df$D_unstd)

# temp stuff
library(imager)
mona <- load.image("output/hires_mona.jpg")
dff <- resize(mona, 200,200, interpolation_type = 2)
bigmona <- as.data.frame(dff) %>% 
  mutate(y = (max(y) - y)/(max(y)-min(y)), x = (x - min(x))/(max(x)-min(x))) %>%
  mutate(y = y * 49 + 1, x = x*49+1)

dff <- resize(mona, 50, 50, interpolation_type = 2)
smallmona <- as.data.frame(dff) %>% 
  mutate(y = (max(y) - y)/(max(y)-min(y)), x = (x - min(x))/(max(x)-min(x))) %>%
  mutate(y = y * 49 + 1, x = x*49+1)

true_density <- mona_df %>% dplyr::select(x,y,value = D_smallD) %>% mutate(value = value/10000) %>% mutate(grid = "True Density")
cov_density <- mona_df %>% dplyr::select(x,y,value = Dblur_smallD) %>% mutate(value = value/10000) %>% mutate(grid = "True Density")
realisation_1 <- simulated_densities_df %>% dplyr::select(x,y,value = dvalue) %>% mutate(value = value * 7451) %>% mutate(grid = "Realisation 1")
maxD <- max(true_density$value, realisation_1$value)
realisation_2 <- simulated_densities_small_df %>% dplyr::select(x,y,value = value) %>% mutate(value = ifelse(value==0,0,maxD)) %>% mutate(grid = "Realisation 2")

#df <- rbind(true_density, realisation_1, realisation_2)
df <- true_density

p1 <- df %>% 
  # filter(grid == "True Density") %>%
  # mutate(grid = factor(grid, levels = c("True Density", "Realisation 1", "Realisation 2"))) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  # geom_point(data = simulated_points_few, inherit.aes = F, aes(x=x,y=y),
  #             colour = "darkorange", pch = 16, size = 1) +
  # scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2, mid = 0.225) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues")) +
  #scale_fill_viridis(limits = c(0, maxD)) +
  #facet_wrap(~ grid, nrow = 1) +
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p1

p2 <- df %>% 
  # filter(grid == "True Density") %>%
  # mutate(grid = factor(grid, levels = c("True Density", "Realisation 1", "Realisation 2"))) %>%
  ggplot(aes(x, y)) + 
  geom_raster(fill = "gray80") +
  geom_point(data = simulated_points_few, inherit.aes = F, aes(x=x,y=y),
             colour = "darkorange", pch = 16, size = 1) +
  # scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2, mid = 0.225) +
  # scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues"), limits = c(0,0.22)) +
  #scale_fill_viridis(limits = c(0, maxD)) +
  #facet_wrap(~ grid, nrow = 1) +
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p2

df2 <- cov_density

p3 <- df2 %>% 
  # filter(grid == "True Density") %>%
  # mutate(grid = factor(grid, levels = c("True Density", "Realisation 1", "Realisation 2"))) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  # geom_point(data = simulated_points_few, inherit.aes = F, aes(x=x,y=y),
  #            colour = "darkorange", pch = 16, size = 1) +
  # scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2, mid = 0.225) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues"), limits = c(0,0.12)) +
  #scale_fill_viridis(limits = c(0, maxD)) +
  #facet_wrap(~ grid, nrow = 1) +
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p3

p4 <- bigmona %>% 
  # filter(grid == "True Density") %>%
  # mutate(grid = factor(grid, levels = c("True Density", "Realisation 1", "Realisation 2"))) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  # geom_point(data = simulated_points_few, inherit.aes = F, aes(x=x,y=y),
  #            colour = "darkorange", pch = 16, size = 1) +
  # scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2, mid = 0.225) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues")) +
  #scale_fill_viridis(limits = c(0, maxD)) +
  #facet_wrap(~ grid, nrow = 1) +
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p4

p5 <- smallmona %>% 
  # filter(grid == "True Density") %>%
  # mutate(grid = factor(grid, levels = c("True Density", "Realisation 1", "Realisation 2"))) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  # geom_point(data = simulated_points_few, inherit.aes = F, aes(x=x,y=y),
  #            colour = "darkorange", pch = 16, size = 1) +
  # scale_fill_continuous_diverging(palette = "Blue-Red 3", l1 = 30, l2 = 100, p1 = .9, p2 = 1.2, mid = 0.225) +
  scale_fill_gradientn(colours = hcl.colors(16, palette = "Blues")) +
  #scale_fill_viridis(limits = c(0, maxD)) +
  #facet_wrap(~ grid, nrow = 1) +
  coord_equal() 
+
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

p5

p1 + p2 + p3 + p4

ggsave("paper/revision/mona_inputdata.png", p1+p2+p3, width=10, height=4, dpi = 300)
