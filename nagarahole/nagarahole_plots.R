library(secr)
library(dplyr)
library(ggplot2)
library(stringr)
library(rgdal)
library(viridis)
library(maptools)
library(colorspace)
library(pals)
library(gridExtra)
library(purrr)

#source("nagarahole/utilityFuncs.R")
#source("code/predicted_densities_for_D0.R")
source("code/add_movement_to_acs.R")

load("nagarahole/output/nagarole_modelruns.RData")

### some pre-processing before plotting

# _0 results have realised AC density per cell, and a cell is 25ha = 0.25km2. Convert this to ACs/1km2
predicted_densities_0 <- predicted_densities_0 %>% dplyr::select(x,y,value=expnumber_ac, traps) %>% mutate(value = value * 4, surface = "Realised AC")
predicted_densities_0_rt <- predicted_densities_0_rt %>% dplyr::select(x,y,value=expnumber_ac, traps) %>% mutate(value = value * 4, surface = "Realised AC")
predicted_densities_0_rt2 <- predicted_densities_0_rt2 %>% dplyr::select(x,y,value=expnumber_ac, traps) %>% mutate(value = value * 4, surface = "Realised AC")

# _y results are from secr and so in per ha. Convert to /1km2.
predicted_densities_y <- predicted_densities_y %>% mutate(value = value * 100, surface = "Expected AC")
predicted_densities_y_rt <- predicted_densities_y_rt %>% mutate(value = value * 100, surface = "Expected AC")
predicted_densities_y_rt2 <- predicted_densities_y_rt2 %>% mutate(value = value * 100, surface = "Expected AC")

### calculate realised usage densities from realised AC densities

# all traps
usage_density_0 <- predicted_densities_0
this_sigma <- exp(fit0$fit$par[3])
densities_per_cell <- map2(usage_density_0$x, usage_density_0$y, .f = add_movement_to_acs, 
                           ac_densities = usage_density_0, sigma = this_sigma, named_density = "value")
usage_density_0 <- usage_density_0 %>% mutate(value = densities_per_cell %>% purrr::reduce(`+`),
                                              traps = "All traps", surface = "Realised usage")

# subset 1
usage_density_0_rt <- predicted_densities_0_rt
this_sigma <- exp(fit0_rt$fit$par[3])
densities_per_cell <- map2(usage_density_0_rt$x, usage_density_0_rt$y, .f = add_movement_to_acs, 
                           ac_densities = usage_density_0_rt, sigma = this_sigma, named_density = "value")
usage_density_0_rt <- usage_density_0_rt %>% mutate(value = densities_per_cell %>% purrr::reduce(`+`),
                                                    traps = "Subset #1", surface = "Realised usage")

# subset 2
usage_density_0_rt2 <- predicted_densities_0_rt2
this_sigma <- exp(fit0_rt2$fit$par[3])
densities_per_cell <- map2(usage_density_0_rt2$x, usage_density_0_rt2$y, .f = add_movement_to_acs, 
                           ac_densities = usage_density_0_rt2, sigma = this_sigma, named_density = "value")
usage_density_0_rt2 <- usage_density_0_rt2 %>% mutate(value = densities_per_cell %>% purrr::reduce(`+`),
                                                      traps = "Subset #2", surface = "Realised usage")


# concatenate all the predicted densities
predicted_densities <- rbind(predicted_densities_0, predicted_densities_0_rt2, predicted_densities_0_rt, 
                             predicted_densities_y,  predicted_densities_y_rt2, predicted_densities_y_rt,
                             usage_density_0, usage_density_0_rt, usage_density_0_rt2)

predicted_densities <- predicted_densities %>% mutate(traps = str_remove(traps, ",.*"))

# concatenate the detectors
detectors_0 <- data.frame(x = rep(cams$x, 3), y = rep(cams$y, 3), traps = "All traps", surface = rep(c("Realised AC", "Expected AC", "Realised usage"), each = nrow(cams)))
detectors_rt <- data.frame(x = rep(reduced_cams$x, 3), y = rep(reduced_cams$y, 3), traps = "Subset #1", surface = rep(c("Realised AC", "Expected AC", "Realised usage"), each = nrow(reduced_cams)))
detectors_rt2 <- data.frame(x = rep(reduced_cams2$x, 3), y = rep(reduced_cams2$y, 3), traps = "Subset #2", surface = rep(c("Realised AC", "Expected AC", "Realised usage"), each = nrow(reduced_cams2)))
detectors <- rbind(detectors_0, detectors_rt, detectors_rt2)

# calculate differences between densities with all traps and densities with reduced traps

# for constant density model
reduced_traps_differences_rac <- predicted_densities %>% 
  filter(str_detect(surface, "Realised AC")) %>%
  left_join(predicted_densities %>% 
              filter(str_detect(surface, "Realised AC")) %>%
              filter(str_detect(traps, "All traps")) %>% 
              select(-traps,-surface), by = c("x", "y")) %>%
  mutate(valuediff = value.x - value.y)

# for covariate model
reduced_traps_differences_eac <- predicted_densities %>% 
  filter(str_detect(surface, "Expected AC")) %>%
  left_join(predicted_densities %>% 
              filter(str_detect(surface, "Expected AC")) %>%
              filter(str_detect(traps, "All traps")) %>% 
              select(-traps,-surface), by = c("x", "y")) %>%
  mutate(valuediff = value.x - value.y)

# for usage surface 
reduced_traps_differences_ru <- predicted_densities %>% 
  filter(str_detect(surface, "Realised usage")) %>%
  left_join(predicted_densities %>% 
              filter(str_detect(surface, "Realised usage")) %>%
              filter(str_detect(traps, "All traps")) %>% 
              select(-traps,-surface), by = c("x", "y")) %>%
  mutate(valuediff = value.x - value.y)

# join them
reduced_traps_differences <- rbind(reduced_traps_differences_rac, reduced_traps_differences_eac, reduced_traps_differences_ru) 

# make a data frame containing the two high density points (activity centers), for plotting
highD_cams_pts <- highD_cams %>% mutate(x = x + c(250, 750), y = y + c(500,0), traps = "Subset #2") %>% as.data.frame()
highD_cams_pts <- rbind(highD_cams_pts,highD_cams_pts,highD_cams_pts) %>%
  mutate(traps = "Subset #2",
         surface = rep(c("Realised AC", "Expected AC", "Realised usage"), each = 2))

###### end of pre-processing for plotting

### plotting

labs_p1 <- data.frame(x = 600000, y = 1302000, 
                      traps = rep(c("All traps", "Subset #1", "Subset #2"), 3),
                      surface = rep(c("Realised AC", "Expected AC", "Realised usage"), each = 3),
                      label = c("(a)", "(b)", "(c)", "(f)", "(g)", "(h)", "(k)", "(l)", "(m)"))
labs_p2 <- data.frame(x = 600000, y = 1302000, 
                      traps = rep(c("Subset #1", "Subset #2"), 3),
                      surface = rep(c("Realised AC", "Expected AC", "Realised usage"), each = 2),
                      label = c("(d)", "(e)", "(i)", "(j)", "(n)", "(o)"))

predicted_densities <- predicted_densities %>% mutate(surface = factor(surface, levels = c("Realised AC", "Expected AC", "Realised usage")))
detectors <- detectors %>% mutate(surface = factor(surface, levels = c("Realised AC", "Expected AC", "Realised usage")))
highD_cams_pts <- highD_cams_pts %>% mutate(surface = factor(surface, levels = c("Realised AC", "Expected AC", "Realised usage")))
reduced_traps_differences <- reduced_traps_differences %>% mutate(surface = factor(surface, levels = c("Realised AC", "Expected AC", "Realised usage")))
labs_p1 <- labs_p1 %>% mutate(surface = factor(surface, levels = c("Realised AC", "Expected AC", "Realised usage")))
labs_p2 <- labs_p2 %>% mutate(surface = factor(surface, levels = c("Realised AC", "Expected AC", "Realised usage")))

# densities for constant density model, different arrays
p1 <- predicted_densities %>% 
  #filter(str_detect(traps, "no cov.")) %>% 
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = value)) +
  geom_point(data = detectors, 
             inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
  geom_point(data = highD_cams_pts, 
             inherit.aes = T, colour = "red", pch = 6, size = 0.8) +
  geom_text(data = labs_p1, aes(label = label)) +
  scale_fill_viridis(name = bquote("Tigers or Tiger ACs/km"^2), direction = -1, limits = c(0,1.1), breaks = c(0,.25,.50,.75,1)) +
  # scale_fill_viridis(name = "Density", direction = -1, limits = c(0.00005,0.00225), 
  #                    breaks = c(2e-4, 7e-4, 1.2e-3, 1.7e-3, 2.2e-3), 
  #                    labels = c("0.0002", "0.0007", "0.0012", "0.0017", "0.0022")) + 
  facet_grid(surface ~ traps) + 
  theme_bw() +
  coord_equal() +
  #xlab("Easting (km)") + ylab("Northing (km)") + coord_equal() +
  theme(legend.position="bottom") + 
  theme(strip.background.y = element_blank(),
        strip.text.y = element_blank(),
        panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.background=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
  guides(fill = guide_colorbar(barwidth = 10, title.vjust=0.75))

p2 <- reduced_traps_differences %>% 
  filter(!str_detect(traps, "All traps")) %>% 
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = valuediff)) +
  geom_point(data = detectors %>% filter(!str_detect(traps, "All traps")), 
             inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
  geom_text(data = labs_p2, aes(label = label)) +
  #scale_fill_continuous_divergingx(name = "% Change", palette = "RdYlBu", rev = TRUE, mid = 0, l3 = -100, c3 = -200) +
  # scale_fill_continuous_divergingx(name = "% Change", palette = "RdBu", rev = TRUE, mid = 0, l3 = -100, p3 = 0.5, p4 = 1.3) +
  scale_fill_gradientn(name = "Change", colours = coolwarm(22), limits = c(-1,1), breaks = c(-1,-.50,0,.50,1)) +
  # scale_fill_gradientn(name = "% Change", colours = coolwarm(22), limits = c(-100,100), 
  #                    breaks = c(-100,-50,0,50,100), labels = c("-100", "-50", "0", "500", "1000")) +  
  #scale_fill_distiller(palette = 'RdBu') +  
  facet_grid(surface ~ traps) + 
  theme_bw() +
  coord_equal() +
  #xlab("Easting (km)") + ylab("Northing (km)") + coord_equal() +
  theme(legend.position="bottom") + 
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        plot.background=element_blank(),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) + 
  guides(fill = guide_colorbar(barwidth = 10, title.vjust=0.75))

grid.arrange(p1,p2,ncol=2,widths=c(5,4))

# this is the one we use in the paper
ggsave("paper/tiger_surfaces.png", grid.arrange(p1,p2,ncol=2,widths=c(6,4.5)),
       width=8, height=8, dpi=600)
