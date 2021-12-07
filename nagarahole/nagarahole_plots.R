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

# concatenate all the predicted densities
predicted_densities <- rbind(predicted_densities_0, predicted_densities_0_rt2, predicted_densities_0_rt, 
                             predicted_densities_y,  predicted_densities_y_rt2, predicted_densities_y_rt)

# concatenate the detectors
detectors <- rbind(data.frame(as.data.frame(cams), traps = "All traps, no cov."), 
                   data.frame(as.data.frame(cams), traps = "All traps, northing"),
                   data.frame(as.data.frame(reduced_cams2), traps = "Subset #2, no cov."),
                   data.frame(as.data.frame(reduced_cams2), traps = "Subset #2, northing"),
                   data.frame(as.data.frame(reduced_cams), traps = "Subset #1, no cov."),
                   data.frame(as.data.frame(reduced_cams), traps = "Subset #1, northing"))

# standardize covariate densities to sum to one (as the no-covariate densities do)
predicted_densities_std_1 <- predicted_densities %>%
  filter(str_detect(traps, "northing")) %>%
  group_by(traps) %>%
  mutate(value = value / sum(value)) %>% # put on desired scale
  ungroup()

predicted_densities_std_2 <- predicted_densities %>%
  filter(str_detect(traps, "no cov."))

predicted_densities_std <- rbind(predicted_densities_std_1, predicted_densities_std_2)

# calculate differences between densities with all traps and densities with reduced traps

# for constant density model
reduced_traps_differences_nocov <- predicted_densities_std %>% 
  filter(str_detect(traps, "no cov.")) %>%
  left_join(predicted_densities_std %>% 
              filter(str_detect(traps, "no cov.")) %>%
              filter(str_detect(traps, "All traps")) %>% 
              select(-traps), by = c("x", "y")) %>%
  mutate(valuediff = value.x - value.y)

# for covariate model
reduced_traps_differences_cov <- predicted_densities_std %>% 
  filter(str_detect(traps, "northing")) %>%
  left_join(predicted_densities_std %>% 
              filter(str_detect(traps, "northing")) %>%
              filter(str_detect(traps, "All traps")) %>% 
              select(-traps), by = c("x", "y")) %>%
  mutate(valuediff = value.x - value.y)

# join the two
reduced_traps_differences <- rbind(reduced_traps_differences_cov, reduced_traps_differences_nocov) %>%
  mutate(valuediff_rel_all = 100 * valuediff / value.y,
         valuediff_rel_red = 100 * valuediff / value.x,
         valuediff_rel = ifelse(valuediff < 0, valuediff_rel_all, valuediff_rel_red)) %>%
  mutate(valuediff_rel_all2 = ifelse(valuediff_rel_all < 0, valuediff_rel_all, valuediff_rel_all/10))

# biggest abs value of differences in all_trap - reduced_trap densities, for setting legend limits
# relative difference
maxabsdiff <- max(abs(reduced_traps_differences$valuediff_rel))
maxabsdiff_nocov <- max(abs(reduced_traps_differences$valuediff_rel[!str_detect(reduced_traps_differences$traps, "no cov.")]))

# # absolute difference
# maxabsdiff <- max(abs(reduced_traps_differences$valuediff))
# maxabsdiff_nocov <- max(abs(reduced_traps_differences$valuediff[!str_detect(reduced_traps_differences$traps, "no cov.")]))

# make a data frame containing the two high density points (activity centers), for plotting
highD_cams_pts <- rbind(highD_cams %>% mutate(x = x + c(250, 750), 
                                              y = y + c(500,0), 
                                              traps = "Subset #2, no cov."),
                        highD_cams %>% mutate(x = x + c(250, 750), 
                                              y = y + c(500,0), 
                                              traps = "Subset #2, northing"))

###### end of pre-processing for plotting

### plotting

# densities for constant density model, different arrays
p1 <- predicted_densities_std %>% 
  filter(str_detect(traps, "no cov.")) %>% 
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = value)) +
  geom_point(data = detectors %>% filter(str_detect(traps, "no cov.")), 
             inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
  geom_point(data = highD_cams_pts %>% filter(str_detect(traps, "no cov.")), 
             inherit.aes = T, colour = "red", pch = 6, size = 0.8) +
  scale_fill_viridis(name = "Density", direction = -1, limits = c(0.00005,0.00225), 
                     breaks = c(2e-4, 7e-4, 1.2e-3, 1.7e-3, 2.2e-3), 
                     labels = c("0.0002", "0.0007", "0.0012", "0.0017", "0.0022")) + 
  facet_grid(. ~ traps) + 
  theme_bw() +
  xlab("Easting (km)") + ylab("Northing (km)") + coord_equal() +
  theme(legend.position="bottom") + 
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) + 
  guides(fill = guide_colorbar(barwidth = 12, title.vjust=0.75))

# density differences for constant density model, different arrays
# differences are relative to the alltraps density
# p2 <- reduced_traps_differences %>% 
#   filter(str_detect(traps, "no cov."), !str_detect(traps, "All traps")) %>% 
#   ggplot(aes(x/1000, y/1000)) + 
#   geom_raster(aes(fill = valuediff_rel_all)) +
#   geom_point(data = detectors %>% filter(str_detect(traps, "no cov."), !str_detect(traps, "All traps")), 
#              inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
#   scale_fill_continuous_divergingx(name = "% Change", palette = "RdYlBu", rev = TRUE, mid = 0, l3 = -100, c3 = -200) +
#  # scale_fill_continuous_divergingx(name = "% Change", palette = "RdBu", rev = TRUE, mid = 0, l3 = -100, p3 = 0.5, p4 = 1.3) +
#  # scale_fill_gradientn(name = "% Change", colours = coolwarm(22), 
#   #                     limits = c(-100, maxabsdiff)) +  
#   #scale_fill_distiller(palette = 'RdBu') +  
#   facet_wrap(~ traps, ncol = 2) + 
#   theme_bw() +
#   xlab("Easting (km)") + ylab("Northing (km)") +
#   theme(legend.position="bottom") + 
#   theme(panel.background=element_blank(),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(), axis.text.y=element_blank(),
#         axis.title.y=element_blank()) + 
#   guides(fill = guide_colorbar(barwidth = 12, title.vjust=0.75))

p2 <- reduced_traps_differences %>% 
  filter(str_detect(traps, "no cov."), !str_detect(traps, "All traps")) %>% 
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = valuediff_rel_all2)) +
  geom_point(data = detectors %>% filter(str_detect(traps, "no cov."), !str_detect(traps, "All traps")), 
             inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
  #scale_fill_continuous_divergingx(name = "% Change", palette = "RdYlBu", rev = TRUE, mid = 0, l3 = -100, c3 = -200) +
  # scale_fill_continuous_divergingx(name = "% Change", palette = "RdBu", rev = TRUE, mid = 0, l3 = -100, p3 = 0.5, p4 = 1.3) +
  scale_fill_gradientn(name = "% Change", colours = coolwarm(22), limits = c(-100,100), 
                       breaks = c(-100,-50,0,50,100), labels = c("-100", "-50", "0", "500", "1000")) +  
  #scale_fill_distiller(palette = 'RdBu') +  
  facet_wrap(~ traps, ncol = 2) + 
  theme_bw() +
  xlab("Easting (km)") + ylab("Northing (km)") + coord_equal() +
  theme(legend.position="bottom") + 
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), axis.text.y=element_blank(),
        axis.title.y=element_blank()) + 
  guides(fill = guide_colorbar(barwidth = 12, title.vjust=0.75))

# # differences are relative to the alltraps density if negative (so min = -100)
# # and relative to the subset of traps density if positive (so max = 100)
# # avoids very big increases where the baseline density is small
# p2r <- reduced_traps_differences %>% 
#   filter(str_detect(traps, "no cov."), !str_detect(traps, "All traps")) %>% 
#   ggplot(aes(x/1000, y/1000)) + 
#   geom_raster(aes(fill = valuediff_rel)) +
#   geom_point(data = detectors %>% filter(str_detect(traps, "no cov."), !str_detect(traps, "All traps")), 
#              inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
#   #scale_fill_continuous_divergingx(name = "% Change", palette = 'RdBu', mid = 0, l1 = -100, p1 = 0.4, p2 = 0.9) +
#   scale_fill_gradientn(name = "% Change", colours = coolwarm(22), limits = c(-100, 100)) +  
#   #scale_fill_distiller(palette = 'RdBu') +  
#   facet_wrap(~ traps, ncol = 2) + 
#   theme_bw() +
#   xlab("Easting (km)") + ylab("Northing (km)") + coord_equal() +
#   theme(legend.position="bottom") + 
#   theme(panel.background=element_blank(),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(), axis.text.y=element_blank(),
#         axis.title.y=element_blank()) + 
#   guides(fill = guide_colorbar(barwidth = 12, title.vjust=0.75))

grid.arrange(p1,p2,ncol=2,widths=c(6,4))
# grid.arrange(p1,p2r,ncol=2,widths=c(6,4))

# this is the one we use in the paper
ggsave("paper/tiger_surfaces_nocovs_v1.png", grid.arrange(p1,p2,ncol=2,widths=c(6,3.7)),
       width=9, height=4.5, dpi=600)

# ggsave("nagarahole/output/tiger_surfaces_nocovs_v2.png", grid.arrange(p1,p2r,ncol=2,widths=c(6,4)),
#        width=9, height=4.5, dpi=600)

# densities for covariate model, different arrays
p3 <- predicted_densities_std %>% 
  filter(str_detect(traps, "northing")) %>% 
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = value)) +
  geom_point(data = detectors %>% filter(str_detect(traps, "northing")), 
             inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
  # scale_fill_viridis(name = "Density", direction = -1, breaks = c(2e-4, 3e-4, 4e-4, 5e-4), 
  #                    labels = c("0.0002", "0.0003", "0.0004", "0.0005")) + 
  scale_fill_viridis(name = "Density", direction = -1, limits = c(0.00005,0.00225), 
                     breaks = c(2e-4, 7e-4, 1.2e-3, 1.7e-3, 2.2e-3), 
                     labels = c("0.0002", "0.0007", "0.0012", "0.0017", "0.0022")) + 
  facet_grid(. ~ traps) + 
  theme_bw() +
  xlab("Easting (km)") + ylab("Northing (km)") + coord_equal() +
  theme(legend.position="bottom") + 
  theme(panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) + 
  guides(fill = guide_colorbar(barwidth = 12, title.vjust=0.75))

# density differences for covariate model, different arrays
p4 <- reduced_traps_differences %>% 
  filter(str_detect(traps, "northing"), !str_detect(traps, "All traps")) %>% 
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = valuediff_rel_all2)) +
  geom_point(data = detectors %>% filter(str_detect(traps, "northing"), !str_detect(traps, "All traps")), 
             inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
  # scale_fill_gradientn(name = "% Change", colours = coolwarm(22),
  #                      limits = c(-maxabsdiff_nocov, maxabsdiff_nocov),
  #                      breaks = c(-15,-10,-5,0,5,10,15)) +
  scale_fill_gradientn(name = "% Change", colours = coolwarm(22), limits = c(-100,100), 
                       breaks = c(-100,-50,0,50,100), labels = c("-100", "-50", "0", "500", "1000")) + 
  # scale_fill_gradientn(name = "% Change", colours = coolwarm(22), limits = c(-100,100), 
  #                      breaks = c(-100,-50,0,50,100), labels = c("-100", "-50", "0", "500", "1000")) +  
  #scale_fill_continuous_divergingx(name = "% Change", palette = 'RdBu', mid = 0, l1 = -100, p1 = 0.4, p2 = 0.9) +
  #scale_fill_distiller(name = "% Change", palette = 'RdBu') +  
  facet_grid(. ~ traps) +
  theme_bw() +
  xlab("Easting (km)") + ylab("Northing (km)") + coord_equal() +
  theme(legend.position="bottom") + 
  theme(panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), axis.text.y=element_blank(),
        axis.title.y=element_blank()) + 
  guides(fill = guide_colorbar(barwidth = 12, title.vjust=0.75))

grid.arrange(p3,p4,ncol=2,widths=c(6,3.7))

ggsave("paper/tiger_surfaces_covs.png", grid.arrange(p3,p4,ncol=2,widths=c(6,3.7)),
       width=9, height=4.5, dpi=600)


###### space use

this_ac_density <- predicted_densities_0 
this_sigma <- exp(fit0$fit$par[3])
densities_per_cell <- map2(this_ac_density$x, this_ac_density$y, .f = add_movement_to_acs, 
                           ac_densities = this_ac_density, sigma = this_sigma, named_density = "value")
this_ac_density <- this_ac_density %>% mutate(value = densities_per_cell %>% purrr::reduce(`+`))

# rbind realised AC density (from before) with this realised usage density and plot togethe
dd <- rbind(predicted_densities_0 %>% mutate(covtype = "Activity centers"),
            this_ac_density %>% mutate(covtype = "Space use"))

su1 <- dd %>% 
  filter(covtype == "Activity centers") %>%
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(name = "Density", direction = -1, breaks = c(2e-4, 1.2e-3, 2.2e-3), 
                     labels = c("0.0002", "0.0012", "0.0022")) + 
  #scale_fill_gradient(low = "black", high = "white") +
  facet_wrap(~ covtype) +
  geom_point(data = detectors, aes(x/1000,y/1000), 
             colour = "red", pch = 4, alpha = 0.2, size = 1) + 
  theme_bw() +
  xlab("Easting (km)") + ylab("Northing (km)") + coord_equal() +
  theme(legend.position="bottom") + 
  theme(panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) + 
  guides(fill = guide_colorbar(barwidth = 8, title.vjust=0.75))

su2 <- dd %>% 
  filter(covtype != "Activity centers") %>%
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = value)) +
  scale_fill_viridis(name = "Density", direction = -1, begin = 0, end = 1,
                     breaks = c(1e-4, 3e-4, 5e-4), 
                     labels = c(expression(1%*%10^{-4}),expression(3%*%10^{-4}),expression(5%*%10^{-4}))) + 
  #scale_fill_gradient(low = "black", high = "white") +
  facet_wrap(~ covtype) +
  geom_point(data = detectors, aes(x/1000,y/1000), 
             colour = "red", pch = 4, alpha = 0.2, size = 1) + 
  theme_bw() +
  xlab("Easting (km)") + ylab("Northing (km)") + coord_equal() +
  theme(legend.position="bottom") + 
  theme(panel.background=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(), axis.text.y=element_blank(),
        axis.title.y=element_blank()) + 
  guides(fill = guide_colorbar(barwidth = 8, title.vjust=0.75, label.hjust = 0.5))

grid.arrange(su1,su2,ncol=2,widths=c(9,7.5))
ggsave("paper/tiger_spaceuse.png", grid.arrange(su1,su2,ncol=2,widths=c(4.85,4)),
       width=5.4, height=4.5, dpi=600)
