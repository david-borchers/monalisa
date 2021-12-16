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
predicted_densities_0 <- predicted_densities_0 %>% dplyr::select(x,y,value=expnumber_ac, traps) %>% mutate(value = value * 4)
predicted_densities_0_rt <- predicted_densities_0_rt %>% dplyr::select(x,y,value=expnumber_ac, traps) %>% mutate(value = value * 4)
predicted_densities_0_rt2 <- predicted_densities_0_rt2 %>% dplyr::select(x,y,value=expnumber_ac, traps) %>% mutate(value = value * 4)

# _y results are from secr and so in per ha. Convert to /1km2.
predicted_densities_y <- predicted_densities_y %>% mutate(value = value * 100)
predicted_densities_y_rt <- predicted_densities_y_rt %>% mutate(value = value * 100)
predicted_densities_y_rt2 <- predicted_densities_y_rt2 %>% mutate(value = value * 100)

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

# calculate differences between densities with all traps and densities with reduced traps

# for constant density model
reduced_traps_differences_nocov <- predicted_densities %>% 
  filter(str_detect(traps, "no cov.")) %>%
  left_join(predicted_densities %>% 
              filter(str_detect(traps, "no cov.")) %>%
              filter(str_detect(traps, "All traps")) %>% 
              select(-traps), by = c("x", "y")) %>%
  mutate(valuediff = value.x - value.y)

# for covariate model
reduced_traps_differences_cov <- predicted_densities %>% 
  filter(str_detect(traps, "northing")) %>%
  left_join(predicted_densities %>% 
              filter(str_detect(traps, "northing")) %>%
              filter(str_detect(traps, "All traps")) %>% 
              select(-traps), by = c("x", "y")) %>%
  mutate(valuediff = value.x - value.y)

# join the two
reduced_traps_differences <- rbind(reduced_traps_differences_cov, reduced_traps_differences_nocov) 

# make a data frame containing the two high density points (activity centers), for plotting
highD_cams_pts <- rbind(highD_cams %>% mutate(x = x + c(250, 750), 
                                              y = y + c(500,0), 
                                              traps = "Subset #2, no cov."),
                        highD_cams %>% mutate(x = x + c(250, 750), 
                                              y = y + c(500,0), 
                                              traps = "Subset #2, northing"))

###### end of pre-processing for plotting

### plotting

labs_p1 <- data.frame(x = 600000, y = 1302000, traps = c("All traps, no cov.", "Subset #1, no cov.", "Subset #2, no cov."), label = c("(a)", "(b)", "(c)"))
labs_p2 <- data.frame(x = 600000, y = 1302000, traps = c("Subset #1, no cov.", "Subset #2, no cov."), label = c("(d)", "(e)"))
labs_p3 <- data.frame(x = 600000, y = 1302000, traps = c("All traps, northing", "Subset #1, northing", "Subset #2, northing"), label = c("(a)", "(b)", "(c)"))
labs_p4 <- data.frame(x = 600000, y = 1302000, traps = c("Subset #1, northing", "Subset #2, northing"), label = c("(d)", "(e)"))

# densities for constant density model, different arrays
p1 <- predicted_densities %>% 
  filter(str_detect(traps, "no cov.")) %>% 
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = value)) +
  geom_point(data = detectors %>% filter(str_detect(traps, "no cov.")), 
             inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
  geom_point(data = highD_cams_pts %>% filter(str_detect(traps, "no cov.")), 
             inherit.aes = T, colour = "red", pch = 6, size = 0.8) +
  geom_text(data = labs_p1, aes(label = label)) +
  scale_fill_viridis(name = bquote("Tiger ACs/km"^2), direction = -1, limits = c(0,1.1), breaks = c(0,.25,.50,.75,1)) +
  # scale_fill_viridis(name = "Density", direction = -1, limits = c(0.00005,0.00225), 
  #                    breaks = c(2e-4, 7e-4, 1.2e-3, 1.7e-3, 2.2e-3), 
  #                    labels = c("0.0002", "0.0007", "0.0012", "0.0017", "0.0022")) + 
  facet_grid(. ~ traps) + 
  theme_bw() +
  xlab("Easting (km)") + ylab("Northing (km)") + coord_equal() +
  theme(legend.position="bottom") + 
  theme(panel.background=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()) + 
  guides(fill = guide_colorbar(barwidth = 12, title.vjust=0.75))

p2 <- reduced_traps_differences %>% 
  filter(str_detect(traps, "no cov."), !str_detect(traps, "All traps")) %>% 
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = valuediff)) +
  geom_point(data = detectors %>% filter(str_detect(traps, "no cov."), !str_detect(traps, "All traps")), 
             inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
  geom_text(data = labs_p2, aes(label = label)) +
  #scale_fill_continuous_divergingx(name = "% Change", palette = "RdYlBu", rev = TRUE, mid = 0, l3 = -100, c3 = -200) +
  # scale_fill_continuous_divergingx(name = "% Change", palette = "RdBu", rev = TRUE, mid = 0, l3 = -100, p3 = 0.5, p4 = 1.3) +
  scale_fill_gradientn(name = "Change", colours = coolwarm(22), limits = c(-1,1), breaks = c(-1,-.50,0,.50,1)) +
    # scale_fill_gradientn(name = "% Change", colours = coolwarm(22), limits = c(-100,100), 
    #                    breaks = c(-100,-50,0,50,100), labels = c("-100", "-50", "0", "500", "1000")) +  
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

grid.arrange(p1,p2,ncol=2,widths=c(6,4))

# this is the one we use in the paper
ggsave("paper/tiger_surfaces_nocovs.png", grid.arrange(p1,p2,ncol=2,widths=c(6,3.7)),
       width=9, height=4.5, dpi=600)

# densities for covariate model, different arrays
p3 <- predicted_densities %>% 
  filter(str_detect(traps, "northing")) %>% 
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = value)) +
  geom_point(data = detectors %>% filter(str_detect(traps, "northing")), 
             inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
  geom_text(data = labs_p3, aes(label = label)) +
  scale_fill_viridis(name = bquote("Tiger ACs/km"^2), direction = -1, limits = c(0,1.10), breaks = c(0,.25,.50,.75,1)) +
  # scale_fill_viridis(name = "Density", direction = -1, breaks = c(2e-4, 3e-4, 4e-4, 5e-4), 
  #                    labels = c("0.0002", "0.0003", "0.0004", "0.0005")) + 
  # scale_fill_viridis(name = "Density", direction = -1, limits = c(0.00005,0.00225), 
  #                    breaks = c(2e-4, 7e-4, 1.2e-3, 1.7e-3, 2.2e-3), 
  #                    labels = c("0.0002", "0.0007", "0.0012", "0.0017", "0.0022")) + 
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
  geom_raster(aes(fill = valuediff)) +
  geom_point(data = detectors %>% filter(str_detect(traps, "northing"), !str_detect(traps, "All traps")), 
             inherit.aes = T, colour = "black", pch = 4, alpha = 0.5, size = 0.5) +
  geom_text(data = labs_p4, aes(label = label)) +
  # scale_fill_gradientn(name = "% Change", colours = coolwarm(22),
  #                      limits = c(-maxabsdiff_nocov, maxabsdiff_nocov),
  #                      breaks = c(-15,-10,-5,0,5,10,15)) +
  scale_fill_gradientn(name = "Change", colours = coolwarm(22), limits = c(-1,1), breaks = c(-1,-.5,0,.5,1)) +
  # scale_fill_gradientn(name = "% Change", colours = coolwarm(22), limits = c(-100,100), 
  #                      breaks = c(-100,-50,0,50,100), labels = c("-100", "-50", "0", "500", "1000")) + 
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
dd <- rbind(predicted_densities_0 %>% mutate(covtype = "Realised AC"),
            this_ac_density %>% mutate(covtype = "Realised usage"))

su1 <- dd %>% 
  filter(covtype == "Realised AC") %>%
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = value)) +
  annotate("text", x = 600, y = 1302, label = "(a)") +
  scale_fill_viridis(name = bquote("Tiger ACs/km"^2), direction = -1, limits = c(0,1.1), breaks = c(0,.25,.50,.75,1)) +
  # scale_fill_viridis(name = "Density", direction = -1, breaks = c(2e-4, 1.2e-3, 2.2e-3), 
  #                    labels = c("0.0002", "0.0012", "0.0022")) + 
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
  filter(covtype != "Realised AC") %>%
  ggplot(aes(x/1000, y/1000)) + 
  geom_raster(aes(fill = value)) +
  annotate("text", x = 600, y = 1302, label = "(b)") +
  #scale_fill_viridis(name = bquote("Tigers/km"^2), direction = -1, limits = c(0,0.5)) +
  scale_fill_viridis(name = bquote("Tigers/km"^2), direction = -1, limits = c(0,1.1), breaks = c(0,.25,.50,.75,1)) +
  # scale_fill_viridis(name = "Density", direction = -1, begin = 0, end = 1,
  #                    breaks = c(1e-4, 3e-4, 5e-4), 
  #                    labels = c(expression(1%*%10^{-4}),expression(3%*%10^{-4}),expression(5%*%10^{-4}))) + 
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
       width=6.5, height=5.4, dpi=600)
