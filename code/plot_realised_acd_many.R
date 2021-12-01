## Plot (a) shows the true AC densities. Plots (b), (c), (d) (e) show the estimated realised AC surfaces 
## (averaged over 100 simulations) estimated using a 4×3 array placed at four different locations. The orange, 
## red, grey and yellow corner marks in plots (b) to (e) indicate the location of each of (b) to (e) in plot (a). 
## The white dashed box located in the centre of the Mona Lisa’s face in plot (a) is also shown in plots (b) to (d) 
## so that one can easily compare the predictions of the centre of the face from each arra 
## Currently Figure 6 (21/11/2021)

library(dplyr)
library(ggplot2)
library(stringr)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(purrr)

load("output/mona_raw_outputs.RData")
load("output/simulated_densities.RData")
load("output/capthist_summaries_100sim.RData")

# add a variable to the simulated data so that can rbind to results a bit later
simulated_densities_df <- simulated_densities_df %>% 
  mutate(occasions = 20, array_origin = "none", value = dvalue * 7451) %>%
  dplyr::select(-dvalue)

# process the outputs
# # # for 1 sim
# predicted_densities_all <- fig34_results %>% purrr::map("predicted_densities") %>% map_df(bind_rows)
# detectors_df_all <- fig34_results %>% purrr::map("detectors_df") %>% map_df(bind_rows)
# average over 100 sim
predicted_densities_all <- res_realised_acd_many %>% purrr::map_depth(1, "predicted_densities") %>% map_df(bind_rows)
# # just one sim (first one)
# predicted_densities_onesim <- predicted_densities_all %>% 
#   group_by(x, y, covtype, occasions, array_size, array_spacing, array_origin, lambda0, sigma, n_pts) %>%
#   summarize(value = first(value)) %>% ungroup()
# average over sims
predicted_densities_all <- predicted_densities_all %>% 
  group_by(x, y, covtype, occasions, array_size, array_spacing, array_origin, lambda0, sigma, n_pts) %>%
  summarize(value = mean(expnumber_ac)) %>% ungroup()

detectors_df_all <- res_realised_acd_many %>% purrr::map_depth(1, "detectors_df") %>% map_df(bind_rows)
detectors_df_all <- detectors_df_all %>% distinct()

sigma <- detectors_df_all$sigma[1]

# change covariate variable names to be consistent with mona_inputs
predicted_densities_all <- predicted_densities_all %>% mutate(covtype = str_remove(covtype, "D~"))
detectors_df_all <- detectors_df_all %>% mutate(covtype = str_remove(covtype, "D~"))

# select only the variables we need
predicted_densities_all <- predicted_densities_all %>% select(x, y, value, occasions, array_origin)

# combine with simulated data
predicted_densities_all <- rbind(predicted_densities_all,simulated_densities_df)
#predicted_densities_all <- rbind(predicted_densities_all,simulated_densities_df %>% select(-dvalue))

brew.cols <- brewer.pal(6, "Accent")[-c(1,5)]
#state.cols <- c(mycols, "#B8DE29FF") 

buffersigmas <- 4
orgn <- data.frame(x = c(15, 15, 27, 27), y = c(15, 31, 15, 31))
yy <- data.frame(array_origin = paste0(orgn$x,"_",orgn$y), xmin = orgn$x - buffersigmas*sigma, xmax = orgn$x + (4+buffersigmas)*sigma, ymin = orgn$y-buffersigmas*sigma, ymax = pmin(50, orgn$y+(6+buffersigmas)*sigma))
all_segs <- data.frame(array_origin = as.character(), xmin = as.numeric(), xmax = as.numeric(), ymin = as.numeric(), ymax = as.numeric())
segl = 2
for(i in 1:nrow(yy)){
  array_origin = yy$array_origin[i]
  xmin = yy$xmin[i]
  xmax = yy$xmax[i]
  ymin = yy$ymin[i]
  ymax = yy$ymax[i]
  segs <- data.frame(array_origin = array_origin,
                     xmin = c(xmin, xmax-segl, xmin, xmax-segl, xmin, xmin, xmax, xmax),
                     xmax = c(xmin + segl, xmax, xmin + segl, xmax, xmin, xmin, xmax, xmax),
                     ymin = c(ymin, ymin, ymax, ymax, ymin, ymax-segl, ymin, ymax-segl),
                     ymax = c(ymin, ymin, ymax, ymax, ymin + segl, ymax, ymin + segl, ymax))
  all_segs <- rbind(all_segs, segs)
}
common_area <- data.frame(xmin = max(orgn$x) - buffersigmas*sigma, 
                          xmax = min(orgn$x) + (4+buffersigmas)*sigma, 
                          ymin = max(orgn$y) - buffersigmas*sigma, 
                          ymax = min(orgn$y) + (6+buffersigmas)*sigma)

trap_labels <- detectors_df_all %>% group_by(array_origin) %>% 
  summarize(x = mean(x), y = mean(y), array_origin = first(array_origin)) %>%
  mutate(label = c("(d)", "(b)", "(e)", "(c)"))

pbig <- predicted_densities_all %>% 
  filter(array_origin == "none", occasions == 20) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  geom_text(data = trap_labels, 
            inherit.aes = T, aes(colour = array_origin, label = label), size = 6) +
  # geom_point(data = detectors_df_all %>% filter(occasions == 20),
  #            inherit.aes = T, aes(colour = array_origin, shape = array_origin), size = 2) +
  geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 1, linetype = 2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_segment(data = all_segs, inherit.aes = F, fill = NA, size = 1,
               aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin)) + 
  scale_fill_viridis(direction = 1, limits = c(0,15)) +  
  scale_colour_manual(name="", 
                      values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                      breaks=c("15_15", "15_31", "27_15", "27_31")) + 
  scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

pbig

orgn <- c(15, 15)
orgn_str <- paste0(orgn[1],"_",orgn[2])

predicted_densities_all %>% 
  filter(array_origin == orgn_str, occasions == 20) %>%
  filter(x >= orgn[1] - buffersigmas*sigma, x <= orgn[1] + (4+buffersigmas)*sigma, y >= orgn[2]-buffersigmas*sigma, y <= orgn[2]+(6+buffersigmas)*sigma) %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  geom_point(data = detectors_df_all %>% filter(array_origin == orgn_str, occasions == 20), 
             inherit.aes = T, aes(colour = array_origin, shape = array_origin), size = 3) +
  geom_rect(data = common_area, inherit.aes = F, colour = "blue", fill = NA, size = 1, linetype = 2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_segment(data = all_segs %>% filter(array_origin == orgn_str), inherit.aes = F, fill = NA, size = 2,
               aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin)) + 
  scale_fill_viridis(limits = c(0,15)) +  
  scale_colour_manual(name = "", 
                      values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                      breaks=c("15_15", "15_31", "27_15", "27_31")) + 
  scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
  coord_equal() +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())


plot_mona_orig <- function(orgn){
  orgn_str <- paste0(orgn[1],"_",orgn[2])
  p <- predicted_densities_all %>% 
    filter(array_origin == "none", occasions == 20) %>%
    filter(x >= orgn[1] - buffersigmas*sigma, x <= orgn[1] + (4+buffersigmas)*sigma, y >= orgn[2]-buffersigmas*sigma, y <= orgn[2]+(6+buffersigmas)*sigma) %>%
    ggplot(aes(x, y)) + 
    geom_raster(aes(fill = value)) +
    geom_point(data = detectors_df_all %>% filter(array_origin == orgn_str, occasions == 20), 
               inherit.aes = T, aes(colour = array_origin, shape = array_origin), size = 3) +
    geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 1, linetype = 2,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    geom_segment(data = all_segs %>% filter(array_origin == orgn_str), inherit.aes = F, fill = NA, size = 2,
                 aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin)) + 
    scale_fill_viridis(direction = 1, limits = c(0,15)) +  
    scale_colour_manual(name = "", 
                        values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                        breaks=c("15_15", "15_31", "27_15", "27_31")) + 
    scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
    coord_equal() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) 
  return(p)
}
p3_true <- plot_mona_orig(orgn = c(15,15))
p3_true

plot_mona <- function(orgn, densities = predicted_densities_all){
  orgn_str <- paste0(orgn[1],"_",orgn[2])
  p <- densities %>% 
    filter(array_origin == orgn_str, occasions == 20) %>%
    filter(x >= orgn[1] - buffersigmas*sigma, x <= orgn[1] + (4+buffersigmas)*sigma, y >= orgn[2]-buffersigmas*sigma, y <= orgn[2]+(6+buffersigmas)*sigma) %>%
    ggplot(aes(x, y)) + 
    geom_raster(aes(fill = value)) +
    geom_point(data = detectors_df_all %>% filter(array_origin == orgn_str, occasions == 20), 
               inherit.aes = T, aes(colour = array_origin), pch = 4, size = 3) +
    geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 1, linetype = 2,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
    geom_segment(data = all_segs %>% filter(array_origin == orgn_str), inherit.aes = F, fill = NA, size = 2,
                 aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin)) + 
    scale_fill_viridis(direction = 1, limits = c(0,15)) +   
    scale_colour_manual(name = "", 
                        values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                        breaks=c("15_15", "15_31", "27_15", "27_31")) + 
    scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
    coord_equal() +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank()) 
  return(p)
}

p3_true <- plot_mona_orig(orgn = c(15,15))
p1_true <- plot_mona_orig(orgn = c(15,31))
p4_true <- plot_mona_orig(orgn = c(27,15))
p2_true <- plot_mona_orig(orgn = c(27,31))

p3 <- plot_mona(orgn = c(15,15))
p1 <- plot_mona(orgn = c(15,31))
p4 <- plot_mona(orgn = c(27,15))
p2 <- plot_mona(orgn = c(27,31))

pp <- pbig + p1 + p2 + p3 + p4 + plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(',
                  tag_sep = '', tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        plot.tag.position = "bottom")

pp

ggsave("paper/mona_torch_higheffort.png", pp, width=12, height=3, dpi = 600)
