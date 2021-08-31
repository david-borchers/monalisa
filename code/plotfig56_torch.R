## Plot (a) shows the true AC densities. Plots (b), (c), (d) (e) show the estimated realised AC surfaces 
## (averaged over 100 simulations) estimated using a 4×3 array placed at four different locations. The orange, 
## red, grey and yellow corner marks in plots (b) to (e) indicate the location of each of (b) to (e) in plot (a). 
## The white dashed box located in the centre of the Mona Lisa’s face in plot (a) is also shown in plots (b) to (d) 
## so that one can easily compare the predictions of the centre of the face from each arra 
## Currently Figure 6 (31/8/2021)

library(tidyverse)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(purrr)

load("output/mona_raw_outputs_100sim.RData")
load("output/simulated_densities_for_paper.RData")
load("output/capthist_summaries_100sim.RData")

# add a variable to the simulated data so that can rbind to results a bit later
simulated_densities_df <- simulated_densities_df %>% mutate(occasions = 1, array_origin = "none") %>%
  rbind(simulated_densities_df %>% mutate(occasions = 20, array_origin = "none"))

simulated_densities_df$value <- simulated_densities_df$dvalue * 7451
simulated_densities_df <- simulated_densities_df %>% dplyr::select(-dvalue)

# process the outputs
# # # for 1 sim
# predicted_densities_all <- fig34_results %>% purrr::map("predicted_densities") %>% map_df(bind_rows)
# detectors_df_all <- fig34_results %>% purrr::map("detectors_df") %>% map_df(bind_rows)
# average over 100 sim
predicted_densities_all <- fig34_results_100sim %>% purrr::map_depth(2, "predicted_densities") %>% map_df(bind_rows)
# # just one sim (first one)
# predicted_densities_onesim <- predicted_densities_all %>% 
#   group_by(x, y, covtype, occasions, array_size, array_spacing, array_origin, lambda0, sigma, n_pts) %>%
#   summarize(value = first(value)) %>% ungroup()
# average over sims
predicted_densities_all <- predicted_densities_all %>% 
  group_by(x, y, covtype, occasions, array_size, array_spacing, array_origin, lambda0, sigma, n_pts) %>%
  summarize(value = mean(countvalue)) %>% ungroup()

detectors_df_all <- fig34_results_100sim %>% purrr::map_depth(2, "detectors_df") %>% map_df(bind_rows)
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

###### if you want to do some standardization of the densities, do it here

# scale the plots to all sum to one
# predicted_densities_all <- predicted_densities_all %>% 
#   select(x, y, value, occasions, array_origin) %>%
#   group_by(occasions, array_origin) %>% 
#   mutate(value = value / sum(value)) %>%
#   ungroup()

# scale the plots to min 0 max 1
# minvalue = predicted_densities_all %>% filter(array_origin == "none", occasions == 20) %>% select(value) %>% min()
# maxvalue = predicted_densities_all %>% filter(array_origin == "none", occasions == 20) %>% select(value) %>% max()
# 
# predicted_densities_all <- predicted_densities_all %>% 
#   select(x, y, value, occasions, array_origin) %>%
#   group_by(occasions, array_origin) %>% 
#   mutate(value = (value - minvalue) / (maxvalue - minvalue)) %>%
#   ungroup()

# predicted_densities_all <- predicted_densities_all %>% 
#   select(x, y, value, occasions, array_origin) %>%
#   group_by(occasions, array_origin) %>% 
#   mutate(value = (value - min(value)) / (max(value) - min(value))) %>%
#   ungroup()

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
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

pbig

orgn <- c(15, 15)

orgn_str <- paste0(orgn[1],"_",orgn[2])

x1 <- predicted_densities_all %>% filter(array_origin == "none", occasions == 20)
x2 <- predicted_densities_all %>% filter(array_origin == "15_15", occasions == 20)
summary(x1$value) 
summary(x2$value) 
sum(x1$value) 
sum(x2$value) 

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
  scale_fill_viridis(limits = c(0,0.00202)) +  
  scale_colour_manual(name = "", 
                      values = c("15_15" = brew.cols[1], "15_31" = brew.cols[2], "27_15" = brew.cols[3], "27_31" = brew.cols[4]),
                      breaks=c("15_15", "15_31", "27_15", "27_31")) + 
  scale_shape_manual(name = "", values = c("15_15" = 1, "15_31" = 2, "27_15" = 3, "27_31" = 4)) +
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

# p3 <- plot_mona(orgn = c(15,15), densities = predicted_densities_onesim)
# p1 <- plot_mona(orgn = c(15,31), densities = predicted_densities_onesim)
# p4 <- plot_mona(orgn = c(27,15), densities = predicted_densities_onesim)
# p2 <- plot_mona(orgn = c(27,31), densities = predicted_densities_onesim)

pbig | ((p1_true | p2_true | p3_true | p4_true) / (p1 | p2 | p3 | p4)) + plot_layout(nrow = 2, byrow = FALSE, width = c(2,1))

pbig | ((p1 + p2) / (p3 + p4)) + plot_layout(nrow = 2, byrow = TRUE)

pp <- pbig + p1 + p2 + p3 + p4 + plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(',
                  tag_sep = '', tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        plot.tag.position = "bottom")

pp

ggsave("paper/mona_torch_higheffort.png", pp, width=12, height=3, dpi = 600)

# END PLOT
# 
# 
# ###### creating the insets
# 
# predicted_densities_inset <- predicted_densities_all %>% 
#   filter(x >= 23, x <= 27, y >= 27, y <= 31) %>%
#   mutate(x = (x - 23) * 40 + ifelse(x == 23, 1, 0),
#          y = (y - 27) * 40 + ifelse(y == 27, 1, 0))
# 
# top5_insets <- predicted_densities_inset %>% 
#   group_by(array_origin, occasions) %>% 
#   slice_max(n = 5, order_by = value, with_ties = FALSE) %>%
#   mutate(rank = dense_rank(desc(value)),
#          x = (5 + 0.25 * x) - ifelse(x == 1, 0.25, 0),
#          y = (5 + 0.25 * y) - ifelse(y == 1, 0.25, 0)) %>%
#   ungroup() %>% arrange(array_origin, occasions, rank)
# 
# top5_insets <- top5_insets %>% mutate(array_origin = str_c("inset_",array_origin))
# 
# predicted_densities_inset <- predicted_densities_inset %>%
#   group_by(array_origin, occasions) %>%
#   expand(x = 1:200, y = 1:200) %>%
#   left_join(predicted_densities_inset)
# 
# predicted_densities_inset <- predicted_densities_inset %>%
#   arrange(array_origin, occasions, x, y) %>%
#   group_by(array_origin, occasions, x) %>%
#   fill(value) %>%
#   ungroup() %>%
#   group_by(array_origin, occasions, y) %>%
#   fill(value) %>%
#   ungroup()
# 
# predicted_densities_inset <- predicted_densities_inset %>%
#   mutate(x = 1 + 49 * (x - 1) / 199,
#          y = 1 + 49 * (y - 1) / 199) 
# 
# predicted_densities_inset <- predicted_densities_inset %>% 
#   mutate(array_origin = str_c("inset_",array_origin))
# 
# ## top left inset (origin = c(15, 31))
# 
# orgn <- c(15, 31)
# predicted_densities_inset_tl <- predicted_densities_all %>% 
#   filter(x >= orgn - 4*sigma, x <= orgn + 7*sigma, y >= 31-4*sigma, y <= 31+8*sigma) 
# 
# predicted_densities_inset_tl %>% filter(occasions == 20) %>%
#   ggplot(aes(x, y)) + 
#   geom_raster(aes(fill = value)) +
#   scale_fill_viridis() 
# 
# 
# +
#   geom_raster(data = predicted_densities_inset %>% filter(occasions == 1), aes(fill = value)) +
#   scale_fill_viridis() +
#   facet_wrap(~ array_origin, nrow = 2) +
#   geom_point(data = detectors_df_all, inherit.aes = T, 
#              colour = "red", pch = 4, alpha = 0.5, size = 0.8) + 
#   geom_rect(data = focus_area, inherit.aes = F, colour = "blue", fill = NA,
#             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
#   geom_point(data = detectors_insets, inherit.aes = T, 
#              colour = "red", pch = 4, size = 4) + 
#   geom_rect(data = focus_area_insets, inherit.aes = F, colour = "blue", fill = NA,
#             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
#   geom_point(data = top5_insets %>% filter(occasions == 1), inherit.aes = T, 
#              colour = "white", pch = 2, size = 1.3) + 
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title=element_blank(),legend.position="none",
#         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),plot.background=element_blank())
# 
# 
# %>%
#   mutate(x = (x - 23) * 40 + ifelse(x == 23, 1, 0),
#          y = (y - 27) * 40 + ifelse(y == 27, 1, 0))
# 
# 
# 
# 
# 
# # combine insets into main dfs
# predicted_densities_all <- rbind(predicted_densities_all, predicted_densities_inset)
# 
# 
# predicted_densities_all %>% group_by(occasions, array_origin) %>% summarize(meanv=mean(value),
#                                                                             sumv = sum(value),
#                                                                             minv = min(value),
#                                                                             maxv = max(value))
# 
# 
# # make the boxes highlighting the inset area
# 
# all_grids <- unique(predicted_densities_all$array_origin)
# is_inset <- str_detect(all_grids, "inset")
# inset_grids <- all_grids[is_inset]
# main_grids <- all_grids[!is_inset]
# 
# focus_area <- data.frame(xmin = 23, xmax = 27, 
#                          ymin = 27, ymax = 31,
#                          array_origin = main_grids)
# 
# focus_area_insets <- data.frame(xmin = 5, xmax = 45, 
#                                 ymin = 5, ymax = 45, 
#                                 array_origin = inset_grids)
# 
# 
# # make the points indicating where detectors are in the inset plots
# 
# detectors_insets <- data.frame(x = c(5,5,45,45,NA), 
#                                y = c(5,45,5,45,NA),
#                                array_origin = inset_grids)
# 
# 
# # recode and reorder factor levels
# 
# array_levels <- c("none", "15_31", "27_31", "15_15","27_15",
#                   "inset_none", "inset_15_31", "inset_27_31", 
#                   "inset_15_15","inset_27_15")
# 
# # capture histories from get_capthist_summaries.r
# paster <- function(d,nd,na){
#   paste0(d,"\n",nd," detections\n(",na, " individuals)")
# }
# 
# capthist_labels <- pmap(list(d = ch_fig5_sim100$detector, nd = ch_fig5_sim100$n_detections, na = ch_fig5_sim100$n_animals), .f = paster) %>% unlist() 
# capthist_labels <- capthist_labels[c(3,4,1,2)]
# 
# array_labels <- c("Simulated", capthist_labels,
#                   "Simulated Inset", "Top Left Inset", "Top Right Inset", "Bottom Left Inset", 
#                   "Bottom Right Inset") 
# 
# predicted_densities_all$array_origin <- factor(predicted_densities_all$array_origin, 
#                                                levels = array_levels,
#                                                labels = array_labels)
# 
# detectors_df_all$array_origin <- factor(detectors_df_all$array_origin, 
#                                         levels = array_levels,
#                                         labels = array_labels)
# 
# detectors_insets$array_origin <- factor(detectors_insets$array_origin, 
#                                         levels = array_levels,
#                                         labels = array_labels)
# 
# focus_area$array_origin <- factor(focus_area$array_origin, 
#                                   levels = array_levels,
#                                   labels = array_labels)
# 
# focus_area_insets$array_origin <- factor(focus_area_insets$array_origin, 
#                                          levels = array_levels,
#                                          labels = array_labels)
# 
# top5_insets$array_origin <- factor(top5_insets$array_origin, 
#                                    levels = array_levels,
#                                    labels = array_labels)
# 
# # split data frame into one for inset and one for the rest, for plotting
# predicted_densities_rest <- predicted_densities_all %>% filter(!str_detect(array_origin, "inset")) %>% droplevels()
# predicted_densities_inset <- predicted_densities_all %>% filter(str_detect(array_origin, "inset")) %>% droplevels()
# 
# # figure 3, 1 occasion
# 
# p1 <- predicted_densities_all %>% filter(occasions == 1) %>%
#   ggplot(aes(x, y)) + 
#   geom_raster(data = predicted_densities_rest %>% filter(occasions == 1), aes(fill = value)) +
#   geom_raster(data = predicted_densities_inset %>% filter(occasions == 1), aes(fill = value)) +
#   scale_fill_viridis() +
#   facet_wrap(~ array_origin, nrow = 2) +
#   geom_point(data = detectors_df_all, inherit.aes = T, 
#              colour = "red", pch = 4, alpha = 0.5, size = 0.8) + 
#   geom_rect(data = focus_area, inherit.aes = F, colour = "blue", fill = NA,
#             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
#   geom_point(data = detectors_insets, inherit.aes = T, 
#              colour = "red", pch = 4, size = 4) + 
#   geom_rect(data = focus_area_insets, inherit.aes = F, colour = "blue", fill = NA,
#             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
#   geom_point(data = top5_insets %>% filter(occasions == 1), inherit.aes = T, 
#              colour = "white", pch = 2, size = 1.3) + 
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title=element_blank(),legend.position="none",
#         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),plot.background=element_blank())
# 
# p1
# 
# ggsave("paper/mona_torch_loweffort.png", p1, width=8, height=4.5, dpi = 600)
# 
# # figure 4, 20 occasions
# 
# capthist_labels <- pmap(list(d = ch_fig6_sim100$detector, nd = ch_fig6_sim100$n_detections, na = ch_fig6_sim100$n_animals), .f = paster) %>% unlist() 
# capthist_labels <- capthist_labels[c(3,4,1,2)]
# 
# array_labels_p2 <- c("Simulated", capthist_labels,
#                      "Simulated Inset", "Top Left Inset", "Top Right Inset", "Bottom Left Inset", 
#                      "Bottom Right Inset") 
# 
# levels(predicted_densities_all$array_origin) <- array_labels_p2
# levels(detectors_df_all$array_origin) <- array_labels_p2
# levels(detectors_insets$array_origin) <- array_labels_p2
# levels(focus_area$array_origin) <- array_labels_p2
# levels(focus_area_insets$array_origin) <- array_labels_p2
# levels(top5_insets$array_origin) <- array_labels_p2
# 
# predicted_densities_all$array_origin <- droplevels(predicted_densities_all$array_origin)
# detectors_df_all$array_origin <- droplevels(detectors_df_all$array_origin)
# detectors_insets$array_origin <- droplevels(detectors_insets$array_origin)
# focus_area$array_origin <- droplevels(focus_area$array_origin)
# focus_area_insets$array_origin <- droplevels(focus_area_insets$array_origin)
# top5_insets$array_origin <- droplevels(top5_insets$array_origin)
# 
# # split data frame into one for inset and one for the rest, for plotting
# predicted_densities_rest <- predicted_densities_all %>% filter(!str_detect(array_origin, "inset")) %>% droplevels()
# predicted_densities_inset <- predicted_densities_all %>% filter(str_detect(array_origin, "inset")) %>% droplevels()
# 
# p2 <- predicted_densities_all %>% filter(occasions == 20) %>%
#   ggplot(aes(x, y)) + 
#   geom_raster(data = predicted_densities_rest %>% filter(occasions == 20), aes(fill = value)) +
#   geom_raster(data = predicted_densities_inset %>% filter(occasions == 20), aes(fill = value)) +
#   scale_fill_viridis() +
#   facet_wrap(~ array_origin, nrow = 2) +
#   geom_point(data = detectors_df_all, inherit.aes = T, 
#              colour = "red", pch = 4, alpha = 0.5, size = 0.8) + 
#   geom_rect(data = focus_area, inherit.aes = F, colour = "blue", fill = NA,
#             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
#   geom_point(data = detectors_insets, inherit.aes = T, 
#              colour = "red", pch = 4, size = 4) + 
#   geom_rect(data = focus_area_insets, inherit.aes = F, colour = "blue", fill = NA,
#             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
#   geom_point(data = top5_insets %>% filter(occasions == 20), inherit.aes = T, 
#              colour = "white", pch = 2, size = 1.3) + 
#   theme(axis.line=element_blank(),axis.text.x=element_blank(),
#         axis.text.y=element_blank(),axis.ticks=element_blank(),
#         axis.title=element_blank(),legend.position="none",
#         panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),plot.background=element_blank())
# 
# ggsave("paper/mona_torch_higheffort.png", p2, width=8, height=4.5, dpi = 600)
