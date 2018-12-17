## plotting figure 3 & 4: activity center surfaces estimated 
# using different arrays after 1 (fig 3) and 20 (fig 4) sampling occasions.

library(tidyverse)
library(viridis)

load("output/mona_raw_outputs.RData")
load("output/simulated_densities_for_paper.RData")
load("output/capthist_summaries.RData")

# add a variable to the simulated data so that can rbind to results a bit later
simulated_densities_df <- simulated_densities_df %>% mutate(occasions = 1, array_origin = "none") %>%
  rbind(simulated_densities_df %>% mutate(occasions = 20, array_origin = "none"))

# process the outputs
predicted_densities_all <- fig34_results %>% purrr::map("predicted_densities") %>% map_df(bind_rows)
detectors_df_all <- fig34_results %>% purrr::map("detectors_df") %>% map_df(bind_rows)

# change covariate variable names to be consistent with mona_inputs
predicted_densities_all <- predicted_densities_all %>% mutate(covtype = str_remove(covtype, "D~"))
detectors_df_all <- detectors_df_all %>% mutate(covtype = str_remove(covtype, "D~"))

# select only the variables we need
predicted_densities_all <- predicted_densities_all %>% select(x, y, value, occasions, array_origin)

# combine with simulated data
predicted_densities_all <- rbind(predicted_densities_all, simulated_densities_df)


###### if you want to do some standardization of the densities, do it here

# scale the plots to all sum to one
# predicted_densities_all <- predicted_densities_all %>% 
#   select(x, y, value, occasions, array_origin) %>%
#   group_by(occasions, array_origin) %>% 
#   mutate(value = value / sum(value)) %>%
#   ungroup()

# scale the plots to min 0 max 1
predicted_densities_all <- predicted_densities_all %>% 
  select(x, y, value, occasions, array_origin) %>%
  group_by(occasions, array_origin) %>% 
  mutate(value = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()


###### creating the insets

predicted_densities_inset <- predicted_densities_all %>% 
  filter(x >= 23, x <= 27, y >= 27, y <= 31) %>%
  mutate(x = (x - 23) * 40 + ifelse(x == 23, 1, 0),
         y = (y - 27) * 40 + ifelse(y == 27, 1, 0))

top5_insets <- predicted_densities_inset %>% 
  group_by(array_origin, occasions) %>% 
  top_n(5, wt = value) %>%
  mutate(rank = dense_rank(desc(value)),
         x = (5 + 0.25 * x) - ifelse(x == 1, 0.25, 0),
         y = (5 + 0.25 * y) - ifelse(y == 1, 0.25, 0)) %>%
  ungroup()

top5_insets <- top5_insets %>% mutate(array_origin = str_c("inset_",array_origin))

predicted_densities_inset <- predicted_densities_inset %>%
  group_by(array_origin, occasions) %>%
  expand(x = 1:200, y = 1:200) %>%
  left_join(predicted_densities_inset)

predicted_densities_inset <- predicted_densities_inset %>%
  arrange(array_origin, occasions, x, y) %>%
  group_by(array_origin, occasions, x) %>%
  fill(value) %>%
  ungroup() %>%
  group_by(array_origin, occasions, y) %>%
  fill(value) %>%
  ungroup()

predicted_densities_inset <- predicted_densities_inset %>%
  mutate(x = 1 + 49 * (x - 1) / 199,
         y = 1 + 49 * (y - 1) / 199) 

predicted_densities_inset <- predicted_densities_inset %>% 
  mutate(array_origin = str_c("inset_",array_origin))

# combine insets into main dfs
predicted_densities_all <- rbind(predicted_densities_all, predicted_densities_inset)


predicted_densities_all %>% group_by(occasions, array_origin) %>% summarize(meanv=mean(value),
                                                                            sumv = sum(value),
                                                                            minv = min(value),
                                                                            maxv = max(value))


# make the boxes highlighting the inset area

all_grids <- unique(predicted_densities_all$array_origin)
is_inset <- str_detect(all_grids, "inset")
inset_grids <- all_grids[is_inset]
main_grids <- all_grids[!is_inset]

focus_area <- data.frame(xmin = 23, xmax = 27, 
                         ymin = 27, ymax = 31,
                         array_origin = main_grids)

focus_area_insets <- data.frame(xmin = 5, xmax = 45, 
                                ymin = 5, ymax = 45, 
                                array_origin = inset_grids)


# make the points indicating where detectors are in the inset plots

detectors_insets <- data.frame(x = c(5,5,45,45,NA), 
                               y = c(5,45,5,45,NA),
                               array_origin = inset_grids)


# recode and reorder factor levels

array_levels <- c("none", "15_31", "27_31", "15_15","27_15",
                  "inset_none", "inset_15_31", "inset_27_31", 
                  "inset_15_15","inset_27_15")

# capture histories from get_capthist_summaries.r
paster <- function(d,nd,na){
  paste0(d,"\n",nd," detections\n(",na, " individuals)")
}

capthist_labels <- pmap(list(d = ch_fig5$detector, nd = ch_fig5$n_detections, na = ch_fig5$n_animals), .f = paster) %>% unlist() 
capthist_labels <- capthist_labels[c(3,4,1,2)]

array_labels <- c("Simulated", capthist_labels,
                  "Simulated Inset", "Top Left Inset", "Top Right Inset", "Bottom Left Inset", 
                  "Bottom Right Inset") 

predicted_densities_all$array_origin <- factor(predicted_densities_all$array_origin, 
                                          levels = array_levels,
                                          labels = array_labels)

detectors_df_all$array_origin <- factor(detectors_df_all$array_origin, 
                                   levels = array_levels,
                                   labels = array_labels)

detectors_insets$array_origin <- factor(detectors_insets$array_origin, 
                                        levels = array_levels,
                                        labels = array_labels)

focus_area$array_origin <- factor(focus_area$array_origin, 
                                        levels = array_levels,
                                        labels = array_labels)

focus_area_insets$array_origin <- factor(focus_area_insets$array_origin, 
                                        levels = array_levels,
                                        labels = array_labels)

top5_insets$array_origin <- factor(top5_insets$array_origin, 
                                         levels = array_levels,
                                         labels = array_labels)

# split data frame into one for inset and one for the rest, for plotting
predicted_densities_rest <- predicted_densities_all %>% filter(!str_detect(array_origin, "inset")) %>% droplevels()
predicted_densities_inset <- predicted_densities_all %>% filter(str_detect(array_origin, "inset")) %>% droplevels()

# figure 3, 1 occasion

p1 <- predicted_densities_all %>% filter(occasions == 1) %>%
  ggplot(aes(x, y)) + 
  geom_raster(data = predicted_densities_rest %>% filter(occasions == 1), aes(fill = value)) +
  geom_raster(data = predicted_densities_inset %>% filter(occasions == 1), aes(fill = value)) +
  scale_fill_viridis() +
  facet_wrap(~ array_origin, nrow = 2) +
  geom_point(data = detectors_df_all, inherit.aes = T, 
             colour = "red", pch = 4, alpha = 0.5, size = 0.8) + 
  geom_rect(data = focus_area, inherit.aes = F, colour = "blue", fill = NA,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_point(data = detectors_insets, inherit.aes = T, 
             colour = "red", pch = 4, size = 4) + 
  geom_rect(data = focus_area_insets, inherit.aes = F, colour = "blue", fill = NA,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_point(data = top5_insets %>% filter(occasions == 1), inherit.aes = T, 
             colour = "white", pch = 2, size = 1.3) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("paper/mona_torch_loweffort.png", p1, width=8, height=4.5, dpi = 600)

# figure 4, 20 occasions

capthist_labels <- pmap(list(d = ch_fig6$detector, nd = ch_fig6$n_detections, na = ch_fig6$n_animals), .f = paster) %>% unlist() 
capthist_labels <- capthist_labels[c(3,4,1,2)]

array_labels_p2 <- c("Simulated", capthist_labels,
                  "Simulated Inset", "Top Left Inset", "Top Right Inset", "Bottom Left Inset", 
                  "Bottom Right Inset") 

levels(predicted_densities_all$array_origin) <- array_labels_p2
levels(detectors_df_all$array_origin) <- array_labels_p2
levels(detectors_insets$array_origin) <- array_labels_p2
levels(focus_area$array_origin) <- array_labels_p2
levels(focus_area_insets$array_origin) <- array_labels_p2
levels(top5_insets$array_origin) <- array_labels_p2

predicted_densities_all$array_origin <- droplevels(predicted_densities_all$array_origin)
detectors_df_all$array_origin <- droplevels(detectors_df_all$array_origin)
detectors_insets$array_origin <- droplevels(detectors_insets$array_origin)
focus_area$array_origin <- droplevels(focus_area$array_origin)
focus_area_insets$array_origin <- droplevels(focus_area_insets$array_origin)
top5_insets$array_origin <- droplevels(top5_insets$array_origin)

# split data frame into one for inset and one for the rest, for plotting
predicted_densities_rest <- predicted_densities_all %>% filter(!str_detect(array_origin, "inset")) %>% droplevels()
predicted_densities_inset <- predicted_densities_all %>% filter(str_detect(array_origin, "inset")) %>% droplevels()

p2 <- predicted_densities_all %>% filter(occasions == 20) %>%
  ggplot(aes(x, y)) + 
  geom_raster(data = predicted_densities_rest %>% filter(occasions == 20), aes(fill = value)) +
  geom_raster(data = predicted_densities_inset %>% filter(occasions == 20), aes(fill = value)) +
  scale_fill_viridis() +
  facet_wrap(~ array_origin, nrow = 2) +
  geom_point(data = detectors_df_all, inherit.aes = T, 
             colour = "red", pch = 4, alpha = 0.5, size = 0.8) + 
  geom_rect(data = focus_area, inherit.aes = F, colour = "blue", fill = NA,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_point(data = detectors_insets, inherit.aes = T, 
             colour = "red", pch = 4, size = 4) + 
  geom_rect(data = focus_area_insets, inherit.aes = F, colour = "blue", fill = NA,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_point(data = top5_insets %>% filter(occasions == 20), inherit.aes = T, 
             colour = "white", pch = 2, size = 1.3) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),
        axis.text.y=element_blank(),axis.ticks=element_blank(),
        axis.title=element_blank(),legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

ggsave("paper/mona_torch_higheffort.png", p2, width=8, height=4.5, dpi = 600)
