## Expected activity center surfaces estimated from a single survey using a model with density a function of one of 
## two simulated spatially-varying covariates.
## Currently Figure 7 (21/11/2021)

library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)
library(viridis)
library(RColorBrewer)
library(patchwork)
library(purrr)

load("output/mona_raw_outputs.RData")
load("output/mona_inputs.RData")

# process the covariates
predicted_densities_covs <- mona_df %>% dplyr::select(x,y,Dgood_bigD,Dblur_bigD) %>% 
  pivot_longer(cols=c(Dgood_bigD,Dblur_bigD), names_to = "covtype") %>% arrange(x,y) %>%
  mutate(value = value/10000,
         array_origin = "none") %>%
  mutate(covtype = str_remove(covtype, "_bigD"))

# process the outputs
predicted_densities_all <- res_expected_acd_many %>% purrr::map("predicted_densities") %>% map_df(bind_rows)
detectors_df_all <- res_expected_acd_many %>% purrr::map("detectors_df") %>% map_df(bind_rows)

# change covariate variable names to be consistent with mona_inputs
predicted_densities_all <- predicted_densities_all %>% mutate(covtype = str_remove(covtype, "D~log\\(")) %>% mutate(covtype = str_remove(covtype, "_bigD\\)"))
detectors_df_all <- detectors_df_all %>% mutate(covtype = str_remove(covtype, "D~log\\(")) %>% mutate(covtype = str_remove(covtype, "_bigD\\)"))

# choose covariates we want 
predicted_densities_all <- predicted_densities_all %>% filter(covtype %in% c("Dgood", "Dblur"))
detectors_df_all <- detectors_df_all %>% filter(covtype %in% c("Dgood", "Dblur"))

# choose variables we need
predicted_densities_all <- predicted_densities_all %>% 
  select(x, y, value = prob_ac, covtype, occasions, array_origin) %>%
  mutate(value = value / 10000) 

predicted_densities_all %>% group_by(covtype, array_origin) %>% summarize(mv = mean(value))

# scale the covariate plots to have the same mean as the density plots
predicted_densities_covs <- predicted_densities_covs %>%
  group_by(covtype, array_origin) %>%
  mutate(value = value * mean(predicted_densities_all$value) / mean(predicted_densities_covs$value)) %>%
  ungroup()

# # scale the plots to all sum to one
# predicted_densities_all <- predicted_densities_all %>%
#   group_by(covtype, array_origin) %>%
#   mutate(value = value / sum(value)) %>%
#   ungroup()

sigma <- detectors_df_all$sigma[1]
buffersigmas <- 4
orgn <- data.frame(x = c(15, 15, 27, 27), y = c(15, 31, 15, 31))
yy <- data.frame(array_origin = paste0(orgn$x,"_",orgn$y), xmin = orgn$x - buffersigmas*sigma, xmax = orgn$x + (4+buffersigmas)*sigma, ymin = orgn$y-buffersigmas*sigma, ymax = pmin(50, orgn$y+(6+buffersigmas)*sigma))
yy <- yy %>% filter(array_origin %in% c("15_15", "27_31"))
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

brew.cols <- brewer.pal(6, "Accent")[-c(1,5)]
#state.cols <- c(mycols, "#B8DE29FF") 

orgn <- c(15, 15)
orgn_str <- paste0(orgn[1],"_",orgn[2])

fill_max <- 15 # 9.1 is max in exp data, or 15 to make fill colour scale same as realised_ac_many plot

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
    scale_fill_viridis(direction = 1, limits = c(0,fill_max)) +   
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

p1good <- plot_mona(orgn = c(15,15), densities = predicted_densities_all %>% filter(covtype == "Dgood"))
p2good <- plot_mona(orgn = c(27,31), densities = predicted_densities_all %>% filter(covtype == "Dgood"))

p1blur <- plot_mona(orgn = c(15,15), densities = predicted_densities_all %>% filter(covtype == "Dblur"))
p2blur <- plot_mona(orgn = c(27,31), densities = predicted_densities_all %>% filter(covtype == "Dblur"))

trap_labels <- detectors_df_all %>% group_by(array_origin) %>% 
  summarize(x = mean(x), y = mean(y), array_origin = first(array_origin)) %>%
  mutate(label = c("(c)", "(b)"))

pgoodcov <- predicted_densities_covs %>% 
  filter(array_origin == "none", covtype == "Dgood") %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  # geom_point(data = detectors_df_all %>% filter(occasions == 20),
  #            inherit.aes = T, aes(colour = array_origin, shape = array_origin), size = 2) +
  geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 1, linetype = 2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_segment(data = all_segs, inherit.aes = F, fill = NA, size = 1,
               aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin)) + 
  geom_text(data = trap_labels, 
            inherit.aes = T, aes(colour = array_origin, label = label), size = 6) +
  scale_fill_viridis(direction = 1, limits = c(0,fill_max)) +  
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

pgoodcov

pblurcov <- predicted_densities_covs %>% 
  filter(array_origin == "none", covtype == "Dblur") %>%
  ggplot(aes(x, y)) + 
  geom_raster(aes(fill = value)) +
  # geom_point(data = detectors_df_all %>% filter(occasions == 20),
  #            inherit.aes = T, aes(colour = array_origin, shape = array_origin), size = 2) +
  geom_rect(data = common_area, inherit.aes = F, colour = "white", fill = NA, size = 1, linetype = 2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) +
  geom_segment(data = all_segs, inherit.aes = F, fill = NA, size = 1,
               aes(x = xmin, xend = xmax, y = ymin, yend = ymax, colour = array_origin)) + 
  geom_text(data = trap_labels %>% mutate(label = c("(f)", "(e)")), 
            inherit.aes = T, aes(colour = array_origin, label = label), size = 6) +
  scale_fill_viridis(direction = 1, limits = c(0,fill_max)) +  
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


pp <- (pgoodcov | p2good | p1good) / (pblurcov | p2blur | p1blur) + plot_layout(nrow = 2) + 
  plot_annotation(tag_levels = 'a',
                  tag_prefix = '(',
                  tag_sep = '', tag_suffix = ')') & 
  theme(plot.tag = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        plot.tag.position = "bottom")

pp

ggsave("paper/mona_covariates.png", pp, width=7.5, height=5.8, dpi = 600)


