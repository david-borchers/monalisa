library(dplyr)
library(secr)

#load("output/simulated_capturehistories.RData")

load("output/mona_raw_outputs_100sim.RData")

ret_ch <- function(x, nocc){
  
  chs <- summary(x)
  
  data.frame(det_x = chs$trapsum$xrange[1],
             det_y = chs$trapsum$yrange[1],
             occasion = nocc,
             n_animals = chs$counts[4,nocc],
             n_detections = sum(chs$counts[6,1:nocc]))
}

ch_fig34_all <- data.frame(det_x = as.numeric(),
                           det_y = as.numeric(),
                           occasion = as.numeric(),
                           n_animals = as.numeric(),
                           n_detections = as.numeric(),
                           n_moves = as.numeric(),
                           sim_id = as.numeric())

for(i in 1:100){
  for(j in 1:8){
    chs <- summary(fig34_results_100sim[[i]][[j]]$capture_history, moves = TRUE)
    nocc <- ncol(chs$counts) - 1
    ch_fig34_all <- rbind(ch_fig34_all, 
                          data.frame(det_x = chs$trapsum$xrange[1],
                                     det_y = chs$trapsum$yrange[1],
                                     occasion = nocc,
                                     n_animals = chs$counts[4,nocc],
                                     n_detections = sum(chs$counts[6,1:nocc]),
                                     n_moves = sum(chs$moves$peranimal * 0:(length(chs$moves$peranimal)-1)),
                                     sim_id = i))
  }
}

# capture histories for Figure 3
ch_fig3_sim100 <- ch_fig34_all %>% filter(occasion == 20) %>% group_by(det_x, det_y, occasion) %>%
  summarize(n_animals = round(mean(n_animals), 0), n_detections = round(mean(n_detections), 0), n_moves = round(mean(n_moves), 0)) %>%
  ungroup() %>% 
  mutate(detector = paste(ifelse(det_y == 15, "Bottom", "Top"),
                          ifelse(det_x == 15, "Left", "Right")))

ch_fig3_firstsim <- ch_fig34_all %>% filter(occasion == 20) %>% group_by(det_x, det_y, occasion) %>%
  summarize(n_animals = first(n_animals), n_detections = first(n_detections), n_moves = first(n_moves)) %>%
  ungroup() %>% 
  mutate(detector = paste(ifelse(det_y == 15, "Bottom", "Top"),
                          ifelse(det_x == 15, "Left", "Right")))

# capture histories for Figure 4
ch_fig4_sim100 <- ch_fig34_all %>% filter(occasion == 1) %>% group_by(det_x, det_y, occasion) %>%
  summarize(n_animals = round(mean(n_animals), 0), n_detections = round(mean(n_detections), 0), n_moves = round(mean(n_moves), 0)) %>%
  ungroup() %>% 
  mutate(detector = paste(ifelse(det_y == 15, "Bottom", "Top"),
                          ifelse(det_x == 15, "Left", "Right")))

ch_fig4_firstsim <- ch_fig34_all %>% filter(occasion == 1) %>% group_by(det_x, det_y, occasion) %>%
  summarize(n_animals = first(n_animals), n_detections = first(n_detections), n_moves = first(n_moves)) %>%
  ungroup() %>% 
  mutate(detector = paste(ifelse(det_y == 15, "Bottom", "Top"),
                          ifelse(det_x == 15, "Left", "Right")))

# # capture histories for Figure 5
# ch_fig5 <- map_df(.x = capturehistory_lots, .f = ret_ch, nocc = 1)
# ch_fig5 <- ch_fig5 %>% mutate(detector = paste(ifelse(det_y == 15, "Bottom", "Top"),
#                                                ifelse(det_x == 15, "Left", "Right")))
# 
# # capture histories for Figure 6
# ch_fig6 <- map_df(.x = capturehistory_lots, .f = ret_ch, nocc = 20)
# ch_fig6 <- ch_fig6 %>% mutate(detector = paste(ifelse(det_y == 15, "Bottom", "Top"),
#                                                ifelse(det_x == 15, "Left", "Right")))

# capture history for Figure 7 is as for Figure 5 (1 occasion, many animals)

# capture histories for Figure 6 and 7
ch_fig67_all <- data.frame(det_x = as.numeric(),
                           det_y = as.numeric(),
                           occasion = as.numeric(),
                           covtype = as.character(),
                           n_animals = as.numeric(),
                           n_detections = as.numeric(),
                           n_moves = as.numeric(),
                           sim_id = as.numeric())

for(i in 1:100){
  for(j in 1:12){
    chs <- summary(fig67_results_100sim[[i]][[j]]$capture_history, moves = TRUE)
    nocc <- ncol(chs$counts) - 1
    covtype <- fig67_results_100sim[[i]][[j]]$predicted_densities$covtype[1]
    ch_fig67_all <- rbind(ch_fig67_all, 
                          data.frame(det_x = chs$trapsum$xrange[1],
                                     det_y = chs$trapsum$yrange[1],
                                     occasion = nocc,
                                     covtype = covtype,
                                     n_animals = chs$counts[4,nocc],
                                     n_detections = sum(chs$counts[6,1:nocc]),
                                     n_moves = sum(chs$moves$peranimal * 0:(length(chs$moves$peranimal)-1)),
                                     sim_id = i))
  }
}

# capture histories for Figure 6 & 7
ch_fig6_sim100 <- ch_fig67_all %>% group_by(det_x, det_y, occasion) %>%
  summarize(n_animals = round(mean(n_animals), 0), n_detections = round(mean(n_detections), 0), n_moves = round(mean(n_moves), 0)) %>%
  ungroup() %>% 
  mutate(detector = paste(ifelse(det_y == 15, "Bottom", "Top"),
                          ifelse(det_x == 15, "Left", "Right")))

ch_fig6_firstsim <- ch_fig67_all %>% group_by(det_x, det_y, occasion) %>%
  summarize(n_animals = first(n_animals), n_detections = first(n_detections), n_moves = first(n_moves)) %>%
  ungroup() %>% 
  mutate(detector = paste(ifelse(det_y == 15, "Bottom", "Top"),
                          ifelse(det_x == 15, "Left", "Right")))

# 
# 
# 
# ch_fig8a <- ret_ch(x = capturehistory_few, nocc = 1)
# ch_fig8b <- ret_ch(x = capturehistory_few, nocc = 3)
# ch_fig8c <- ret_ch(x = capturehistory_few, nocc = 10)
# ch_fig8d <- ret_ch(x = capturehistory_few, nocc = 20)
# ch_fig8 <- rbind(ch_fig8a, ch_fig8b, ch_fig8c, ch_fig8d)
# rm(ch_fig8a, ch_fig8b, ch_fig8c, ch_fig8d)

# capture history for Figure 9 is as for Figure 8 

save(ch_fig3_sim100, ch_fig4_sim100, ch_fig6_sim100, ch_fig3_firstsim, ch_fig4_firstsim, ch_fig6_firstsim, file = "output/capthist_summaries_100sim.RData")
